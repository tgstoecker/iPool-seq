#!/usr/bin/env python
# trim.tag.py, Copyright 2016, 2017 Florian G. Pflug
# 
# This file is part of the iPool-Seq Analysis Pipeline
#
# The iPool-Seq Analysis Pipeline is free software: you can redistribute it
# and/or modify it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# The iPool-Seq Analysis Pipeline is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with the iPool-Seq Analysis Pipeline.  If not, see
# <http://www.gnu.org/licenses/

# *****************************************************************************
# Implements "UMI extraction & technical sequence removal"
# ./trim.tag.py <cpu cores> <in read1> <in read2> <out read1> <out read2>
# Input and output files are read/written as gzip-compressed FastQ files
# *****************************************************************************
import sys
import regex
import itertools
import io
import gzip
import string
import distance
import pprint
import signal
import multiprocessing
from copy import copy
from recordtype import recordtype
from Bio.SeqIO.QualityIO import FastqGeneralIterator

DNA_REVCOMP = string.maketrans("ACGT", "TGCA")
def revcomp(s):
  return s.translate(DNA_REVCOMP)[::-1]

ALLOW_N_PATTERN = regex.compile('([ACGT])', regex.V1)
def allow_N(s):
  return ALLOW_N_PATTERN.sub('[\\1N]', s)

def allow_partial(s):
  return ('($|'.join(list(s)) + (')' * (len(s) - 1)))

# *** Patterns used to trim read 1 and read2.
#
# Read1 is assumed to have the following form
#   <filler> <12bp barcode> <one of TRIM_R1> || <rest>
# Read2 is assumed to have the form
#   <filler> <one of TRIM_R2> || <rest>
# Everthing before || is removed from the reads. The length of the filler
# can be 0 up to TRIM_MAX_5PFILLER bases. The fixed parts TRIM_R1 resp. TRIM_R2
# are allowed to contain up to ALLOWED_MISSMATCHES missmatches. N is treated as a
# *match* within these regions. This is necessary because these fixed regions cause
# problems for basecalling, are often contain manY Ns.
#
# Within these patterns, the following match group names are used:
#   trim: While trimmed part, i.e. everthing before ||
#    i5p: <filler>, i.e. 0 up to TRIM_MAX_5PFILLER bases at the 5' end
#     bc: Barcode (currently only on read 1)
#   t<i>: i-th sequence in TRIM_R1 resp. TRIM_R2. Can be used to decide which of the
#         sequences was found in a read.
#    seq: The surviving part of the read, i.e. everthing after ||
# 
TRIM_R1 = [ 'AGATGTGTATAAGAGACAG' ]
TRIM_R2 = [ 'CTGTGGTATCCTGTGGCGATC', 'CTGTGGTATCCTGTGGCGTGAGTGGC' ]
TRIM_MAX_MM = 4
TRIM_MAX_5PFILLER = 1
R1_PATTERN = regex.compile('^(?P<trim>(?P<i5p>[ACGTN]{0,%d})(?P<bc>[AGCT]{12})(%s))(?P<seq>[ACGTN]*)$' % (
  TRIM_MAX_5PFILLER,
  "|".join("(?P<t%d>(%s){s<=%d})" % (i+1, allow_N(t), TRIM_MAX_MM) for i, t in enumerate(TRIM_R1))
), regex.BESTMATCH | regex.V1)
R2_PATTERN = regex.compile('^(?P<trim>(?P<i5p>[ACGTN]{0,%d})(%s))(?P<seq>[ACGTN]*)$' % (
  TRIM_MAX_5PFILLER,
  "|".join("(?P<t%d>(%s){s<=%d})" % (i+1, allow_N(t), TRIM_MAX_MM) for i, t in enumerate(TRIM_R2))
), regex.BESTMATCH | regex.V1)
Ri_PATTERN = [ R1_PATTERN, R2_PATTERN ]
R2_PATTERN_RC = regex.compile("%s" % (
  "|".join('(?P<t%d>(%s){s<=%d})' % (i, revcomp(t), TRIM_MAX_MM) for i, t in enumerate(TRIM_R2))
), regex.BESTMATCH | regex.V1)
#R1_PATTERN = regex.compile('^(?P<trim>(?P<bc>[AGCT]{12})(%s){s<=%d})(?P<seq>[ACGTN]*)$' % (allow_N(TRIM_R1), TRIM_MAX_MM), regex.V1)
#R2_PATTERN = regex.compile('^(?P<trim>(?P<t1>%s){s<=%d}|(?P<t2>%s){s<=%d})(?P<seq>[ACGTN]*)$' % (allow_N(TRIM_R2[0]), R2_MISSMATCHES, allow_N(TRIM_R2[1]), R2_MISSMATCHES), regex.V1)
#R2_PATTERN_RC = regex.compile('(?P<t1>%s){s<=%d}|(?P<t2>%s){s<=%d}' % (revcomp(TRIM_R2[0]), TRIM_MAX_MM, revcomp(TRIM_R2[1]), TRIM_MAX_MM), regex.BESTMATCH | regex.V1)

# *** Overlap detection.
# The last OVERLAP_SEED_LENTH bases from each read and searched for in the mate, allowing
# for OVERLAP_SEED_MISSMATCHES missmatches. This seed is extended to a gapless alignment,
# and if that alignment contains no more than OVERLAP_MISSMATCHES, it is accepted. Otherwisel
# the reads are assumed not to overlap
OVERLAP_SEED_LENGTH=10
OVERLAP_SEED_MISSMATCHES=2
OVERLAP_IDENTITY=0.9

# *** Readname pattern. Group "bn" is the part that is the same for read1 and read2
NAME_PATTERN = regex.compile('^(?P<bn>.*)/[12]$', regex.V1)

# Holds the statistics collected while processings pairs
Stats = recordtype('Stats', 'invalid_r1 invalid_r2 emptymate overlap', default=0)

# Default
DEBUG = False

# Stored a pair of two (Illumina) reads, i.e. two reverse-complemented reads
# generally assumed to come from the same fragment of DNA (but from different
# strands, hence the different directions). Allows such a pair to be trimmed
# in various ways, and to be either aligned (i.e. the relativ position of
# the two reads is known) or unaligned. Alignment information is kept across
# trimming operations. Also tracks per-base qualities, and ensures that the
# qualities and bases are not 'shiften' relative' to each other when the
# sequence is trimmed.
class FwdRevPair:
  # Creates a pair of forward and reverse reads from sequences and qualities.
  # If idx1_seq2start is specified, the pair is aligned, i.e. it is assumed
  # that the reads came from the same DNA fragment, and that their relative
  # position is known. idx1_seq2start is the (0-base) index of of the first
  # base on read1 that is also contained in read2 (where it is the last base)
  def __init__(self, seq1, qual1, seq2, qual2, idx1_seq2start = None):
    self.__seq = [seq1, seq2]
    self.__qual = [qual1, qual2]
    if (len(self.__seq[0]) != len(self.__qual[0])) or (len(self.__seq[1]) != len(self.__qual[1])):
      raise ValueError('sequences and qualities must have the same length')
    self.__startof1_on_0 = None

  def __copy__(self):
    return FwdRevPair(self.seq(0), self.qual(0), self.seq(1), self.qual(1), self.pos_otherstart(0))

  # Returns the (0-based) index of the first base in sequence i (0 or 1)
  # that is also contained in ther other sequence (1-i).  
  def pos_otherstart(self, i):
    if self.__startof1_on_0 == None:
      return None
    return self.len(i) - self.len(0) + self.__startof1_on_0

  # Marks the pair as aligned. The (0-based) position <pos_on_i>
  # of the <i>-th sequence (0 or 1) is aligned to the position
  # <pos_on_other> on the other sequence (i.e. sequence <1-i>).
  def align(self, i, pos_on_i, pos_on_other):
    if i == 0:
      #  startof1_on_0   pos_on_i
      #       v              v
      #     |----------------:------------------>     (0)
      #       <--------------:---------------------|  (1)
      #       ^              ^
      #    last(1)     pos_on_other
      self.__startof1_on_0 = pos_on_i - (self.last(1) - pos_on_other)
    elif i == 1:
      #  startof1_on_0  pos_on_other
      #       v              v
      #     |----------------:------------------>     (0)
      #       <--------------:---------------------|  (1)
      #       ^              ^
      #    last(1)       pos_on_i
      self.__startof1_on_0 = pos_on_other - (self.last(1) - pos_on_i)
    else: raise IndexError('sequence index must be 0 or 1')
    # DEBUGGING CHECK
    if pos_on_other != self.aligned_to(i, pos_on_i):
      raise RuntimeError('internal inconsistency')

  # Find the (0-based) position on the sequence <1-i> that is
  # aligned to the (0-based) position <pos_on_i> of sequence <i>.
  def aligned_to(self, i, pos_on_i):
    return self.last(1-i) + self.pos_otherstart(i) - pos_on_i
  
  # Checks if the pair is aligned, if so returns True.
  def is_aligned(self):
    return self.__startof1_on_0 != None
  
  # Returns the lenth of sequence <i> (0 or 1)
  def len(self, i):
    return len(self.__seq[i])
  
  # Returns the last (0-based) position of sequence <i> (0 or 1)
  def last(self, i):
    return len(self.__seq[i]) - 1
  
  # Returns the sequence <i> (0 or 1)
  def seq(self, i):
    return self.__seq[i]

  # Returns the qualities for sequence <i> (0 or 1)
  def qual(self, i):
    return self.__qual[i]
    
  # Returns the part of sequence <i> (0 or 1) that overlaps the
  # other sequence (i.e. sequence <1-i>).
  def overlap(self, i):
    return self.seq(i)[max(0, self.aligned_to(1-i, self.last(1-i))):min(self.len(i), self.aligned_to(1-i, -1))]

  # Removes the first <n> bases at the 5' end (i.e. beginning)
  # of sequence <i> (0 or 1)
  def trim_5p(self, i, n):
    n = min(n, self.len(i))
    self.__seq[i] = self.seq(i)[n:]
    self.__qual[i] = self.qual(i)[n:]
    if (i == 0) and self.is_aligned():
      self.__startof1_on_0 -= n  

  # Removes bases at the 3' end (i.e. ending) of sequence <i>
  # (0 or 1) such that n bases remain after trimming.
  def trim_3p(self, i, n):
    n = min(m, self.len(i))
    self.__seq[i] = self.seq(i)[:n]
    self.__qual[i] = self.qual(i)[:n]
    if (i == 1) and self.is_aligned():
      self.__startof1_on_0 += n  

  # Removes the bases at the 3' ends of both sequences that
  # extend beyond the start of the mate, i.e. for alignment
  #            |---------------------->  (1)
  #      <------------------------|      (2)
  # returns:
  #            |------------------|      (1)
  #            |------------------|      (2)
  def trim_3p_overhangs(self):
    if self.__startof1_on_0 != None:
      for i in [0, 1]:
        n = -self.pos_otherstart(1-i)
        if n > 0:
          self.__seq[i] = self.seq(i)[:-n]
          self.__qual[i] = self.qual(i)[:-n]
          if i == 1:
            self.__startof1_on_0 = 0

def find_overlap_seed(seq1, seq2):
  p = regex.compile("(%s){s<=%d}" % (revcomp(seq1[-OVERLAP_SEED_LENGTH:]), OVERLAP_SEED_MISSMATCHES), regex.BESTMATCH | regex.V1)
  return(p.search(seq2))

def align_if_similar(r, r_putative, stats):
  if r.is_aligned():
    return
  # Check if seed extends to a sensible alignment
  if r_putative.is_aligned():
    # Extract overlapping region from both reads, and compare their hamming distance.
    r_putative.trim_3p_overhangs()
    if distance.hamming(r_putative.overlap(0), revcomp(r_putative.overlap(1)), normalized=True) <= OVERLAP_IDENTITY:
      # Actual overlap, update r's alignment information
      r.align(0, 0, r_putative.aligned_to(0, 0))
      stats.overlap += 1

def align_using_overlap(r, stats):
  if r.is_aligned():
    return
  r_putative = copy(r)
  # Figure out if reads overlap by using the last 10 bases of each read as a seed
  pa = find_overlap_seed(r_putative.seq(1), r_putative.seq(0))
  if pa != None:
    # Found last OVERLAP_SEED_LENGTH bases of read 2 within read 1, with up to OVERLAP_SEED_MISSMATCHES missmatches
    r_putative.align(0, pa.start(), r_putative.last(1))
  else:
    pa = find_overlap_seed(r_putative.seq(0), r_putative.seq(1))
    if pa != None:
      # Found last OVERLAP_SEED_LENGTH bases of read 1 within read 2, with up to OVERLAP_SEED_MISSMATCHES missmatches
      r_putative.align(1, pa.start(), r_putative.last(0))
  # Align if overlapping region is similar enough
  align_if_similar(r, r_putative, stats)

def align_using_pattern(r, stats):
  if r.is_aligned():
    return
  # Figure out if reads overlap by searching for the (reverse-complemented)
  # read2 pattern in read1. WE ASSUME THAT THE PATTERN WAS ALREADY REMOVED
  # FROM READ2!!
  r_putative = copy(r)
  palign = R2_PATTERN_RC.search(r.seq(0), partial=True)
  if ((palign != None) and (palign.end() - palign.start() >= OVERLAP_SEED_LENGTH)):
    # Found read2 pattern within read1.
    r_putative.align(0, palign.start() - 1, 0)
  # Align if overlapping region is similar enough
  align_if_similar(r, r_putative, stats)
  
def process(input):
  # Stats
  stats = Stats()  

  # Extract and check read names
  if input[0][0] != input[1][0]:
    pn = [ None, None ]
    for i in [0, 1]:
      pn[i] = NAME_PATTERN.match(input[i][0])
      if pn[i] == None:
        raise ValueError('invalid query name %s', input[i][0])
    if pn[0].group('bn') != pn[1].group('bn'):
      raise ValueError('non-matching query names of mates, "%s" and "%s"' % (pn[1].group('bn'), pn[2].group('bn')))
    r_name = pn[1].group('bn')
  else:
    r_name = input[0][0]

  # Construct read pair
  r = FwdRevPair(input[0][1], input[0][2], input[1][1], input[1][2])

  # Attempt to align overlapping reads
  align_using_overlap(r, stats)

  # Match sequences against patterns and trim, and collect barcodes
  p = [ None, None]
  bc = []
  valid = True
  for i in [0, 1]:
    p[i] = Ri_PATTERN[i].match(r.seq(i)) if Ri_PATTERN[i] else None
    if p[i] != None:
      if use[i]:
        bc.append(p[i].group('bc'))
      r.trim_5p(i, len(p[i].group('trim')))
    else:
      if i == 0: stats.invalid_r1 += 1
      if i == 1: stats.invalid_r2 += 1
      valid = False
  if not valid:
    return (stats, None, None)
  
  # Attempt to align overlapping reads again we we didn't succeed earlier
  if not r.is_aligned():
    align_using_overlap(r, stats)
  
  # Attempt to align overlapping reads by searching for read 2's (reverse-
  # complemented) pattern in read1
  if not r.is_aligned():
    align_using_pattern(r, stats)

  # Remove 3' overhangs
  r.trim_3p_overhangs()

  # XXX: Deal with empty reads!!! Replace with a single N
  # Return transformed read pair 
  return (stats, ("%s|%s/1" % (r_name, ":".join(bc)), r.seq(0), r.qual(0)),
                 ("%s|%s/2" % (r_name, ":".join(bc)), r.seq(1), r.qual(1)))

def init_worker():
  signal.signal(signal.SIGINT, signal.SIG_IGN)
    
if __name__ == '__main__':
  # Command-line arguments
  cores = int(sys.argv[1])
  ifile1 = sys.argv[2]
  ifile2 = sys.argv[3]
  ofile1 = sys.argv[4]
  ofile2 = sys.argv[5]
  use = [ True, False ]
  print >> sys.stderr, "reading from %s and %s, writing to %s and %s, using %d cores" % (ifile1, ifile2, ofile1, ofile2, cores)

  # Creater workers
  pool = multiprocessing.Pool(cores, init_worker)
  try:
    # Open input files
    input1 = FastqGeneralIterator(gzip.open(ifile1, 'rb'))
    input2= FastqGeneralIterator(gzip.open(ifile2, 'rb'))

    # Open output files (separate ones for "first" and "second" reads).
    output1 = io.BufferedWriter(gzip.open(ofile1, 'wb', compresslevel=1), buffer_size=100) #1024*1024)
    output2 = io.BufferedWriter(gzip.open(ofile2, 'wb', compresslevel=1), buffer_size=100) #1024*1024)

    # Scan input files (separate ones for "first" and "second" reads).
    # The loop iterates over both in parallel. Each read is represented
    # by a triple (name, bases, qualities).
    written = 0
    skipped = 0
    index = 0
    stats = Stats()
    def log():
      print >> sys.stderr, "%d pairs read" % (written+skipped)
      print >> sys.stderr, "%d (%.2f%%) of pairs did overlap" % (stats.overlap, float(100 * stats.overlap) / float(written+skipped))
      if R1_PATTERN:
        print >> sys.stderr, "%d (%.2f%%) of pairs contains invalid 1st read" % (stats.invalid_r1, float(100 * stats.invalid_r1) / float(written+skipped))
      if R2_PATTERN:
        print >> sys.stderr, "%d (%.2f%%) of pairs contains invalid 2nd read" % (stats.invalid_r2, float(100 * stats.invalid_r2) / float(written+skipped))
      print >> sys.stderr, "%d (%.2f%%) of pairs ended up with an empty mate after trimming" % (stats.emptymate, float(100 * stats.emptymate) / float(written+skipped))
      print >> sys.stderr, "%d (%.2f%%) of pairs written to %s and %s" % (written, float(100 * written) / float(written+skipped), ofile1, ofile2)
      print >> sys.stderr, "%d (%.2f%%) of pairs skipped" % (skipped, float(100 * skipped) / float(written+skipped))
    for output in pool.imap(process, itertools.izip(input1, input2), chunksize=8192):
    #for output in map(process, itertools.izip(input1, input2)):
      # Split
      stats_delta, r1, r2 = output
      
      # Update stats
      stats.emptymate += stats_delta.emptymate
      stats.overlap += stats_delta.overlap
      stats.invalid_r1 += stats_delta.invalid_r1
      stats.invalid_r2 += stats_delta.invalid_r2
      
      # Write trimmed reads to output files
      if r1 != None and r2 != None:
        written += 1
        output1.write("@%s\n%s\n+\n%s\n" % r1)
        output2.write("@%s\n%s\n+\n%s\n" % r2)
      else:
        skipped += 1

      # Update statistics and report progress    
      if (written + skipped) % 100000 == 0:
        log()

    # Close output files and print statistics
    output1.close()
    output2.close()
    log()
  except KeyboardInterrupt:
    pool.terminate()
    pool.join()
