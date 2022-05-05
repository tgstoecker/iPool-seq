#!/usr/bin/env python
# umicounts.tag.py, Copyright 2016, 2017 Florian G. Pflug
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
# Implements "Correcting for sequencing errors"
# ./umicounts.tag.py <mapq threshold> <input bam> <output UMI table>
# Tabel is written as gzipped tab-separated file
# *****************************************************************************
import sys
import re
import pysam
import io
import gzip
from pprint import pprint
from collections import defaultdict
from itertools import chain
from umi_error_correction import Umi, UmiCount, ErrorCorrector
from pysamutils import read_pair_iterator

TAG_LENGTH = 12
EC_MAX_HAMMINGDISTANCE = 1
EC_MAX_MAPPINGDISTANCE = 3
EC_REPORT_LIMIT = 100
EC_REPORT_STEP = 500
EC_REPORT_BP = 50

# Handle reads belonging to a particular transcript fragment,
def process(seq, umi_counts_per_rg, output, do_report = True):
  for rg, umi_counts in umi_counts_per_rg.iteritems():
    sys.stderr.write('*** PROCESSING READ GROUP %s\n' % rg)
    
    # Create error-corrector and process UMIs
    ec = ErrorCorrector(TAG_LENGTH, EC_MAX_MAPPINGDISTANCE, EC_MAX_HAMMINGDISTANCE, True)
    ecumi_counts = ec.process(umi_counts)

    # Output read-count per tag-cluster
  #  sys.stderr.write('rg\tseq\tpos\tend\ttag\trawumis\tcount\tplus\tminus\n')
    for u, c in sorted(ecumi_counts.iteritems(), key=lambda u: (u[0].start, u[0].end, u[0].tag)):
  #    sys.stderr.write('%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n' % (rg, seq, u.start, u.end, u.tag, c.rawumis, c.total(), c.plus, c.minus))
      output.write('%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n' % (rg, seq, u.start, u.end, u.tag, c.rawumis, c.total(), c.plus, c.minus))

    sys.stderr.write('=== DONE PROCESSING READ GROUP %s\n' % rg)
    sys.stderr.flush()

# Command-line arguments
MIN_MQ = int(sys.argv[1])
ifile = sys.argv[2]
ofile = sys.argv[3]
sys.stderr.write('Skipping reads with mapping quality below %d\n' % MIN_MQ)

# Open input file
input = pysam.AlignmentFile(ifile, mode="rb")
input_pairs = read_pair_iterator(input.fetch(until_eof=True), skip_unpaired=False)
if input.header['HD']['SO'] != 'coordinate':
  raise RuntimeError("input file must be sorted by coordinate")
  
# Open output file
output = io.BufferedWriter(gzip.open(ofile, 'wb', compresslevel=1), buffer_size=64*1024)

# Statistics
seq = None
umi_counts_per_rg = defaultdict(lambda : defaultdict(lambda: UmiCount(0,0,1)))
processed = 0
skipped = defaultdict(int)

# Scan reads
output.write('rg\tseq\tpos\tend\ttag\trawumis\tcount\tplus\tminus\n')
for read_fwd, read_rev in input_pairs:
  # One of the reads in a pair can be missing
  read_m = read_fwd if read_fwd is not None else read_rev
  read_o = read_fwd if read_fwd is None else read_rev

  # If the sequence changed, process accumulated per-tag read counts
  # XXX: Don't wait for sequence to change, keep rolling deque of
  # UMIs instead and process once they're EC_MAX_MAPPINGDISTANCE
  # to the left of the current fragments..
  p_seq = input.getrname(read_m.reference_id)
  seq_changed = False
  if seq != p_seq:
    if seq is not None:
      process(seq, umi_counts_per_rg, output)
      seq_changed = True
    seq = p_seq
    umi_counts_per_rg = defaultdict(lambda : defaultdict(lambda: UmiCount(0,0,1)))

  # Progress report
  if seq_changed or (input_pairs.total + 1 % 100000 == 0):
    sys.stderr.write('PROGRESS: processed %d of %d pairs (%.2f%%)\n' % (processed, input_pairs.total, 100*float(processed)/input_pairs.total))
    for reason, count in chain(input_pairs.skipped.iteritems(), skipped.iteritems()):
      sys.stderr.write('PROGRESS:  skipped %d pairs (%.2f%%) because \'%s\'\n' % (count, 100*float(count)/input_pairs.total, reason))
    sys.stderr.flush()

  # Skip pairs with low avg mapping quality.
  mapq = (read_m.mapping_quality + read_o.mapping_quality) / 2.0 if read_o is not None else read_m.mapping_quality
  if mapq < MIN_MQ:
    skipped['low mapping quality'] += 1
    continue

  # Pair is accepted
  processed += 1

  # Get barcode and read group
  tag = read_m.get_tag('XT')
  if (read_o is not None) and (read_o.get_tag('XT') != tag):
    raise RuntimeError('inconsistent tags (%s and %s) in pair %s' % (tag, read_o.get_tag('XT'), read_m.queryname))
  rg = read_m.get_tag('RG')
  if (read_o is not None) and (read_o.get_tag('RG') != rg):
    raise RuntimeError('inconsistent readgroups (%s and %s) in pair %s' % (rg, read_o.get_tag('RG'), read_m.queryname))

  # Extract UMI and increment count
  if read_o is not None:
    tag = Umi(tag=tag, start=read_fwd.reference_start, end=read_rev.reference_end)
    umi_counts_per_rg[rg][tag].plus += 1
  else:
    tag = Umi(tag=tag, start=read_m.reference_start, end=read_m.reference_end)
    umi_counts_per_rg[rg][tag].plus += 1

# Process last position and close output
if seq is not None:
  process(seq, umi_counts_per_rg, output)
output.close()
