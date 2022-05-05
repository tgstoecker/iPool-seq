#!/usr/bin/env python
# assign_to_features.py, Copyright 2016, 2017 Florian G. Pflug
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
# Implements "Assignment to mutants"
# ./assign_to_features.py <gff-file> <input bam> <output bam>
# *****************************************************************************
import pysam
import re
import sys
import itertools
import copy
import pysamutils
import pprint
from BCBio import GFF
from collections import namedtuple
from collections import defaultdict

READNAME_PATTERN = re.compile('^(.*)\|([ACGTN]{12})(:([ACGTN]{12}))?$')

FEATURE_BOUNDARIES_FUZZ = 10
FEATURE_FLANK_LENGTH = 1000

# Command-line arguments
ifile_features = sys.argv[1]
ifile = sys.argv[2]
ofile = sys.argv[3]

# Read features from input GFF file
Feature = namedtuple('Feature', ['name', 'start', 'end', 'strand', 'left', 'right', 'hits'])
features_bychr = defaultdict(list)
with open(ifile_features) as input_features:
  for r in GFF.parse(input_features):
    features_chr = features_bychr[r.id]
    for f in r.features:
      features_chr.append(Feature(name=f.qualifiers['Name'][0],
                                  start=int(f.location.start),
                                  end=int(f.location.end),
                                  strand=int(f.strand),
                                  left='5p' if f.strand == +1 else '3p',
                                  right='3p' if f.strand == +1 else '5p',
                                  hits=[0, 0]))

# Sort features within chromosomes by their start position, and
# Build index which allows looking up features by their ID.
features_byid = dict()
for chr, features_chr in features_bychr.iteritems():
  features_chr.sort(key=lambda x: x.start)
  p = None
  for f in features_chr:
    if (p is not None) and (p.end >= f.start):
      raise ValueError('features must not overlap', p.name, f.name)
    p = f
    features_byid[id(f)] = f

# Open input file
input = pysam.AlignmentFile(ifile, mode="rb")
if input.header['HD']['SO'] != 'coordinate':
  raise RuntimeError("input file must be sorted by coordinate")

# Open output file
header = copy.deepcopy(dict(input.header))
header['RG'] = ([ { 'ID': '%s:%s' % (f.name, s), 'SM' : '1' }
                  for f in features_byid.itervalues()
                  for s in ['3p', '5p'] ] +
                [ { 'ID': 'ambiguous', 'SM' : '1' },
                  { 'ID': 'unmatched', 'SM' : '1' } ])
output = pysamutils.reorder_buffer(pysam.AlignmentFile(ofile, mode="wb", header=header))

# Information about previously seen mate
Mate = namedtuple('Mate', ['start', 'end'])
mates = {}

fragments = 0
discordant = 0
unmatched = 0
ambiguous = 0
for read in input.fetch(until_eof=True):
  if read.is_unmapped:
    # Skip unmapped reads
    continue
  elif not read.is_proper_pair and not read.mate_is_unmapped:
    # Skip discordant pairs where both mates are mapped
    if read.is_read1:
      discordant = discordant + 1
    continue
  elif read.is_proper_pair and read.query_name not in mates:
    # Wait for mate before processing read
    mates[read.query_name] = (read, output.placeholder())
    continue

  # Switch to correct feature list if chromosome changed
  r_chr = input.getrname(read.rname)
  if (chr != r_chr):
    features_chr = features_bychr[r_chr]
    chr = r_chr

  # Initialize
  fragments = fragments + 1
  hits = set()

  # Fetch mate (must have been seen earlier)
  start, end = read.reference_start, read.reference_end
  mate, mate_placeholder = None, None
  if read.is_proper_pair:
    mate, mate_placeholder = mates[read.query_name]
    del mates[read.query_name]
    start, end = min(mate.reference_start, start), max(mate.reference_end, end)

  # Move RMID from read name to XT tag.
  m = READNAME_PATTERN.match(read.query_name)
  if m is None:
    read.query_name = m
    rmid1, rmid2 = "%d-%d" % (start, end), None
  else:
    read.query_name = m.group(1)
    if mate is not None:
      mate.query_name = m.group(1)
    rmid1, rmid2 = m.group(2), m.group(4)
  swapped = 0 if (rmid2 is None) or (rmid1 < rmid2) else 1 if rmid1 > rmid2 else None
  rmid = rmid1 if rmid2 is None else ":".join(sorted([rmid1, rmid2]))
  read.set_tag(tag='XT', value=rmid, value_type='Z', replace=True)
  read.set_tag(tag='XS', value=swapped, value_type='i', replace=True)
  if mate is not None:
    mate.set_tag(tag='XT', value=rmid, value_type='Z', replace=True)
    mate.set_tag(tag='XS', value=swapped, value_type='i', replace=True)

  # Find features whose left flank coincides with the fragment's end,
  # or whose right flank coincides with the fragment's start. Feature boundaries
  # are assumed to be only be determined up to FEATURE_BOUNDARIES_FUZZ bases.
  for f in features_chr:
    if end >= f.start - FEATURE_BOUNDARIES_FUZZ and end < f.start + FEATURE_BOUNDARIES_FUZZ:
      hits.add((id(f), f.left))
    if start >= f.end - FEATURE_BOUNDARIES_FUZZ and start <= f.end + FEATURE_BOUNDARIES_FUZZ:
      hits.add((id(f), f.right))

  # For single-ended reads, additionally consider features a match if the read points
  # toward the feature and starts no further than FEATURE_FLANK_LENGTH away. Feature
  # boundaries are again assumed to only be determined up to FEATURE_BOUNDARIES_FUZZ. 
  if not read.is_proper_pair:
    start = read.reference_start
    end = read.reference_end
    if not read.is_reverse:
      for f in features_chr:
        if start >= f.start - FEATURE_BOUNDARIES_FUZZ - FEATURE_FLANK_LENGTH and end <= f.start + FEATURE_BOUNDARIES_FUZZ:
          hits.add((id(f), f.left))
    else:
      for f in features_chr:
        if start >= f.end - FEATURE_BOUNDARIES_FUZZ and end <= f.end + FEATURE_BOUNDARIES_FUZZ + FEATURE_FLANK_LENGTH:
          hits.add((id(f), f.right))
        
  # Determine whether the fragments matches a feature uniquely
  rg = None
  if len(hits) == 0:
    unmatched = unmatched + 1
    rg = 'unmatched'
  elif len(hits) > 1:
    ambiguous = ambiguous + 1
    rg = 'ambiguous'
  else:
    f_id, s = hits.pop()
    rg = '%s:%s' % (features_byid[f_id].name, s)

  # Output current read's mate (we process fragments after seeing the second read)
  if mate is not None:
    mate.set_tag(tag='RG', value_type='Z', value=rg, replace=True)
    output.fillin(mate, placeholder=mate_placeholder)

  # And output this read.
  read.set_tag(tag='RG', value_type='Z', value=rg, replace=True)
  output.write(read)

output.close()    
print >> sys.stderr, 'skipped %d discordant fragments' % discordant
print >> sys.stderr, 'skipped %d fragments matching none of the features' % unmatched
print >> sys.stderr, 'skipped %d fragments matching multiple features' % ambiguous
