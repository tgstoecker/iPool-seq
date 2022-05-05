# umi_error_correction.py, Copyright 2016, 2017 Florian G. Pflug
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

import sys
import distance
from math import floor
from collections import defaultdict
from collections import namedtuple
from recordtype import recordtype

# Immutable representation of an UMI
Umi = namedtuple('Umi', ['start', 'end', 'tag'])

# Mutable representation of a node in the UMI forest. The field 'umi' contains
# an Umi, and 'parents' a set of parent nodes, which the Umi is similar too,
# *and* which have a higher total read count. If the total read counts of two
# similar Umis are equal, the less-than (<) operator of the Umi type is used to
# break the tie, and the "greater" Umi is made the parent. This ensures
# non-circularity of the graph in all cases.
UmiEC = recordtype('UmiEC', ['umi', 'parents'])

# Represents Umi read counts, wich reads of the sense (+) and anti-sense (-)
# strands tracked separetely. Note that UmiCounts are not ordered!
# Also tracks the number of raw umis that contributed to the read counts
class UmiCount:
  def __init__(self, plus, minus, rawumis):
    self.rawumis = rawumis
    self.plus = plus
    self.minus = minus
    
  def total(self):
    return self.plus + self.minus
    
  def __add__(self, other):
    return UmiCount(self.plus + other.plus, self.minus + other.minus, self.rawumis + other.rawumis)
    
  def __eq__(self, other):
    return (self.plus == other.plus) and (self.minus == other.minus)

  def __ne__(self, other):
    return not (self == other)
    
  def __str__(self):
    return '%d+,%d-' % (self.plus, self.minus)

# Error-Correction implementation.
# Two Umis whose start and end-positions each differ by at most max_mapping_distance
# bases, and whose tags have hamming distance at most max_hamming_distance are merged
# into one.
class ErrorCorrector:
  def __init__(self, tag_length, max_mapping_distance, max_hamming_distance, report=False):
    self.tag_length = tag_length
    self.max_mapping_distance = max_mapping_distance
    self.max_hamming_distance = max_hamming_distance
    self.tsw = int(floor(tag_length / (self.max_hamming_distance+1)))
    self.report = report
    self.report_step_umis = 500
    self.report_step_basepairs = 50

  def process(self, umi_counts):
    if self.report:
      sys.stderr.write('%d UMIs with %s reads before error-correction\n' % (len(umi_counts), sum(umi_counts.itervalues(), UmiCount(0,0,0))))

    # Build forest 
    ecinfo = self.build_forest(umi_counts)
    # Find roots and sum counts
    ecumi_counts = self.sum_trees(umi_counts, ecinfo)
    
    if self.report:
      sys.stderr.write('%d UMIs with %s reads after error-correction\n' % (len(ecumi_counts), sum(ecumi_counts.itervalues(), UmiCount(0,0,0))))

    return ecumi_counts

  def link_if_similar(self, umi, umi2, umi_counts, ecinfo):
    log = False
#      if (umi2.tag == 'CTTCGCGTTCACGAACTAAG:GCAGTGAATTAGTTGTCCAT') and (umi.tag == 'CTTCGCGTTCACGAACTAAG:GCAGTGAATTAGTTGTCCAT'):
#        log = True
    if log:
      print >> sys.stderr, 'DEBUG: UMI1: %s x %d' % (umi, umi_counts[umi])
      print >> sys.stderr, 'DEBUG: UMI2: %s x %d' % (umi2, umi_counts[umi2])
    # Check mapping distance (start)
    if abs(umi2.start - umi.start) > self.max_mapping_distance:
      # Mapping positions lie too far apart to unify UMIs.
      if log:
        print >> sys.stderr, 'DEBUG: NOT MERGED, start distance >  %d' % self.max_mapping_distance
      return
    elif umi2.start < umi.start:
      # Impossible, list is sorted!
      raise RuntimError('internal error, UMI list corrupted')
    # Check mapping distance (end)
    if abs(umi2.end - umi.end) > self.max_mapping_distance:
      # Mapping positions lie too far apart to unify UMIs
      if log:
        print >> sys.stderr, 'DEBUG: NOT MERGED, end distance >  %d' % self.max_mapping_distance
      return
    # Check hamming distance
    hd = False
    if self.max_hamming_distance > 0:
      hd = (distance.hamming(umi.tag, umi2.tag) > self.max_hamming_distance)
    else:
      hd = (umi.tag != umi2.tag)
    if hd:
      # Hamming distance too large to allow unification of UMIs
      if log:
        print >> sys.stderr, 'DEBUG: NOT MERGED, hamming distance > %d' % self.max_hamming_distance
      return
    # Determine which to merge into which
    umi_cnt = umi_counts[umi].total()
    umi2_cnt = umi_counts[umi2].total()
    if (umi_cnt, umi) < (umi2_cnt, umi2):
      p, c = umi2, umi
    elif (umi2_cnt, umi2) < (umi_cnt, umi):
      p, c = umi, umi2
    else:
      raise RuntimeError("Tertium non datur fails for %s and %s" % (umi, umi2))
    if log:
      print >> sys.stderr, 'DEBUG: WILL MERGE %s x %d --> %s x %d' % (c, umi_counts[c], p, umi_counts[p])
    # Merge UMIs as described above 
    ecinfo[c].parents.add(ecinfo[p])

  # Splits a tag into enough disjoint parts such that whenever the hamming
  # distance of two tags is at most max_hamming_distance, then at least one
  # part is identical in both tags.
  def splittag(self, tag):
    if len(tag) != self.tag_length:
      raise ValueError("tag %s has wrong length, %d instead of %d" % (tag, len(tag), self.tag_length))
    return [ tag[ i*self.tsw : ( (i+1)*self.tsw if i < self.max_hamming_distance else self.tag_length ) ]
             for i in range(self.max_hamming_distance+1) ]
  
  def build_forest(self, umi_counts):
    if self.report:
      sys.stderr.write('Finding similar UMIs...\n')
      sys.stderr.flush()

    # Data-structure which represents the forest. ecinfo associates each Umi 
    # with an UmIEC instance, which contains the Umi itself, and a set of
    # parents. Note that since parents is a set, the "forest" does not actually
    # contain trees, but DAGs.
    ecinfo = dict()
    for u in umi_counts.iterkeys():
      ecinfo[u] = UmiEC(u, set())

    # Build Umi index. The first level splits Umis by their starting position.
    # Below that, there are max_hamming_distance + 1 associative maps (dicts),
    # M_k, one for each of the parts returned by splittag(). The map M_k maps a
    # string s to the set of Umis whose k-th tag part equals s. Therefore, given
    # a tag x, all tags y whose hamming distance is at most max_hamming_distance,
    # are contained in at least one of the sets M_k(x_k), where x_k is the k-th
    # tag part of x. 
    umi_index = defaultdict(lambda: [ defaultdict(set) for i in range(0, self.max_hamming_distance+1) ])
    for umi in umi_counts.iterkeys():
      t = umi_index[umi.start]
      for i, p in enumerate(self.splittag(umi.tag)):
        t[i][p].add(umi)
    
    # Populate ecinfo
    last_pos = 0
    last_report_pos = 0
    # Scan Umis
    for i, umi in enumerate(sorted(umi_counts.iterkeys(), key=lambda u: u.start)):
      if umi.start < last_pos:
        raise RuntimeError('UMIs not sorted by start position (%d after %d)!' % (umi.start, last_pos))
      if self.report and (i % self.report_step_umis == 0):
        sys.stderr.write('  processed %d UMIs\n' % i)
        sys.stderr.flush()
      if self.report and ((umi.start - last_report_pos) >= self.report_step_basepairs):
        last_report_pos = umi.start
        sys.stderr.write('  processed %d basepairs\n' % umi.start)
        sys.stderr.flush()
      last_pos = umi.start
      # Scan part of first index level corresponding to a similar start position.
      for s2 in range(umi.start, umi.start + self.max_mapping_distance + 1):
        # Scan tag parts of 1st Umi
        for k2, p2 in enumerate(self.splittag(umi.tag)):
          # Scan Umis having the same k2-th tag part p2 as the 1st Umi
          for umi2 in umi_index[s2][k2][p2]:
            # Avoid calling link_if_similar twice for every pair, using the total
            # order imbued on the Umis to decide which call to suppress. Note that
            # whatever condition we use here *must* ensure that if umi has a lower
            # start position than umi2, link_if_similar *will* be called, since we
            # only let s2 range within [ umi.start, umi.start + max_dist ] above.
            if umi < umi2:
              self.link_if_similar(umi, umi2, umi_counts, ecinfo)
    
    # Return forest
    return ecinfo
              
  def sum_trees(self, umi_counts, ecinfo):
    if self.report:
     sys.stderr.write('Merging similar UMIs...\n')
     sys.stderr.flush()
     
    # All counts within an UMI tree contribute to the roots count.
    # If the tree is actually a net, i.e. if there is more than one "root"
    # reachable from a certain leaf, that leaf's assignment is ambiguous,
    # and its count is not attributed to any root.
    processed = 0
    ambiguous = 0
    ambiguous_count = UmiCount(0,0,0)
    ecumi_counts = defaultdict(lambda: UmiCount(0,0,0))
    for umi, umi_ec in ecinfo.iteritems():
      if self.report and (processed % self.report_step_umis == 0):
        sys.stderr.write('  processed %d UMIs\n' % processed)
        sys.stderr.flush()
      processed += 1
      # Start with current node as "active set"
      nodes = set([umi_ec])
      root = None
      while nodes:
        # Pick node in active set,
        n = nodes.pop()
        if n.parents:
          # and replace it by its parents
          nodes.update(set(n.parents))
        elif root == None:
          # or mark as root if it has none,
          root = n.umi
        elif root != n.umi:
          # but if we already found a root, it's ambiguous.
          root = None
      if root != None:
        # Found a unique root, so add node's counts to root's counts
        if umi_counts[umi].rawumis != 1:
          raise RuntimeError('before error correction, umis should have rawumis=1')
        ecumi_counts[root] += umi_counts[umi]
      else:
        # Found more than one root, so throw node away
        if umi_counts[umi].rawumis != 1:
          raise RuntimeError('before error correction, umis should have rawumis=1')
        ambiguous_count += umi_counts[umi]
    # Check that we haven't missed any UMI
    if sum(umi_counts.itervalues(), UmiCount(0,0,0)) != (sum(ecumi_counts.itervalues(), UmiCount(0,0,0)) + ambiguous_count):
      raise RuntimeError('error correction buggy!')
    if self.report:
      sys.stderr.write('%d UMIs with %s reads were ambiguous\n' % (ambiguous_count.rawumis, ambiguous_count))
      sys.stderr.flush()
    return ecumi_counts
