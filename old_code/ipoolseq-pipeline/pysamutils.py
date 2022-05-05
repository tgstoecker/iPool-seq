# pysamutils.py, Copyright 2016, 2017 Florian G. Pflug
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

from collections import deque
from collections import defaultdict
from recordtype import recordtype

class read_pair_iterator:
  read_pair = recordtype("read_pair", ['left', 'right', 'complete'])

  def __init__(self, read_iterator, skip_unpaired = True):
    self.skip_unpaired = skip_unpaired
    self.iterator = read_iterator
    self.queue = deque()
    self.left_names = dict()
    self.skipped = defaultdict(int)
    self.total = 0

  def __iter__(self):
    return self
    
  def next(self):
    # Repeat until a valid pair is found
    while True:
      # Scan ahead until the front of the queue contains a complete pair.
      while not self.queue or not self.queue[0].complete:
        # Queue is empty or first entry is not yet complete. Fetch read.
        r = self.iterator.next()
        # (A) If both mates are unmapped, the pair is ignored. (B) If exactly
        # one mate is mapped, the pair is accepted, and the other mate
        # set to "None". If skip_unpaired was set, the pair is dropped instead.
        # (C) If both mates are mapped but discordantely, the pair is ignored.
        # (D) If both are mapped concordantly, the pair is (naturally) accepted.
        if r.is_unmapped and r.mate_is_unmapped:
          # A. Ignore.
          if r.is_read1:
            self.total += 1
            self.skipped['unmapped'] += 1
        elif r.is_unmapped and not r.mate_is_unmapped:
          # B, unmapped mate. Ignore
          pass
        elif r.mate_is_unmapped and not r.is_unmapped:
          self.total += 1
          if not self.skip_unpaired:
            # B, mapped mate. Create pair entry with missing mate, enqueue.
            p = read_pair_iterator.read_pair(r, None, True)
            self.queue.append(p)
          else:
            # B, mapped mate. Ignore because skip_unpaired is set
            self.skipped['unpaired'] += 1
        elif not r.is_proper_pair:
          # C. Ignore.
          if r.is_read1:
            self.total += 1
            self.skipped['discordant'] += 1
        elif r.query_name in self.left_names:
          # D, and mate was already seen, i.e. pair already has an entry.
          p = self.left_names[r.query_name]
          del self.left_names[r.query_name]
          p.right = r
          p.complete = True
        else:
          # D, and pair has no entry yet. Create one and flag as incomplete.
          p = read_pair_iterator.read_pair(r, None, False)
          self.queue.append(p)
          self.left_names[r.query_name] = p
          self.total += 1

      # Dequeue and return frontmost pair, which is now known to be complete
      # We assume that exactly one of the reads maps to the forward strand,
      # and return that one as the first read. If the forward read lies fully
      # to the right of the reverse read (i.e., if they look like this:
      # <--| |-->), the pair is skipped.
      p = self.queue.popleft()
      if (p.right is not None) and (p.left.is_reverse == p.right.is_reverse):
        raise RuntimeException("concordant pairs expected to have orientation FR or RF")
      if not p.left.is_reverse:
        fwd, rev = p.left, p.right
      else:
        fwd, rev = p.right, p.left
      if (fwd is not None) and (rev is not None) and (fwd.reference_start >= rev.reference_end):
        self.skipped['outwards'] += 1
        continue
      return (fwd, rev)

# Iterator wrapper which allows objects to be dynamically
# inserted into the iteration. Used to postpone handling
# of dovetailed reads below.
class injectable_iterator:
  def __init__(self, iterator):
    self.iterator = iterator
    self.queue = list()
    
  def __iter__(self):
    return self
    
  def inject(self, object, flag = None):
    self.queue.append((object, flag))
    
  def next(self):
    if self.queue:
      return self.queue.pop(0)
    else:
      return (self.iterator.next(), None)

# AlignmentFile writter which allows placeholders to be written, which
# are later replaced by actual reads. Once a placeholder was written,
# further reads are buffered in memory until the placeholder is filled in.
class reorder_buffer:
  object_placeholder = recordtype('reorder_buffer_object_placeholder', ['reorder_buffer', 'slot'])
  vacant = object()
  dont_fill = object()

  def __init__(self, writer):
    self.postponed = list()
    self.writer = writer
    
  def write(self, object):
    # If there are postponed objects, enqueue new object at the end
    # Otherwise, output immediately.
    if self.postponed:
      self.postponed.append(object)
    else:
      self.writer.write(object)
  
  def placeholder(self):
    # Create placeholder object and put onto postponed list
    ph = reorder_buffer.object_placeholder(reorder_buffer=self, slot=reorder_buffer.vacant)
    self.postponed.append(ph)
    return ph
  
  def fillin(self, object, placeholder = None):
    # If no placeholder is specified, assume it goes at the end
    if placeholder == None:
      return self.write(object)
    # Validate placeholder
    if placeholder.__class__ != reorder_buffer.object_placeholder:
      raise TypeError('placeholder must be of type reorder_buffer.object_placeholder')
    if id(placeholder.reorder_buffer) != id(self):
      raise ValueError('placeholder belongs to different reorder_buffer')
    if id(placeholder.slot) != id(reorder_buffer.vacant):
      raise ValueError('placeholder was filled in already')
    # Put onto postponed list at correct position
    placeholder.slot = object
    # Write postponed objects up to the first postponed placeholder
    while self.postponed:
      head = self.postponed[0]
      if head.__class__ != reorder_buffer.object_placeholder:
        self.postponed.pop(0)
        self.writer.write(head)
      elif id(head.slot) == id(reorder_buffer.vacant):
        break
      elif id(head.slot) == id(reorder_buffer.dont_fill):
        self.postponed.pop(0)
      else:
        self.postponed.pop(0)
        self.writer.write(head.slot)

  def dontfill(self, placeholder):
    self.fillin(reorder_buffer.dont_fill, placeholder=placeholder)

  def close(self):
    if self.postponed:
      raise RuntimeError('cannot close reorder_buffer containing placeholders')
    self.writer.close()
