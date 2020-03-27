# -*- coding: utf-8 -*-


"""
Evolving sequence.

Classes in this module:
    - EvoSeqElement     (site element in an evolving sequence)
    - EvoSeq            (evolving sequence)
"""

from enum import Enum, auto

from asymmetree.tools.DoublyLinkedList import DLList, DLListElement, DLListIterator


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2020, David Schaller"


class State(Enum):
    
    ROOT       = auto()
    INSERTION  = auto()
    INHERITED  = auto()


class EvoSeqElement(DLListElement):
    """Site element in an evolving sequence."""
    
    __slots__ = ('status', 'site_id', 'parent_el')
    
    def __init__(self, value,
                 status, site_id,
                 prev_el=None, next_el=None,
                 parent_el=None):
        
        super().__init__(value, prev_el=prev_el, next_el=next_el)
        
        self.status = status
        self.site_id = site_id
        self.parent_el = parent_el


class EvoSeq(DLList):
    """Evolving sequence."""
    
    __slots__ = ('_first', '_last', 'count', 'current')
    
    def __init__(self, *args):
        
        super().__init__(args)
        
        
    def __iter__(self):
        
        return EvoSeqIterator(self)


    def __next__(self):
        
        pass
    
    
    def append(self, value, status, site_id, parent_el=None):
        
        new_end = EvoSeqElement(value, status, site_id,
                                prev_el=self._last,
                                parent_el=parent_el)
        if self._last:
            self._last._next_el = new_end
        self._last = new_end
        if not self._first:
            self._first = new_end
        self._count += 1
        return new_end
    
    
    def append_left(self, value, status, site_id, parent_el=None):
        
        new_start = EvoSeqElement(value, status, site_id,
                                  next_el=self._first,
                                  parent_el=parent_el)
        if self._first:
            self._first._prev_el = new_start
        self._first = new_start
        if not self._last:
            self._last = new_start
        self._count += 1
        return new_start
        
    
    def clone(self):
        
        child_seq = EvoSeq()
        
        parent_site = self._first
        while parent_site:
            child_seq.append(parent_site._value,
                             State.INHERITED, parent_site.site_id,
                             parent_el=parent_site)
            parent_site = parent_site._next_el
            
        return child_seq
    
    
    def count_inherited(self):
        
        counter = 0
        
        site = self._first
        while site:
            if site.status == State.INHERITED:
                counter += 1
            site = site._next_el
                
        return counter
        
          
            
class EvoSeqIterator(DLListIterator):
    """Iterator class for evolving sequence."""
    
    def __init__(self, evoseq):
        
        super().__init__(evoseq)
        
    
    def __next__(self):
        """Next site element in the evolving sequence.
        
        Overrides the super class methods and returns the element instead
        of the value only."""
        
        if self._current:
            x = self._current
            self._current = self._current._next_el
            return x
        else:
            raise StopIteration