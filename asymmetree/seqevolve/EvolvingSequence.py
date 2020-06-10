# -*- coding: utf-8 -*-

"""
Evolving sequence.

Classes in this module:
    - EvoSeqElement     (site element in an evolving sequence)
    - EvoSeq            (evolving sequence)
"""

from enum import Enum, auto

from asymmetree.tools.DoublyLinkedList import DLList, DLListElement, DLListIterator


__author__ = 'David Schaller'


class State(Enum):
    
    ROOT       = auto()
    INSERTION  = auto()
    INHERITED  = auto()


class EvoSeqElement(DLListElement):
    """Site element in an evolving sequence."""
    
    __slots__ = ('status', 'site_id', 'parent_el', 'rate_class', 'rate_factor')
    
    def __init__(self, value,
                 status, site_id,
                 prev_el=None, next_el=None,
                 parent_el=None,
                 rate_class=0, rate_factor=1.0):
        
        super().__init__(value, prev_el=prev_el, next_el=next_el)
        
        self.status = status
        self.site_id = site_id
        self.parent_el = parent_el
        self.rate_class = rate_class
        self.rate_factor = rate_factor


class EvoSeq(DLList):
    """Evolving sequence."""
    
    def __init__(self):
        
        super().__init__()
        
        
    def __iter__(self):
        
        return EvoSeqIterator(self)


    def __next__(self):
        
        pass
    
    
    def append(self, value, status, site_id, parent_el=None,
               rate_class=0, rate_factor=1.0):
        
        new_end = EvoSeqElement(value, status, site_id,
                                prev_el=self._last,
                                parent_el=parent_el,
                                rate_class=rate_class,
                                rate_factor=rate_factor)
        if self._last:
            self._last._next = new_end
        self._last = new_end
        if not self._first:
            self._first = new_end
        self._count += 1
        return new_end
    
    
    def append_left(self, value, status, site_id, parent_el=None,
                    rate_class=0, rate_factor=1.0):
        
        new_start = EvoSeqElement(value, status, site_id,
                                  next_el=self._first,
                                  parent_el=parent_el,
                                  rate_class=rate_class,
                                  rate_factor=rate_factor)
        if self._first:
            self._first._prev = new_start
        self._first = new_start
        if not self._last:
            self._last = new_start
        self._count += 1
        return new_start
    
    
    def insert_right_of(self, element, value, status, site_id, parent_el=None,
                        rate_class=0, rate_factor=1.0):
        
        if element is self._last:
            new_element = self.append(value, status, site_id,
                                      parent_el=parent_el,
                                      rate_class=rate_class,
                                      rate_factor=rate_factor)
            
        else:
            new_element = EvoSeqElement(value, status, site_id,
                                        prev_el=element,
                                        next_el=element._next,
                                        parent_el=parent_el,
                                        rate_class=rate_class,
                                        rate_factor=rate_factor)
            new_element._next._prev = new_element
            element._next = new_element
            self._count += 1
        
        return new_element    
        
    
    def clone(self):
        
        child_seq = EvoSeq()
        
        parent_site = self._first
        while parent_site:
            child_seq.append(parent_site._value,
                             State.INHERITED, parent_site.site_id,
                             parent_el=parent_site,
                             rate_class=parent_site.rate_class,
                             rate_factor=parent_site.rate_factor)
            parent_site = parent_site._next
            
        return child_seq
    
    
    def element_pairs(self):
        """Generator for subsequent element pairs."""
        
        if self._count >= 2:
            
            first = self._first
            second = self._first._next
            
            while second:
                yield (first, second)
                first = second
                second = second._next
    
    
    def count_status(self, status):
        
        counter = 0
        
        site = self._first
        while site:
            if site.status == status:
                counter += 1
            site = site._next
                
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
            self._current = self._current._next
            return x
        else:
            raise StopIteration