# -*- coding: utf-8 -*-

"""
Doubly-linked list.

Implementation of a doubly-linked list in Python 3. The list enables access
to single list elements in order to modify/delete values in constant time.

Classes in this module:
    - DLListElement     (doubly-linked list element)
    - DLList            (doubly-linked list)
"""

import collections


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2020, David Schaller"


class DLListElement:
    """Doubly-linked list element."""
    
    __slots__ = ['value', 'prev_el', 'next_el']
    
    def __init__(self, value, prev_el=None, next_el=None):
        
        self.value = value
        self.prev_el = prev_el
        self.next_el = next_el


class DLList:
    """Doubly-linked list."""
    
    __slots__ = ['_first', '_last', 'count', 'current']
    
    def __init__(self, *args):
        
        self._first = None
        self._last = None
        self.count = 0
        for arg in args:
            if isinstance(arg, collections.Iterable):
                for item in arg:
                    self.append(item)
            else:
                self.append(item)
        
    
    def __len__(self):
        
        return self.count
    
    
    def __nonzero__(self):
        
        return True if self.count > 0 else False
    
    
    def __iter__(self):
        
        return DLListIterator(self)


    def __next__(self):
        
        pass
    
    
    def __getitem__(self, index):
        
        if not isinstance(index, int):
            raise TypeError("Index must be of type 'int'!")
        if index >= 0:
            if index >= self.count:
                raise IndexError("Index " + str(index) + " is out of bounds!")
            element = self._first
            for i in range(index):
                element = element.next_el
            return element.value
        else:
            if index < (-self.count):
                raise IndexError("Index " + str(index) + " is out of bounds!")
            element = self._last
            for i in range(-index - 1):
                element = element.prev_el
            return element.value
        raise IndexError("Index " + str(index) + " is out of bounds!")
    
    
    def first(self):
        
        return self._first.value
    
    
    def last(self):
        
        return self._last.value
    
    
    def first_element(self):
        
        return self._first
    
    
    def append(self, value):
        
        new_end = DLListElement(value, prev_el=self._last)
        if self._last:
            self._last.next_el = new_end
        self._last = new_end
        if not self._first:
            self._first = new_end
        self.count += 1
        return new_end
    
    
    def extend(self, iterable):
        
        for value in iterable:
            self.append(value)
    
    
    def append_left(self, value):
        
        new_start = DLListElement(value, next_el=self._first)
        if self._first:
            self._first.prev_el = new_start
        self._first = new_start
        if not self._last:
            self._last = new_start
        self.count += 1
        return new_start
    
    
    def remove_element(self, element):
        """Remove an item by reference to the 'DLListElement' instance in O(1)."""
        
        if element.prev_el:
            element.prev_el.next_el = element.next_el
        if element.next_el:
            element.next_el.prev_el = element.prev_el
        if self._first is element:
            self._first = element.next_el
        if self._last is element:
            self._last = element.prev_el
        element.prev_el, element.next_el = None, None
        self.count -= 1
    
    
    def remove(self, value):
        """Remove an item by value in O(n)."""
        
        element = self._first
        while element:
            if element.value == value:
                self.remove_element(element)
                return
            element = element.next_el
        raise KeyError("Value", value, "is not in the doubly-linked list!")
    
    
    def popright(self):
        """Removes the last element of the list and returns its value."""
        
        if self._last:
            value = self._last.value
            self.remove_element(self._last)
            return value
        else:
            return None
    
    
    def popleft(self):
        """Removes the first element of the list and returns its value."""
        
        if self._first:
            value = self._first.value
            self.remove_element(self._first)
            return value
        else:
            return None
    
    
    def clear(self):
        
        self._first = None
        self._last = None
        self.count = 0
        

class DLListIterator:
    """Iterator class for doubly-linked list."""
    
    def __init__(self, dllist):
        
        self.dllist = dllist
        self._current = dllist._first
        
    
    def __next__(self):
        
        if self._current:
            x = self._current
            self._current = self._current.next_el
            return x.value
        else:
            raise StopIteration
        

if __name__ == "__main__":
    dllist = DLList([7, 8, 9])
    dllist.append(1)
    dllist.append(2)
    dllist.append(3)
    dllist.append(4)
    dllist.append(5)
    dllist.remove(3)
    dllist.append_left(30)
    for x in dllist:
        print(x, end=" ")
    print("\n", dllist[0], dllist[1], dllist[2], dllist[3], dllist[-4])