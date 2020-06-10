# -*- coding: utf-8 -*-

"""
Doubly-linked list.

Implementation of a doubly-linked list in Python 3. The list enables access
to single list elements in order to modify/delete values in constant time.
"""

import collections


__author__ = 'David Schaller'


class DLListElement:
    """Doubly-linked list element."""
    
    __slots__ = ('_value', '_prev', '_next')
    
    def __init__(self, value, prev_el=None, next_el=None):
        
        self._value = value
        self._prev = prev_el
        self._next = next_el


class DLList:
    """Doubly-linked list."""
    
    __slots__ = ('_first', '_last', '_count')
    
    def __init__(self, *args):
        
        self._first = None
        self._last = None
        self._count = 0
        for arg in args:
            if isinstance(arg, collections.abc.Iterable):
                for item in arg:
                    self.append(item)
            else:
                self.append(item)
        
    
    def __len__(self):
        
        return self._count
    
    
    def __nonzero__(self):
        
        return True if self._count > 0 else False
    
    
    def __iter__(self):
        
        return DLListIterator(self)


    def __next__(self):
        
        pass
    
    
    def __getitem__(self, index):
        
        return self.element_at(index)._value
    
    
    def first(self):
        
        return self._first._value
    
    
    def last(self):
        
        return self._last._value
    
    
    def first_element(self):
        
        return self._first
    
    
    def element_at(self, index):
        
        if not isinstance(index, int):
            raise TypeError("index must be of type 'int'")
            
        if index < (-self._count) or index >= self._count:
            raise IndexError('index {} is out of bounds'.format(index))
            
        if index < 0:
            index += self._count
        
        if index <= self._count // 2:
            # start from beginning
            element = self._first
            for _ in range(index):
                element = element._next
        
        else:
            # start from end
            element = self._last
            for _ in range(self._count - index - 1):
                element = element._prev
            
        return element
    
    
    def append(self, value):
        
        new_end = DLListElement(value, prev_el=self._last)
        if self._last:
            self._last._next = new_end
        self._last = new_end
        if not self._first:
            self._first = new_end
        self._count += 1
        return new_end
    
    
    def extend(self, iterable):
        
        for value in iterable:
            self.append(value)
    
    
    def append_left(self, value):
        
        new_start = DLListElement(value, next_el=self._first)
        if self._first:
            self._first._prev = new_start
        self._first = new_start
        if not self._last:
            self._last = new_start
        self._count += 1
        return new_start
    
    
    def remove_element(self, element):
        """Remove an item by reference to the 'DLListElement' instance in O(1)."""
        
        if element._prev:
            element._prev._next = element._next
        if element._next:
            element._next._prev = element._prev
        if self._first is element:
            self._first = element._next
        if self._last is element:
            self._last = element._prev
        element._prev, element._next = None, None
        self._count -= 1
    
    
    def remove(self, value):
        """Remove an item by value in O(n)."""
        
        element = self._first
        while element:
            if element._value == value:
                self.remove_element(element)
                return
            element = element._next
        raise KeyError('value {} is not in the doubly-linked list'.format(value))
        
        
    def remove_range(self, index, length=None):
        """Removes a range from the index of the specified length.
        
        Removes the range [index, index+length) from the sequence. If no length
        is specified or index+length is out of bounds, the list gets truncated."""
        
        if not isinstance(index, int):
            raise TypeError("index must be of type 'int'")
        elif index < 0:
            index = self._count + index
            
        if length is not None and (not isinstance(length, int) or length < 1):
            raise TypeError("length must be of type 'int' and >0")
        
        if length is None or index + length >= self._count:
            self.truncate(index)
            
        elif index == 0:
            self.truncate_left(length)
            
        else:
            cut_start = self.element_at(index)
            cut_end = self.element_at(index + length - 1)
            
            cut_start._prev._next = cut_end._next
            cut_end._next._prev = cut_start._prev
            cut_start._prev, cut_end._next = None, None
            
            self._count -= length
    
    
    def insert_right_of(self, element, value):
        """Insert a new item right of an element of the list in O(1)."""
        
        if element is self._last:
            new_element = self.append(value)
            
        else:
            new_element = DLListElement(value,
                                        prev_el=element,
                                        next_el=element._next)
            new_element._next._prev = new_element
            element._next = new_element
            self._count += 1
        
        return new_element
            
    
    def truncate(self, index):
        """Truncate all elements starting from the specified index."""
        
        if index <= 0:
            self.clear()
        else:
            new_end = self.element_at(index-1)
            self._last = new_end
            if new_end._next:
                new_end._next._prev = None
                new_end._next = None
            self._count = index
            
    
    def truncate_left(self, n):
        """Truncate n elements on the left side."""
        
        if n >= self._count:
            self.clear()
        else:
            new_start = self.element_at(n)
            self._first = new_start
            if new_start._prev:
                new_start._prev._next = None
                new_start._prev = None
            self._count -= n
        
    
    def popright(self):
        """Removes the last element of the list and returns its value."""
        
        if self._last:
            value = self._last._value
            self.remove_element(self._last)
            return value
        else:
            return None
    
    
    def popleft(self):
        """Removes the first element of the list and returns its value."""
        
        if self._first:
            value = self._first._value
            self.remove_element(self._first)
            return value
        else:
            return None
    
    
    def clear(self):
        
        self._first = None
        self._last = None
        self._count = 0
        
    
    def _count_actual(self):
        """Counts the actual number of elements."""
        
        current = self._first
        counter = 0
        
        while current:
            counter += 1
            current = current._next
            
        return counter
        

class DLListIterator:
    """Iterator class for doubly-linked list."""
    
    def __init__(self, dllist):
        
        self.dllist = dllist
        self._current = dllist._first
        
    
    def __next__(self):
        
        if self._current:
            x = self._current
            self._current = self._current._next
            return x._value
        else:
            raise StopIteration