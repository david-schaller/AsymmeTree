"""Data structure for evolving sequences.

The EvoSeq class represents a sequence that evolves along a tree. It is implemented as a doubly
linked list, which allows for efficient insertions and deletions of sites, which are necessary to
simulate indels. Each site in the sequence is represented by an EvoSeqNode, which contains
information about the state of the site (whether it is inherited from the parent sequence, inserted,
or part of the root sequence), its site ID, and its parent element in the case of inherited sites.
"""

from __future__ import annotations

from enum import auto
from enum import Enum
from typing import Iterator

from tralda.datastructures.doubly_linked import DLList
from tralda.datastructures.doubly_linked import DLListIterator
from tralda.datastructures.doubly_linked import DLListNode


class State(Enum):
    """State of the sites in evolving sequences."""

    ROOT = auto()
    INSERTION = auto()
    INHERITED = auto()


class EvoSeqNode(DLListNode):
    """Site in an evolving sequence."""

    __slots__ = ("status", "site_id", "parent_el", "rate_class", "rate_factor")

    def __init__(
        self,
        value: int,
        status: State,
        site_id: int,
        prev_node: EvoSeqNode | None = None,
        next_node: EvoSeqNode | None = None,
        parent_el: EvoSeqNode | None = None,
        rate_class: int = 0,
        rate_factor: float = 1.0,
    ) -> None:
        """Constructor for the EvoSeqNode class.

        Args:
            value: The character at the site represented as an integer (index in the alphabet).
            status: The state of the site (whether it is inherited from the parent sequence,
                inserted, or part of the root sequence).
            site_id: The site ID, which is a unique identifier for the site across the entire tree.
            prev_node: The previous site in the sequence.
            next_node: The next site in the sequence.
        """
        super().__init__(value, prev_node=prev_node, next_node=next_node)

        self.status = status
        self.site_id = site_id
        self.parent_el = parent_el
        self.rate_class = rate_class
        self.rate_factor = rate_factor


class EvoSeq(DLList):
    """Evolving sequence.

    The EvoSeq class represents a sequence that evolves along a tree. It is implemented as a doubly
    linked list, which allows for efficient insertions and deletions of sites, which are necessary
    to simulate indels.
    """

    def __init__(self) -> None:
        """Constructor for the EvoSeq class."""
        super().__init__()

    def __iter__(self):
        """Iterator for the sites in the evolving sequence."""
        return EvoSeqIterator(self)

    def __next__(self):
        pass

    def append(
        self,
        value: int,
        status: State,
        site_id: int,
        parent_el: EvoSeqNode | None = None,
        rate_class: int = 0,
        rate_factor: float = 1.0,
    ) -> EvoSeqNode:
        """Append a new site to the end of the evolving sequence.

        Args:
            value: The character at the site represented as an integer (index in the alphabet).
            status: The state of the site (whether it is inherited from the parent sequence,
                inserted, or part of the root sequence).
            site_id: The site ID, which is a unique identifier for the site across the entire tree.
            parent_el: The parent element in the case of inherited sites.
            rate_class: The rate class of the site.
            rate_factor: The rate factor of the site.

        Returns:
            The newly appended EvoSeqNode.
        """
        new_end = EvoSeqNode(
            value,
            status,
            site_id,
            prev_node=self._last,
            parent_el=parent_el,
            rate_class=rate_class,
            rate_factor=rate_factor,
        )
        if self._last:
            self._last._next = new_end
        self._last = new_end
        if not self._first:
            self._first = new_end
        self._count += 1

        return new_end

    def append_left(
        self,
        value: int,
        status: State,
        site_id: int,
        parent_el: EvoSeqNode | None = None,
        rate_class: int = 0,
        rate_factor: float = 1.0,
    ) -> EvoSeqNode:
        """Append a new site to the beginning of the evolving sequence.

        Args:
            value: The character at the site represented as an integer (index in the alphabet).
            status: The state of the site (whether it is inherited from the parent sequence,
                inserted, or part of the root sequence).
            site_id: The site ID, which is a unique identifier for the site across the entire tree.
            parent_el: The parent element in the case of inherited sites.
            rate_class: The rate class of the site.
            rate_factor: The rate factor of the site.

        Returns:
            The newly appended EvoSeqNode (at the beginning of the sequence).
        """
        new_start = EvoSeqNode(
            value,
            status,
            site_id,
            next_node=self._first,
            parent_el=parent_el,
            rate_class=rate_class,
            rate_factor=rate_factor,
        )
        if self._first:
            self._first._prev = new_start
        self._first = new_start
        if not self._last:
            self._last = new_start
        self._count += 1

        return new_start

    def insert_right_of(
        self,
        element: EvoSeqNode,
        value: int,
        status: State,
        site_id: int,
        parent_el: EvoSeqNode | None = None,
        rate_class: int = 0,
        rate_factor: float = 1.0,
    ) -> EvoSeqNode:
        """Insert a new site to the right of the given element in the evolving sequence.

        Args:
            element: The element to the left of the new site.
            value: The character at the site represented as an integer (index in the alphabet).
            status: The state of the site (whether it is inherited from the parent sequence,
                inserted, or part of the root sequence).
            site_id: The site ID, which is a unique identifier for the site across the entire tree.
            parent_el: The parent element in the case of inherited sites.
            rate_class: The rate class of the site.
            rate_factor: The rate factor of the site.

        Returns:
            The newly inserted EvoSeqNode.
        """
        if element is self._last:
            new_element = self.append(
                value,
                status,
                site_id,
                parent_el=parent_el,
                rate_class=rate_class,
                rate_factor=rate_factor,
            )

        else:
            new_element = EvoSeqNode(
                value,
                status,
                site_id,
                prev_node=element,
                next_node=element._next,
                parent_el=parent_el,
                rate_class=rate_class,
                rate_factor=rate_factor,
            )
            new_element._next._prev = new_element
            element._next = new_element
            self._count += 1

        return new_element

    def clone(self) -> EvoSeq:
        """Clone the evolving sequence to create a child sequence from a parent sequence.

        All sites in the child sequence are marked as inherited, and they reference the
        corresponding sites in the parent sequence.

        Returns:
            A new EvoSeq that is a clone of the current sequence.
        """
        child_seq = EvoSeq()

        parent_site = self._first
        while parent_site:
            child_seq.append(
                parent_site._value,
                State.INHERITED,
                parent_site.site_id,
                parent_el=parent_site,
                rate_class=parent_site.rate_class,
                rate_factor=parent_site.rate_factor,
            )
            parent_site = parent_site._next

        return child_seq

    def element_pairs(self) -> Iterator[tuple[EvoSeqNode, EvoSeqNode]]:
        """Generator for subsequent element pairs.

        Yields:
            Tuples of subsequent elements in the evolving sequence.
        """
        if self._count >= 2:
            first = self._first
            second = self._first._next

            while second:
                yield (first, second)
                first = second
                second = second._next

    def count_status(self, status: State) -> int:
        """Count the number of sites with a given status.

        Args:
            status: The status to count.

        Returns:
            The number of sites with the given status.
        """
        counter = 0

        site = self._first
        while site:
            if site.status == status:
                counter += 1
            site = site._next

        return counter


class EvoSeqIterator(DLListIterator):
    """Iterator class for evolving sequence."""

    def __init__(self, evoseq: EvoSeq) -> None:
        """Constructor for the EvoSeqIterator class.

        Args:
            evoseq: The evolving sequence to iterate over.
        """
        super().__init__(evoseq)

    def __next__(self) -> EvoSeqNode:
        """Next site element in the evolving sequence.

        Overrides the super class methods and returns the element instead of the value only.

        Returns:
            The next EvoSeqNode in the evolving sequence.
        """
        if self._current:
            x = self._current
            self._current = self._current._next
            return x
        else:
            raise StopIteration
