# -*- coding: utf-8 -*-

## Copyright 2015 Kevin B Jacobs
##
## Licensed under the Apache License, Version 2.0 (the "License"); you may
## not use this file except in compliance with the License.  You may obtain
## a copy of the License at
##
##        http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
## WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
## License for the specific language governing permissions and limitations
## under the License.


import random
import collections

from operator import itemgetter, mul
from itertools import chain, combinations, count, cycle, groupby, islice, repeat, starmap, tee
from itertools import filterfalse, zip_longest


# Python itertools recipes
#   http://docs.python.org/2/library/itertools.html#recipes

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


def tabulate(function, start=0):
    "Return function(0), function(1), ..."
    return map(function, count(start))


def consume(iterator, n):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(islice(iterator, n, n), None)


def nth(iterable, n, default=None):
    "Returns the nth item or a default value"
    return next(islice(iterable, n, None), default)


def quantify(iterable, pred=bool):
    "Count how many times the predicate is true"
    return sum(map(pred, iterable))


def padnone(iterable):
    """Returns the sequence elements and then returns None indefinitely.

    Useful for emulating the behavior of the built-in map() function.
    """
    return chain(iterable, repeat(None))


def ncycles(iterable, n):
    "Returns the sequence elements n times"
    return chain.from_iterable(repeat(tuple(iterable), n))


def dotproduct(vec1, vec2):
    return sum(map(mul, vec1, vec2))


def flatten(listOfLists):
    "Flatten one level of nesting"
    return chain.from_iterable(listOfLists)


def repeatfunc(func, times=None, *args):
    """Repeat calls to func with specified arguments.

    Example:  repeatfunc(random.random)
    """
    if times is None:
        return starmap(func, repeat(args))
    return starmap(func, repeat(args, times))


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks"

    >>> [''.join(g) for g in grouper('ABCDEFG', 3, 'x')]
    ['ABC', 'DEF', 'Gxx']
    """
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def roundrobin(*iterables):
    "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    # Recipe credited to George Sakkis
    pending = len(iterables)
    nexts = cycle(iter(it).__next__ for it in iterables)
    while pending:
        try:
            for item in nexts:
                yield item()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))


def roundrobin2(*iterables):
    "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    sentinel = object()
    return (x for x in chain(*zip_longest(fillvalue=sentinel, *iterables)) if x is not sentinel)


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def unique_everseen(iterable, key=None):
    """List unique elements, preserving order. Remember all elements ever seen.

    >>> ''.join(unique_everseen('AAAABBBCCDAABBB'))
    'ABCD'
    >>> ''.join(unique_everseen('ABBCcAD', str.lower))
    'ABCD'
    """
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in filterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element


def unique_justseen(iterable, key=None):
    """
    List unique elements, preserving order. Remember only the element just seen."

    >>> ''.join(unique_justseen('AAAABBBCCDAABBB'))
    'ABCDAB'
    >>> ''.join(unique_justseen('ABBCcAD', str.lower))
    'ABCAD'
    """
    return map(next, map(itemgetter(1), groupby(iterable, key)))


def iter_except(func, exception, first=None):
    """Call a function repeatedly until an exception is raised.

    Converts a call-until-exception interface to an iterator interface.
    Like __builtin__.iter(func, sentinel) but uses an exception instead
    of a sentinel to end the loop.

    Examples:
        bsddbiter = iter_except(db.next, bsddb.error, db.first)
        heapiter = iter_except(functools.partial(heappop, h), IndexError)
        dictiter = iter_except(d.popitem, KeyError)
        dequeiter = iter_except(d.popleft, IndexError)
        queueiter = iter_except(q.get_nowait, Queue.Empty)
        setiter = iter_except(s.pop, KeyError)

    """
    try:
        if first is not None:
            yield first()
        while 1:
            yield func()
    except exception:
        pass


def random_product(*args, **kwds):
    "Random selection from itertools.product(*args, **kwds)"
    pools = list(map(tuple, args)) * kwds.get('repeat', 1)
    return tuple(random.choice(pool) for pool in pools)


def random_permutation(iterable, r=None):
    "Random selection from itertools.permutations(iterable, r)"
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(random.sample(pool, r))


def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)


def random_combination_with_replacement(iterable, r):
    "Random selection from itertools.combinations_with_replacement(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.randrange(n) for i in range(r))
    return tuple(pool[i] for i in indices)


def tee_lookahead(t, i):
    """Inspect the i-th upcomping value from a tee object
       while leaving the tee object at its current position.

       Raise an IndexError if the underlying iterator doesn't
       have enough values.

    """
    for value in islice(t.__copy__(), i, None):
        return value
    raise IndexError(i)

# Contributed itertools recipes

_nothing = object()


def is_empty_iter(it):
    return next(it, _nothing) is _nothing


def first(iterable, default=_nothing):
    """Return the first item of an iterable, ``default`` if there is none."""
    try:
        return next(iter(iterable))
    except StopIteration:
        if default is _nothing:
            raise ValueError('first element does not exist')
        return default


def only_one(iterable, default=_nothing, sentinel=_nothing):
    """
    Return the first item from iterable, if and only if iterable contains a
    single element.  Raises ValueError if iterable contains more than a
    single element.  If iterable is empty, then return default value, if
    provided.  Otherwise raises ValueError.
    """
    it = iter(iterable)

    try:
        item = next(it)
    except StopIteration:
        if default is not _nothing:
            return default
        raise ValueError('zero length sequence')

    try:
        next(it)
        if sentinel is not _nothing:
            return sentinel
        raise ValueError('there can be only one')
    except StopIteration:
        return item


def chunked(iterable, n):
    """Break an iterable into lists of a given length::

        >>> list(chunked([1, 2, 3, 4, 5, 6, 7], 3))
        [(1, 2, 3), (4, 5, 6), (7,)]

    If the length of ``iterable`` is not evenly divisible by ``n``, the last
    returned list will be shorter.
    """
    for group in zip_longest(*[iter(iterable)] * n, fillvalue=_nothing):
        if group[-1] is _nothing:
            # If this is the last group, shuck off the padding:
            group = group[:group.index(_nothing)]
        yield group


def ilen(iterable):
    """Return the number of items in ``iterable``."""
    return sum(1 for _ in iterable)


# Helpers for zip_exact

class LengthMismatch(Exception):
    pass


def _zip_exact_thow():
    raise LengthMismatch
    yield None  # unreachable


def _zip_exact_check(rest):
    for i in rest:
        try:
            next(i)
        except LengthMismatch:
            pass
        else:
            raise LengthMismatch
    return
    yield None  # unreachable


def zip_exact(*iterables):
    """
    zip_exact(iter1 [,iter2 [...]]) --> iterator object

    Return an iterator whose .next() method returns a tuple where the i-th
    element comes from the i-th iterable argument.  The .next() method
    continues until the shortest iterable in the argument sequence is
    exhausted.  If all sequences are exhausted a StopIteration exception is
    raised, otherwise a LengthMismatch exception is raised.

    Works like itertools.zip(), but throws a LengthMismatch exception if any
    iterable's length differs.

    Zero length iterable lists return iter(iterables).  Unit length iterables
    return iter(iterables[0]).  Otherwise an zip object is returned.

    If the len() of all iterables is available and match, then this function
    returns zip(*iterables).  In the length of any iterable cannot be
    determined, then sentinel objects are appended to each iterable and an
    zip object of these augmented iterables is returned.  If the length a
    proper subset of iterables can be determined, the second strategy is
    employed, even if known lengths do not match.  This is in anticipation of
    iterable side-effects that may ultimately balance the lengths.

    Inspired largely by Peter Otten's zip_exc, modified to check lengths.
    (http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/497006)

    >>> list(zip_exact())
    []
    >>> list(zip_exact([]))
    []
    >>> list(zip_exact((), (), ()))
    []

    >>> list(zip_exact("abc", range(3)))
    [('a', 0), ('b', 1), ('c', 2)]

    >>> list(zip_exact("", range(3)))
    Traceback (most recent call last):
         ...
    vgraph.iterstuff.LengthMismatch

    >>> list(zip_exact(range(3), ()))
    Traceback (most recent call last):
         ...
    vgraph.iterstuff.LengthMismatch

    >>> list(zip_exact(range(3), range(2), range(4)))
    Traceback (most recent call last):
         ...
    vgraph.iterstuff.LengthMismatch

    >>> items = zip_exact(iter(range(3)), range(2), range(4))
    >>> next(items)
    (0, 0, 0)
    >>> next(items)
    (1, 1, 1)
    >>> next(items)
    Traceback (most recent call last):
         ...
    vgraph.iterstuff.LengthMismatch
    """
    if not iterables:
        return iter(iterables)
    elif len(iterables) == 1:
        return iter(iterables[0])

    first = iterables[0]
    rest  = iterables[1:]

    try:
        n = len(first)
        if all(len(i) == n for i in rest):
            return zip(*iterables)
    except (TypeError, AttributeError):
        pass

    # Must use sentinel objects to enforce length equality
    rest  = [chain(i, _zip_exact_thow()) for i in rest]
    first = chain(first, _zip_exact_check(rest))
    return zip(*[first] + rest)


class OrderError(ValueError):
    pass


def sort_almost_sorted(iterable, key=None, windowsize=1000, stable=True):
    """
    sort_almost_sorted(iterable, key=None, windowsize=1000, stable=True)

    Sorts an almost sorted iterable of items provided that all misordered
    items are within windowsize elements of the correct location in the final
    iterable.  Returns a generator that yields the correctly sorted iterable
    or raises OrderError if the correct iterable cannot be constructed with
    the given window size.

    If a key function is provided or the stable argument is True, then the
    resulting sort order is stable, with items returned in proper sort order
    and equal keys returned in the same relative order.  Otherwise, the
    resulting ordering is not guaranteed to be stable.

    Test window sizes:

        >>> list(sort_almost_sorted([1,2,4,5,3],windowsize=3))
        [1, 2, 3, 4, 5]

        >>> list(sort_almost_sorted([1,2,5,6,4],windowsize=1))
        Traceback (most recent call last):
             ...
        vgraph.iterstuff.OrderError: Misordered keys beyond window size

        >>> list(sort_almost_sorted([2,3,4,5,6,7,9,10,1],windowsize=10))
        [1, 2, 3, 4, 5, 6, 7, 9, 10]

        >>> list(sort_almost_sorted([2,3,4,5,6,7,9,10,1],windowsize=8))
        Traceback (most recent call last):
        ...
        vgraph.iterstuff.OrderError: Misordered keys beyond window size

    Test stability:

        # key=int(x), so all keys will compare equal and values are not in the
        # natural ordering.
        >>> list(sort_almost_sorted([1.7,1.5,1.6,1.4], key=lambda x: int(x)))
        [1.7, 1.5, 1.6, 1.4]

    Test key function:

        >>> list(sort_almost_sorted([1,2,3,4], key=lambda x: -x))
        [4, 3, 2, 1]
    """
    from operator import itemgetter

    # Key sorts are always stable, since we already pay the price for a
    # decorated iterable.
    if key is not None:
        decorated = ((key(item), i, item) for i, item in enumerate(iterable))
        ordered   = _sort_almost_sorted(decorated, windowsize)
        return map(itemgetter(2), ordered)

    # Otherwise, use a similar method as above to ensure stability
    elif stable:
        decorated = ((item, i) for i, item in enumerate(iterable))
        ordered   = _sort_almost_sorted(decorated, windowsize)
        return map(itemgetter(0), ordered)

    # Unstable, undecorated sort
    else:
        return _sort_almost_sorted(iterable, windowsize)


def _sort_almost_sorted(iterable, windowsize):
    """
    Internal function.  See sort_almost_sorted
    """
    from heapq import heapify, heappushpop, heappop

    # STAGE 1: Fill initial window and heapify
    it   = iter(iterable)
    heap = list(islice(it, windowsize))

    # Establish invariant len(heap)>0
    if not heap:
        return

    heapify(heap)

    # STAGE 2: Slide window until end of iterable
    last = heap[0]

    # Loop invariants:
    #   len(heap)==c, where c is a constant 0 < c <= windowsize
    #   all(e>=last for e in heap)
    for item in it:
        item = heappushpop(heap, item)

        if item < last:
            raise OrderError('Misordered keys beyond window size')

        last = item

        yield item

    # Release reference to last item, since it is no longer needed
    del last

    # STAGE 3: Drain window, no need to check sort order of remaining elements
    #          since remaining elements must be within windowsize distance
    while heap:
        yield heappop(heap)


def ensure_ordered(iterable, key=None):
    """
    ensure_ordered(iterable, key=None) -> iterable

    Returns a generator that yields all elements of iterable, provided that
    the elements are sorted in non-descending order.  Otherwise an OrderError
    exception is raised.

    >>> list(ensure_ordered([1, 2, 3]))
    [1, 2, 3]

    >>> list(ensure_ordered([3, 2, 1]))
    Traceback (most recent call last):
         ...
    vgraph.iterstuff.OrderError: Invalid sort order

    >>> list(ensure_ordered([1.7, 1.5, 1.6, 1.4], key=lambda x: int(x)))
    [1.7, 1.5, 1.6, 1.4]
    """
    it = iter(iterable)

    try:
        last = next(it)
    except StopIteration:
        return

    yield last

    if key is None:
        for item in it:
            if item < last:
                raise OrderError('Invalid sort order')
            last = item
            yield item
    else:
        last = key(last)

        for item in it:
            current = key(item)
            if current < last:
                raise OrderError('Invalid sort order')
            last = current
            yield item


def ensure_unique_everseen(iterable, key=None):
    """
    ensure_unique_everseen(iterable, key=None) -> iterable

    Returns a generator that yields all elements of iterable, provided that
    the elements are unique based on hashability and equality.  Otherwise a
    ValueError exception is raised.

    >>> list(ensure_unique_everseen([1, 2, 3]))
    [1, 2, 3]

    >>> list(ensure_unique_everseen([1, 2, 1]))
    Traceback (most recent call last):
         ...
    ValueError: non-unique element detected

    >>> list(ensure_unique_everseen([1.7, 1.5, 1.6, 1.4]))
    [1.7, 1.5, 1.6, 1.4]

    >>> list(ensure_unique_everseen([1.7, 1.5, 1.6, 1.4], key=lambda x: int(x)))
    Traceback (most recent call last):
         ...
    ValueError: non-unique element detected
    """
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in iterable:
            if element in seen:
                raise ValueError('non-unique element detected')
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k in seen:
                raise ValueError('non-unique element detected')
            seen_add(k)
            yield element
