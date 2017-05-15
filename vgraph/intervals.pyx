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


from heapq    import heappop, heapreplace, heapify, merge as heapq_merge
from operator import attrgetter


def augment_intervals(it, i, key):
    last = None
    for value in it:
        k = key(value)
        if last is not None and k < last:
            raise ValueError('invalid order')
        yield k, i, value


def interval_merge(iterables, key=None):
    if key is None:
        return heapq_merge(*iterables)

    def _augmented_merge(iterables, key):
        augmented = [_augment_items(it, i, key) for i, it in enumerate(iterables)]
        for min_key, index, value in heapq_merge(*augmented):
            yield index, value

    return _augmented_merge(iterables, key)


def _augment_items(it, i, key):
    """Help merge function."""
    last = None
    for value in it:
        k = key(value)
        if last is not None and k < last:
            raise ValueError('invalid order')
        yield k, i, value


def union(record_list, min_distance=0, interval_func=attrgetter('start', 'stop')):
    '''
    >>> from operator import itemgetter
    >>> ifunc = itemgetter(0,1)

    >>> l1 = [(0,1),(1,2),(2,3),(3,4),(4,5)]
    >>> l2 = [(1,2),(3,4)]
    >>> for start, stop, vals in union([l1,l2], 0, ifunc):
    ...     print(start, stop, vals)
    0 1 [[(0, 1)], []]
    1 2 [[(1, 2)], [(1, 2)]]
    2 3 [[(2, 3)], []]
    3 4 [[(3, 4)], [(3, 4)]]
    4 5 [[(4, 5)], []]

    >>> for start, stop, vals in union([l1,l2], 1, ifunc):
    ...     print(start, stop, vals)
    0 5 [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)], [(1, 2), (3, 4)]]

    >>> l1 = [(0,5),(10,15),(20,25)]
    >>> l2 = [(5,11),(21,21)]
    >>> for start, stop, vals in union([l1,l2], 0, ifunc):
    ...     print(start, stop, vals)
    0 5 [[(0, 5)], []]
    5 15 [[(10, 15)], [(5, 11)]]
    20 25 [[(20, 25)], [(21, 21)]]

    >>> l1 = [(0,5),(10,15),(20,25)]
    >>> for start, stop, vals in union([l1, []], 0, ifunc):
    ...     print(start, stop, vals)
    0 5 [[(0, 5)], []]
    10 15 [[(10, 15)], []]
    20 25 [[(20, 25)], []]

    >>> for start, stop, vals in union([[], l1], 0, ifunc):
    ...     print(start, stop, vals)
    0 5 [[], [(0, 5)]]
    10 15 [[], [(10, 15)]]
    20 25 [[], [(20, 25)]]
    '''
    n = len(record_list)
    start = stop = 0
    items = 0
    records = [[] for _ in range(n)]

    for i, rec in interval_merge(record_list, key=interval_func):
        rec_start, rec_stop = interval_func(rec)

        if not records or min_distance <= rec_start - stop:
            if items:
                yield start, stop, records
                items = 0
                records = [[] for _ in range(n)]
            start, stop = rec_start, rec_stop

        assert rec_start >= start
        items += 1
        records[i].append(rec)
        stop = max(stop, rec_stop)

    if items:
        yield start, stop, records
