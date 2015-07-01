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


'''
Normalization of alleles relative to a reference sequence
'''


from __future__       import division, print_function

import sys

from collections import namedtuple
from itertools   import chain


normalized_alleles = namedtuple('shuffled_alleles', 'start stop alleles')


cpdef trim_common_suffixes(strs, int min_len=0, int max_trim=0):
    '''trim common suffixes'''

    if len(strs) < 2:
        return 0, strs

    rev_strs = [ s[::-1] for s in strs ]

    trimmed, rev_strs = trim_common_prefixes(rev_strs, min_len, max_trim)

    if trimmed:
        strs = [ s[::-1] for s in rev_strs ]

    return trimmed, strs


cpdef trim_common_prefixes(strs, int min_len=0, int max_trim=0):
    '''trim common prefixes'''

    cdef int i, trimmed = 0
    cdef bytes s1, s2

    if len(strs) > 1:
        s1 = min(strs)
        s2 = max(strs)

        for i in range(len(s1) - min_len):
            if s1[i] != s2[i]:
                break
            trimmed = i + 1

    if 0 < max_trim < trimmed:
        trimmed = max_trim

    if trimmed > 0:
        strs = [ s[trimmed:] for s in strs ]

    return trimmed, strs


cpdef normalize_alleles(bytes ref, int start, int stop, alleles, int bound=-1, int ref_step=24, left=True, bint shuffle=True):
    if left:
        if bound < 0:
            bound = 0
        return normalize_alleles_left(ref, start, stop, alleles, bound, ref_step, shuffle)
    else:
        if bound < 0:
            bound = len(ref)
        return normalize_alleles_right(ref, start, stop, alleles, bound, ref_step, shuffle)


cdef normalize_alleles_left(bytes ref, int start, int stop, alleles, int bound, int ref_step=24, bint shuffle=True):
    '''Normalize loci by removing extraneous reference padding'''
    cdef int trimmed

    if alleles[0] != ref[start:stop]:
        raise ValueError('Reference alleles does not match reference sequence: {} != {}'.format(alleles[0], ref[start:stop]))

    if len(alleles) < 2 or start <= 0 or start <= 0:
        return normalized_alleles(start, stop, alleles)

    # STEP 0: Trim prefixes if needed to clear bound
    if start < bound:
        trimmed, alleles = trim_common_prefixes(alleles, max_trim=bound - start)
        start += trimmed

    # STEP 1: Trim common suffix
    trimmed, alleles = trim_common_suffixes(alleles)
    stop -= trimmed

    # STEP 2: Trim common prefix
    trimmed, alleles = trim_common_prefixes(alleles)
    start += trimmed

    #assert bound <= start,'start={:d}, left bound={:d} alleles={}'.format(start, bound, alleles)

    # STEP 3: While a null allele exists, left shuffle by prepending alleles
    #         with reference and trimming common suffixes
    while shuffle and '' in alleles and start > bound:
        step = min(ref_step, start - bound)

        r = ref[start - step:start].upper()
        new_alleles = [ r+a for a in alleles ]

        trimmed, new_alleles = trim_common_suffixes(new_alleles)

        if not trimmed:
            break

        start -= trimmed
        stop  -= trimmed

        if trimmed == step:
            alleles = new_alleles
        else:
            left    = step - trimmed
            alleles = [ a[left:] for a in new_alleles ]
            break

    return normalized_alleles(start, stop, tuple(alleles))


cdef normalize_alleles_right(bytes ref, int start, int stop, alleles, int bound, int ref_step=24, bint shuffle=True):
    '''Normalize loci by removing extraneous reference padding'''
    cdef int trimmed, chrom_stop = len(ref)

    if alleles[0] != ref[start:stop]:
        raise ValueError('Reference alleles does not match reference sequence: {} != {}'.format(alleles[0], ref[start:stop]))

    if len(alleles) < 2 or stop >= chrom_stop:
        return normalized_alleles(start, stop, alleles)

    # STEP 0: Trim suffixes if needed to clear bound
    if stop > bound:
        trimmed, alleles = trim_common_suffixes(alleles, max_trim=stop - bound)
        stop -= trimmed

    # STEP 1: Trim common prefix
    trimmed, alleles = trim_common_prefixes(alleles)
    start += trimmed

    # STEP 2: Trim common suffix
    trimmed, alleles = trim_common_suffixes(alleles)
    stop -= trimmed

    #assert bound >= stop,'stop={:d}, right bound={:d}'.format(stop, bound)

    # STEP 3: While a null allele exists, right shuffle by appending alleles
    #         with reference and trimming common prefixes
    while shuffle and '' in alleles and stop < bound:
        step = min(ref_step, bound - stop)

        r = ref[stop:stop+step].upper()
        new_alleles = [ a+r for a in alleles ]

        trimmed, new_alleles = trim_common_prefixes(new_alleles)

        if not trimmed:
            break

        start += trimmed
        stop  += trimmed

        if trimmed == step:
            alleles = new_alleles
        else:
            left    = step - trimmed
            alleles = [ a[:-left] for a in new_alleles ]
            break

    return normalized_alleles(start, stop, tuple(alleles))


def prefixes(s):
    if not s:
        yield ''

    for i in range(1, len(s) + 1):
        yield s[:i]


def suffixes(s):
    if not s:
        yield ''

    for i in range(1, len(s) + 1):
        yield s[-i:]


class NormalizedLocus(object):
    '''Normalization data for a single VCF record'''
    __slots__ = ('recnum', 'record', 'start', 'stop', 'left', 'right')

    def __init__(self, recnum, record, ref, left_bound=0):
        self.recnum = recnum
        self.record = record

        start, stop, alleles = normalize_alleles(ref, record.start, record.stop, record.alleles, left=True, shuffle=False)
        refa, alts = alleles[0], alleles[1:]

        # Left shuffle locus with all alt alleles considered simultaneously and left bound of previous locus
        self.left = normalize_alleles(ref, start, stop, alleles, bound=left_bound, left=True)

        # Right shuffle locus with all alt alleles considered simultaneously
        self.right = normalize_alleles(ref, start, stop, alleles, left=False)

        # Minimum start and stop coordinates over each alt allele
        # n.b. may be broader than with all alleles or with bounds
        lefts = [[start, self.left.start],
                 (normalize_alleles(ref, start, stop,           (refa, alt),    left=True).start for alt in alts),
                 (normalize_alleles(ref, start, start + len(r), (r,    ''),     left=True).start for r in prefixes(refa) if r),
                 (normalize_alleles(ref, start, start,          ('',   prealt), left=True).start for alt in alts for prealt in prefixes(alt) if prealt)]

        rights = [[stop, self.right.stop],
                  (normalize_alleles(ref, start,         stop, (refa, alt),    left=False).stop for alt in alts),
                  (normalize_alleles(ref, stop - len(r), stop, (r,    ''),     left=False).stop for r in suffixes(refa) if r),
                  (normalize_alleles(ref, stop,          stop, ('',   sufalt), left=False).stop for alt in alts for sufalt in suffixes(alt) if sufalt)]

        self.start = min(chain.from_iterable(lefts))
        self.stop  = max(chain.from_iterable(rights))

    def extreme_order_key(self):
        return self.start, self.stop

    def left_order_key(self):
        return self.left.start, self.recnum
