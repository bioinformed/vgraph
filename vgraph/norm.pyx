# -*- coding: utf-8 -*-
# cython: language_level=3
# cython: embedsignature=True
# cython: profile=True
###############################################################################

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

import sys

from collections import namedtuple
from itertools   import chain


class NormalizationError(ValueError):
    pass


class ReferenceMismatch(NormalizationError):
    pass


NormalizedAlleles = namedtuple('NormalizedAlleles', 'start stop alleles')


cpdef trim_common_suffixes(strs, int min_len=0, int max_trim=0):
    '''trim common suffixes

    >>> trim_common_prefixes([])
    (0, [])
    >>> trim_common_prefixes(['AA'])
    (0, ['AA'])
    >>> trim_common_prefixes(['A','AB'])
    (1, ['', 'B'])
    >>> trim_common_prefixes(['A','BA'])
    (0, ['A', 'BA'])
    >>> trim_common_prefixes(['A','AB','AAAB'])
    (1, ['', 'B', 'AAB'])
    >>> trim_common_prefixes(['AB','ABA','ABAA'])
    (2, ['', 'A', 'AA'])
    >>> trim_common_prefixes(['AB','ABA','ABAA','C'])
    (0, ['AB', 'ABA', 'ABAA', 'C'])
    >>> trim_common_prefixes(['AAAAAAAAA','AAAAAAAAA','AAAAAAAAA','AAAAAAAAA'])
    (9, ['', '', '', ''])
    '''

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
    cdef str s1, s2

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


cdef shuffle_left(str ref, int *start, int *stop, alleles, int bound, int ref_step):
    cdef int trimmed, step, left, n = len(alleles)

    while 0 < alleles.count('') < n and start[0] > bound:
        step = min(ref_step, start[0] - bound)

        r = ref[start[0] - step:start[0]].upper()
        new_alleles = [ r+a for a in alleles ]

        trimmed, new_alleles = trim_common_suffixes(new_alleles)

        if not trimmed:
            break

        start[0] -= trimmed
        stop[0]  -= trimmed

        if trimmed == step:
            alleles = new_alleles
        else:
            left    = step - trimmed
            alleles = [ a[left:] for a in new_alleles ]
            break

    return alleles


cdef shuffle_right(str ref, int *start, int *stop, alleles, int bound, int ref_step):
    cdef int trimmed, step, left, n = len(alleles)

    while 0 < alleles.count('') < n and stop[0] < bound:
        step = min(ref_step, bound - stop[0])
        r = ref[stop[0]:stop[0]+step].upper()
        new_alleles = [ a+r for a in alleles ]

        trimmed, new_alleles = trim_common_prefixes(new_alleles)

        if not trimmed:
            break

        start[0] += trimmed
        stop[0]  += trimmed

        if trimmed == step:
            alleles = new_alleles
        else:
            left    = step - trimmed
            alleles = [ a[:-left] for a in new_alleles ]
            break

    return alleles


cpdef normalize_alleles(str ref, int start, int stop, alleles, int bound=-1, int ref_step=24, left=True, bint shuffle=True):
    if left:
        if bound < 0:
            bound = 0
        return normalize_alleles_left(ref, start, stop, alleles, bound, ref_step, shuffle)
    else:
        if bound < 0:
            bound = len(ref)
        return normalize_alleles_right(ref, start, stop, alleles, bound, ref_step, shuffle)


cdef normalize_alleles_left(str ref, int start, int stop, alleles, int bound, int ref_step=24, bint shuffle=True):
    '''Normalize loci by removing extraneous reference padding'''
    cdef int trimmed

    if alleles[0] != ref[start:stop]:
        raise ReferenceMismatch('Reference alleles does not match reference sequence: {} != {}'.format(alleles[0], ref[start:stop]))

    if len(alleles) < 2 or start <= 0:
        return NormalizedAlleles(start, stop, alleles)

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

    # STEP 3: Force shuffle right if start doesn't clear bound
    if start < bound:
        alleles = shuffle_right(ref, &start, &stop, alleles, stop + bound - start, ref_step)

    #assert bound <= start,'start={:d}, left bound={:d} alleles={}'.format(start, bound, alleles)

    # STEP 4: While a null allele exists, left shuffle by prepending alleles
    #         with reference and trimming common suffixes
    if shuffle:
        alleles = shuffle_left(ref, &start, &stop, alleles, bound, ref_step)

    return NormalizedAlleles(start, stop, tuple(alleles))


cdef normalize_alleles_right(str ref, int start, int stop, alleles, int bound, int ref_step=24, bint shuffle=True):
    '''Normalize loci by removing extraneous reference padding'''
    cdef int trimmed, chrom_stop = len(ref)

    if alleles[0] != ref[start:stop]:
        raise ReferenceMismatch('Reference alleles does not match reference sequence: {} != {}'.format(alleles[0], ref[start:stop]))

    if len(alleles) < 2 or stop >= chrom_stop:
        return NormalizedAlleles(start, stop, alleles)

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

    # STEP 3: Force shuffle left if stop doesn't clear bound
    if stop > bound:
        alleles = shuffle_left(ref, &start, &stop, alleles, start - stop - bound, ref_step)

    #assert bound >= stop,'stop={:d}, right bound={:d}'.format(stop, bound)

    # STEP 4: While a null allele exists, right shuffle by appending alleles
    #         with reference and trimming common prefixes
    if shuffle:
        alleles = shuffle_right(ref, &start, &stop, alleles, bound, ref_step)

    return NormalizedAlleles(start, stop, tuple(alleles))


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
    '''Normalization data for a single VCF record and genotype'''
    __slots__ = ('recnum', 'record', 'start', 'stop', 'alleles', 'allele_indices', 'phased',
                 'min_start', 'max_stop', 'left', 'right')

    def __init__(self, recnum, record, ref, name=None):
        self.recnum = recnum
        self.record = record

        refa = ref[record.start:record.stop]
        if record.alleles[0] != refa:
            raise ReferenceMismatch('Reference mismatch at {}:{}-{}, found={}, expected={}'
                      .format(record.contig, record.start + 1, record.stop, record.alleles[0], refa))

        if name is not None:
            sample = record.samples[name]
            self.allele_indices = sample.allele_indices
            self.phased = sample.phased

            geno_alleles = (record.alleles[0],) + tuple(a for i, a in enumerate(record.alleles[1:], 1) if i in self.allele_indices)
            self.allele_indices = tuple(geno_alleles.index(record.alleles[i]) for i in self.allele_indices)
            self.alleles = geno_alleles
        else:
            self.alleles = record.alleles
            self.allele_indices = self.phased = None

        # Left shuffle locus with all alt alleles considered simultaneously
        self.left = normalize_alleles(ref, record.start, record.stop, self.alleles, left=True)

        # Right shuffle locus with all alt alleles considered simultaneously
        self.right = normalize_alleles(ref, record.start, record.stop, self.alleles, left=False)

        start, stop, alleles = normalize_alleles(ref, record.start, record.stop, self.alleles, left=True, shuffle=False)
        refa, alts = alleles[0], alleles[1:]
        self.start, self.stop, self.alleles = start, stop, alleles

        # Minimum start and stop coordinates over each alt allele
        # n.b. may be broader than with all alleles or with bounds
        lefts = [[start, record.start, self.left.start],
                 (normalize_alleles(ref, start, stop,           (refa, alt),    left=True).start for alt in alts if refa),
                 (normalize_alleles(ref, start, start + len(r), (r,    ''),     left=True).start for r in prefixes(refa) if r),
                 (normalize_alleles(ref, start, start,          ('',   prealt), left=True).start for alt in alts for prealt in prefixes(alt) if prealt)]

        rights = [[stop, record.stop, self.right.stop],
                  (normalize_alleles(ref, start,         stop, (refa, alt),    left=False).stop for alt in alts if refa),
                  (normalize_alleles(ref, stop - len(r), stop, (r,    ''),     left=False).stop for r in suffixes(refa) if r),
                  (normalize_alleles(ref, stop,          stop, ('',   sufalt), left=False).stop for alt in alts for sufalt in suffixes(alt) if sufalt)]

        self.min_start = min(chain.from_iterable(lefts))
        self.max_stop  = max(chain.from_iterable(rights))

    def extreme_order_key(self):
        return self.min_start, self.max_stop

    def left_order_key(self):
        return self.left.start, self.recnum

    def natural_order_key(self):
        return self.start, self.recnum

    def record_order_key(self):
        return self.recnum
