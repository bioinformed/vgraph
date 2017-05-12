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


cdef dict seq_normalization_table = str.maketrans('ACGTNacgtNRYSWKMBDHV', 'ACGTNACGTNNNNNNNNNNN')


cpdef str normalize_seq(str seq):
    """Normalize DNA sequence for case and IUPAC ambiguity codes.

    Mapping table::

        Un-normalized: ACGTacgtNnRYSWKMBDHV
        Normalized:    ACGTACGTNNNNNNNNNNNN

    Args:
        seq (str): DNA sequence containing only ACGTNacgtNRYSWKMBDHV characters

    Returns:
        str: Normalized sequence

    """
    return seq.translate(seq_normalization_table)


cdef inline bint intersects(int a_start, int a_stop, int b_start, int b_stop):
    '''
    Test if interval [start0, stop0) overlaps [start1, stop1) where [a, b) repesents a zero-based
    half-open interval where a <= b.

    Args:
        a_start (int): start coordinate of first interval (int, 0-based, inclusive)
        a_stop  (int): stop coordinate of first interval (int, 0-based, exclusive)
        b_start (int): start coordinate of second interval (int, 0-based, inclusive)
        b_stop  (int): stop coordinate of second interval (int, 0-based, exclusive)

    Returns:
        bool: True if intersect
    '''
    assert a_start <= a_stop
    assert b_start <= b_stop

    return (a_start == a_stop == b_start == b_stop) or a_stop > b_start and a_start < b_stop


cpdef fancy_match(str s1, str s2):
    cdef bytes b1 = s1.encode()
    cdef bytes b2 = s2.encode()
    cdef int n1 = len(b1)
    cdef int n2 = len(b2)
    cdef int i

    if n1 != n2:
        return False

    cdef char *c1 = <char*>b1
    cdef char *c2 = <char*>b2

    for i in range(n1):
        if c1[i] == b'*' or c2[i] == b'*':
            pass
        elif c1[i] == b'.' or c2[i] == b'.':
            return None
        elif c1[i] != c2[i]:
            return False
    return True


cpdef trim_common_suffixes(strs, int max_trim=0):
    '''trim common suffixes

    >>> trim_common_suffixes([])
    (0, [])
    >>> trim_common_suffixes(['AA'])
    (0, ['AA'])
    >>> trim_common_suffixes(['A','BA'])
    (1, ['', 'B'])
    >>> trim_common_suffixes(['A','AB'])
    (0, ['A', 'AB'])
    >>> trim_common_suffixes(['A','BA','BAAA'])
    (1, ['', 'B', 'BAA'])
    >>> trim_common_suffixes(['BA','ABA','AABA'])
    (2, ['', 'A', 'AA'])
    >>> trim_common_suffixes(['AB','ABA','ABAA','C'])
    (0, ['AB', 'ABA', 'ABAA', 'C'])
    >>> trim_common_suffixes(['AAAAAAAAA','AAAAAAAAA','AAAAAAAAA','AAAAAAAAA'])
    (9, ['', '', '', ''])

    '''
    if len(strs) <= 1:
        return 0, strs

    cdef int i, trim = 0

    cdef str s
    cdef bytes str0 = strs[0].encode()
    cdef bytes str1 = strs[1].encode()
    cdef char *s0 = <char*>str0
    cdef char *s1 = <char*>str1
    cdef int n0 = len(str0) - 1
    cdef int n1 = len(str1) - 1

    while n0 >= 0 and n1 >= 0 and s0[n0] == s1[n1]:
        trim += 1
        n0 -= 1
        n1 -= 1

    for s in strs[2:]:
        str1 = s.encode()
        s1 = <char*>str1
        n0 = len(str0) - 1
        n1 = len(str1) - 1
        i = 0
        while n0 >= 0 and n1 >= 0 and i < trim and s0[n0] == s1[n1]:
            i  += 1
            n0 -= 1
            n1 -= 1

        if i < trim:
            trim = i

    if 0 < max_trim < trim:
        trim = max_trim

    if trim > 0:
        strs = [s[:-trim] for s in strs]

    return trim, strs


cpdef trim_common_prefixes(strs, int max_trim=0):
    '''trim common prefixes

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
    if len(strs) <= 1:
        return 0, strs

    cdef int i, trim = 0
    cdef str s
    cdef bytes str0 = strs[0].encode()
    cdef bytes str1 = strs[1].encode()
    cdef char *s0 = <char*>str0
    cdef char *s1 = <char*>str1

    while s0[trim] != 0 and s1[trim] != 0 and s0[trim] == s1[trim]:
        trim += 1

    for s in strs[2:]:
        str1 = s.encode()
        s1 = <char*>str1
        i = 0
        while i < trim and s1[i] != 0 and s0[i] == s1[i]:
            i += 1
        if i < trim:
            trim = i

    if 0 < max_trim < trim:
        trim = max_trim

    if trim > 0:
        strs = [s[trim:] for s in strs]

    return trim, strs


cdef shuffle_left(str ref, int *start, int *stop, alleles, int bound, int ref_step):
    cdef int trimmed, step, left, n = len(alleles)

    while 0 < alleles.count('') < n and start[0] > bound:
        step = min(ref_step, start[0] - bound)

        r = normalize_seq(ref[start[0] - step:start[0]])
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
        r = normalize_seq(ref[stop[0]:stop[0]+step])
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

    alleles = list(map(normalize_seq, alleles))

    if alleles[0] != normalize_seq(ref[start:stop]):
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

    alleles = list(map(normalize_seq, alleles))

    if alleles[0] != normalize_seq(ref[start:stop]):
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
    __slots__ = ('recnum', 'record', 'contig', 'start', 'stop', 'alleles', 'allele_indices', 'phased',
                 'min_start', 'max_stop', 'left', 'right')

    def __init__(self, recnum, record, ref, name=None, variant_padding=0):
        self.recnum = recnum
        self.record = record
        self.contig = record.contig

        refa = normalize_seq(ref[record.start:record.stop])
        rec_ref = normalize_seq(record.ref)

        if rec_ref != refa[0] and rec_ref != refa:
            raise ReferenceMismatch('Reference mismatch at {}:{}-{}, found={}, expected={}'
                      .format(record.contig, record.start + 1, record.stop, rec_ref, refa))

        if name is not None:
            sample = record.samples[name]
            self.allele_indices = sample.allele_indices
            self.phased = sample.phased
        else:
            self.phased = False
            self.allele_indices = (0, 1)

        geno_alleles = (record.alleles[0],) + tuple(a for i, a in enumerate(record.alleles[1:], 1) if i in self.allele_indices)
        self.allele_indices = tuple(geno_alleles.index(record.alleles[i]) if i is not None else None for i in self.allele_indices)
        self.alleles = geno_alleles

        # Detect and handle refcalls and no-calls
        ploidy         = len(self.allele_indices)
        ref_alleles    = self.allele_indices.count(0)
        nocall_alleles = self.allele_indices.count(None)

        #else:
        #    self.alleles = record.alleles
        #    self.allele_indices = self.phased = None
        #    ploidy = ref_alleles = nocall_alleles = 0

        if ref_alleles + nocall_alleles == ploidy:
            self.left = self.right = self
            self.start, self.stop = record.start, record.stop
            self.min_start = record.start
            self.max_stop  = record.stop
            return

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

        self.min_start = max(0, min(chain.from_iterable(lefts)) - variant_padding)
        self.max_stop  = max(chain.from_iterable(rights)) + variant_padding

    def extreme_order_key(self):
        return self.min_start, self.max_stop

    def left_order_key(self):
        return self.left.start, self.recnum

    def natural_order_key(self):
        return self.start, self.recnum

    def record_order_key(self):
        return self.recnum

    def intersects(self, other):
        return intersects(self.start, self.stop, other.start, other.stop)

    def extremes_intersect(self, other):
        return intersects(self.min_start, self.max_stop, other.min_start, other.max_stop)

    def is_ref(self):
        return len(self.allele_indices) == self.allele_indices.count(0)
