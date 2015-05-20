'''
Normalization of alleles relative to a reference sequence
'''
import sys
from collections import namedtuple


normalized_alleles = namedtuple('shuffled_alleles', 'start stop alleles')


cpdef trim_common_suffixes(strs, int min_len=0):
    '''trim common suffixes'''

    if len(strs) < 2:
        return 0, strs

    rev_strs = [ s[::-1] for s in strs ]

    trimmed, rev_strs = trim_common_prefixes(rev_strs, min_len)

    if trimmed:
        strs = [ s[::-1] for s in rev_strs ]

    return trimmed, strs


cpdef trim_common_prefixes(strs, int min_len=0):
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

    assert alleles[0] == ref[start:stop]

    if len(alleles) < 2 or start <= 0 or start <= 0:
        return normalized_alleles(start, stop, alleles)

    # STEP 1: Trim common suffix
    trimmed, alleles = trim_common_suffixes(alleles)
    stop -= trimmed

    # STEP 2: Trim common prefix
    trimmed, alleles = trim_common_prefixes(alleles)
    start += trimmed

    #assert bound <= start,'start={:d}, left bound={:d}'.format(start, bound)

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

    assert alleles[0] == ref[start:stop]

    if len(alleles) < 2 or stop >= chrom_stop:
        return normalized_alleles(start, stop, alleles)

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
