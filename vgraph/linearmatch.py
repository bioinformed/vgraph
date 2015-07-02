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


from __future__  import division, print_function

from itertools   import combinations, combinations_with_replacement
from collections import defaultdict


class RefAllele(object):
    __slots__ = ('seq',)
    phase = None

    def __init__(self, seq):
        self.seq = seq

    def __repr__(self):
        return 'RefAllele({})'.format(self.seq or '-')


class AltAllele(object):
    __slots__ = ('seq', 'phase')

    def __init__(self, seq, phase):
        self.seq = seq
        self.phase = phase

    def __repr__(self):
        if self.phase is None:
            return 'AltAllele({})'.format(self.seq or '-')
        else:
            return 'AltAllele({}, phase={})'.format(self.seq or '-', self.phase)


def is_valid_geno(zygosity_constraints, alts1, alts2):
    observed_zygosity = defaultdict(int)
    for allele in alts1:
        if isinstance(allele, AltAllele):
            observed_zygosity[allele] += 1
    for allele in alts2:
        if isinstance(allele, AltAllele):
            observed_zygosity[allele] += 1
    return zygosity_constraints == observed_zygosity


def generate_graph(ref, start, stop, loci, name, debug=False):
    zygosity_constraints = defaultdict(int)
    graph = []
    pos = start

    for locus in loci:
        sample = locus.record.samples[name]
        indices = sample.allele_indices
        phased = sample.phased
        left = locus.left
        assert ref[left.start:left.stop] == left.alleles[0]

        if pos < left.start:
            graph.append([RefAllele(ref[pos:left.start])])
        elif pos > left.start:
            raise ValueError('unordered locus previous start={}, current start={}'.format(pos, left.start))

        alleles = _make_alleles(left.alleles, indices, phased, zygosity_constraints)
        graph.append(list(alleles))
        pos = left.stop

    if pos < stop:
        graph.append([RefAllele(ref[pos:stop])])

    if debug:
        print('-' * 80)
        print('linear VG [{:d}, {:d})'.format(start, stop))
        for i, alleles in enumerate(graph, 1):
            print('  {}: {}'.format(i, alleles))
        print()
        print('ZYGOSITY CONSTRAINTS:', zygosity_constraints)
        print()

    return graph, zygosity_constraints


def _make_alleles(alleles, indices, phased, zygosity_constraints):
    index_set = set(indices)
    het = len(index_set) > 1

    for i in index_set:
        # Alt allele at phased and heterozygous locus
        if i and phased and het:
            # each alt allele is distinct
            for phasename in (j for j, idx in enumerate(indices) if i == idx):
                allele = AltAllele(alleles[i], phasename)
                zygosity_constraints[allele] = 1
                yield allele

        # Alt allele at unphased or homozygous locus
        elif i:
            # single alt allele
            allele = AltAllele(alleles[i], None)
            zygosity_constraints[allele] = indices.count(i)
            yield allele

        # Ref allele
        else:
            yield RefAllele(alleles[0])


def generate_paths(graph, debug=False):
    # Initial path of (seq, alts, phasesets, antiphasesets)
    paths = [('', [], set(), set())]

    for alleles in graph:
        paths = _extend_paths(paths, alleles)

    paths = list(tuple(p[:2]) for p in paths)

    return paths


def _extend_paths(inpaths, alleles):
    for seq, alts, phasesets, antiphasesets in inpaths:
        # Set of new phase sets being added from this allele
        # nb: must occur prior to pruning
        add_phasesets = set(allele.phase for allele in alleles if allele.phase and allele.phase not in phasesets)

        # prune adjacent alleles based on phase and anti-phase constraints
        pruned_alleles = _apply_phase_constrants(alleles, phasesets, antiphasesets)

        for allele in pruned_alleles:
            new_phasesets = _update_phasesets(phasesets, allele.phase)
            new_antiphasesets = _update_antiphasesets(antiphasesets, add_phasesets, allele.phase)
            new_alts = alts + [allele] if isinstance(allele, AltAllele) else alts
            yield seq + allele.seq, new_alts, new_phasesets, new_antiphasesets


def _apply_phase_constrants(alleles, phasesets, antiphasesets):
    if phasesets:
        # If any adjacent allele belongs to an already selected phase set,
        # then allow only those alleles.  Otherwise, all alleles are valid.
        alleles = [allele for allele in alleles if allele.phase in phasesets] or alleles

    if antiphasesets:
        # Remove any allele belonging to an anti-phaseset
        alleles = [allele for allele in alleles if allele.phase not in antiphasesets]

    return alleles


def _update_phasesets(phasesets, phaseset):
    if phaseset is not None and phaseset not in phasesets:
        # only copy when adding a new phaseset
        phasesets = phasesets.copy()
        phasesets.add(phaseset)
    return phasesets


def _update_antiphasesets(antiphasesets, add_phasesets, phaseset):
    if add_phasesets:
        # Copy if adding a new anti-phaseset, one not equal to the current allele's phaseset
        add_anti = set(p for p in add_phasesets if p != phaseset and p not in antiphasesets)
        if add_anti:
            antiphasesets = antiphasesets | add_anti
    return antiphasesets


def intersect_paths(paths1, paths2):
    if not isinstance(paths1, list):
        paths1 = list(paths1)
    if not isinstance(paths2, list):
        paths2 = list(paths2)

    index1 = set(seq for seq, alts in paths1)
    index  = set(seq for seq, alts in paths2 if seq in index1)

    paths1 = (p for p in paths1 if p[0] in index)
    paths2 = (p for p in paths2 if p[0] in index)

    return paths1, paths2


def generate_genotypes(paths, zygosity_constraints, debug=False):
    # Any het constraint remove the need to consider homozygous genotypes
    # FIXME: diploid assumption
    if debug and not isinstance(paths, list):
        paths = list(paths)

    if 1 in zygosity_constraints.values():
        pairs = combinations(paths, 2)
    else:
        pairs = combinations_with_replacement(paths, 2)

    genos = ((seq1, seq2) if seq1 <= seq2 else (seq2, seq1)
             for (seq1, alts1), (seq2, alts2) in pairs
             if is_valid_geno(zygosity_constraints, alts1, alts2))

    genos = sorted(set(genos))

    if debug:
        print('PATHS:')
        for i, (seq, path) in enumerate(paths, 1):
            assert len(path) == len(set(path))
            print('{:4d}: {}'.format(i, seq))
        print()

        print('POSSIBLE HAPLOTYPES:')
        for i, (seq, path) in enumerate(paths, 1):
            assert len(path) == len(set(path))
            print('{:4d}: {}'.format(i, seq))
        print()

        print('GENOTYPES:')

        for i, (allele1, allele2) in enumerate(genos, 1):
            print('{:4d}: {}/{}'.format(i, allele1, allele2))
        print()

    return genos
