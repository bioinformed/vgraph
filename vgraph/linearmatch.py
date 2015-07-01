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


from __future__       import division, print_function

from itertools        import combinations, combinations_with_replacement
from collections      import namedtuple, defaultdict, Counter

from vgraph.norm      import NormalizedLocus


RefAllele = namedtuple('RefAllele', 'start stop seq phase')
AltAllele = namedtuple('AltAllele', 'start stop seq phase')


def is_valid_path(phase_constraints, path):
    if phase_constraints:
        pathset = set(p for p in path if isinstance(p, AltAllele))
        for phaseset in phase_constraints.itervalues():
            n = len(pathset & phaseset)
            if n != 0 and n != len(phaseset):
                raise ValueError('should not happen')
                return False
    return True


def is_valid_geno(zygosity_constraints, paths):
    observed_zygosity = Counter(p for path in paths for p in path if isinstance(p, AltAllele))
    return zygosity_constraints == observed_zygosity


def get_path_seq(path):
    return ''.join(p.seq for p in path)


def generate_linear_graph(ref, start, stop, loci, name):
    zygosity_constraints = defaultdict(int)
    phase_constraints = defaultdict(set)
    graph = []
    pos = start

    for locus in sorted(loci, key=NormalizedLocus.left_order_key):
        sample = locus.record.samples[name]
        indices = sample.allele_indices
        phased = sample.phased
        left = locus.left

        if pos < left.start:
            graph.append([RefAllele(pos, left.start, ref[pos:left.start], None)])
        elif pos > left.start:
            raise ValueError('unordered locus previous start={}, current start={}'.format(pos, left.start))

        pos = left.start
        alleles = []
        for i in set(indices):
            if i:
                # Add phase constraint for hets only
                phasename = None
                if phased and len(indices) == 2 and indices[0] != indices[1]:
                    phasename = 'M' if indices[0] == i else 'F'

                allele = AltAllele(left.start, left.stop, left.alleles[i], phasename)
                alleles.append(allele)

                # Add zygosity constraint
                zygosity_constraints[allele] = indices.count(i)

                if phasename:
                    phase_constraints[phasename].add(allele)

            else:
                alleles.append(RefAllele(left.start, left.stop, ref[left.start:left.stop], None))

        graph.append(alleles)
        pos = left.stop

    if pos < stop:
        graph.append([RefAllele(pos, stop, ref[pos:stop], None)])

    return graph, zygosity_constraints, phase_constraints


def generate_paths(graph):
    # queue of (seq, start, stop, alts, phasesets, antiphasesets)
    paths = [('', [], set(), set())]

    #print()
    #print('  PATHS')
    for i, alleles in enumerate(graph, 1):
        paths = _extend_paths(paths, alleles)
        #paths = list(paths)
        #print('    {:d} Alleles: {}'.format(i, alleles))
        #for j, (seq, path, ps, aps) in enumerate(paths, 1):
        #    print('      {:d}.{:d}: seq={} ps={} aps={} alts={}'.format(i, j, seq, ps, aps, [a for a in path if a.phase]))
        #print()
    paths = (p[1] for p in paths)
    return paths


def _extend_paths(inpaths, alleles):
    for seq, path, phasesets, antiphasesets in inpaths:
        # Set of new phase sets being added from this allele
        # nb: must occur prior to pruning
        add_phasesets = set(allele.phase for allele in alleles if allele.phase and allele.phase not in phasesets)

        # prune adjacent alleles based on phase and anti-phase constraints
        pruned_alleles = _apply_phase_constrants(alleles, phasesets, antiphasesets)

        for allele in pruned_alleles:
            assert not path or path[-1].stop == allele.start
            assert allele.phase not in antiphasesets
            new_phasesets = _update_phasesets(phasesets, allele.phase)
            new_antiphasesets = _update_antiphasesets(antiphasesets, add_phasesets, allele.phase)
            yield seq + allele.seq, path + [allele], new_phasesets, new_antiphasesets


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


def generate_genotypes(ref, start, stop, loci, name, debug=False):
    graph, zygosity_constraints, phase_constraints = generate_linear_graph(ref, start, stop, loci, name)

    if debug:
        print('-' * 80)
        print('linear VG [{:d}, {:d})'.format(start, stop))
        for i, alleles in enumerate(graph, 1):
            print('  {}: {}'.format(i, alleles))
        print()
        print('ZYGOSITY CONSTRAINTS:', zygosity_constraints)
        print()
        print('PHASE CONSTRAINTS:', phase_constraints)
        print()

    paths = map(tuple, generate_paths(graph))
    assert all(is_valid_path(phase_constraints, path) for path in paths)

    # Any het constraint remove the need to consider homozygous genotypes
    # FIXME: diploid assumptions
    if 1 in zygosity_constraints.values():
        pairs = combinations(paths, 2)
    else:
        pairs = combinations_with_replacement(paths, 2)

    valid_pairs = [pair for pair in pairs if is_valid_geno(zygosity_constraints, pair)]
    valid_genos = sorted(set(tuple(sorted([get_path_seq(p1), get_path_seq(p2)])) for p1, p2 in valid_pairs))

    if debug:
        print('PATHS:')
        for i, path in enumerate(paths, 1):
            assert len(path) == len(set(path))
            print('{:4d}: {}'.format(i, path))
        print()

        print('POSSIBLE HAPLOTYPES:')
        for i, path in enumerate(valid_paths, 1):
            assert len(path) == len(set(path))
            print('{:4d}: {}'.format(i, get_path_seq(path)))
        print()

        print('VALID GENOTYPES:')

        for i, (allele1, allele2) in enumerate(valid_genos, 1):
            print('{:4d}: {}/{}'.format(i, allele1, allele2))
        print()

    return valid_genos
