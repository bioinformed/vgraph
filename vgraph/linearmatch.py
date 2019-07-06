# Copyright 2015 Kevin B Jacobs
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License.  You may obtain
# a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations
# under the License.


"""Match variant acyclic variant graphs."""


import sys

from itertools   import combinations, combinations_with_replacement
from collections import defaultdict

from vgraph.norm import normalize_seq


EMPTY_PATH = ('', [], set(), set())


TRIM_MIN, TRIM_MARGIN = 15, 5


class OverlapError(ValueError):
    """Exception raised when alleles overlap."""
    pass


def trim_seq(seq):
    """Trim long sequences for convenient display formatting."""
    if len(seq) >= TRIM_MIN:
        return '{}...{}'.format(seq[:TRIM_MARGIN], seq[-TRIM_MARGIN:])
    else:
        return seq


def trim_ref(ref, start, stop):
    """Fill in reference sequence, trimming if excessively long."""
    if stop - start >= TRIM_MIN:
        return '{}...{}'.format(ref[start:start + TRIM_MARGIN], ref[stop - TRIM_MARGIN:stop])
    else:
        return ref[start:stop]


class RefAllele:
    """Graph node for a reference haplotype."""
    __slots__ = ('start', 'stop', 'ref')
    phase = None

    def __init__(self, start, stop, ref):
        """Build new RefAllele."""
        self.start = start
        self.stop = stop
        self.ref = ref

    @property
    def seq(self):
        """Return the sequence of this allele."""
        return normalize_seq(self.ref[self.start:self.stop])

    def __len__(self):
        """Return the length of this allele."""
        return self.stop - self.start

    def __repr__(self):
        """Return a convenient string representation of this allele."""
        seq = trim_ref(self.ref, self.start, self.stop) or '-'
        return 'RefAllele({}, {}, {})'.format(self.start, self.stop, seq)


class HomAltAllele:
    """Graph node for a homozygous alternative haplotype."""
    __slots__ = ('start', 'stop', 'seq')
    phase = None

    def __init__(self, start, stop, seq):
        """Build new HomAltAllele."""
        self.start = start
        self.stop = stop
        self.seq = normalize_seq(seq)

    def __len__(self):
        """Return the length of this allele."""
        return self.stop - self.start

    def __repr__(self):
        """Return a convenient string representation of this allele."""
        seq = trim_seq(self.seq) or '-'
        return 'HomAltAllele({}, {}, {})'.format(self.start, self.stop, seq)


class NocallAllele:
    """Graph node for a no-call haplotype."""
    __slots__ = ('start', 'stop')
    phase = None

    def __init__(self, start, stop):
        """Build new NocallAllele."""
        self.start = start
        self.stop = stop

    def __len__(self):
        """Return the length of this allele."""
        return self.stop - self.start

    @property
    def seq(self):
        """Return the sequence of this allele."""
        return '.' * (self.stop - self.start)

    def __repr__(self):
        """Return a convenient string representation of this allele."""
        n = len(self)
        if n < TRIM_MIN:
            seq = '.' * n or '-'
        else:
            seq = '...'
        return 'NocallAllele({}, {}, {})'.format(self.start, self.stop, seq)


class HetAltAllele:
    """Graph node for a heterozygous alternative haplotype."""
    __slots__ = ('start', 'stop', 'seq', 'phase')

    def __init__(self, start, stop, seq, phase=None):
        """Build new HetAltAllele."""
        self.start = start
        self.stop = stop
        self.seq = normalize_seq(seq)
        self.phase = phase

    def __len__(self):
        """Return the length of this allele."""
        return len(self.seq)

    def __repr__(self):
        """Return a convenient string representation of this allele."""
        seq = trim_seq(self.seq) or '-'
        if self.phase is None:
            return 'HetAltAllele({}, {}, {})'.format(self.start, self.stop, seq)
        else:
            return 'HetAltAllele({}, {}, {}, phase={})'.format(self.start, self.stop, seq, self.phase)


def is_valid_geno(zygosity_constraints, alts1, alts2):
    """Check if a genotype is valid given zygosity constraints."""
    if not zygosity_constraints:
        return True
    observed_zygosity = defaultdict(int)
    for allele in alts1:
        observed_zygosity[allele] += 1
    for allele in alts2:
        observed_zygosity[allele] += 1
    return zygosity_constraints == observed_zygosity


def generate_graph(ref, start, stop, loci, debug=False):
    """Generate an acyclic variant graph."""
    zygosity_constraints = defaultdict(int)

    last = start
    for locus in loci:
        assert(locus.start >= last)
        last = locus.start

    def _generate_graph():
        pos = start

        for locus in loci:
            if pos is not None:
                if pos < locus.start:
                    yield pos, locus.start, [RefAllele(pos, locus.start, ref)]
                elif pos > locus.start:
                    raise OverlapError('overlapping locus: previous stop={}, current start={}'.format(pos, locus.start))

            alleles = _make_alleles(ref, locus, zygosity_constraints)
            yield locus.start, locus.stop, list(alleles)
            pos = locus.stop

        if stop is not None and pos < stop:
            yield pos, stop, [RefAllele(pos, stop, ref)]

    return _generate_graph(), zygosity_constraints


def _make_alleles(ref, locus, zygosity_constraints):
    """Make graph alleles."""
    indices = locus.allele_indices

    if None in indices or 'PASS' not in locus.record.filter:
        yield NocallAllele(locus.start, locus.stop)
        return

    index_set = set(indices)
    alleles = locus.alleles
    phased = locus.phased
    het = len(index_set) > 1

    for i in index_set:
        # Alt allele at phased and heterozygous locus
        if i and phased and het:
            # each alt allele is distinct
            for phasename in (j for j, idx in enumerate(indices) if i == idx):
                allele = HetAltAllele(locus.start, locus.stop, alleles[i], phasename)
                zygosity_constraints[allele] = 1
                yield allele

        # Alt allele at unphased or homozygous locus
        elif i:
            # single alt allele
            if het:
                allele = HetAltAllele(locus.start, locus.stop, alleles[i])
                zygosity_constraints[allele] = indices.count(i)
            else:
                allele = HomAltAllele(locus.start, locus.stop, alleles[i])
            yield allele

        # Ref allele
        else:
            yield RefAllele(locus.start, locus.stop, ref)


def generate_paths(graph, feasible_paths=None, debug=False):
    """Generate paths through a variant graph."""
    # if debug:
    #     print('-' * 80)
    #     print('linear VG [{:d}, {:d})'.format(start, stop))
    #     for i, alleles in enumerate(graph, 1):
    #         print('  {}: {}'.format(i, alleles))
    #     print()
    #     print('ZYGOSITY CONSTRAINTS:', zygosity_constraints)
    #     print()

    # Initial path of (seq, alts, phasesets, antiphasesets)
    paths = [EMPTY_PATH]

    for i, (start, stop, alleles) in enumerate(graph):
        if debug:
            print('GRAPH: step={}, start={}, stop={}, alleles={}'.format(i + 1, start, stop, alleles), file=sys.stderr)

        paths = extend_paths(paths, alleles)

        if feasible_paths is not None:
            paths = prune_paths(paths, feasible_paths)

        if debug:
            paths = list(paths)
            for j, p in enumerate(paths):
                print('     PATH{}: {}'.format(j + 1, p), file=sys.stderr)
            print(file=sys.stderr)

    paths = [tuple(p[:2]) for p in paths]

    if debug:
        print('FINAL PATHS:', file=sys.stderr)
        for j, p in enumerate(paths):
            print('     PATH{}: {}'.format(j + 1, p), file=sys.stderr)
        print(file=sys.stderr)

    return paths


def extend_paths(inpaths, alleles):
    """Expand paths in a variant graph with additional alleles."""
    if not inpaths:
        inpaths = [EMPTY_PATH]

    for seq, alts, phasesets, antiphasesets in inpaths:
        # Set of new phase sets being added from this allele
        # nb: must occur prior to pruning
        add_phasesets = {allele.phase for allele in alleles if allele.phase and allele.phase not in phasesets}

        # prune adjacent alleles based on phase and anti-phase constraints
        pruned_alleles = _apply_phase_constrants(alleles, phasesets, antiphasesets)

        for allele in pruned_alleles:
            new_phasesets = _update_phasesets(phasesets, allele.phase)
            new_antiphasesets = _update_antiphasesets(antiphasesets, add_phasesets, allele.phase)
            new_alts = alts + [allele] if isinstance(allele, HetAltAllele) else alts
            yield seq + allele.seq, new_alts, new_phasesets, new_antiphasesets


def _apply_phase_constrants(alleles, phasesets, antiphasesets):
    """Apply phase constraints to alleles."""
    if phasesets:
        # If any adjacent allele belongs to an already selected phase set,
        # then allow only those alleles.  Otherwise, all alleles are valid.
        alleles = [allele for allele in alleles if allele.phase in phasesets] or alleles

    if antiphasesets:
        # Remove any allele belonging to an anti-phaseset
        alleles = [allele for allele in alleles if allele.phase not in antiphasesets]

    return alleles


def _update_phasesets(phasesets, phaseset):
    """Update phase sets."""
    if phaseset is not None and phaseset not in phasesets:
        # only copy when adding a new phaseset
        phasesets = phasesets.copy()
        phasesets.add(phaseset)
    return phasesets


def _update_antiphasesets(antiphasesets, add_phasesets, phaseset):
    """Update anti-phase sets."""
    if add_phasesets:
        # Copy if adding a new anti-phaseset, one not equal to the current allele's phaseset
        add_anti = {p for p in add_phasesets if p != phaseset and p not in antiphasesets}
        if add_anti:
            antiphasesets = antiphasesets | add_anti
    return antiphasesets


def prune_paths(paths, feasible_paths):
    """Prune graph paths."""
    for path in paths:
        p = path[0]
        for f in feasible_paths:
            if f[0].startswith(p):
                yield path
                break


def intersect_paths(paths1, paths2):
    """Intersect graph paths."""
    if not isinstance(paths1, list):
        paths1 = list(paths1)
    if not isinstance(paths2, list):
        paths2 = list(paths2)

    index1 = {seq for seq, alts in paths1}
    index  = {seq for seq, alts in paths2 if seq in index1}

    paths1 = (p for p in paths1 if p[0] in index)
    paths2 = (p for p in paths2 if p[0] in index)

    return paths1, paths2


def generate_genotypes(paths, zygosity_constraints, debug=False):
    """Generate genotypes from feasible graph paths and constraints."""
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
        print('PATHS:', file=sys.stderr)
        for i, (seq, path) in enumerate(paths, 1):
            assert len(path) == len(set(path))
            print('{:4d}: {}'.format(i, seq), file=sys.stderr)
        print(file=sys.stderr)

        print('POSSIBLE HAPLOTYPES:', file=sys.stderr)
        for i, (seq, path) in enumerate(paths, 1):
            assert len(path) == len(set(path))
            print('{:4d}: {}'.format(i, seq), file=sys.stderr)
        print(file=sys.stderr)

        print('GENOTYPES:', file=sys.stderr)

        for i, (allele1, allele2) in enumerate(genos, 1):
            print('{:4d}: {}/{}'.format(i, allele1, allele2), file=sys.stderr)
        print(file=sys.stderr)

    return genos
