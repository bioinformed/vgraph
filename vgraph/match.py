# Copyright 2016 Kevin B Jacobs
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

"""Utilities to perform variant graph matching."""

import sys

from collections        import Counter
from dataclasses        import dataclass
from itertools          import chain
from typing             import Optional

from vgraph.bed         import load_bedmap
from vgraph.norm        import NormalizedLocus, fancy_match, normalize_seq, ReferenceMismatch
from vgraph.intervals   import union
from vgraph.iterstuff   import sort_almost_sorted, is_empty_iter, unique_everseen
from vgraph.lazy_fasta  import LazyFastaContig
from vgraph.linearmatch import generate_graph, generate_paths, generate_genotypes_with_paths, generate_genotypes, intersect_paths, OverlapError


@dataclass(frozen=True)
class AlleleMatch:
    """Dataclass for allele matching results."""
    allele_ploidy: int
    allele_depth:  int
    ref_ploidy:    int
    ref_depth:     Optional[int]
    other_ploidy:  Optional[int]
    other_depth:   Optional[int]


def valid_alleles(alleles):
    """Return if alleles are valid (i.e. at least one non-symbolic alt)."""
    return len(alleles) > 1 and not any('<' in a or '[' in a or ']' in a for a in alleles)


def is_alt_genotype(record, name):
    """Return if the named sample has a non-reference genotype call."""
    sample = record.samples[name]
    indices = sample.allele_indices
    return not (not indices or None in indices or indices.count(0) == len(indices) or max(indices) >= len(record.alleles))


def records_to_loci(ref, records, name, variant_padding):
    """Convert variant records to NormalizedLocus records."""
    for recnum, record in enumerate(records):
        try:
            yield NormalizedLocus(recnum, record, ref, name, variant_padding)
        except ReferenceMismatch:
            print(f'Reference mismatch: {record.contig}:{record.start}-{record.stop}')


def all_contigs(varfiles):
    """Return all contigs in order seen in the variant file header and index.

    Args:
        varfiles: Input variant file object

    Returns:
        All unique contigs of the file

    """
    contigs = list(varfiles.header.contigs)
    if varfiles.index is not None:
        contigs.extend(varfiles.index)
    return unique_everseen(contigs)


def informative_contigs(varfile):
    """Scan an indexed variant file to determine which contigs have alignments.

    Args:
        varfile: Input variant file object

    Returns:
        All contigs that have data

    """
    if varfile.index is None:
        raise ValueError('Variant file requires index')
    return (contig for contig in varfile.index if not is_empty_iter(varfile.fetch(contig)))


def region_filter_include(records, include):
    """Remove records that do not overlap those provided."""
    for _, _, (rec, inc) in union([records, include]):
        if inc:
            yield from rec


def region_filter_exclude(records, exclude):
    """Remove records that overlap those provided."""
    for _, _, (rec, exc) in union([records, exclude]):
        if not exc:
            yield from rec


def filter_gq(records, name, min_gq):
    """Filter records based on a minimum genotype quality value."""
    for record in records:
        if 'GQ' in record.format:
            gq = record.samples[name]['GQ']
            if gq is not None and gq < min_gq:
                continue
        yield record


def filter_records(records, name, args):
    """Filter records based on multiple criteria."""
    if args.min_gq is not None:
        records = filter_gq(records, name, args.min_gq)

    if args.include_filter:
        include = {f.strip() for fs in args.include_filter for f in fs.split(',')}
        records = (record for record in records if not include.isdisjoint(record.filter))

    if args.exclude_filter:
        exclude = {f.strip() for fs in args.exclude_filter for f in fs.split(',')}
        records = (record for record in records if exclude.isdisjoint(record.filter))

    return records


def records_by_chromosome(refs, varfiles, names, args, get_all=False):
    """Group variant records by chromosome."""
    contigs_all   = unique_everseen(chain.from_iterable(all_contigs(var) for var in varfiles))
    contigs_seen  = set(chain.from_iterable(informative_contigs(var) for var in varfiles))
    contigs_fetch = [contig for contig in contigs_all if contig in contigs_seen]

    if args.include_regions is not None:
        include = load_bedmap(args.include_regions)
        contigs_fetch = [contig for contig in contigs_fetch if contig in include]

    if args.exclude_regions is not None:
        exclude = load_bedmap(args.exclude_regions)

    for contig in contigs_fetch:
        try:
            if args.lazy_ref:
                ref = LazyFastaContig(refs, contig)
            else:
                ref = refs.fetch(contig)
        except KeyError:
            continue

        records = [var.fetch(contig) if contig in var.index else [] for var in varfiles]

        if get_all:
            all_records = records = [list(l) for l in records]

        records = [filter_records(r, name, args) for r, name in zip(records, names)]

        if args.include_regions is not None:
            records = [region_filter_include(r, include[contig]) for r in records]
        if args.exclude_regions is not None:
            records = [region_filter_exclude(r, exclude[contig]) for r in records]

        loci = [records_to_loci(ref, r, name, args.reference_padding) for name, r in zip(names, records)]
        loci = [sort_almost_sorted(l, key=NormalizedLocus.natural_order_key) for l in loci]

        if get_all:
            yield contig, ref, loci, all_records
        else:
            yield contig, ref, loci


def get_superlocus_bounds(superloci):
    """Get the most 5' and 3' boundaries of a superlocus."""
    start = min(locus.min_start for super in superloci for locus in super)
    stop  = max(locus.max_stop  for super in superloci for locus in super)
    return start, stop


def locus_equal_trivial(locus1, locus2):
    """Compare two loci for trivial equality (i.e. no graph-based matching)."""
    left1, left2 = locus1.left, locus2.left

    if left1.start != left2.start:
        return False

    if left1.stop != left2.stop:
        return False

    alleles1, alleles2 = left1.alleles, left2.alleles
    g1 = tuple(alleles1[i] for i in locus1.allele_indices)
    g2 = tuple(alleles2[i] for i in locus2.allele_indices)

    if not locus1.phased or not locus2.phased:
        g1, g2 = tuple(sorted(g1)), tuple(sorted(g2))

    return g1 == g2


def superlocus_equal_trivial(super1, super2):
    """Compare two superloci for trivial equality (i.e. no graph-based matching)."""
    if len(super1) != len(super2):
        return False

    for locus1, locus2 in zip(super1, super2):
        if not locus_equal_trivial(locus1, locus2):
            return False

    return True


def superlocus_equal(ref, start, stop, super1, super2, debug=False):
    """Compare two superloci."""
    if superlocus_equal_trivial(super1, super2):
        return True, 'T'

    # Bounds come from normalized extremes
    start, stop = get_superlocus_bounds([super1, super2])

    # Create genotype sets for each superlocus
    try:
        graph1, constraints1 = generate_graph(ref, start, stop, super1, debug)
        graph2, constraints2 = generate_graph(ref, start, stop, super2, debug)

        paths1 = generate_paths(graph1, debug=debug)
        paths2 = generate_paths(graph2, feasible_paths=paths1, debug=debug)

        paths1, paths2 = intersect_paths(paths1, paths2)

    except OverlapError:
        status = None, 'N'

    else:
        genos1 = set(generate_genotypes(paths1, constraints1, debug))
        genos2 = set(generate_genotypes(paths2, constraints2, debug))

        # Test whether genotype sets intersect
        if genos1.isdisjoint(genos2):
            status = False, 'H'
        else:
            status = True, 'H'

    return status


def find_allele_exact_match(ref, allele, superlocus):
    """Search for allele using a fast exact match criteria."""
    for locus in superlocus:
        a  = allele.left
        ll = locus.left

        if a.start == ll.start and a.stop == ll.stop and 'PASS' in locus.record.filter:
            return sum(
                locus.allele_indices.count(locus.alleles.index(alt))
                for alt in set(ll.alts).intersect(a.alts)
            )

    return 0


def path_allele_counts(paths):
    """Count number ploidy of each allele at each node in paths."""
    for ps in zip(*paths):
        yield Counter(p.index for p in ps if p.locus)


def path_to_ads(path, counts):
    """Convert a path through a variant graph into a sequence of allele depths."""
    for p, c in zip(path, counts):
        # Skip nodes with no VCF record
        if not p.locus:
            continue
        record = p.locus.record
        sample = record.samples[0]
        dp     = sample['AD'][p.index] if 'AD' in sample and p.index is not None else sample.get('MIN_DP', 0)
        yield dp // c[p.index]


def path_to_ref_ads(path):
    """Convert a path through a variant graph into a sequence of reference allele depths."""
    for p in path:
        # Skip nodes with no VCF record
        if not p.locus:
            continue
        record = p.locus.record
        sample = record.samples[0]
        yield sample['AD'][0] if 'AD' in sample else sample.get('MIN_DP', 0)


def build_match_result(geno, matches, super_ref):
    """Build match results."""
    seqs, paths   = zip(*geno)
    allele_counts = list(path_allele_counts(paths))
    allele_depths = [list(path_to_ads(path, allele_counts)) for path in paths]

    found, ref, other = [], [], []
    for m, seq, ad in zip(matches, seqs, allele_depths):
        if m:
            found.append(ad)
        elif fancy_match(super_ref, seq):
            ref.append(ad)
        else:
            other.append(ad)

    ref_ploidy    = len(ref)
    allele_ploidy = len(found)
    other_ploidy  = len(other)

    # If a reference allele was not called, then collect all reference allele depths
    # from an arbitrary path, as AD always contains reference counts.
    if not ref_ploidy and paths:
        ref = [list(path_to_ref_ads(paths[0]))]

    allele_ad     = empty_min(chain.from_iterable(found), default=None)
    ref_ad        = empty_min(chain.from_iterable(ref),   default=None)
    other_ad      = empty_min(chain.from_iterable(other), default=None)

    return AlleleMatch(allele_ploidy, allele_ad, ref_ploidy, ref_ad, other_ploidy, other_ad)


nothing = object()


def empty_min(items, default=nothing):
    items = list(items)
    if not items and default is not nothing:
        return default
    return min(items)


def int_mean(items, default=nothing):
    """Take the rounded mean of a sequence of items."""
    n = len(items)
    if not n and default is not nothing:
        return default
    return round(sum(items) / n)


def build_match_strings(ref, start, stop, allele, mode='sensitive', debug=False):
    """Build allele matching strings."""
    alts = allele.alts

    if debug:
        print('  Allele: start={}, stop={}, size={}, ref={}, seq={}'.format(
            allele.start,
            allele.stop,
            allele.stop - allele.start,
            allele.ref,
            ','.join(alts),
        ), file=sys.stderr)

    # Require reference matches within the wobble zone + padding built into each normalized allele
    if mode == 'specific':
        super_alleles = [normalize_seq(ref[start:allele.start] + alt + ref[allele.stop:stop]) for alt in alts]
        super_ref     = normalize_seq(ref[start:stop])
    elif mode == 'sensitive':
        super_alleles = [
            normalize_seq(
                '*' * (allele.min_start - start)
                + ref[allele.min_start:allele.start]
                + alt
                + ref[allele.stop:allele.max_stop]
                + '*' * (stop - allele.max_stop)
            ) for alt in alts
        ]

        super_ref = normalize_seq(
            '*' * (allele.min_start - start)
            + ref[allele.min_start:allele.max_stop]
            + '*' * (stop - allele.max_stop)
        )
    else:
        raise ValueError(f'invalid match mode specified: {mode}')

    if debug:
        print('                MODE:', mode,          file=sys.stderr)
        print('       SUPER ALLELES:', super_alleles, file=sys.stderr)
        print('        SUPER REF:   ', super_ref,     file=sys.stderr)

    assert all(len(a) == stop - start - len(allele.ref) + len(alt) for a, alt in zip(super_alleles, alts))

    return super_ref, super_alleles


def compare_alleles(alleles, seq):
    """Compare sequence with each potential allele using fancy matching.

    Args:
        alleles: candidate alleles
        seq: sequence with which to compare

    Returns:
        True: if any allele fancy-matches seq
        False: if no alleles match and none are uncertain
        None: if no alleles match and one or more is uncertain

    """
    notfound = False
    for allele in alleles:
        match = fancy_match(allele, seq)
        if match:
            return True
        elif match is None:
            notfound = None
    return notfound


def find_allele_matches(ref, start, stop, allele, genos, ploidy, mode, debug=False):
    """Analyze graph paths to find allele matches."""
    # superlocus contains impossible genotypes and no paths are valid
    super_ref, super_alleles = build_match_strings(ref, start, stop, allele, mode, debug)

    if not genos:
        return None

    all_matches = (
        (
            [compare_alleles(super_alleles, seq) for (seq, _) in geno],
            geno,
        ) for geno in genos
    )

    zygosity, nocalls, _, geno, matches = max(
        (
            sum(m or 0 for m in matches),
            -matches.count(None),
            i,
            geno,
            matches,
        ) for i, (matches, geno) in enumerate(all_matches)
    )

    if not zygosity and nocalls:
        return None

    result = build_match_result(geno, matches, super_ref)

    if debug:
        for super_allele in super_alleles:
            print('   ALLELE:{} {}'.format(len(super_allele), super_allele), file=sys.stderr)
        for i, (g, m) in enumerate(zip(genos, matches)):
            print('   GENO{:02d}:{} {}'.format(i, tuple(map(len, g)),  g), file=sys.stderr)
            print(f'  MATCH{i:02d}: {m}', file=sys.stderr)
        print(file=sys.stderr)
        print(
            f'ALLELE: id={allele.record.id}, allele_ploidy={result.allele_ploidy}, '
            f'ref_ploidy={result.ref_ploidy}, other_ploidy={result.other_ploidy}, ploidy={result.ploidy}',
            file=sys.stderr
        )
        print(f'  ZYGOSITY: {zygosity}', file=sys.stderr)

    return result


def find_allele(ref, allele, superlocus, mode='sensitive', debug=False):
    """Check for the presence of an allele within a superlocus."""
    start, stop = get_superlocus_bounds([[allele], superlocus])

    # Create genotype sets for each superlocus
    try:
        graph, constraints = generate_graph(ref, start, stop, superlocus, debug)

        if debug:
            graph = list(graph)
            for i, (astart, astop, alleles) in enumerate(graph):
                print(f'  GRAPH{i:02d}: start={astart}, stop={astop}, alleles={alleles}', file=sys.stderr)
            print(file=sys.stderr)

        paths = list(generate_paths(graph, debug=debug))

        if debug:
            for i, p in enumerate(paths):
                print(f'  PATH{i:02d}: {p}', file=sys.stderr)
            print(file=sys.stderr)

    except OverlapError:
        return None

    # Generate the set of diploid genotypes (actually haplotypes)
    ploidy = max(len(locus.allele_indices) for locus in superlocus) if superlocus else 2
    genos  = list(generate_genotypes_with_paths(paths, constraints, ploidy))

    return find_allele_matches(ref, start, stop, allele, genos, ploidy, mode, debug)
