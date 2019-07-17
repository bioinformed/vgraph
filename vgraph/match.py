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
from itertools          import chain

from vgraph.bed         import load_bedmap
from vgraph.norm        import NormalizedLocus, fancy_match, normalize_seq, ReferenceMismatch
from vgraph.intervals   import union
from vgraph.iterstuff   import sort_almost_sorted, is_empty_iter, unique_everseen
from vgraph.lazy_fasta  import LazyFastaContig
from vgraph.linearmatch import generate_graph, generate_paths, generate_genotypes, intersect_paths, OverlapError


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
            if len(record.alleles) > 1:
                yield NormalizedLocus(recnum, record, ref, name, variant_padding)
        except ReferenceMismatch:
            print('Reference mismatch: {}:{}-{}'.format(record.contig, record.start, record.stop))


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


def find_allele(ref, allele, superlocus, debug=False):
    """Check for the presence of an allele within a superlocus."""
    # FASTPATH: Avoid constructing the graph match if the allele and the superlocus
    #           match perfectly.
    if (len(superlocus) == 1
            and allele.start == superlocus[0].start
            and allele.stop == superlocus[0].stop
            and allele.alleles[1] in superlocus[0].alleles[1:]
            and 'PASS' in superlocus[0].record.filter):

        i = superlocus[0].alleles.index(allele.alleles[1])
        z = superlocus[0].allele_indices.count(i)

        return z

    # Bounds come from normalized extremes
    start, stop = get_superlocus_bounds([[allele], superlocus])

    if debug:
        print('  Allele: start={}, stop={}, size={}, seq={}'.format(
            allele.start,
            allele.stop,
            allele.stop - allele.start,
            allele.alleles[1]
        ), file=sys.stderr)

    # Require reference matches within the wobble zone + padding built into each normalized allele
    super_allele = ('*' * (allele.min_start - start)
                 +  ref[allele.min_start:allele.start]
                 +  allele.alleles[1]
                 +  ref[allele.stop:allele.max_stop]
                 +  '*' * (stop - allele.max_stop))

    super_allele = normalize_seq(super_allele)

    assert len(super_allele) == stop - start - len(allele.alleles[0]) + len(allele.alleles[1])

    # Create genotype sets for each superlocus
    try:
        graph, constraints = generate_graph(ref, start, stop, superlocus, debug)
        graph = list(graph)

        if debug:
            for i, (start, stop, alleles) in enumerate(graph):
                print('  GRAPH{:02d}: start={}, stop={}, alleles={}'.format(i, start, stop, alleles), file=sys.stderr)
            print(file=sys.stderr)

        paths = generate_paths(graph, debug=debug)

        if debug:
            paths = list(paths)
            for i, p in enumerate(paths):
                print('  PATH{:02d}: {}'.format(i, p), file=sys.stderr)
            print(file=sys.stderr)

    except OverlapError:
        return None

    # Generate the set of diploid genotypes (actually haplotypes)
    genos = set(generate_genotypes(paths, constraints, debug))

    # Apply matcher to each pair of allele
    matches = [(fancy_match(super_allele, a1), fancy_match(super_allele, a2))
            for a1, a2 in genos]

    # Find the highest zygosity
    z = max(((a1 or 0) + (a2 or 0)) for a1, a2 in matches)

    # If no match, check for the presense of any nocalls
    if not z and any(None in m for m in matches):
        z = None

    if debug:
        print('   ALLELE:{} {}'.format(len(super_allele), super_allele), file=sys.stderr)
        for i, (g, m) in enumerate(zip(genos, matches)):
            print('   GENO{:02d}:{} {}'.format(i, tuple(map(len, g)),  g), file=sys.stderr)
            print('  MATCH{:02d}: {}'.format(i, m), file=sys.stderr)
        print(file=sys.stderr)
        print('  ZYGOSITY: {}'.format(z), file=sys.stderr)

    return z
