#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Copyright 2016 Kevin B Jacobs
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

from __future__ import division, print_function


from collections        import defaultdict
from operator           import attrgetter

from vgraph.bed         import BedFile
from vgraph.norm        import NormalizedLocus
from vgraph.intervals   import union
from vgraph.iterstuff   import sort_almost_sorted, is_empty_iter
from vgraph.linearmatch import generate_graph, generate_paths, generate_genotypes, intersect_paths


def valid_alleles(alleles):
    return not any('<' in a or '[' in a or ']' in a for a in alleles)


def is_alt_genotype(record, name):
    sample = record.samples[name]
    indices = sample.allele_indices
    return bool(indices and -1 not in indices and indices.count(0) != len(indices))


def records_to_loci(ref, records, name):
    # Prevents left shuffling the locus through the previous locus
    last_left_stop = 0
    for recnum, record in enumerate(records):
        if valid_alleles(record.alleles) and is_alt_genotype(record, name):
            locus = NormalizedLocus(recnum, record, ref, name, last_left_stop)
            last_left_stop = locus.left.stop
            yield locus


def informative_chromosomes(vars):
    if not vars.index:
        raise ValueError('Variant file requires index')
    return (chrom for chrom in vars.index if not is_empty_iter(vars.fetch(chrom)))


def make_bedmap(bedfile):
    bedmap = defaultdict(list)
    for record in bedfile:
        bedmap[record.contig].append(record)
    return bedmap


def region_filter_include(records, include):
    for start, stop, (rec, inc) in union([records, include]):
        if inc:
            for record in rec:
                yield record


def region_filter_exclude(records, exclude):
    for start, stop, (rec, exc) in union([records, exclude]):
        if not exc:
            for record in rec:
                yield record


def filter_gq(records, name, min_gq):
    for record in records:
        if 'GQ' in record.format:
            gq = record.samples[name]['GQ']
            if gq is not None and gq < min_gq:
                continue
        yield record


def filter_records(records, name, args):
    if args.min_gq is not None:
        records = filter_gq(records, name, args.min_gq)

    if args.include_filter:
        include = set(f.strip() for fs in args.include_filter for f in fs.split(','))
        records = (record for record in records if not include.isdisjoint(record.filter))

    if args.exclude_filter:
        exclude = set(f.strip() for fs in args.exclude_filter for f in fs.split(','))
        records = (record for record in records if exclude.isdisjoint(record.filter))

    return records


def variants_by_chromosome(refs, vars, names, args):
    for var in vars:
        if not var.index:
            raise ValueError('Input variant file `{}` is missing an index'.format(var.filename))

    chroms = [set(informative_chromosomes(var)) for var in vars]
    chroms = set.union(*chroms)

    if args.include_regions is not None:
        include = make_bedmap(BedFile(args.include_regions))
        chroms &= set(include)

    if args.exclude_regions is not None:
        exclude = make_bedmap(BedFile(args.exclude_regions))

    if args.include_file_regions:
        assert len(args.include_file_regions) == len(vars)
        include_files = [make_bedmap(BedFile(fn)) for fn in args.include_file_regions]

    if args.exclude_file_regions:
        assert len(args.exclude_file_regions) == len(vars)
        exclude_files = [make_bedmap(BedFile(fn)) for fn in args.exclude_file_regions]

    for chrom in chroms:
        ref  = refs.fetch(chrom).upper()
        loci = [filter_records(var.fetch(chrom), name, args) for var,name in zip(vars, names)]

        if args.include_file_regions:
            loci = [region_filter_include(l, inc[chrom]) for l,inc in zip(loci, include_files)]
        if args.exclude_file_regions:
            loci = [region_filter_exclude(l, exl[chrom]) for l,exl in zip(loci, exclude_files)]

        loci = [records_to_loci(ref, l, name) for name, l in zip(names, loci)]
        loci = [sort_almost_sorted(l, key=NormalizedLocus.extreme_order_key) for l in loci]

        if args.include_regions is not None:
            loci = [region_filter_include(l, include[chrom]) for l in loci]
        if args.exclude_regions is not None:
            loci = [region_filter_exclude(l, exclude[chrom]) for l in loci]

        yield chrom, ref, loci


def get_superlocus_bounds(superloci):
    start = min(super[ 0].left.start for super in superloci if super)
    stop  = max(super[-1].left.stop  for super in superloci if super)
    return start, stop


def locus_equal_trivial(locus1, locus2):
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
    if len(super1) != len(super2):
        return False

    for locus1, locus2 in zip(super1, super2):
        if not locus_equal_trivial(locus1, locus2):
            return False

    return True


def superlocus_equal(ref, start, stop, super1, super2, debug=False):
    if superlocus_equal_trivial(super1, super2):
       return True, 'T'

    if not super1 or not super2:
        return False, 'T'

    # Bounds come from left normalized extremes
    start, stop = get_superlocus_bounds([super1, super2])

    # Create genotype sets for each superlocus
    try:
        graph1, constraints1 = generate_graph(ref, start, stop, super1, debug)
        graph2, constraints2 = generate_graph(ref, start, stop, super2, debug)

        paths1 = generate_paths(graph1, debug)
        paths2 = generate_paths(graph2, debug)

        paths1, paths2 = intersect_paths(paths1, paths2)

        genos1 = set(generate_genotypes(paths1, constraints1, debug))
        genos2 = set(generate_genotypes(paths2, constraints2, debug))

        # Test whether genotype sets intersect
        if genos1.isdisjoint(genos2):
            status = False, 'H'
        else:
            status = True, 'H'
    except ValueError:
        status = None, 'N'

    return status
