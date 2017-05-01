#!/usr/bin/env python
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

from os.path            import expanduser
from operator           import attrgetter

from pysam              import VariantFile, Fastafile

from vgraph.norm        import NormalizedLocus
from vgraph.intervals   import union
from vgraph.iterstuff   import sort_almost_sorted
from vgraph.match       import variants_by_chromosome, get_superlocus_bounds, superlocus_equal, find_allele


def match_database(args):
    # Load FASTA reference
    refs = Fastafile(expanduser(args.reference))

    # Open input variant files
    db = VariantFile(args.database)
    sample = VariantFile(args.sample)

    format_meta = []
    for fmt, meta in db.header.formats.items():
        if fmt not in sample.header.formats:
            format_meta.append(meta.name)
            sample.header.formats.add(meta.name + '_FOUND',    meta.number, meta.type,
                                      'Allele(s) found: ' + meta.description)
            sample.header.formats.add(meta.name + '_NOTFOUND', meta.number, meta.type,
                                      'Allele(s) not found: ' + meta.description)
            sample.header.formats.add(meta.name + '_NOCALL',   meta.number, meta.type,
                                      'Allele(s) with uncertain presense: ' + meta.description)

    info_meta = []
    for info, meta in db.header.info.items():
        if info not in sample.header.info:
            info_meta.append(meta.name)
            sample.header.info.add(meta.name + '_FOUND',    meta.number, meta.type,
                                   'Allele(s) found: ' + meta.description)
            sample.header.info.add(meta.name + '_NOTFOUND', meta.number, meta.type,
                                   'Allele(s) not found: ' + meta.description)
            sample.header.info.add(meta.name + '_NOCALL',   meta.number, meta.type,
                                   'Allele(s) with uncertain presense: ' + meta.description)

    with VariantFile(args.output, 'w', header=sample.header) as out:
        # Create parallel locus iterator by chromosome
        for chrom, ref, loci in variants_by_chromosome(refs, [sample, db], [args.name, None], args):
            # Create superloci by taking the union of overlapping loci across all of the locus streams
            loci = [sort_almost_sorted(l, key=NormalizedLocus.extreme_order_key) for l in loci]
            superloci = union(loci, interval_func=attrgetter('min_start', 'max_stop'))

            # Proceed by superlocus
            for _, _, (super2, alleles) in superloci:
                alleles.sort(key=NormalizedLocus.natural_order_key)
                super2.sort(key=NormalizedLocus.natural_order_key)

                for allele in alleles:
                    super_all = [locus for locus in super2 if locus.extremes_intersect(allele)]

                    super_trimmed = super_all.copy()
                    while super_trimmed and super_trimmed[-1].is_ref():
                        super_trimmed.pop()
                    while super_trimmed and super_trimmed[0].is_ref():
                        super_trimmed.pop(0)

                    super_start, super_stop = get_superlocus_bounds([[allele], super_trimmed])

                    print('-'*80)
                    print('{}:[{:d}-{:d}):'.format(chrom, super_start, super_stop))
                    print()

                    print('  ALLELE: {} {}:[{}-{}) ref={} alt={}'.format(allele.record.id, allele.contig,
                                                                         allele.start, allele.stop,
                                                                         allele.alleles[0] or '-', allele.alleles[1] or '-'))
                    print()

                    for i, locus in enumerate(super_trimmed, 1):
                        lref = locus.alleles[0] or '-'
                        indices = locus.allele_indices
                        if indices.count(None) == len(indices):
                            geno = 'nocall'
                        elif indices.count(0) == len(indices):
                            geno = 'refcall'
                        else:
                            sep = '|' if locus.phased else '/'
                            geno = sep.join(locus.alleles[a] or '-' if a is not None else '.' for a in indices)
                        print('  VAR{:d}: {}[{:5d}-{:5d}) ref={} geno={}'.format(i, locus.contig, locus.start, locus.stop, lref, geno))
                    print()

                    match_zygosity = find_allele(ref, allele, super_trimmed, debug=args.debug)

                    print('    MATCH={}'.format(match_zygosity))
                    print()

                    if match_zygosity is None:
                        suffix = '_NOCALL'
                    elif match_zygosity == 0:
                        suffix = '_NOTFOUND'
                    else:
                        suffix = '_FOUND'

                    for locus in super_all:
                        if not locus.intersects(allele):
                            continue

                        if suffix == '_FOUND':
                            locus.record.id = allele.record.id

                        for name in info_meta:
                            if name in allele.record.info:
                                sname = name + suffix
                                value = locus.record.info[sname] if sname in locus.record.info else ()
                                locus.record.info[sname] = value + (allele.record.info[name],)

                        sample = locus.record.samples[0]
                        for name in format_meta:
                            if name in allele.record.format:
                                sname = name + suffix
                                value = sample[sname] if sample.get(sname) else ()
                                sample[sname] = value + allele.record.samples[0][name]

                for locus in sorted(super2, key=NormalizedLocus.record_order_key):
                    out.write(locus.record)
    