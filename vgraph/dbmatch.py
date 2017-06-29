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

import sys
from os.path            import expanduser
from operator           import attrgetter

from pysam              import VariantFile, Fastafile

from vgraph.norm        import NormalizedLocus
from vgraph.intervals   import union
from vgraph.iterstuff   import sort_almost_sorted
from vgraph.match       import records_by_chromosome, get_superlocus_bounds, find_allele


def annotate_info(locus, allele, info_meta, suffix, times):
    for name in info_meta:
        if name in allele.record.info:
            sname = name + suffix
            orig_value = locus.record.info.get(sname, ())
            new_value = allele.record.info[name]
            if not isinstance(new_value, tuple):
                new_value = (new_value,)
            locus.record.info[sname] = orig_value + new_value*times


def annotate_format(locus, allele, format_meta, suffix, times):
    sample = locus.record.samples[0]
    for name in format_meta:
        if name in allele.record.format:
            sname = name + suffix
            orig_value = sample.get(sname, ())
            new_value = allele.record.samples[0][name]
            if not isinstance(new_value, tuple):
                new_value = (new_value,)
            sample[sname] = orig_value + new_value*times


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
            sample.header.formats.add(meta.name + '_FOUND',    number='.', type=meta.type,
                                      description='Allele(s) found: ' + meta.description)
            sample.header.formats.add(meta.name + '_NOTFOUND', number='.', type=meta.type,
                                      description='Allele(s) not found: ' + meta.description)
            sample.header.formats.add(meta.name + '_NOCALL',   number='.', type=meta.type,
                                      description='Allele(s) with uncertain presense: ' + meta.description)

    info_meta = []
    for info, meta in db.header.info.items():
        if info not in sample.header.info:
            info_meta.append(meta.name)
            sample.header.info.add(meta.name + '_FOUND',    number='.', type=meta.type,
                                   description='Allele(s) found: ' + meta.description)
            sample.header.info.add(meta.name + '_NOTFOUND', number='.', type=meta.type,
                                   description='Allele(s) not found: ' + meta.description)
            sample.header.info.add(meta.name + '_NOCALL',   number='.', type=meta.type,
                                   description='Allele(s) with uncertain presense: ' + meta.description)

    with VariantFile(args.output, 'w', header=sample.header) as out:
        # Create parallel locus iterator by chromosome
        for chrom, ref, loci in records_by_chromosome(refs, [sample, db], [args.name, None], args):
            # Create superloci by taking the union of overlapping loci across all of the locus streams
            loci = [sort_almost_sorted(l, key=NormalizedLocus.extreme_order_key) for l in loci]
            superloci = union(loci, interval_func=attrgetter('min_start', 'max_stop'))

            # Proceed by superlocus
            for _, _, (superlocus, alleles) in superloci:
                alleles.sort(key=NormalizedLocus.natural_order_key)
                superlocus.sort(key=NormalizedLocus.natural_order_key)

                for allele in alleles:
                    super_allele = [locus for locus in superlocus if locus.extremes_intersect(allele)]

                    # Remove all reference calls from the superlocus.
                    # This is primarily done to remove long leading and trailing reference regions.
                    # Interstitial reference regions will be added back, based on how gaps are handled.
                    super_non_ref = [locus for locus in super_allele if not locus.is_ref()]

                    if args.debug:
                        super_start, super_stop = get_superlocus_bounds([[allele], super_non_ref])
                        print('-'*80, file=sys.stderr)
                        print('{}:[{:d}-{:d}):'.format(chrom, super_start, super_stop), file=sys.stderr)
                        print(file=sys.stderr)

                        print('  ALLELE: {} {}:[{}-{}) ref={} alt={}'.format(allele.record.id, allele.contig,
                                                                             allele.start, allele.stop,
                                                                             allele.alleles[0] or '-', allele.alleles[1] or '-'), file=sys.stderr)
                        print(file=sys.stderr)

                        for i, locus in enumerate(super_non_ref, 1):
                            lref = locus.alleles[0] or '-'
                            indices = locus.allele_indices
                            if indices.count(None) == len(indices):
                                geno = 'nocall'
                            elif indices.count(0) == len(indices):
                                geno = 'refcall'
                            else:
                                sep = '|' if locus.phased else '/'
                                geno = sep.join(locus.alleles[a] or '-' if a is not None else '.' for a in indices)
                            print('  VAR{:d}: {}[{:5d}-{:5d}) ref={} geno={}'.format(i, locus.contig, locus.start, locus.stop, lref, geno), file=sys.stderr)

                    # Search superlocus for allele
                    match_zygosity = find_allele(ref, allele, super_non_ref, debug=args.debug)

                    if args.debug:
                        print(file=sys.stderr)
                        print('    MATCH={}'.format(match_zygosity), file=sys.stderr)
                        print(file=sys.stderr)

                    # Annotate results of search
                    if match_zygosity is None:
                        suffix = '_NOCALL'
                    elif match_zygosity == 0:
                        suffix = '_NOTFOUND'
                    else:
                        suffix = '_FOUND'

                    # Number of times to repeat the copied metadata
                    times = match_zygosity if suffix == '_FOUND' else 1

                    for locus in super_allele:
                        annotate_info(locus, allele, info_meta, suffix, times)
                        annotate_format(locus, allele, format_meta, suffix, times)

                for locus in sorted(superlocus, key=NormalizedLocus.record_order_key):
                    out.write(locus.record)
