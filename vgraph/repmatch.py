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
from vgraph.match       import variants_by_chromosome, get_superlocus_bounds, superlocus_equal


def match_replicates(args):
    # Load FASTA reference
    refs = Fastafile(expanduser(args.reference))

    # Open input variant files
    in_vars = [VariantFile(var) for var in [args.vcf1, args.vcf2]]

    out_vars = [None, None]

    if args.out1:
        in_vars[0].header.formats.add('BD', '1', 'String', 'Match decision for call (match: =, mismatch: X, error: N)')
        in_vars[0].header.formats.add('BK', '1', 'String', 'Sub-type for match decision (trivial: T, haplotype: H, error: N)')
        out_vars[0] = VariantFile(args.out1, 'w', header=in_vars[0].header)

    if args.out2:
        in_vars[1].header.formats.add('BD', '1', 'String', 'Match decision for call (match: =, mismatch: X, error: N)')
        in_vars[1].header.formats.add('BK', '1', 'String', 'Sub-type for match decision (trivial: T, haplotype: H, error: N)')
        out_vars[1] = VariantFile(args.out2, 'w', header=in_vars[1].header)

    match_status_map = {True : '=', False : 'X', None : '.'}

    # Create parallel locus iterator by chromosome
    for chrom, ref, loci in variants_by_chromosome(refs, in_vars, [args.name1, args.name2], args):
        # Create superloci by taking the union of overlapping loci across all of the locus streams
        loci = [sort_almost_sorted(l, key=NormalizedLocus.extreme_order_key) for l in loci]
        superloci = union(loci, min_distance=args.reference_padding,
                                interval_func=attrgetter('min_start', 'max_stop'))

        # Proceed by superlocus
        for _, _, (super1, super2) in superloci:
            super1.sort(key=NormalizedLocus.natural_order_key)
            super2.sort(key=NormalizedLocus.natural_order_key)

            super_start, super_stop = get_superlocus_bounds([super1, super2])

            print('-'*80)
            print('{}:[{:d}-{:d}):'.format(chrom, super_start, super_stop))
            print()

            for i, (name, superlocus) in enumerate([(args.name1, super1), (args.name2, super2)], 1):
                for locus in superlocus:
                    lstart = locus.start
                    lstop = locus.stop
                    lref = locus.alleles[0] or '-'
                    indices = locus.allele_indices
                    sep = '|' if locus.phased else '/'
                    geno = sep.join(locus.alleles[a] or '-' if a is not None else '.' for a in indices)
                    print('  NORM{:d}: [{:5d}-{:5d}) ref={} geno={}'.format(i, lstart, lstop, lref, geno))
            print()

            match, match_type = superlocus_equal(ref, super_start, super_stop, super1, super2, debug=args.debug)
            match_status = match_status_map[match]

            print('    MATCH={} TYPE={}'.format(match_status, match_type))
            print()

            # The hard work is done.  The rest is just output and formatting...

            if out_vars[0]:
                for locus in sorted(super1, key=NormalizedLocus.record_order_key):
                    locus.record.samples[args.name1]['BD'] = match_status
                    locus.record.samples[args.name1]['BK'] = match_type
                    out_vars[0].write(locus.record)

            if out_vars[1]:
                for locus in sorted(super2, key=NormalizedLocus.record_order_key):
                    locus.record.samples[args.name2]['BD'] = match_status
                    locus.record.samples[args.name2]['BK'] = match_type
                    out_vars[1].write(locus.record)

            for i, superlocus in enumerate([super1, super2], 1):
                for locus in superlocus:
                    print('  VCF{:d}: {}'.format(i, locus.record), end='')
            print()

    for out_var in out_vars:
        if out_var is not None:
            out_var.close()
