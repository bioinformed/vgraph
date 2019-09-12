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

"""Match a genome against another presumably identical genome (i.e. replicates)."""

from os.path            import expanduser
from operator           import attrgetter

from pysam              import VariantFile, Fastafile

from vgraph.norm        import NormalizedLocus
from vgraph.intervals   import union
from vgraph.iterstuff   import sort_almost_sorted
from vgraph.match       import records_by_chromosome, get_superlocus_bounds, superlocus_equal


def make_outputs(in_vars, out1, out2):
    """Make output files."""
    out_vars = [None, None]

    if out1:
        in_vars[0].header.formats.add('BD', '1', 'String', 'Match decision for call (match: =, mismatch: X, error: N)')
        in_vars[0].header.formats.add('BK', '1', 'String', 'Sub-type for match decision (trivial: T, haplotype: H, error: N)')
        out_vars[0] = VariantFile(out1, 'w', header=in_vars[0].header)

    if out2:
        in_vars[1].header.formats.add('BD', '1', 'String', 'Match decision for call (match: =, mismatch: X, error: N)')
        in_vars[1].header.formats.add('BK', '1', 'String', 'Sub-type for match decision (trivial: T, haplotype: H, error: N)')
        out_vars[1] = VariantFile(out2, 'w', header=in_vars[1].header)

    return out_vars


def write_match(out, superlocus, name, match_status, match_type):
    """Write match to output file."""
    if not out:
        return

    for locus in sorted(superlocus, key=NormalizedLocus.record_order_key):
        sample       = locus.record.samples[name]
        sample['BD'] = match_status
        sample['BK'] = match_type
        out.write(locus.record)


def match_replicates(args):
    """Match a genome against another presumably identical genome (i.e. replicates)."""
    refs     = Fastafile(expanduser(args.reference))
    in_vars  = [VariantFile(var) for var in [args.vcf1, args.vcf2]]
    out_vars = make_outputs(in_vars, args.out1, args.out2)

    match_status_map = {True: '=', False: 'X', None: '.'}

    # Create parallel locus iterator by chromosome
    for chrom, ref, loci in records_by_chromosome(refs, in_vars, [args.name1, args.name2], args):
        # Create superloci by taking the union of overlapping loci across all of the locus streams
        loci = [sort_almost_sorted(l, key=NormalizedLocus.extreme_order_key) for l in loci]
        superloci = union(loci, interval_func=attrgetter('min_start', 'max_stop'))

        # Proceed by superlocus
        for _, _, (super1, super2) in superloci:
            super1.sort(key=NormalizedLocus.natural_order_key)
            super2.sort(key=NormalizedLocus.natural_order_key)

            super_start, super_stop = get_superlocus_bounds([super1, super2])

            print('-' * 80)
            print('{}:[{:d}-{:d}):'.format(chrom, super_start, super_stop))
            print()

            for i, superlocus in enumerate([super1, super2], 1):
                for locus in superlocus:
                    lstart = locus.start
                    lstop = locus.stop
                    lref = locus.ref or '-'
                    indices = locus.allele_indices
                    sep = '|' if locus.phased else '/'
                    geno = sep.join(locus.alleles[a] or '-' if a is not None else '.' for a in indices)
                    print('  NORM{:d}: [{:5d}-{:5d}) ref={} geno={}'.format(i, lstart, lstop, lref, geno))
            print()

            match, match_type = superlocus_equal(ref, super_start, super_stop, super1, super2, debug=args.debug)
            match_status = match_status_map[match]

            print('    MATCH={} TYPE={}'.format(match_status, match_type))
            print()

            write_match(out_vars[0], super1, args.name1, match_status, match_type)
            write_match(out_vars[1], super2, args.name2, match_status, match_type)

            for i, superlocus in enumerate([super1, super2], 1):
                for locus in superlocus:
                    print('  VCF{:d}: {}'.format(i, locus.record), end='')
            print()

    for out_var in out_vars:
        if out_var is not None:
            out_var.close()
