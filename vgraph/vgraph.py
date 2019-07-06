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

"""Match variants in a representation agnostic manner using graph methods."""

import sys

from os.path            import expanduser
from argparse           import ArgumentParser

from pysam              import Fastafile, VariantFile

from vgraph.norm        import NormalizedLocus
from vgraph.iterstuff   import sort_almost_sorted
from vgraph.match       import records_by_chromosome
from vgraph.repmatch    import match_replicates
from vgraph.dbmatch     import match_database, match_database2


def normalize(args):
    """Normalize variants."""
    refs = Fastafile(expanduser(args.reference))
    variants = VariantFile(args.sample)

    with VariantFile(args.output, 'w', header=variants.header) as out:
        # Create parallel locus iterator by chromosome
        for _, ref, loci in records_by_chromosome(refs, [variants], [None], args):
            loci = sort_almost_sorted(loci[0], key=NormalizedLocus.left_order_key)

            for locus in loci:
                record  = locus.record
                start   = locus.left.start
                stop    = locus.left.stop
                alleles = locus.left.alleles

                if '' in alleles:
                    pad = ref[start - 1:start]
                    start -= 1
                    alleles = [pad + a for a in alleles]

                record.alleles = alleles
                record.start   = start
                record.stop    = stop

                out.write(record)


def tryint(s):
    """Try to convert s into an integer.  Otherwise return s unchanged.

    >>> tryint(1)
    1
    >>> tryint('1')
    1
    >>> tryint('one')
    'one'

    """
    try:
        return int(s)
    except ValueError:
        return s


def add_common_args(parser):
    """Add arguments common to multiple commands."""
    parser.add_argument('--reference', metavar='FASTA', required=True,
                        help='Reference FASTA+FAI (required)')
    parser.add_argument('--lazy-ref', action='store_true', help='Read reference file as needed (lazily), rather than all upfront.')
    parser.add_argument('-p', '--reference-padding', metavar='N', type=int, default=2,
                        help='Pad variants by N bp when forming superloci (default=2)')
    parser.add_argument('--include-regions', metavar='BED', help='BED file of regions to include in comparison')
    parser.add_argument('--exclude-regions', metavar='BED', help='BED file of regions to exclude from comparison')
    parser.add_argument('--include-filter', metavar='F', action='append',
                        help='Include records with filter status F.  Option may be specified multiple times or F can be comma delimited')
    parser.add_argument('--exclude-filter', metavar='F', action='append',
                        help='Exclude records with filter status F.  Option may be specified multiple times or F can be comma delimited')
    parser.add_argument('--min-gq', metavar='N', type=int,
                        help='Exclude records with genotype quality (GQ) < N')


def arg_parser():
    """Build vgraph's argument parser."""
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='Commands', dest='command')

    repmatch_parser = subparsers.add_parser('repmatch', help='compare two replicate samples')
    repmatch_parser.add_argument('vcf1', help='Sample 1 VCF/BCF input (- for stdin)')
    repmatch_parser.add_argument('vcf2', help='Sample 2 VCF/BCF input (- for stdin)')
    repmatch_parser.add_argument('--out1', help='Sample 1 VCF/BCF output (optional)')
    repmatch_parser.add_argument('--out2', help='Sample 2 VCF/BCF output (optional)')
    repmatch_parser.add_argument('--name1', metavar='N', default=0, type=tryint,
                                 help='Name or index of sample in sample 1 file (default=0)')
    repmatch_parser.add_argument('--name2', metavar='N', default=0, type=tryint,
                                 help='Name or index of sample in sample 2 file (default=0)')

    add_common_args(repmatch_parser)
    repmatch_parser.set_defaults(func=match_replicates)

    dbmatch_parser = subparsers.add_parser('dbmatch', help='compare a database of alleles to a sample')
    dbmatch_parser.add_argument('database', help='Database of alleles VCF/BCF input (- for stdin)')
    dbmatch_parser.add_argument('sample',   help='Sample VCF/BCF input (- for stdin)')
    dbmatch_parser.add_argument('--name', metavar='N', default=0, type=tryint, help='Name or index of sample in sample file (default=0)')
    dbmatch_parser.add_argument('-o', '--output', default='-', help='VCF/BCF output')

    add_common_args(dbmatch_parser)
    dbmatch_parser.set_defaults(func=match_database)

    dbmatch2_parser = subparsers.add_parser('dbmatch2', help='compare a database of alleles to a sample')
    dbmatch2_parser.add_argument('database', help='Database of alleles VCF/BCF input (- for stdin)')
    dbmatch2_parser.add_argument('sample',   help='Sample VCF/BCF input (- for stdin)')
    dbmatch2_parser.add_argument('--name', metavar='N', default=0, type=tryint,
                                 help='Name or index of sample in sample file (default=0)')
    dbmatch2_parser.add_argument('-o', '--output', default='-', help='VCF/BCF output')
    dbmatch2_parser.add_argument('-t', '--table', default='-', help='tabular variant output by id')

    add_common_args(dbmatch2_parser)
    dbmatch2_parser.set_defaults(func=match_database2)

    norm_parser = subparsers.add_parser('norm', help='Normalize a VCF/BCF file')
    norm_parser.add_argument('sample',   help='Sample VCF/BCF input (- for stdin)')
    norm_parser.add_argument('-o', '--output', default='-', help='VCF/BCF output')

    add_common_args(norm_parser)
    norm_parser.set_defaults(func=normalize)

    parser.add_argument('--debug', action='store_true', help='Output extremely verbose debugging information')
    parser.add_argument('--profile', action='store_true', help='Profile code performance')

    return parser


def run_vgraph(parser, args):
    """Run vgraph."""
    if args.command and args.func:
        args.func(args)
    else:
        parser.print_help()


def main():
    """Just vgraph's main CLI function."""
    parser = arg_parser()
    args = parser.parse_args()

    if args.profile:
        import yappi
        yappi.start()
        try:
            run_vgraph(parser, args)
        finally:
            yappi.stop()
            stats = yappi.get_func_stats().sort('tsub').strip_dirs()
            stats.print_all(out=sys.stderr, columns={0: ('name', 45), 1: ('ncall', 10), 2: ('tsub', 8), 3: ('ttot', 8), 4: ('tavg', 8)})
    else:
        run_vgraph(parser, args)


if __name__ == '__main__':
    main()
