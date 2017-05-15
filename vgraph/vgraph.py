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

from argparse           import ArgumentParser

from vgraph.repmatch    import match_replicates
from vgraph.dbmatch     import match_database


def tryint(s):
    '''Try to convert s into an integer.  Otherwise return s unchanged.

    >>> tryint(1)
    1
    >>> tryint('1')
    1
    >>> tryint('one')
    'one'
    '''
    try:
        return int(s)
    except ValueError:
        return s


def add_common_args(parser):
    parser.add_argument('--reference', metavar='FASTA', required=True,
                        help='Reference FASTA+FAI (required)')
    parser.add_argument('--lazy-ref', action='store_true', help='Read reference file as needed (lazily), rather than all upfront.')
    parser.add_argument('-p', '--reference-padding', metavar='N', default=2,
                        help='Pad variants by N bp when forming superloci (default=2)')
    parser.add_argument('--include-regions', metavar='BED', help='BED file of regions to include in comparison')
    parser.add_argument('--exclude-regions', metavar='BED', help='BED file of regions to exclude from comparison')
    parser.add_argument('--include-file-regions', metavar='BED', action='append',
                        help='BED file of regions to include for each input file')
    parser.add_argument('--exclude-file-regions', metavar='BED', action='append',
                        help='BED file of regions to exclude from comparison for each input file')
    parser.add_argument('--include-filter', metavar='F', action='append',
                        help='Include records with filter status F.  Option may be specified multiple times or F can be comma delimited')
    parser.add_argument('--exclude-filter', metavar='F', action='append',
                        help='Exclude records with filter status F.  Option may be specified multiple times or F can be comma delimited')
    parser.add_argument('--min-gq', metavar='N', type=int,
                        help='Exclude records with genotype quality (GQ) < N')
    #parser.add_argument('-o', '--out-vcf', default='-', help='Output VCF (- for stdout)')


def arg_parser():
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

    dbmatch_parser = subparsers.add_parser('dbmatch', help='compare a database of alleles to a sample')
    dbmatch_parser.add_argument('database', help='Database of alleles VCF/BCF input (- for stdin)')
    dbmatch_parser.add_argument('sample',   help='Sample VCF/BCF input (- for stdin)')
    dbmatch_parser.add_argument('--name', metavar='N', default=0, type=tryint,
                                 help='Name or index of sample in sample file (default=0)')
    dbmatch_parser.add_argument('-o', '--output', default='-', help='Sample VCF/BCF output')

    add_common_args(dbmatch_parser)

    parser.add_argument('--debug', action='store_true', help='Output extremely verbose debugging information')
    parser.add_argument('--profile', action='store_true', help='Profile code performance')

    return parser


def run_analysis(parser, args):
    if not args.command:
        parser.print_help()
    elif args.command == 'repmatch':
        match_replicates(args)
    elif args.command == 'dbmatch':
        match_database(args)
    else:
        raise ValueError('Unknown command: {}'.format(args.command))


def main():
    parser = arg_parser()
    args = parser.parse_args()

    if args.profile:
        import yappi
        yappi.start()
        try:
            run_analysis(parser, args)
        finally:
            yappi.stop()
            stats = yappi.get_func_stats().sort('tsub').strip_dirs()
            stats.print_all(out=sys.stderr, columns={0: ('name', 45), 1: ('ncall', 10), 2: ('tsub', 8), 3: ('ttot', 8), 4: ('tavg', 8)})
    else:
        run_analysis(parser, args)


if __name__ == '__main__':
    main()
