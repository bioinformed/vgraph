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


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('vcf1', help='VCF/BCF input 1 (- for stdin)')
    parser.add_argument('vcf2', help='VCF/BCF input 2 (- for stdin)')
    parser.add_argument('--out1', help='VCF/BCF output 1')
    parser.add_argument('--out2', help='VCF/BCF output 2')
    parser.add_argument('--name1', metavar='N', default=0, type=tryint,
                        help='Name or index of sample in vcf1 (default=0)')
    parser.add_argument('--name2', metavar='N', default=0, type=tryint,
                        help='Name or index of sample in vcf2 (default=0)')
    parser.add_argument('-p', '--reference-padding', metavar='N', default=3,
                        help='Force loci within N bp into the same super locus (default=3)')
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
    parser.add_argument('--reference', metavar='FASTA', required=True, help='Reference FASTA+FAI')
    parser.add_argument('--debug', action='store_true', help='Output extremely verbose debugging information')
    parser.add_argument('--profile', action='store_true', help='Profile code performance')

    return parser.parse_args()


def main():
    args = parse_args()
    if args.profile:
        import yappi
        yappi.start()
        match_replicates(args)
        yappi.stop()
        stats = yappi.get_func_stats().sort('tsub').strip_dirs()
        stats.print_all(out=sys.stderr, columns={0: ('name', 45), 1: ('ncall', 10), 2: ('tsub', 8), 3: ('ttot', 8), 4: ('tavg', 8)})
    else:
        match_replicates(args)


if __name__ == '__main__':
    main()
