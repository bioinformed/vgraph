from __future__ import division, print_function

from os.path          import expanduser
from argparse         import ArgumentParser
from collections      import defaultdict
from operator         import attrgetter

from pysam            import VariantFile, Fastafile

from vgraph.bed       import BedFile
from vgraph.norm      import NormalizedLocus
from vgraph.intervals import intersect
from vgraph.iterstuff import sort_almost_sorted, is_empty_iter
from vgraph.vargraph  import generate_genotypes


def records_to_loci(ref, records):
    # Prevents left shuffling the locus through the previous locus
    last_left_stop = 0
    for i, record in enumerate(records):
        if all(not alt.startswith('<') for alt in record.alts):
            locus = NormalizedLocus(i, record, ref, last_left_stop)
            last_left_stop = locus.left.stop
            yield locus


def informative_chromosomes(vars):
    return (chrom for chrom in vars.index if not is_empty_iter(vars.fetch(chrom)))


def make_bedmap(bedfile):
    bedmap = defaultdict(list)
    for record in bedfile:
        bedmap[record.contig].append(record)
    return bedmap


def include_filter(records, include):
    for start, stop, (rec, inc) in intersect([records, include]):
        if inc:
            for record in rec:
                yield record


def variants_by_chromosome(refs, vars, include=None):
    vars   = [VariantFile(var) for var in vars]

    chroms = [set(informative_chromosomes(var)) for var in vars]
    chroms = set.union(*chroms)

    if include is not None:
        include = make_bedmap(BedFile(include))
        chroms &= set(include)

    for chrom in chroms:
        ref  = refs.fetch(chrom).upper()
        loci = [records_to_loci(ref, var.fetch(chrom)) for var in vars]
        loci = [sort_almost_sorted(l, key=attrgetter('start', 'stop')) for l in loci]
        if include is not None:
            loci = [include_filter(l, include[chrom]) for l in loci]
        yield chrom, ref, loci


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
    parser.add_argument('vcf1', help='VCF/BCF input 1 (- for stdin).')
    parser.add_argument('vcf2', help='VCF/BCF input 2 (- for stdin).')
    parser.add_argument('--name1', metavar='N', default=0, type=tryint,
                        help='Name or index of sample in vcf1 (default=0).')
    parser.add_argument('--name2', metavar='N', default=0, type=tryint,
                        help='Name or index of sample in vcf2 (default=0).')
    parser.add_argument('-p', '--reference-padding', metavar='N', default=2,
                        help='Force loci within N bp into the same super locus (default=2).')
    parser.add_argument('-i', '--include', help='BED file of high confidence regions to compare')
    parser.add_argument('-o', '--out-vcf', default='-', help='Output VCF (- for stdout).')
    parser.add_argument('--reference', required=True, help='Reference FASTA')

    return parser.parse_args()


def main(args):
    # Load FASTA reference
    refs  = Fastafile(expanduser(args.reference))

    # Create parallel locus iterator by chromosome
    vars  = variants_by_chromosome(refs, [args.vcf1, args.vcf2], args.include)

    # Proceed by chromosome
    for chrom, ref, loci in vars:
        # Create superloci by intersecting locus streams
        superloci = intersect(loci, min_distance=args.reference_padding)

        # Proceed by superlocus
        for start, stop, (super1, super2) in superloci:
            # Resort each superlocus back into original VCF order
            super1 = sorted(super1, key=attrgetter('recnum'))
            super2 = sorted(super2, key=attrgetter('recnum'))

            # Create genotype sets for each superlocus
            genos1 = set(generate_genotypes(ref, start, stop, super1, args.name1))
            genos2 = set(generate_genotypes(ref, start, stop, super2, args.name2))

            # Test whether genotype sets intersect
            match  = bool(genos1 & genos2)

            # The hard work is done.  The rest is just output and formatting...
            status = 'MATCH!' if match else 'MISMATCH!'
            print('{}:[{:d}-{:d}): {}'.format(chrom, start, stop, status))

            for i, superlocus in enumerate([super1, super2], 1):
                for locus in superlocus:
                    print('  VCF{:d}: {}'.format(i, locus.record), end='')
            print()

            for i, (name, superlocus) in enumerate([(args.name1, super1), (args.name2, super2)], 1):
                for locus in superlocus:
                    lstart = locus.left.start
                    lstop = locus.left.stop
                    sample = locus.record.samples[name]
                    indices = sample.allele_indices
                    sep = '|' if sample.phased else '/'
                    geno = sep.join(locus.left.alleles[a] or '-' if a is not None else '.' for a in indices)
                    print('  NORM{:d}: [{:5d}-{:5d}) ref={} geno={}'.format(i, lstart, lstop, locus.left.alleles[0], geno))
            print()


def cli():
    main(parse_args())

if __name__ == '__main__':
    cli()
