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

"""Match a genome to a database of alleles."""

import csv
import sys

from os.path            import expanduser
from operator           import attrgetter

from pysam              import VariantFile, Fastafile

from vgraph.norm        import NormalizedLocus
from vgraph.intervals   import union
from vgraph.iterstuff   import sort_almost_sorted
from vgraph.match       import records_by_chromosome, get_superlocus_bounds, find_allele


def annotate_info(locus, allele, info_meta, suffix, times):
    """Annotate INFO in sample with fields from database."""
    for name in info_meta:
        if name in allele.record.info:
            sname = name + suffix
            orig_value = locus.record.info.get(sname, ())
            new_value = allele.record.info[name]
            if not isinstance(new_value, tuple):
                new_value = (new_value,)
            locus.record.info[sname] = orig_value + new_value * times


def annotate_format(locus, allele, format_meta, suffix, times):
    """Annotate FORMAT in sample with fields from database."""
    sample = locus.record.samples[0]
    for name in format_meta:
        if name in allele.record.format:
            sname = name + suffix
            orig_value = sample.get(sname, ())
            new_value = allele.record.samples[0][name]
            if not isinstance(new_value, tuple):
                new_value = (new_value,)
            sample[sname] = orig_value + new_value * times


def match_sort_key(match):
    """Return number of matching alleles to order matches."""
    if not match:
        return 0
    return match.allele_ploidy


def generate_superlocus_matches(chrom, superlocus, ref, alleles, mode, debug=False):
    """Generate allele matches for a superlocus."""
    for allele in alleles:
        super_allele = [locus for locus in superlocus if locus.extremes_intersect(allele)]

        # Remove all reference calls from the superlocus.
        # This is primarily done to remove long leading and trailing reference regions.
        # Interstitial reference regions will be added back, based on how gaps are handled.
        super_non_ref = [locus for locus in super_allele if not locus.is_ref()]

        if debug:
            super_start, super_stop = get_superlocus_bounds([[allele], super_non_ref])
            print('-' * 80, file=sys.stderr)
            print('{}:[{:d}-{:d}):'.format(chrom, super_start, super_stop), file=sys.stderr)
            print(file=sys.stderr)

            print('  ALLELE: {} {}:[{}-{}) ref={} alt={}'.format(
                allele.record.id,
                allele.contig,
                allele.start,
                allele.stop,
                allele.alleles[0] or '-',
                allele.alleles[1] or '-'
            ), file=sys.stderr)
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

        # Search superlocus for each ALT allele; stop if any are found
        matches = (find_allele(ref, allele, allele_index, super_non_ref, mode=mode, debug=debug) for allele_index in range(1, len(allele.alleles)))
        matches = sorted(matches, key=match_sort_key)

        match = matches[-1] if matches else None

        if debug:
            print(file=sys.stderr)
            print('    MATCH={}'.format(match), file=sys.stderr)
            print(file=sys.stderr)

        yield super_allele, allele, match


def generate_matches(refs, sample, db, args):
    """Generate allele matches over all chromosomes."""
    for chrom, ref, loci in records_by_chromosome(refs, [sample, db], [args.name, None], args):
        # Create superloci by taking the union of overlapping loci across all of the locus streams
        loci = [sort_almost_sorted(l, key=NormalizedLocus.extreme_order_key) for l in loci]
        superloci = union(loci, interval_func=attrgetter('min_start', 'max_stop'))

        for _, _, (superlocus, alleles) in superloci:
            alleles.sort(key=NormalizedLocus.natural_order_key)
            superlocus.sort(key=NormalizedLocus.natural_order_key)

            yield superlocus, generate_superlocus_matches(chrom, superlocus, ref, alleles, args.mode, args.debug)


def build_new_metadata(db, sample):
    """Build new metadata definitions for sample matches."""
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

    return format_meta, info_meta


def translate_match(match):
    """Translate match to STATUS and TIMES."""
    status, times = 'NOTFOUND', 1
    if match is None:
        status = 'NOCALL'
    elif match.allele_ploidy:
        status = 'FOUND'
        times  = match.allele_ploidy

    return status, times


def match_database(args):
    """Match a genome to a database of alleles."""
    refs   = Fastafile(expanduser(args.reference))
    db     = VariantFile(expanduser(args.database))
    sample = VariantFile(expanduser(args.sample))

    format_meta, info_meta = build_new_metadata(db, sample)

    with VariantFile(args.output, 'w', header=sample.header) as out:
        for superlocus, matches in generate_matches(refs, sample, db, args):
            for allele_locus, allele, match in matches:
                # Annotate results of search
                status, times = translate_match(match)
                suffix = '_' + status

                for locus in allele_locus:
                    annotate_info(locus, allele, info_meta, suffix, times)
                    annotate_format(locus, allele, format_meta, suffix, times)

            for locus in sorted(superlocus, key=NormalizedLocus.record_order_key):
                out.write(locus.record)


def update_info_header(header):
    """Add match INFO fields VCF header for dbmatch2."""
    info_header = header.info
    if 'FOUND' not in info_header:
        info_header.add('FOUND',    number='.', type='String', description='Allele(s) found')
    if 'NOTFOUND' not in info_header:
        info_header.add('NOTFOUND', number='.', type='String', description='Allele(s) not found')
    if 'NOCALL' not in info_header:
        info_header.add('NOCALL',   number='.', type='String', description='Allele(s) not called due to uncertainty')


def clear_info_fields(loci):
    """Clear INFO fields, if present, prior to setting them."""
    for locus in loci:
        info = locus.record.info
        for status in ('FOUND', 'NOTFOUND', 'NOCALL'):
            if status in info:
                del info[status]


def write_table_header(out):
    """Write a match header for the tabular output of dbmatch2."""
    out.writerow([
        'SAMPLE_ID',
        'ALLELE_ID',
        'STATUS',
        'CALL_QUALITY',
        'ALLELE_PLOIDY',
        'REFERENCE_PLOIDY',
        'OTHER_PLOIDY',
        'ALLELE_READS',
        'REFERENCE_READS',
        'OTHER_READS',
    ])


def write_table_row(out, sample_name, var_id, superlocus, status, match):
    """Write a match row for the tabular output of dbmatch2."""
    if not out:
        return

    gts  = [locus.record.samples[sample_name] for locus in superlocus]
    qual = min(gt.get('GQ', 0) for gt in gts) if gts else ''

    row = [
        sample_name,
        var_id,
        status,
        qual,
    ]

    if match:
        row += [
            match.allele_ploidy, match.ref_ploidy, match.other_ploidy,
            match.allele_depth if match.allele_depth is not None else 'NOT_CALLED',
            match.ref_depth,   # ref allele depth should always be reported if records are present
            match.other_depth  if match.other_depth  is not None else 'NOT_CALLED',

        ]
    else:
        row += ['NO_CALL'] * 6

    out.writerow(row)


def match_database2(args):
    """Match a genome to a database of alleles."""
    refs   = Fastafile(expanduser(args.reference))
    db     = VariantFile(expanduser(args.database))
    sample = VariantFile(expanduser(args.sample))

    try:
        sample_name = sample.header.samples[args.name]
    except TypeError:
        sample_name = args.name

    if db.index is None:
        raise ValueError('database file must be indexed')
    if sample.index is None:
        raise ValueError('sample file must be indexed')

    # Open tabluar output file, if requested
    table = None
    if args.table:
        tablefile = open(args.table, 'w') if args.table != '-' else sys.stdout
        table = csv.writer(tablefile, delimiter='\t', lineterminator='\n')
        write_table_header(table)

    update_info_header(sample.header)

    with VariantFile(args.output, 'w', header=sample.header) as out:
        for superlocus, matches in generate_matches(refs, sample, db, args):
            clear_info_fields(superlocus)

            for allele_locus, allele, match in matches:
                dbvar  = allele.record
                var_id = dbvar.id or f'{dbvar.chrom}_{dbvar.start+1}_{dbvar.stop}_{dbvar.alts[0]}'

                status, times = translate_match(match)

                for locus in allele_locus:
                    info = locus.record.info
                    info[status] = info.get(status, ()) + (var_id, ) * times

                write_table_row(table, sample_name, var_id, allele_locus, status, match)

            for locus in sorted(superlocus, key=NormalizedLocus.record_order_key):
                out.write(locus.record)
