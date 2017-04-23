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


from collections import defaultdict
from operator    import attrgetter
from itertools   import groupby

from .smartfile  import smartfile


class BedRecord(object):
    """
    Simple class for working with records from BED files.
    """
    __slots__ = ('contig', 'start', 'stop', 'name', 'score', 'strand', 'thick_start', 'thick_end', 'item_rgb')
    field_names = __slots__

    def __init__(self, contig, start, stop, name=None, score=None, strand=None, thick_start=None,
                 thick_end=None, item_rgb=None):

        self.contig = contig
        self.start = start
        self.stop = stop
        self.name = name
        self.score = score
        self.strand = strand
        self.thick_start = thick_start
        self.thick_end = thick_end
        self.item_rgb = item_rgb

    @staticmethod
    def from_line(line):
        line = line.rstrip()
        fields = line.split('\t')
        contig, start, stop = fields[0], int(fields[1]), int(fields[2])

        n = len(fields)

        name        = fields[3] or None if n >= 4 else None
        score       = fields[4] or None if n >= 5 else None
        strand      = fields[5] or None if n >= 6 else None
        thick_start = fields[6] or None if n >= 7 else None
        thick_end   = fields[7] or None if n >= 8 else None
        item_rgb    = fields[8] or None if n >= 9 else None

        return BedRecord(contig, start, stop, name, score, strand, thick_start, thick_end, item_rgb)

    def to_tuple(self):
        return (self.contig, self.start, self.stop, self.name, self.score,
                self.strand, self.thick_start, self.thick_end, self.item_rgb)

    def to_line(self):
        line = '\t'.join(str(f if f is not None else '') for f in self.to_tuple())
        return line.rstrip()

    def __repr__(self):
        fields = ('%s=%r' % (k, v) for k, v in zip(self.field_names, self.to_tuple()) if v not in (None, ''))
        return 'BedRecord(%s)' % ', '.join(fields)


class BedFile(object):
    """
    Simple class for iterating through the records of a BED file.
    """
    def __init__(self, filename):
        self.filename = filename
        self._tabix = None

    def __iter__(self):
        return BedFile.parse_bed_lines(smartfile(self.filename))

    @property
    def tabix(self):
        if self._tabix:
            return self._tabix

        import pysam
        self._tabix = pysam.Tabixfile(self.filename)

        return self._tabix

    def query(self, contig=None, start=None, stop=None):
        records = self.tabix.fetch(contig, start, stop)
        return BedFile.parse_bed_lines(records)

    @property
    def contigs(self):
        return self.tabix.contigs

    @staticmethod
    def parse_bed_lines(lines):
        for line in lines:
            line = line.rstrip()

            if not line or line.startswith(';') or line.startswith('track '):
                continue

            yield BedRecord.from_line(line)


def load_bedmap(filename):
    '''Load BED file as a dictionary mapping contig to list of BedRecords

    Args:
        filename (str): input filename

    Returns:
        dict: dictionary mapping contig to list of BedRecords
    '''
    bed = sorted(BedFile(filename), key=attrgetter('contig', 'start', 'stop'))

    bedmap = defaultdict(list)
    for contig, contig_records in groupby(bed, attrgetter('contig')):
        bedmap[contig] = list(contig_records)

    return dict(bedmap)
