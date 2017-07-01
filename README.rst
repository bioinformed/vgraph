Variant Graph Comparison Tool
=============================

Kevin Jacobs <jacobs@bioinfomed.com>

``vgraph`` is a command line application and Python library to compare
genetic variants using variant graphs.  Conventional methods used to compare
variants apply heuristic normalization rules and then compare variants
individually by matching based on genomic position and allele information.
In contrast, ``vgraph`` utilizes a graph representation of genomic variants
to precisely compare complex variants that are refractory to comparison
using conventional methods.

Input
-----

``vgraph`` currently accepts block gzipped and indexed VCF and BCF files.
Support for Complete Genomics ``var`` and ``masterVar`` formats will be
added in a future version.  ``vgraph`` also requires an indexed reference
genome in FASTA+FAI format.

Output
------

``vgraph`` outputs diagnostic only information to stdout.  In `repmatch`
mode, there are options to output either of the two input files with match
status annotations.  In `dbmatch` mode, the `sample` input file is output
after copying all new INFO and FORMT annotations from the `database` file.

Usage
-----

``vgraph`` takes the following command line parameters::

    usage: vgraph [-h] [--debug] [--profile] {repmatch,dbmatch} ...

    positional arguments:
      {repmatch,dbmatch}  Commands
        repmatch          compare two replicate samples
        dbmatch           compare a database of alleles to a sample

    optional arguments:
      -h, --help          show this help message and exit
      --debug             Output extremely verbose debugging information
      --profile           Profile code performance

The parameters for ``repmatch`` are::

      usage: vgraph repmatch [-h] [--out1 OUT1] [--out2 OUT2] [--name1 N]
                             [--name2 N] --reference FASTA [-p N]
                             [--include-regions BED] [--exclude-regions BED]
                             [--include-file-regions BED]
                             [--exclude-file-regions BED] [--include-filter F]
                             [--exclude-filter F] [--min-gq N]
                             vcf1 vcf2

      positional arguments:
        vcf1                  Sample 1 VCF/BCF input (- for stdin)
        vcf2                  Sample 2 VCF/BCF input (- for stdin)

      optional arguments:
        -h, --help            show this help message and exit
        --out1 OUT1           Sample 1 VCF/BCF output (optional)
        --out2 OUT2           Sample 2 VCF/BCF output (optional)
        --name1 N             Name or index of sample in sample 1 file (default=0)
        --name2 N             Name or index of sample in sample 2 file (default=0)
        --reference FASTA     Reference FASTA+FAI (required)
        -p N, --reference-padding N
                              Pad variants by N bp when forming superloci
                              (default=2)
        --include-regions BED
                              BED file of regions to include in comparison
        --exclude-regions BED
                              BED file of regions to exclude from comparison
        --include-file-regions BED
                              BED file of regions to include for each input file
        --exclude-file-regions BED
                              BED file of regions to exclude from comparison for
                              each input file
        --include-filter F    Include records with filter status F. Option may be
                              specified multiple times or F can be comma delimited
        --exclude-filter F    Exclude records with filter status F. Option may be
                              specified multiple times or F can be comma delimited
        --min-gq N            Exclude records with genotype quality (GQ) < N


The parameters for ``dbmatch`` are::

      usage: vgraph dbmatch [-h] [--name N] [-o OUTPUT] --reference FASTA [-p N]
                            [--include-regions BED] [--exclude-regions BED]
                            [--include-file-regions BED]
                            [--exclude-file-regions BED] [--include-filter F]
                            [--exclude-filter F] [--min-gq N]
                            database sample

      positional arguments:
        database              Database of alleles VCF/BCF input (- for stdin)
        sample                Sample VCF/BCF input (- for stdin)

      optional arguments:
        -h, --help            show this help message and exit
        --name N              Name or index of sample in sample file (default=0)
        -o OUTPUT, --output OUTPUT
                              Sample VCF/BCF output
        --reference FASTA     Reference FASTA+FAI (required)
        -p N, --reference-padding N
                              Pad variants by N bp when forming superloci
                              (default=2)
        --include-regions BED
                              BED file of regions to include in comparison
        --exclude-regions BED
                              BED file of regions to exclude from comparison
        --include-file-regions BED
                              BED file of regions to include for each input file
        --exclude-file-regions BED
                              BED file of regions to exclude from comparison for
                              each input file
        --include-filter F    Include records with filter status F. Option may be
                              specified multiple times or F can be comma delimited
        --exclude-filter F    Exclude records with filter status F. Option may be
                              specified multiple times or F can be comma delimited
        --min-gq N            Exclude records with genotype quality (GQ) < N


Installation
------------

Before ``vgraph`` may be installed, your systems requires a C compiler, a
functioning version of Python 3.5 or newer with development libraries
installed, and the ``pip`` installer.  The steps to install and ensuring
these tools are functional depend on your operating system and personal
configuration.  Proceed only once these pre-requisites are available.

First install the latest version of the Cython and pysam packages::

    pip install -U Cython
    pip install -U pysam

If all these steps have succeeded, then install ``vgraph``::

    pip install -U git+https://github.com/bioinformed/vgraph.git
