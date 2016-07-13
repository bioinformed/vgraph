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

**NOTICE:** ``vgraph`` is under heavy development and is not yet ready for
production use.  The graph model and comparison logic is becoming mature and
continues to be validated and fine tuned on data from diverse variant
callers.  However, output is still very verbose, diagnostic heavy and
doesn't have built-in aggregation, but is usable with basic technical
knowledge.  Please treat these versions as "alpha" quality that are suitable
for evaluation by other tool developers and skilled technical users.  The
first beta version will be ready shortly and will be more "end user"
friendly and will have pretty output and nice aggregation and summary output
statistics.  A production-ready release for the masses is planned in the
near future.

Input
-----

``vgraph`` currently accepts block gzipped and indexed VCF and BCF files.
Support for Complete Genomics ``var`` and ``masterVar`` formats will be
added in a future version.  ``vgraph`` also requires an indexed reference
genome in FASTA+FAI format.

Output
------

``vgraph`` is current in an early development state and currently outputs
diagnostic only information.  The output will be improved in future versions
to generate::

    1. Annotated VCF/BCF files with per-record match status
    2. Detailed output and diagnostic information for mismatching loci
    3. Summary statistics (overall and stratified by variant type and context)

Usage
-----

``vgraph`` takes the following command line parameters::

    usage: vgraph [-h] [--name1 N] [--name2 N] [-p N] [-i INCLUDE] [-o OUT_VCF]
                  --reference FASTA [--debug] [--profile]
                  vcf1 vcf2

    positional arguments:
      vcf1                  VCF/BCF input 1 (- for stdin).
      vcf2                  VCF/BCF input 2 (- for stdin).

    optional arguments:
      -h, --help            show this help message and exit
      --name1 N             Name or index of sample in vcf1 (default=0).
      --name2 N             Name or index of sample in vcf2 (default=0).
      -p N, --reference-padding N
                            Force loci within N bp into the same super locus
                            (default=2).
      -i INCLUDE, --include INCLUDE
                            BED file of high confidence regions to compare
      -o OUT_VCF, --out-vcf OUT_VCF
                            Output VCF (- for stdout).
      --reference FASTA     Reference FASTA+FAI
      --debug               Output extremely verbose debugging information
      --profile             Profile code performance

Installation
------------

Before ``vgraph`` may be installed, your systems requires a C compiler, a
functioning version of Python 2.7 with development libraries installed, and
the ``pip`` installer.  The steps to install and ensuring these tools are
functional depend on your operating system and personal configuration. 
Proceed only once these pre-requisites are available.

First install the latest version of the Cython package::

    pip install -U Cython

``vgraph`` currently requires a pre-release version of the ``pysam``
package.  This requirement will be lifted soon and an official release of
``pysam`` will be required instead.  In the mean time, install this
version::

    pip install -U git+https://github.com/pysam-developers/pysam

If all these steps have succeeded, then finally install ``vgraph``::

    pip install -U git+https://github.com/bioinformed/vgraph.git
