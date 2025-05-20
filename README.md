# Haplotype Comparison Tools (hap.py)

Peter Krusche <pkrusche@illumina.com>

This is a set of programs based on [htslib](https://github.com/samtools/htslib)
to benchmark variant calls against gold standard truth datasets.

> **Note:** This project has been migrated to Python 3.
> See the [Python 3 migration guide](doc/python3_migration.md)
> for details about the migration and compatibility.

To compare a VCF against a gold standard dataset, use the following commmand line
to perform genotype-level haplotype comparison.

```bash
# Python 3 version
python3 /path/to/hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa

# Or if installed via pip
hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa
```

We also have a script to perform comparisons only based on chromosome, position,
and allele identity. This comparison will not resolve haplotypes and only verify
that the same alleles were observed at the same positions (e.g. for comparison
of somatic callsets).

```bash
som.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa
```

More information can be found below in the [usage section](#usage).

## Contents

* [Motivation](#motivation)
* [Complex variant comparison](#complex-variant-comparison)
* [Variant preprocessing](#variant-preprocessing)
* [Variant counting](#variant-counting)
* [Usage](#usage)
  * [hap.py](#happy)
  * [som.py](#sompy)
* [Installation](#installation)
  * [Building from Source (Advanced)](#building-from-source-advanced)
* [Quick Start](#quick-start)
* [Legacy Installation (Deprecated)](#legacy-installation-deprecated)
* [System requirements](#system-requirements)
  * [Hardware](#hardware)
  * [Linux](#linux)
  * [OS X](#os-x)
  * [Windows](#windows)
  * [Other requirements](#other-requirements)
* [Python 3 Migration](#python-3-migration)
* [Key Features](#key-features)

## Motivation

### Complex variant comparison

A major challenge when comparing VCF files for diploid samples is the handling
of complex variant representations. In a VCF file, we describe two haplotype
sequences by means of REF-ALT pairs and genotypes. These variant calls do not
uniquely represent the haplotype sequences: since alignments are always not unique
even when using a fixed set of gap and substitution scores,
different variant calling methods may produce different variant representations.
While some of these representational differences can be handled using
pre-processing of VCF files (e.g. variant trimming and left-shifting), others
cannot be fixed easily.

In addition to comparing VCF records individually, we produce a graph-based representation
of the VCF alleles, create all possible haplotype sequences, and compare
these by alignment / exact matching. Here is an example where this is needed:

*Variant representation 1 (shown in purple in the image below):*

```
CHROM POS   REF  ALT             GT
chrQ  10    G    GTGTGTGCATGCT   0/1
```

*Variant representation 2 (shown in green in the image below):*

```
CHROM POS   REF  ALT             GT
chrQ  16    G    GCATGCT         0/1
chrQ  19    T    TGTGTG          0/1
```

```bash
# Command to run hap.py for complex comparison (example)
./hap.py truth.vcf query.vcf -o output/prefix -r ref.fa --engine=vcfeval --eval-outside-conf
```

```bash
# Another command example
eval_out=output/prefix
./hap.py ${truth_vcf} ${query_vcf} \
    -f ${conf_bed} \
    -r ${ref_fa} -o ${eval_out} \
    -l chr1,chr2 --no-json --engine=xcmp \
    --force-interactive # for testing / demo purposes only
```

![Example representation of variants](doc/rep_ex.PNG "Example Variant Representation")

Both representations in this example are able to produce the same alt sequences,
but we are not able to match them up with standard VCF tools. In particular,
we can see from this example that the second representation actually may allow us
to create two different sets of alt sequences if they are part of unphased
heterozygous variant calls. When we don't know the phasing
of our variants, the insertions could have occurred on different haplotypes when using
representation 2.

With this tool, we can produce all haplotypes sequences by enumerating paths
through a reference graph. By finding the paths / alt alleles that are
consistent between two VCFs files we can produce accurate benchmarking
numbers for comparing a VCF to a gold standard truth set.
See [doc/spec.md](doc/spec.md) for more information.

An alternative method to compare  complex variant calls is implemented in
[RTG vcfeval](https://github.com/RealTimeGenomics/rtg-tools). It is possible
to use vcfeval with hap.py, and to use hap.py only for pre-processing,
stratification and counting.

The comparison method in vcfeval is more sophisticated than ours and can
resolve some corner cases more accurately.
For whole-genome comparisons, the difference between the two benchmarking
methods is small, but when focusing on difficult subsets of the genome or
when using variant calling methods that produce many complex variant calls,
these corner cases can become relevant. Moreover, when benchmarking against
gold-standard datasets that cover difficult regions of the genome (e.g.
[Platinum Genomes](http://www.illumina.com/platinumgenomes/)), the more complicated
subsets of the genome will be responsible for most of the difference between
methods.

### Variant preprocessing

Another component of hap.py is a variant pre-processing method which
deals with complex variant representations and MNPs. When different callers
may represent variants using a different number of VCF records, we should
attempt to count these in a consistent fashion between methods. One example
is the representation of MNVs as individual SNPs vs. as complex variants.

Consider the following case:

*Complex variant representation*:

```
chrQ  16    GGG    TTT         0/1
```

vs.

*Atomized representation*:

```
chrQ  16    G      T         0/1
chrQ  17    G      T         0/1
chrQ  18    G      T         0/1
```

If this variant is a false-positive, the first representation would naively
contribute a single FP record. A variant caller that outputs the second
representation would instead receive a penalty of three FPs for making
the same variant call. Overall, the difference between the two representations
might show significantly when looking at precision levels or false-positive
rates (since these are relative to the total number of query counts, which
use the same representations), but become important when we need to compare
absolute numbers of false-positives. For this case, hap.py can perform a re-alignment
of REF and ALT alleles on the query VCF, and splits the records into atomic
variant alleles to produce more granular counts using [pre.py](doc/normalisation.md).
Left-shifting and trimming are also supported.

```bash
# Example of preprocessing command
./pre.py in.vcf -o out.vcf -r ref.fa
```

```bash
# Another preprocessing example
./pre.py ${in_vcf} -r ${ref_fa} -o ${norm_vcf} --verbose --logfile ${norm_log} --profile
```

### Variant counting

Hap.py includes a module to produce stratified variant counts. Variant types
are determined using a re-alignment of REF and ALT alleles. This is more reliable
than only using allele lengths. Consider the following complex deletion.

```
chr1    201586350       .       CTCTCTCTCT      CA
```

This complex variant call is equivalent to a deletion, followed by a SNP. Our
quantification code will recognize this variant as a deletion and a SNP, and will
count it in both categories (so a TP call for this variant will contribute a
SNP and an INDEL). This effectively deals with variant calling methods that
prefer to combine local haplotypes in the same variant records
(e.g. Freebayes / Platypus), which would otherwise fall into a hard-to-assess
"COMPLEX" variant call category that varies substantially between
different variant calling methods.

```
chr1    201586350       .       CTCTCTCTC       C
chr1    201586359       .       T               A
```

Another feature of the quantification module in hap.py is stratification into
variant sub-types and into genomic regions. For example, precision and recall
can be computed at the same time for all
[GA4GH stratification regions](https://github.com/ga4gh/benchmarking-tools/tree/master/resources/stratification-bed-files),
and for different INDEL lengths (\<5, 7-15, 16+). Hap.py also calculates
het-hom and Ti/Tv ratios for all subsets of benchmarked variants.
Note that all region matching in hap.py is based on reference coordinates
only. One case where this can lead to counterintuitive results is when considering
hompolymer insertions:

```
Reference:

>chrQ
CAAAAA

VCF:
chrQ    1   C   CA  0/1

BED for homopolymers:
1   6
```

In this example, the variant call given above would not be captured by the bed region for the
homopolymers because it is associated with the reference base just before. To account for this,
the bed intervals need to be expanded to include the padding base just before the regions.

Finally, we produce input data for ROC and precision/recall curves. An
[example](doc/microbench.md) is included.

```bash
# Example of qfy.py command
./qfy.py truth.vcf query.vcf -o output/prefix -r ref.fa
```

```bash
# Another qfy.py example
./qfy.py ${truth_vcf} ${query_vcf} -o ${eval_out} -r ${ref_fa} -f ${conf_bed} -T ${target_bed} --threads 2
```

```bash
# Example of xcmp.py command
./xcmp.py truth.vcf query.vcf -o output/prefix -r ref.fa
```

## Usage

The main two tools are hap.py (diploid precision/recall evaluation) and som.py
(somatic precision/recall evaluation -- this ignores the GT and just checks for
presence of alleles). Other tools are qfy.py (which just executes the quantification
step of the analysis pipeline, this requires a
[GA4GH-intermediate](https://github.com/ga4gh/benchmarking-tools/) VCF file), and
[pre.py](doc/normalisation.md), which is hap.py's input cleaning and
variant normalisation step.

Here are some small example command lines. Advanced features like confident call
 / ambiguity / FP regions are also available, see the documentation for each
 tool for these.

Below, we assume that the code has been installed to the directory `${HAPPY}`.

### hap.py

See also [doc/happy.md](doc/happy.md).

```bash
$ ${HAPPY}/bin/hap.py  \
      example/happy/PG_NA12878_chr21.vcf.gz \
      example/happy/NA12878_chr21.vcf.gz \
      -f example/happy/PG_Conf_chr21.bed.gz \
      -o test
$ ls test.*
test.metrics.json  test.summary.csv
```

This example compares an example run of GATK 1.6 on NA12878 agains the Platinum
Genomes reference dataset (***Note: this is a fairly old version of GATK, so
don't rely on these particular numbers for competitive comparisons!***).

The summary CSV file contains all high-level metrics:

| Type          |  TRUTH.TOTAL|  QUERY.TOTAL | METRIC.Recall | METRIC.Precision | METRIC.Frac\_NA | TRUTH.TOTAL.TiTv\_ratio | QUERY.TOTAL.TiTv\_ratio | TRUTH.TOTAL.het\_hom\_ratio | QUERY.TOTAL.het\_hom\_ratio|
|---------------|-------------|--------------|---------------|------------------|-----------------|-------------------------|-------------------------|-----------------------------|----------------------------|
|INDEL          |         9124|         9905 |      0.869406 |         0.978441 |        0.194548 |                     NaN |                     NaN |                    1.463852 |                    1.209105|
|SNP            |        52520|        48078 |      0.894478 |         0.998258 |        0.021070 |                2.081002 |                2.082603 |                    1.595621 |                    1.487599|

These numbers tell us the SNP and indel recall of our query VCF against the
truth dataset. See [doc/happy.md](doc/happy.md) for more documentation and some
advice for their interpretation.

### som.py

Som.py is a simple comparison tool based on bcftools. It does not perform genotype or haplotype matching.

See [doc/sompy.md](doc/sompy.md) for more documentation.

```bash
# Example of som.py command
./som.py truth.vcf query.vcf -o output/prefix -r ref.fa
```

## Installation

### Using pip (Recommended)

hap.py can be installed using pip:

```bash
# Install from PyPI
pip install hap.py

# Or install from source directory
git clone https://github.com/Illumina/hap.py.git
cd hap.py
pip install .
```

This will build all necessary components including the C++ parts and install the Python package with command-line entry points.

To install with optional dependencies for C++/Cython extensions (recommended for performance) or development tools:

```bash
pip install .[cpp]      # For C++/Cython accelerated features
pip install .[dev]      # For development tools (testing, linting)
pip install .[cpp,dev]  # For both
```

### Building from Source (Advanced)

If you need to build from source and `pip install .` does not meet your needs (e.g., you need to customize the C++ build process extensively or are working in an environment without pip):

1. **Prerequisites**:

   * A C++14 compatible compiler (e.g., GCC, Clang, MSVC)
   * CMake (version 3.10 or newer)
   * Python (version 3.7 or newer, including development headers)
   * Boost libraries (version 1.55.0 or newer - iostreams, regex, filesystem, system, program_options). These can be automatically built by our scripts if not found system-wide.
   * Zlib development libraries.

2. **Configure and Build**:
   The `pyproject.toml` and CMake setup are designed to be handled by `pip`. For manual control, you would typically invoke CMake directly, but this is now an advanced use case. The `install.py` script is being deprecated.

   For developers, the standard Python build frontends should be used:

   ```bash
   python -m build
   ```

   This will produce a wheel in the `dist/` directory, which can then be installed with `pip install dist/hap.py-*.whl`.

## Quick Start

After installation, the `hap.py` command-line tool will be available.

```bash
hap.py --help # Show help message

# Example: Compare a VCF file against a truth VCF
hap.py truth.vcf.gz query.vcf.gz -r reference.fa -o output_prefix
```

(Further examples and detailed usage can be found in the documentation.)

## Legacy Installation (Deprecated)

The old `install.py` script is deprecated and will be removed in v1.0.0. Please migrate to using `pip install .` as described above.

## System requirements

### Hardware

Compiling and testing can be done on a standard desktop system with 8GB of RAM. Whole-genome
comparisons (e.g. comparing a gVCF file against the [Platinum Genomes truth dataset](http://www.illumina.com/platinumgenomes/))
can use up to 64GB of RAM (20GB typical, depending on the input VCF) and about 4-12 minutes
using 40 processor cores. Whole exome comparison (using an exome bed mask and the `-T` switch)
can be carried out on a desktop system.

### Linux

Tested on:

```text
Ubuntu 12.04,14.04,16.04,18.04
CentOS 6.x, 7.x
```

Hap.py must be compiled with g++ version 4.9.x or later, or with a recent version of Clang (testing is performed
with g++).

### OS X

Hap.py builds and passes basic tests on OS X 10.9+, but full WGS analyses are not tested for this platform.

### Windows

Hap.py is not tested on Windows. The main dependency that fails compilation is htslib. Given a build
of htslib and pysam, using hap.py on Windows should be possible.

### Other requirements

Hap.py requires a human genome reference sequence which contains at least
chromosomes 1-22,X,Y,M. The chromosomes should be named chr1-chr22, chrX, chrY,
chrM. there is a script  in [src/sh/make_hg19.sh](src/sh/make_hg19.sh) to create
such a sequence, but you can also  specify your own. In order for the
integration tests to run successfully, it is necessary  to point hap.py to the
reference sequence using

```bash
export HGREF=<path-to-hg19.fa>
```

Note that, while the test cases are based on hg19, other reference sequences are
usable as well  once the tool is installed.

Hap.py also requires a copy of the [Boost libraries](http://www.boost.org) to
work, with version >=  1.55. If compilation should fail using the included version
of boost, you can compile a subset of boost like this:

```bash
cd ~
wget http://downloads.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.bz2
tar xjf boost_1_55_0.tar.bz2
cd boost_1_55_0
./bootstrap.sh --with-libraries=filesystem,chrono,thread,iostreams,system,regex,test,program_options
./b2 --prefix=$HOME/boost_1_55_0_install install
```

You can point Cmake to your version of boost as follows:

```bash
export BOOST_ROOT=$HOME/boost_1_55_0_install
```

The complete list of dependencies / packages to install beforehand can be found
in the [Dockerfile](Dockerfile).

## Python 3 Migration

This project is undergoing a migration to Python 3. Key goals include:

* Full Python 3.7+ compatibility.
* Modernized build system using `pyproject.toml` (PEP 517/518).
* Improved packaging and installation via `pip`.
* Adoption of modern Python development practices (type hinting, linting, automated testing).

For more details, see:

* [Migration Status](PYTHON3_MIGRATION_FINAL.md) - Current state and test plan.
* [Migration Tools](PYTHON3_MIGRATION_TOOLS.md) - Tools for fixing remaining issues (if applicable, link may be outdated).
* [Core Documentation](PYTHON3_CORE.md) - Technical details of the Python 3 implementation (if applicable, link may be outdated).

## Key Features

* Python 3.7+ compatibility
* Focus on vcfeval as the primary comparison engine
* Stratified performance metrics using BED files
* Improved build system and dependency management
* Better string handling and error reporting

## Quick Start

After installation, the command-line tools will be available through entry points:

```bash
# Show help message
hap -h

# Example: Compare a VCF file against a truth VCF
hap truth.vcf.gz query.vcf.gz -r reference.fa -o output_prefix

# Run preprocessing on a VCF file
pre input.vcf -o output.vcf -r reference.fa

# Run somatic comparison
som truth.vcf.gz query.vcf.gz -r reference.fa -o output_prefix
```

Other tools (`qfy`, `ftx`, etc.) are also available as entry points after installation.
