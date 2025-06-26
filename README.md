# Haplotype Comparison Tools (hap.py)

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![PyPI Version](https://img.shields.io/pypi/v/hap_py.svg)](https://pypi.org/project/hap_py/)
[![License](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](LICENSE.txt)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A modern Python implementation of tools for benchmarking variant calls against gold standard truth datasets.

**Key Features**:
- 🐍 Pure Python implementation with optional C++ extensions for performance
- 📊 Comprehensive variant comparison metrics and reporting
- 🔄 Support for complex variant representations and haplotypes
- 🧪 Extensive test suite with high coverage
- 📦 Easy installation via pip
- 🚀 Optimized for modern Python (3.8+)

> **Note**: This is a modernized fork of the original hap.py project, focusing on maintainability and Python 3 support.
> Legacy helpers `hapcmp` and `hapenum` are no longer included in this fork.

## Quick Start

### Installation

```bash
# Install from PyPI
pip install hap_py

# Or install from source
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py
pip install -e .
```

### Basic Usage

To compare a VCF against a gold standard dataset, use the following command to perform genotype-level haplotype comparison:

```bash
# Basic usage
hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa

# For detailed help
hap.py --help
```

### Somatic Variant Comparison

For comparing somatic variant calls, we provide `som.py` which performs position-based comparison without resolving haplotypes. This is particularly useful for somatic variant calling where phasing information is less critical.

```bash
# Compare somatic variant calls
som.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa
```

## Features

### Core Functionality

- **Haplotype-Aware Comparison**: Accurate benchmarking of complex variant representations
- **Somatic Variant Analysis**: Specialized tools for somatic variant comparison
- **Variant Normalization**: Consistent representation of equivalent variants
- **Comprehensive Metrics**: Detailed performance statistics and quality metrics

### Advanced Features

- **Region-based Analysis**: Focus on specific genomic regions using BED files
- **Stratification**: Performance breakdown by variant type and genomic context
- **Parallel Processing**: Efficient handling of whole-genome datasets
- **Extensible Architecture**: Plugin system for custom comparison methods

## Documentation

For detailed documentation, please refer to:

- [User Guide](doc/happy.md) - Comprehensive guide to using hap.py
- [Somatic Variant Analysis](doc/sompy.md) - Special considerations for somatic variants
- [Variant Normalization](doc/normalisation.md) - Details on variant representation
- [Performance Tuning](doc/microbench.md) - Optimization guide for large datasets

Additional documentation, including an architecture overview, is available in the `docs/` directory and can be built locally using **MkDocs**.

## Table of Contents

- [Quick Start](#quick-start)
  - [Installation](#installation)
  - [Basic Usage](#basic-usage)
- [Features](#features)
  - [Core Functionality](#core-functionality)
  - [Advanced Features](#advanced-features)
- [Documentation](#documentation)
- [Motivation](#motivation)
  - [Complex Variant Comparison](#complex-variant-comparison)
  - [Variant Preprocessing](#variant-preprocessing)
  - [Variant Counting](#variant-counting)
  - [Enhanced ROC Analysis](#enhanced-roc-analysis-and-statistical-confidence)
- [System Requirements](#system-requirements)
  - [Hardware](#hardware)
  - [Software Dependencies](#software-dependencies)
- [Contributing](#contributing)
- [License](#license)
- [Citing hap.py](#citing-happy)

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
export RTG=./rtg-core-<ver>/rtg  # path to bundled RTG Tools if available
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

### RTG vcfeval Integration

hap.py relies on the external [RTG vcfeval](https://github.com/RealTimeGenomics/rtg-tools)
binary when using the ``--engine=vcfeval`` option. The executable is located by
checking the ``RTG`` or ``RTGTOOLS_PATH`` environment variables and then falling
back to ``rtg`` on the ``PATH``. If the executable cannot be found, hap.py will
raise an error.

Set ``RTG`` or ``RTGTOOLS_PATH`` to the full path of the ``rtg`` binary if it is
not available globally. If this repository includes a directory like
``rtg-core-<ver>``, you can point the environment variable to the bundled
``rtg`` script:

```bash
export RTG=/path/to/hap.py/rtg-core-<ver>/rtg
```

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

### Enhanced ROC Analysis and Statistical Confidence

**✅ Phase 2 Complete (January 2025)** - The modernized quantify module provides sophisticated ROC (Receiver Operating Characteristic) analysis
with statistical rigor through bootstrap confidence intervals. This enables robust evaluation of
variant caller performance across different quality score thresholds.

**Key ROC Analysis Features:**

* **✅ Statistical Confidence Intervals**: Bootstrap sampling with Jeffreys method for reliable uncertainty estimates
* **✅ Quality Score Stratification**: Performance analysis across quality bins (Q1-10, Q10-20, Q20-30, Q30-40, Q40+)
* **✅ Multi-threshold Analysis**: Standardized evaluation at Q10, Q20, Q30, Q40, Q50 thresholds
* **✅ Variant Type Stratification**: Separate ROC curves for SNPs, INDELs, and combined analysis
* **✅ Precision-Recall Curves**: Comprehensive performance characterization across all quality thresholds

**ROC Analysis Output Files:**

* `.roc.tsv` - ROC curve data with confidence intervals for each variant type
* `.quality_stratification.tsv` - Performance metrics within each quality score bin
* `.multi_threshold.tsv` - Standardized threshold analysis for consistent benchmarking

**Documentation:**
* [Quantify Module Overview](doc/quantify.md) - Comprehensive documentation and configuration
* [ROC Analysis User Guide](doc/roc_analysis_guide.md) - Practical examples and interpretation
* [QuantifyEngine API Reference](doc/api/quantify_engine.md) - Technical API documentation
* [Integration Testing Guide](doc/testing/roc_analysis_integration_tests.md) - Testing framework documentation
* [ROC Analysis Migration Guide](doc/migration/roc_analysis_migration_guide.md) - Upgrade guide for existing users

```bash
# ROC analysis is enabled by default in hap.py
hap.py truth.vcf query.vcf -r reference.fa -o benchmark_results

# Direct quantify usage with ROC analysis
qfy.py truth.vcf query.vcf -o output/prefix -r ref.fa
```

```bash
# Example showing ROC analysis output files
ls benchmark_results.*
# benchmark_results.summary.csv
# benchmark_results.roc.tsv
# benchmark_results.quality_stratification.tsv
# benchmark_results.multi_threshold.tsv
```

**Legacy ROC Output:** We also produce input data for ROC and precision/recall curves compatible
with external plotting tools. An [example](doc/microbench.md) is included.

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

### Prerequisite Tools

The `bgzip` and `tabix` utilities from **htslib** are optional. By default
hap.py uses `pysam` for compression and indexing, but it can fall back to
these command line tools if they are available. When running with
`--engine=vcfeval`, the `rtg` executable from **rtg-tools** must also be
installed.

```bash
# Debian/Ubuntu
sudo apt-get install -y tabix

# Conda
conda install -c bioconda htslib rtg-tools
```



Installing with `pip` builds the Python package and bundled C++ components.
For most users, running `pip install hap.py` or creating the provided
`environment.yml` is sufficient. A compiler is only required when developing the
Cython extensions.

Prebuilt wheels for Linux and macOS are available on PyPI, so installing on a
supported platform does not require a full C++ toolchain.

### Using Conda

Alternatively, you can create a fully configured conda environment using
[mamba](https://github.com/mamba-org/mamba):

```bash
mamba env create -f environment.yml
conda activate hap-py
```

The provided `environment.yml` installs hap.py with the optional C++ extras and
includes `rtg-tools` from the Bioconda channel.

For a development environment with additional tools (e.g., for ROC analysis and testing),
use `environment-dev.yml`:

```bash
micromamba env create -f environment-dev.yml
# or
conda env create -f environment-dev.yml

micromamba activate happy-dev
```

To install with optional dependencies for C++/Cython extensions (recommended for performance) or development tools:

```bash
pip install .[cpp]      # For optional Cython accelerated features
pip install .[dev]      # For development tools (testing, linting)
pip install .[cpp,dev]  # For both
pip install .[rtgtools] # For using the vcfeval engine via RTG Tools
```

For an editable installation that includes both development and Cython extras, you
can run:

```bash
pip install -e .[dev,cpp]
```

Before installing, ensure that build tools and Python headers are available. On
Debian/Ubuntu systems the required packages can be installed with:

```bash
sudo apt-get install -y build-essential python3-dev cmake zlib1g-dev libbz2-dev
```

Before running the tests, consult
[ENVIRONMENT_SETUP.md](ENVIRONMENT_SETUP.md) for the full list of
testing dependencies and activation commands.

### Building from Source (Advanced)

If you need to build from source and `pip install .` does not meet your needs (e.g., you want to customize the build process or are working in an environment without pip):

1. **Prerequisites**:

   * CMake (version 3.10 or newer)
   * Python (version 3.7 or newer, including development headers)
   * Zlib development libraries

2. **Configure and Build**:
   The `pyproject.toml` and CMake setup are designed to be handled by `pip`. For manual control, you can invoke CMake directly, but this is now an advanced use case. The `install.py` script is deprecated.

   For developers, the standard Python build frontends should be used:

   ```bash
   python -m build
   ```

   This will produce a wheel in the `dist/` directory, which can then be installed with `pip install dist/hap.py-*.whl`.

## Quick Start

After installation, the `hap.py` command-line tool will be available.

```bash
hap.py --help # Show help message

# Verify that external dependencies are available
hap.py --check-deps

# Example: Compare a VCF file against a truth VCF
hap.py truth.vcf.gz query.vcf.gz -r reference.fa -o output_prefix
```

(Further examples and detailed usage can be found in the documentation.)



## System requirements

### Hardware

Compiling and testing can be done on a standard desktop system with 8GB of RAM. Whole-genome
comparisons (e.g. comparing a gVCF file against the [Platinum Genomes truth dataset](http://www.illumina.com/platinumgenomes/))
can use up to 64GB of RAM (20GB typical, depending on the input VCF) and about 4-12 minutes
using 40 processor cores. Whole exome comparison (using an exome bed mask and the `-T` switch)
can be carried out on a desktop system.


### Linux

Tested on Ubuntu 18.04+ and CentOS 7 or newer.

If you plan to build the optional C++/Cython extensions yourself, a C++14
compiler such as a recent g++ or Clang is required.

### macOS

Hap.py builds and passes basic tests on macOS 10.15+. Full WGS analyses are not routinely tested on this platform.

### Windows

Hap.py is not regularly tested on Windows. The main dependency that fails compilation is htslib. Given a build
of htslib and pysam, using hap.py on Windows should be possible.

### Other requirements

Hap.py requires a human genome reference sequence which contains at least
chromosomes `1-22`, `X`, `Y`, and `M`. The chromosomes should be named
`chr1`-`chr22`, `chrX`, `chrY`, `chrM`. Point the tests to your reference with

```bash
export HGREF=<path-to-reference.fa>
```

All other dependencies can be installed via `pip` or `conda`. Use the
[environment.yml](environment.yml) file to create a development environment.
Manual Boost builds are rarely needed—only developers working on the Cython
extensions may need to provide a custom Boost installation.

### Required system packages

The integration tests rely on several external tools being available on your
`$PATH`:

- `bcftools`
- `htslib` (provides `bgzip` and `tabix`)
- `bgzip`
- `tabix`
- `rtg-tools`

Example installation commands:

```bash
# Debian/Ubuntu
sudo apt-get install -y bcftools tabix

# Conda
conda install -c bioconda bcftools htslib rtg-tools
```

Ensure these executables are discoverable before running the integration tests.

## Python 3 Migration

This project has been migrated to Python 3. Please see the [Python 3 migration guide](doc/python3_migration.md) for details.

## Citing hap.py

If you use hap.py in your research, please cite the original publication:

```
Krusche, P., Trigg, L., Boutros, P.C. et al.
Best practices for benchmarking germline small-variant calls in human genomes.
Nat Biotechnol 37, 555–560 (2019).
https://doi.org/10.1038/s41587-019-0054-x
```

BibTeX entry:
```bibtex
@article{krusche2019best,
  title={Best practices for benchmarking germline small-variant calls in human genomes},
  author={Krusche, Peter and Trigg, Len and Boutros, Paul C and Mason, Christopher E and De La Vega, Francisco M and Moore, Barry L and Gonzalez-Porta, Mar and Eberle, Michael A and Tezak, Ziv and Lababidi, Samir and others},
  journal={Nature biotechnology},
  volume={37},
  number={5},
  pages={555--560},
  year={2019},
  publisher={Nature Publishing Group US New York}
}
```

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE.txt](LICENSE.txt) file for details.



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
