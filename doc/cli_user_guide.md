# User Guide: Modernized hap.py CLI Tools

This guide covers the command-line interface for the modernized hap.py tools.

## Table of Contents

1. [Overview](#overview)
2. [Command Line Tools](#command-line-tools)
3. [hap.py - Main Comparison Tool](#happy---main-comparison-tool)
4. [som.py - Somatic Comparison](#sompy---somatic-comparison)
5. [pre.py - Preprocessing](#prepy---preprocessing)
6. [qfy.py - Quantification](#qfypy---quantification)
7. [Common Options](#common-options)
8. [Output Files](#output-files)
9. [Examples](#examples)
10. [Best Practices](#best-practices)

## Overview

The modernized hap.py provides four main command-line tools:

* **hap.py** - Main diploid variant benchmarking tool
* **som.py** - Somatic variant comparison (ignores genotype information)
* **pre.py** - Variant preprocessing and normalization
* **qfy.py** - Quantification-only analysis (requires GA4GH intermediate VCF)

All tools are available as entry points after pip installation and follow consistent command-line patterns.

## Command Line Tools

### Installation Verification

After installation, verify all tools are available:

```bash
# Check tool availability
hap.py --help
som.py --help
pre.py --help
qfy.py --help

# Check version information
hap.py --version
```

## hap.py - Main Comparison Tool

The primary tool for diploid variant benchmarking with haplotype-aware comparison.

### Basic Syntax

```bash
hap.py [truth.vcf] [query.vcf] [options]
```

### Required Arguments

* **truth.vcf** - Ground truth VCF file (reference standard)
* **query.vcf** - Query VCF file to be evaluated
* **-r/--reference** - Reference genome FASTA file
* **-o/--output** - Output prefix for result files

### Essential Options

```bash
# Basic comparison
hap.py truth.vcf query.vcf -r reference.fa -o benchmark_results

# With confident regions
hap.py truth.vcf query.vcf -r reference.fa -o results -f confident.bed

# Specify engine
hap.py truth.vcf query.vcf -r reference.fa -o results --engine=vcfeval

# Limit to specific chromosomes
hap.py truth.vcf query.vcf -r reference.fa -o results -l chr1,chr2,chr3
```

### Advanced Options

**Engine Selection:**
```bash
--engine=xcmp        # Default hap.py engine (graph-based)
--engine=vcfeval     # RTG vcfeval engine (more sophisticated)
--engine=scmp-somatic # Somatic comparison engine
```

**Output Control:**
```bash
--no-json           # Disable JSON metrics output
--no-csv            # Disable CSV summary output
--verbose           # Enable verbose logging
--logfile=log.txt   # Write logs to file
```

**Analysis Options:**
```bash
--stratification=regions.bed  # Stratify results by genomic regions
--roc-filter=QUAL            # Generate ROC curves using QUAL field
--threads=8                  # Set number of parallel threads
```

**GA4GH Compliance:**
```bash
--quantify-method=ga4gh             # Enable GA4GH-compliant output
--ga4gh-stratification=regions.bed  # GA4GH stratification regions
```

### Example Commands

**Basic whole-genome comparison:**
```bash
hap.py PlatinumGenomes_truth.vcf.gz GATK_calls.vcf.gz \
    -r GRCh38.fa \
    -f confident_regions.bed.gz \
    -o GATK_vs_PG \
    --threads=8
```

**Exome comparison with RTG engine:**
```bash
hap.py truth_exome.vcf.gz query_exome.vcf.gz \
    -r GRCh38.fa \
    -f exome_confident.bed \
    -o exome_benchmark \
    --engine=vcfeval \
    --threads=4
```

**ROC analysis with quality stratification:**
```bash
hap.py truth.vcf.gz query.vcf.gz \
    -r reference.fa \
    -o roc_analysis \
    --roc-filter=QUAL \
    --quantify-method=ga4gh
```

## som.py - Somatic Comparison

> **⚠️ COMPONENT NOT AVAILABLE**: The som.py tool was not included in this modernized Python 3 version of hap.py. Please see the [removed components](removed_components.md) documentation for more information and alternatives.

Tool for somatic variant comparison that ignores genotype information (unavailable in this version).

### Original Functionality

The original som.py tool provided:

* Allele-based comparison (ignoring genotype)
* Suitable for tumor-normal comparisons
* Simple presence/absence evaluation

### Example Usage (For Reference Only)

These examples are provided for reference but **will not work** in this version:

```bash
# UNAVAILABLE IN THIS VERSION
som.py somatic_truth.vcf tumor_calls.vcf \
    -r reference.fa \
    -o somatic_benchmark \
    -f callable_regions.bed
```

If you need som.py functionality, please see [removed components](removed_components.md) for alternatives or to request its restoration.

## pre.py - Preprocessing

Variant preprocessing and normalization tool.

### Basic Syntax

```bash
pre.py [input.vcf] [options]
```

### Preprocessing Features

* **Left-shifting** - Move indels to leftmost position
* **Trimming** - Remove redundant bases from alleles
* **Decomposition** - Split multi-allelic sites
* **Normalization** - Standardize variant representation

### Common Options

```bash
-r/--reference     # Reference genome (required)
-o/--output        # Output VCF file
--window-size      # Window size for complex variant processing
--left-shift       # Enable left-shifting (default: true)
--trim-alts        # Trim ALT alleles (default: true)
--split-complex    # Split complex variants
```

### Example Commands

**Basic normalization:**
```bash
pre.py raw_variants.vcf \
    -r reference.fa \
    -o normalized_variants.vcf \
    --left-shift --trim-alts
```

**Complex variant decomposition:**
```bash
pre.py complex_variants.vcf \
    -r reference.fa \
    -o decomposed_variants.vcf \
    --split-complex \
    --window-size=1000
```

## qfy.py - Quantification

Quantification-only tool for analyzing pre-processed GA4GH intermediate VCF files.

### Basic Syntax

```bash
qfy.py [ga4gh_intermediate.vcf] [options]
```

### Use Cases

* Reanalyze existing benchmarking results
* Apply different stratification regions
* Generate additional metrics
* Custom quantification workflows

### Example Commands

**Basic quantification:**
```bash
qfy.py intermediate.vcf \
    -o quantification_results \
    --reference=reference.fa
```

**With stratification:**
```bash
qfy.py intermediate.vcf \
    -o stratified_results \
    --stratification=genomic_regions.bed \
    --reference=reference.fa
```

## Common Options

### Input/Output Options

```bash
-r/--reference FASTA    # Reference genome file (required)
-o/--output PREFIX      # Output file prefix
-f/--false-regions BED  # Confident/callable regions
-T/--target-regions BED # Restrict analysis to target regions
```

### Analysis Control

```bash
-l/--location REGIONS   # Chromosomes/regions to analyze (chr1,chr2)
--threads NUMBER        # Parallel processing threads
--engine ENGINE         # Comparison engine (xcmp, vcfeval)
--verbose              # Verbose logging
--logfile FILE         # Log file path
```

### Quality Control

```bash
--roc-filter FIELD     # ROC analysis using specified field
--no-haprec            # Disable haplotype-based record splitting
--force-interactive    # Enable interactive mode for debugging
```

## Output Files

### Standard Output Files

All hap.py runs produce standard output files:

**Summary Files:**
* `{prefix}.summary.csv` - High-level precision/recall metrics
* `{prefix}.metrics.json` - Detailed metrics in JSON format

**VCF Files:**
* `{prefix}.vcf.gz` - Annotated VCF with benchmarking decisions
* `{prefix}.vcf.gz.tbi` - Tabix index for VCF file

**Region Files:**
* `{prefix}.false-positives.bed` - False positive regions
* `{prefix}.false-negatives.bed` - False negative regions

### Enhanced Output Files

**ROC Analysis (when enabled):**
* `{prefix}.roc.tsv` - ROC curve data with confidence intervals
* `{prefix}.quality_stratification.tsv` - Quality score bin analysis
* `{prefix}.multi_threshold.tsv` - Standard threshold analysis

**GA4GH Compliance (when enabled):**
* `{prefix}.ga4gh.vcf` - GA4GH-compliant annotated VCF
* `{prefix}.ga4gh.json` - GA4GH metrics in JSON format
* `{prefix}.ga4gh.tsv` - GA4GH metrics in TSV format

**Stratification (when enabled):**
* `{prefix}.stratified.tsv` - Results by genomic region
* `{prefix}.regions.bed` - Analyzed regions

## Examples

### Example 1: Basic Whole Genome Analysis

```bash
# Download example data (if needed)
wget https://example.com/NA12878_truth.vcf.gz
wget https://example.com/NA12878_query.vcf.gz
wget https://example.com/GRCh38.fa

# Run analysis
hap.py NA12878_truth.vcf.gz NA12878_query.vcf.gz \
    -r GRCh38.fa \
    -o NA12878_benchmark \
    --threads=8 \
    --verbose

# View results
cat NA12878_benchmark.summary.csv
```

### Example 2: Exome Analysis with Preprocessing

```bash
# Preprocess variants
pre.py raw_exome_calls.vcf \
    -r GRCh38.fa \
    -o normalized_exome_calls.vcf \
    --left-shift --trim-alts

# Run benchmarking
hap.py exome_truth.vcf.gz normalized_exome_calls.vcf \
    -r GRCh38.fa \
    -f exome_targets.bed \
    -o exome_benchmark \
    --engine=vcfeval
```

### Example 3: ROC Analysis with GA4GH Compliance

```bash
# Run with ROC analysis and GA4GH output
hap.py truth.vcf.gz query.vcf.gz \
    -r reference.fa \
    -o comprehensive_analysis \
    --roc-filter=QUAL \
    --quantify-method=ga4gh \
    --ga4gh-stratification=stratification_regions.bed \
    --verbose

# Results include ROC curves and GA4GH-compliant output
ls comprehensive_analysis.*
```

### Example 4: Multi-Sample Analysis

```bash
# Analyze multiple samples
for sample in sample1 sample2 sample3; do
    hap.py truth_${sample}.vcf.gz query_${sample}.vcf.gz \
        -r reference.fa \
        -o ${sample}_benchmark \
        --threads=4 &
done
wait

# Combine results (custom script)
python combine_results.py *_benchmark.summary.csv > combined_results.csv
```

## Best Practices

### Input Preparation

1. **Index all VCF files:**
   ```bash
   bgzip input.vcf
   tabix -p vcf input.vcf.gz
   ```

2. **Validate VCF files:**
   ```bash
   bcftools view -h input.vcf.gz | head -20  # Check header
   bcftools stats input.vcf.gz               # Generate statistics
   ```

3. **Prepare reference genome:**
   ```bash
   samtools faidx reference.fa  # Create .fai index
   ```

### Performance Optimization

1. **Use appropriate threads:**
   ```bash
   --threads=$(nproc)  # Use all available cores
   ```

2. **Limit analysis regions for testing:**
   ```bash
   -l chr22  # Test on smaller chromosome first
   ```

3. **Use vcfeval for complex regions:**
   ```bash
   --engine=vcfeval  # More sophisticated for difficult regions
   ```

### Quality Control

1. **Always use confident regions:**
   ```bash
   -f confident_callable_regions.bed
   ```

2. **Enable verbose logging:**
   ```bash
   --verbose --logfile analysis.log
   ```

3. **Validate results:**
   ```bash
   # Check summary metrics are reasonable
   cat results.summary.csv
   
   # Verify output files exist
   ls results.*
   ```

### Troubleshooting

1. **Memory issues:**
   ```bash
   # Use smaller regions
   -l chr22
   
   # Reduce threads
   --threads=2
   ```

2. **RTG engine issues:**
   ```bash
   # Check RTG availability
   hap.py --list-engines
   
   # Use alternative engine
   --engine=xcmp
   ```

3. **VCF format issues:**
   ```bash
   # Preprocess problematic VCFs
   pre.py input.vcf -r reference.fa -o clean.vcf
   ```

This guide provides comprehensive coverage of the modernized hap.py CLI tools. For additional information, see the tool-specific help (`--help`) and the detailed documentation in the `doc/` directory.
