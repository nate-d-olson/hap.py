# hap.py - Comprehensive Usage Guide

## Table of Contents
- [Basic Usage](#basic-usage)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Common Use Cases](#common-use-cases)
- [Advanced Options](#advanced-options)
- [Performance Tuning](#performance-tuning)
- [Troubleshooting](#troubleshooting)

## Basic Usage

### Basic Command Structure

```bash
hap.py <truth_vcf> <query_vcf> -f <confidence_regions> -r <reference_fasta> -o <output_prefix>
```

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `truth_vcf` | Path to truth VCF/BCF file | `truth.vcf.gz` |
| `query_vcf` | Path to query VCF/BCF file | `query.vcf.gz` |
| `-f, --confidence-regions` | BED file with confident regions | `confident_regions.bed` |
| `-r, --reference` | Reference FASTA file | `hg38.fa` |
| `-o, --output-prefix` | Output file prefix | `output/result` |

### Example: Basic Comparison

```bash
hap.py \
  truth.vcf.gz \
  query.vcf.gz \
  -f confident_regions.bed \
  -r hg38.fa \
  -o output/result
```

## Input Files

### VCF/BCF Files

- Must be sorted and indexed
- Can be gzipped (.vcf.gz, .bcf) or uncompressed (.vcf)
- Should follow VCF 4.2+ specification

### Confidence Regions BED File

- Defines regions where variants are considered confident
- Should be sorted and indexed
- Can be generated using tools like GATK's `CallableLoci`

### Reference Genome

- Must be indexed (`samtools faidx`)
- Should match the reference used for variant calling
- Can be compressed with bgzip

## Output Files

| File | Description |
|------|-------------|
| `{prefix}.summary.csv` | Summary metrics in CSV format |
| `{prefix}.extended.csv` | Extended metrics and counts |
| `{prefix}.metrics.json` | Detailed metrics in JSON format |
| `{prefix}.roc.*.csv` | ROC curve data for different variant types |
| `{prefix}.vcf.gz` | Annotated VCF with comparison results |

## Common Use Cases

### Germline Variant Comparison

```bash
hap.py \
  truth.vcf.gz \
  query.vcf.gz \
  -f confident_regions.bed \
  -r hg38.fa \
  -o germline_comparison/result \
  --threads 8
```

### Somatic Variant Comparison

```bash
som.py \
  truth.vcf.gz \
  query.vcf.gz \
  -f confident_regions.bed \
  -r hg38.fa \
  -o somatic_comparison/result \
  --threads 8
```

### Stratified Analysis

```bash
hap.py \
  truth.vcf.gz \
  query.vcf.gz \
  -f confident_regions.bed \
  -r hg38.fa \
  -o stratified/result \
  --stratification stratification.tsv \
  --threads 8
```

## Advanced Options

### Performance Optimization

```bash
hap.py \
  truth.vcf.gz \
  query.vcf.gz \
  -f confident_regions.bed \
  -r hg38.fa \
  -o output/result \
  --threads 16 \
  --engine vcfeval \
  --no-json
```

### Region-based Analysis

```bash
hap.py \
  truth.vcf.gz \
  query.vcf.gz \
  -f confident_regions.bed \
  -r hg38.fa \
  -o output/result \
  -T target_regions.bed \
  --threads 8
```

### Customizing Comparison Parameters

```bash
hap.py \
  truth.vcf.gz \
  query.vcf.gz \
  -f confident_regions.bed \
  -r hg38.fa \
  -o output/result \
  --min-gq 20 \
  --min-hom-alt 5 \
  --adjust-conf-regions \
  --threads 8
```

## Performance Tuning

### Memory Usage

- Use `--chunk-size` to control memory usage
- Default: 1000000 (1Mbp)
- Decrease for lower memory usage, increase for better performance

### Threading

- Use `--threads` to specify number of worker threads
- Default: 1
- Optimal: number of CPU cores - 1

### Engine Selection

| Engine | Description | Memory | Speed |
|--------|-------------|--------|-------|
| `xcmp` | Default engine | Medium | Fast |
| `vcfeval` | RTG vcfeval | High | Medium |
| `cmp_ordered` | Simple ordered comparison | Low | Slow |

## Troubleshooting

### Common Issues

1. **Missing Index Files**
   - Ensure all VCF/BCF files are indexed with tabix
   - Ensure reference FASTA is indexed with samtools faidx

2. **Memory Issues**
   - Reduce chunk size with `--chunk-size`
   - Use `--no-json` to reduce memory usage
   - Increase system swap space if needed

3. **Performance Problems**
   - Use `--threads` to enable parallel processing
   - Consider using `--engine vcfeval` for large datasets
   - Ensure input files are on fast storage (SSD recommended)

### Getting Help

For additional help, please:
1. Check the [documentation](doc/)
2. Search the [issue tracker](https://github.com/nate-d-olson/hap.py/issues)
3. Open a new issue if your problem isn't addressed
