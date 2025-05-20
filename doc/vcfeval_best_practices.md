# Best Practices for Running hap.py with vcfeval Engine

This document provides guidance on how to effectively use hap.py with the vcfeval execution engine for variant calling benchmarking.

## Overview

hap.py is a benchmarking tool for variant calls that supports multiple comparison engines. The vcfeval engine, developed by Real Time Genomics (RTG), offers advantages in certain scenarios, particularly for complex variant comparisons.

## Basic Usage

```bash
hap.py \
    --engine vcfeval \
    -r <reference.fa> \
    -f <truth_regions.bed> \
    -o <output_prefix> \
    <truth.vcf.gz> \
    <query.vcf.gz>
```

## Key Parameters

| Parameter | Description |
|-----------|-------------|
| `--engine vcfeval` | Specifies vcfeval as the comparison engine (default is xcmp) |
| `-r, --reference` | Reference FASTA file (required) |
| `-f, --false-positives` | BED file with regions where variants are considered |
| `-T, --target-regions` | BED file with regions to restrict analysis to |
| `--threads` | Number of threads to use (vcfeval scales well with multiple threads) |
| `-o, --output-prefix` | Prefix for output files |
| `--stratification` | TSV file defining stratification regions for detailed analysis |
| `--gender` | Sex of the sample (determines how chrX/Y are handled) |

## Stratification Analysis

Stratification allows analysis of variant performance in different genomic contexts:

```bash
hap.py \
    --engine vcfeval \
    -r <reference.fa> \
    --stratification <stratification.tsv> \
    -o <output_prefix> \
    <truth.vcf.gz> \
    <query.vcf.gz>
```

The stratification TSV file defines different regions to analyze separately, such as:

- High/low GC content regions
- Repeat regions
- Segmental duplications
- Specific genes or functional elements

## Best Practices

1. **Always specify the correct reference genome** that was used to create both VCF files.

2. **Use appropriate region BED files**:
   - Truth regions (`-f`) should contain high-confidence regions from your truth set
   - Target regions (`-T`) should contain areas you intend to evaluate

3. **Set the proper gender/sex parameter** when analyzing sex chromosomes.

4. **Use stratification** for more detailed analysis of variant caller performance in different genomic contexts.

5. **Consider vcfeval-specific options** when working with complex variants:
   - `--engine-vcfeval-path`: Path to RTG Tools installation
   - `--engine-vcfeval-template`: Pre-built RTG reference SDF

6. **Adjust resources based on dataset size**:
   - Increase threads for whole genome comparisons
   - Allocate sufficient memory for large VCF files

## Example Pipeline Configuration

When using hap.py within a workflow (like Snakemake), configure your pipeline with parameters similar to:

```yaml
references:
  GRCh38:
    stratifications:
      id: "GRCh38"
      file: "stratifications/GRCh38.tar.gz"

params:
  happy:
    threads: 16
    engine: "vcfeval"
    engine_extra: "--engine-vcfeval-path /path/to/rtg-tools"
```

## Output Files

With output prefix `-o prefix`, hap.py with vcfeval produces:

- `prefix.summary.csv`: Summary statistics (precision/recall)
- `prefix.vcf.gz`: Annotated VCF with TP/FP/FN annotations
- `prefix.extended.csv`: Detailed metrics per variant type
- `prefix.roc.all.csv.gz`: Data for ROC curve generation
- `prefix.runtime_metrics.json`: Performance metrics for the comparison

## Troubleshooting

- **Memory issues**: Reduce threads or partition your analysis by chromosome
- **Long runtime**: Check if your regions files are appropriate and not too expansive
- **Missing variants**: Ensure normalizing settings match between truth and query VCFs

## References

- [hap.py GitHub repository](https://github.com/Illumina/hap.py)
- [RTG vcfeval documentation](https://github.com/RealTimeGenomics/rtg-tools)
