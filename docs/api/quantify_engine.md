# QuantifyEngine API Documentation

## Overview

The `QuantifyEngine` class provides comprehensive variant analysis capabilities with enhanced ROC (Receiver Operating Characteristic) analysis features implemented in Phase 2. This class replaces the original C++ quantify component with a pure Python implementation using pandas and pysam.

## Class: QuantifyEngine

### Location
```python
from hap_py.haplo.python_quantify import QuantifyEngine
```

### Constructor

```python
def __init__(
    self,
    truth_vcf: str,
    query_vcf: str,
    reference: Optional[str] = None,
    regions: Optional[str] = None,
    apply_filters: bool = False,
    output_vtc: bool = False,
    quantify_method: str = "xcmp",
    enable_roc_analysis: bool = True,
    roc_bootstrap_samples: int = 1000,
    quality_stratification: bool = True,
)
```

#### Parameters

- **truth_vcf** (*str*): Path to truth VCF file containing known variants
- **query_vcf** (*str*): Path to query/test VCF file to be evaluated
- **reference** (*Optional[str]*): Path to reference FASTA file (optional)
- **regions** (*Optional[str]*): BED file with regions to quantify (optional). The engine indexes regions per chromosome for efficient queries.
- **apply_filters** (*bool*): Whether to apply filters from VCF FILTER field
- **output_vtc** (*bool*): Whether to output variant truth categories
- **quantify_method** (*str*): Quantification method - 'xcmp' or 'ga4gh' (default: 'xcmp')
- **enable_roc_analysis** (*bool*): Enable Phase 2 ROC analysis with confidence intervals (default: True)
- **roc_bootstrap_samples** (*int*): Number of bootstrap samples for confidence intervals (default: 1000)
- **quality_stratification** (*bool*): Enable quality score-based stratification (default: True)

#### Raises

- **ValueError**: If `quantify_method` is not 'xcmp' or 'ga4gh'

### Core Analysis Methods

#### quantify()

```python
def quantify(self) -> Dict[str, Any]:
```

Execute the complete variant quantification workflow including Phase 2 ROC analysis.

**Returns:**
- **Dict[str, Any]**: Dictionary containing quantification results including:
  - Basic metrics (TP, FP, FN, precision, recall, F1)
  - ROC analysis data (if enabled)
  - Quality stratification metrics (if enabled)
  - Multi-threshold analysis results

### Phase 2 ROC Analysis Methods

#### _perform_roc_analysis()

```python
def _perform_roc_analysis(self) -> None:
```

Orchestrates the complete ROC analysis workflow for different variant types.

**Workflow:**
1. Generate ROC curves for SNPs, INDELs, and all variants
2. Calculate bootstrap confidence intervals using Jeffreys method
3. Perform quality score stratification
4. Analyze standard quality thresholds

**Results Storage:**
- `self.roc_data`: ROC curve data for each variant type
- `self.bootstrap_confidence_intervals`: Confidence intervals for precision/recall
- `self.quality_metrics`: Quality bin performance metrics

**Dependencies:**
- Requires SciPy for confidence interval calculations
- Requires scikit-learn for AUC calculations (optional)

#### _generate_roc_curve()

```python
def _generate_roc_curve(self, variant_type: str, bvt_values: list) -> None:
```

Generates ROC curve data for a specific variant type by calculating precision and recall at different quality thresholds.

**Parameters:**
- **variant_type** (*str*): Type of variant for ROC curve ('snp', 'indel', 'all')
- **bvt_values** (*list*): List of BVT (Biallelic Variant Type) values to include

**Algorithm:**
1. Filter variants by specified BVT values
2. Sort variants by quality score (descending)
3. Calculate TP, FP, FN at each quality threshold
4. Compute precision and recall at each threshold
5. Calculate AUC if scikit-learn is available

**Output Structure:**
```python
self.roc_data[variant_type] = {
    "thresholds": [float, ...],      # Quality thresholds
    "tp": [int, ...],                # True positives at each threshold
    "fp": [int, ...],                # False positives at each threshold
    "fn": [int, ...],                # False negatives at each threshold
    "precision": [float, ...],       # Precision at each threshold
    "recall": [float, ...],          # Recall at each threshold
    "auc": float                     # Area under curve (if sklearn available)
}
```

#### _calculate_bootstrap_confidence_intervals()

```python
def _calculate_bootstrap_confidence_intervals(self) -> None:
```

Calculates bootstrap confidence intervals for precision and recall using the Jeffreys confidence interval method.

**Statistical Method:**
- Uses `jeffreysCI` from `hap_py.tools.ci`
- Provides asymptotically correct confidence intervals
- Handles edge cases (zero denominators) gracefully

**Requirements:**
- SciPy must be available for statistical calculations
- Logs warning if SciPy is not available

**Output Structure:**
```python
self.bootstrap_confidence_intervals[variant_type] = {
    "precision_ci": [{"lower": float, "upper": float}, ...],
    "recall_ci": [{"lower": float, "upper": float}, ...],
    "auc": float  # AUC value (if available)
}
```

#### _perform_quality_stratification()

```python
def _perform_quality_stratification(self) -> None:
```

Stratifies variants by quality score into predefined bins and calculates performance metrics for each bin.

**Quality Bins:**
- Q1-10: Quality scores 1-10
- Q10-20: Quality scores 10-20
- Q20-30: Quality scores 20-30
- Q30-40: Quality scores 30-40
- Q40+: Quality scores 40 and above

**Metrics Calculated:**
- True Positives (TP)
- False Positives (FP)
- False Negatives (FN)
- Precision
- Recall
- F1 Score
- Variant count per bin

**Output Structure:**
```python
self.quality_metrics = {
    "bin_metrics": {
        "Q1-10": {
            "quality_range": "1-10",
            "TP": int, "FP": int, "FN": int,
            "PRECISION": float, "RECALL": float, "F1": float,
            "variant_count": int
        },
        # ... other bins
    },
    "total_variants": int
}
```

#### _perform_multi_threshold_analysis()

```python
def _perform_multi_threshold_analysis(self) -> None:
```

Analyzes variant calling performance at standard quality thresholds commonly used in genomics.

**Standard Thresholds:**
- Q10: Quality ≥ 10
- Q20: Quality ≥ 20
- Q30: Quality ≥ 30
- Q40: Quality ≥ 40
- Q50: Quality ≥ 50

**Variant Types Analyzed:**
- SNPs: Single nucleotide polymorphisms
- INDELs: Insertions and deletions
- All: Combined analysis

**Output Structure:**
```python
self.roc_data["multi_threshold"] = {
    "snp": {
        "Q10": {
            "threshold": 10,
            "TP": int, "FP": int, "FN": int,
            "PRECISION": float, "RECALL": float, "F1": float,
            "variants_above_threshold": int
        },
        # ... other thresholds
    },
    # ... other variant types
}
```

#### _write_roc_results()

```python
def _write_roc_results(self, output_prefix: str) -> None:
```

Writes comprehensive ROC analysis results to multiple output files.

**Parameters:**
- **output_prefix** (*str*): Prefix for all output files

**Generated Files:**

1. **`{prefix}.roc.tsv`**: ROC curve data with confidence intervals
   - Columns: Type, Threshold, TP, FP, FN, Precision, Recall, Precision_Lower, Precision_Upper, Recall_Lower, Recall_Upper

2. **`{prefix}.quality_stratification.tsv`**: Quality bin performance metrics
   - Columns: Quality_Bin, Quality_Range, TP, FP, FN, Precision, Recall, F1, Variant_Count

3. **`{prefix}.multi_threshold.tsv`**: Standard threshold analysis
   - Columns: Type, Threshold_Name, Threshold, TP, FP, FN, Precision, Recall, F1, Variants_Above_Threshold

4. **`{prefix}.roc_plot.png`**: Precision-recall curves visualization (if matplotlib available)
   - Separate curves for SNPs, INDELs, and all variants
   - Includes AUC values in legend

**Dependencies:**
- matplotlib (optional): For plot generation
- Results are written even if matplotlib is not available

### Data Structures

#### ROC Data Structure
```python
roc_data = {
    "snp": {
        "thresholds": List[float],
        "tp": List[int],
        "fp": List[int],
        "fn": List[int],
        "precision": List[float],
        "recall": List[float],
        "auc": float
    },
    "indel": { ... },  # Same structure
    "all": { ... },    # Same structure
    "multi_threshold": {
        "snp": {
            "Q10": {"threshold": 10, "TP": int, ...},
            "Q20": {"threshold": 20, "TP": int, ...},
            # ...
        },
        # ... other variant types
    }
}
```

#### Bootstrap Confidence Intervals
```python
bootstrap_confidence_intervals = {
    "snp": {
        "precision_ci": [{"lower": float, "upper": float}, ...],
        "recall_ci": [{"lower": float, "upper": float}, ...],
        "auc": float
    },
    # ... other variant types
}
```

#### Quality Metrics
```python
quality_metrics = {
    "bin_metrics": {
        "Q1-10": {
            "quality_range": str,
            "TP": int, "FP": int, "FN": int,
            "PRECISION": float, "RECALL": float, "F1": float,
            "variant_count": int
        },
        # ... other bins
    },
    "total_variants": int
}
```

## Dependencies

### Required
- **pandas**: Data manipulation and analysis
- **pysam**: VCF file handling
- **numpy**: Numerical computations

### Optional (with graceful degradation)
- **scipy**: Statistical functions for confidence intervals
- **scikit-learn**: AUC calculations
- **matplotlib**: ROC curve visualization

### Internal Dependencies
- **hap_py.tools.ci**: Confidence interval calculations (jeffreysCI)

## Usage Examples

### Basic ROC Analysis
```python
from hap_py.haplo.python_quantify import QuantifyEngine

# Initialize engine with ROC analysis enabled
engine = QuantifyEngine(
    truth_vcf="truth.vcf",
    query_vcf="query.vcf",
    enable_roc_analysis=True,
    quality_stratification=True
)

# Run complete analysis
results = engine.quantify()

# Access ROC data
roc_curves = engine.roc_data
confidence_intervals = engine.bootstrap_confidence_intervals
quality_metrics = engine.quality_metrics
```

### Customized ROC Configuration
```python
# Configure ROC analysis parameters
engine = QuantifyEngine(
    truth_vcf="truth.vcf",
    query_vcf="query.vcf",
    enable_roc_analysis=True,
    roc_bootstrap_samples=2000,  # More bootstrap samples
    quality_stratification=False  # Skip quality stratification
)

results = engine.quantify()
```

### Writing Results
```python
# Generate all ROC output files
engine._write_roc_results("analysis_output")

# This creates:
# - analysis_output.roc.tsv
# - analysis_output.quality_stratification.tsv
# - analysis_output.multi_threshold.tsv
# - analysis_output.roc_plot.png (if matplotlib available)
```

## Performance Considerations

### Memory Usage
- ROC analysis stores data points for all quality thresholds
- Memory usage scales with number of variants and unique quality values
- Quality stratification reduces memory footprint by binning

### Computational Complexity
- ROC curve generation: O(n log n) due to sorting by quality
- Bootstrap confidence intervals: O(n) per threshold
- Quality stratification: O(n) for binning variants
- Multi-threshold analysis: O(n) per threshold

### Optimization Recommendations
- For large datasets (>1M variants), consider:
  - Disabling quality stratification if not needed
  - Reducing bootstrap samples for faster computation
  - Using regional analysis to reduce dataset size

## Error Handling

### Common Issues
1. **Missing Dependencies**: Graceful degradation when optional packages unavailable
2. **Empty Variant Sets**: Handles cases with no variants of specific types
3. **Zero Denominators**: Safe handling in precision/recall calculations
4. **File I/O Errors**: Proper error reporting for output file creation

### Logging
- Comprehensive logging at INFO level for workflow progress
- DEBUG level for detailed metrics
- WARNING level for missing dependencies or edge cases

## Version History

### Phase 2 (Current)
- Enhanced ROC analysis with confidence intervals
- Quality score stratification
- Multi-threshold analysis
- Comprehensive output formats
- Statistical rigor with bootstrap methods

### Phase 1
- Basic variant matching and classification
- Core precision/recall metrics
- Foundation for ROC analysis
