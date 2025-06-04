# Integration Test Documentation for ROC Analysis

## Overview

This document describes the integration testing approach for the Phase 2 ROC analysis functionality in the hap.py quantify module. The tests validate the end-to-end functionality of enhanced ROC analysis features including confidence intervals, quality stratification, and multi-threshold analysis.

## Test Structure

### Test Files
- **`tests/unit/test_roc_analysis.py`**: Unit tests for individual ROC methods
- **`tests/integration/test_quantify_roc.py`**: Integration tests with real VCF data
- **`test_roc_functionality.py`**: Standalone ROC functionality demonstration

### Test Categories

#### 1. Unit Tests
**Purpose**: Validate individual ROC analysis methods in isolation

**Coverage:**
- ROC analysis initialization
- ROC curve generation algorithms
- Bootstrap confidence interval calculations
- Quality stratification logic
- Multi-threshold analysis
- Result output formatting

#### 2. Integration Tests
**Purpose**: Validate complete ROC analysis workflow with real data

**Coverage:**
- End-to-end ROC analysis with VCF files
- Output file generation and validation
- Compatibility with different VCF formats
- Performance with large datasets
- Error handling and edge cases

#### 3. Functional Tests
**Purpose**: Demonstrate practical usage scenarios

**Coverage:**
- Real-world variant caller evaluation
- Comparison with known benchmarks
- Visualization output validation

## Test Data Requirements

### VCF Test Files

#### Small Test Dataset
- **File**: `tests/data/small_test.vcf`
- **Description**: 100-500 variants with quality scores
- **Purpose**: Fast unit tests and basic functionality validation

#### Medium Test Dataset
- **File**: `tests/data/medium_test.vcf`
- **Description**: 1000-5000 variants across multiple chromosomes
- **Purpose**: Integration testing and performance validation

#### Large Test Dataset
- **File**: `tests/data/large_test.vcf`
- **Description**: 10000+ variants with diverse quality distributions
- **Purpose**: Stress testing and performance benchmarks

### Test Data Characteristics

#### Quality Score Distribution
```
Q1-10:   20% of variants
Q10-20:  25% of variants
Q20-30:  25% of variants
Q30-40:  20% of variants
Q40+:    10% of variants
```

#### Variant Type Distribution
```
SNPs:    70% of variants
INDELs:  25% of variants
Other:   5% of variants
```

#### Truth/Query Matching
```
True Positives:  70-80% of query variants
False Positives: 15-20% of query variants
False Negatives: 5-10% of truth variants
```

## Key Integration Test Cases

### 1. Basic ROC Analysis Workflow

```python
def test_basic_roc_analysis():
    """Test complete ROC analysis workflow."""
    engine = QuantifyEngine(
        truth_vcf="tests/data/truth_small.vcf",
        query_vcf="tests/data/query_small.vcf",
        enable_roc_analysis=True
    )

    # Run analysis
    results = engine.quantify()

    # Validate ROC data structure
    assert "roc_data" in results
    assert "snp" in engine.roc_data
    assert "indel" in engine.roc_data
    assert "all" in engine.roc_data

    # Validate confidence intervals
    assert hasattr(engine, "bootstrap_confidence_intervals")

    # Validate quality metrics
    assert hasattr(engine, "quality_metrics")
```

### 2. Output File Generation

```python
def test_roc_output_files():
    """Test ROC analysis output file generation."""
    engine = QuantifyEngine(
        truth_vcf="tests/data/truth_medium.vcf",
        query_vcf="tests/data/query_medium.vcf"
    )

    results = engine.quantify()

    with tempfile.TemporaryDirectory() as tmpdir:
        output_prefix = os.path.join(tmpdir, "test_output")
        engine._write_roc_results(output_prefix)

        # Validate output files exist
        assert os.path.exists(f"{output_prefix}.roc.tsv")
        assert os.path.exists(f"{output_prefix}.quality_stratification.tsv")
        assert os.path.exists(f"{output_prefix}.multi_threshold.tsv")

        # Validate file contents
        roc_df = pd.read_csv(f"{output_prefix}.roc.tsv", sep="\t")
        assert len(roc_df) > 0
        assert "Type" in roc_df.columns
        assert "Precision" in roc_df.columns
        assert "Recall" in roc_df.columns
```

### 3. Quality Stratification Validation

```python
def test_quality_stratification():
    """Test quality score stratification functionality."""
    engine = QuantifyEngine(
        truth_vcf="tests/data/truth_medium.vcf",
        query_vcf="tests/data/query_medium.vcf",
        quality_stratification=True
    )

    results = engine.quantify()

    # Validate quality metrics structure
    assert hasattr(engine, "quality_metrics")
    assert "bin_metrics" in engine.quality_metrics

    # Validate expected quality bins
    expected_bins = ["Q1-10", "Q10-20", "Q20-30", "Q30-40", "Q40+"]
    for bin_name in expected_bins:
        assert bin_name in engine.quality_metrics["bin_metrics"]

    # Validate metrics for each bin
    for bin_name, metrics in engine.quality_metrics["bin_metrics"].items():
        assert "TP" in metrics
        assert "FP" in metrics
        assert "FN" in metrics
        assert "PRECISION" in metrics
        assert "RECALL" in metrics
        assert "F1" in metrics
```

### 4. Multi-Threshold Analysis

```python
def test_multi_threshold_analysis():
    """Test standard quality threshold analysis."""
    engine = QuantifyEngine(
        truth_vcf="tests/data/truth_medium.vcf",
        query_vcf="tests/data/query_medium.vcf"
    )

    results = engine.quantify()

    # Validate multi-threshold data
    assert "multi_threshold" in engine.roc_data

    # Validate variant types
    expected_types = ["snp", "indel", "all"]
    for variant_type in expected_types:
        assert variant_type in engine.roc_data["multi_threshold"]

    # Validate thresholds
    expected_thresholds = ["Q10", "Q20", "Q30", "Q40", "Q50"]
    for threshold in expected_thresholds:
        assert threshold in engine.roc_data["multi_threshold"]["snp"]
```

### 5. Confidence Interval Validation

```python
def test_confidence_intervals():
    """Test bootstrap confidence interval calculations."""
    # Skip if SciPy not available
    pytest.importorskip("scipy")

    engine = QuantifyEngine(
        truth_vcf="tests/data/truth_small.vcf",
        query_vcf="tests/data/query_small.vcf",
        roc_bootstrap_samples=100  # Smaller for faster testing
    )

    results = engine.quantify()

    # Validate confidence intervals exist
    assert hasattr(engine, "bootstrap_confidence_intervals")

    for variant_type in ["snp", "indel", "all"]:
        if variant_type in engine.bootstrap_confidence_intervals:
            ci_data = engine.bootstrap_confidence_intervals[variant_type]
            assert "precision_ci" in ci_data
            assert "recall_ci" in ci_data

            # Validate CI structure
            for ci_list in [ci_data["precision_ci"], ci_data["recall_ci"]]:
                for ci in ci_list:
                    assert "lower" in ci
                    assert "upper" in ci
                    assert 0.0 <= ci["lower"] <= 1.0
                    assert 0.0 <= ci["upper"] <= 1.0
                    assert ci["lower"] <= ci["upper"]
```

### 6. Large Dataset Performance

```python
@pytest.mark.slow
def test_large_dataset_performance():
    """Test ROC analysis performance with large datasets."""
    import time

    engine = QuantifyEngine(
        truth_vcf="tests/data/truth_large.vcf",
        query_vcf="tests/data/query_large.vcf"
    )

    start_time = time.time()
    results = engine.quantify()
    execution_time = time.time() - start_time

    # Performance assertions
    assert execution_time < 300  # Should complete in under 5 minutes

    # Validate results completeness
    assert len(engine.roc_data["all"]["thresholds"]) > 0
    assert hasattr(engine, "quality_metrics")
```

### 7. Error Handling and Edge Cases

```python
def test_empty_vcf_handling():
    """Test handling of empty or minimal VCF files."""
    engine = QuantifyEngine(
        truth_vcf="tests/data/empty_truth.vcf",
        query_vcf="tests/data/empty_query.vcf"
    )

    # Should not raise exception
    results = engine.quantify()

    # Should create empty but valid data structures
    assert isinstance(engine.roc_data, dict)

def test_missing_quality_scores():
    """Test handling of VCF files without quality scores."""
    engine = QuantifyEngine(
        truth_vcf="tests/data/truth_no_qual.vcf",
        query_vcf="tests/data/query_no_qual.vcf"
    )

    results = engine.quantify()

    # Should handle gracefully with default quality values
    assert "all" in engine.roc_data
    assert len(engine.roc_data["all"]["thresholds"]) >= 1

def test_missing_dependencies():
    """Test graceful degradation when optional dependencies unavailable."""
    # Mock missing scipy
    with patch.dict('sys.modules', {'scipy': None}):
        engine = QuantifyEngine(
            truth_vcf="tests/data/truth_small.vcf",
            query_vcf="tests/data/query_small.vcf"
        )

        results = engine.quantify()

        # Should complete but skip confidence intervals
        assert not hasattr(engine, "bootstrap_confidence_intervals") or \
               len(engine.bootstrap_confidence_intervals) == 0
```

## Test Data Generation

### Creating Test VCF Files

```python
def create_test_vcf(filename, num_variants=1000, quality_distribution=None):
    """Create synthetic VCF file for testing."""
    if quality_distribution is None:
        quality_distribution = {
            "low": (1, 10, 0.2),      # min, max, fraction
            "medium": (10, 30, 0.5),
            "high": (30, 60, 0.3)
        }

    with open(filename, 'w') as f:
        # VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##INFO=<ID=BVT,Number=1,Type=String,Description=\"Variant Type\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        # Generate variants
        for i in range(num_variants):
            chrom = f"chr{random.randint(1, 22)}"
            pos = random.randint(1000, 100000000)
            ref = random.choice(["A", "T", "G", "C"])
            alt = random.choice(["A", "T", "G", "C"])

            # Generate quality based on distribution
            quality = generate_quality_score(quality_distribution)

            # Determine variant type
            variant_type = "SNP" if len(ref) == len(alt) == 1 else "INDEL"

            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{quality}\tPASS\t")
            f.write(f"BVT={variant_type}\tGT\t1/1\n")
```

## Continuous Integration

### GitHub Actions Configuration

```yaml
name: ROC Analysis Integration Tests

on: [push, pull_request]

jobs:
  test-roc-analysis:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, 3.10, 3.11]

    steps:
    - uses: actions/checkout@v3

    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}

    - name: Install dependencies
      run: |
        pip install -e .[dev]

    - name: Run ROC analysis unit tests
      run: |
        pytest tests/unit/test_roc_analysis.py -v

    - name: Run ROC analysis integration tests
      run: |
        pytest tests/integration/test_quantify_roc.py -v

    - name: Run performance tests
      run: |
        pytest tests/integration/test_quantify_roc.py::test_large_dataset_performance -v -m slow
```

## Performance Benchmarks

### Expected Performance Metrics

| Dataset Size | Variants | Expected Time | Memory Usage |
|-------------|----------|---------------|--------------|
| Small       | 100      | < 1 second    | < 50 MB      |
| Medium      | 1,000    | < 10 seconds  | < 100 MB     |
| Large       | 10,000   | < 2 minutes   | < 500 MB     |
| Extra Large | 100,000  | < 10 minutes  | < 2 GB       |

### Performance Test Implementation

```python
@pytest.mark.performance
def test_roc_analysis_performance():
    """Benchmark ROC analysis performance."""
    import psutil
    import time

    # Monitor memory usage
    process = psutil.Process()
    initial_memory = process.memory_info().rss

    # Time execution
    start_time = time.time()

    engine = QuantifyEngine(
        truth_vcf="tests/data/benchmark_truth.vcf",
        query_vcf="tests/data/benchmark_query.vcf"
    )
    results = engine.quantify()

    execution_time = time.time() - start_time
    peak_memory = process.memory_info().rss
    memory_increase = peak_memory - initial_memory

    # Log performance metrics
    logger.info(f"Execution time: {execution_time:.2f} seconds")
    logger.info(f"Memory increase: {memory_increase / 1024 / 1024:.2f} MB")

    # Performance assertions
    assert execution_time < 120  # 2 minutes max
    assert memory_increase < 1024 * 1024 * 1024  # 1 GB max
```

## Test Result Validation

### Output File Validation

```python
def validate_roc_output_file(filepath):
    """Validate structure and content of ROC output files."""
    df = pd.read_csv(filepath, sep="\t")

    # Validate required columns
    required_columns = ["Type", "Threshold", "TP", "FP", "FN", "Precision", "Recall"]
    for col in required_columns:
        assert col in df.columns, f"Missing required column: {col}"

    # Validate data types and ranges
    assert df["TP"].dtype == int
    assert df["FP"].dtype == int
    assert df["FN"].dtype == int
    assert (df["Precision"] >= 0.0).all() and (df["Precision"] <= 1.0).all()
    assert (df["Recall"] >= 0.0).all() and (df["Recall"] <= 1.0).all()

    # Validate variant types
    expected_types = {"SNP", "INDEL", "ALL"}
    assert set(df["Type"].unique()).issubset(expected_types)

    return True
```

### Statistical Validation

```python
def validate_roc_statistics(roc_data):
    """Validate statistical properties of ROC data."""
    for variant_type, data in roc_data.items():
        if variant_type == "multi_threshold":
            continue

        # Validate monotonicity (recall should generally increase as threshold decreases)
        thresholds = data["thresholds"]
        recalls = data["recall"]

        # Check for generally increasing recall trend
        increasing_trend = sum(recalls[i] <= recalls[i+1] for i in range(len(recalls)-1))
        total_transitions = len(recalls) - 1

        if total_transitions > 0:
            # Allow some tolerance for noise in real data
            assert increasing_trend / total_transitions >= 0.8, \
                f"Recall should generally increase for {variant_type}"

        # Validate precision and recall are in valid range
        assert all(0.0 <= p <= 1.0 for p in data["precision"])
        assert all(0.0 <= r <= 1.0 for r in data["recall"])
```

## Documentation Examples

### Basic Usage Example

```python
def example_basic_roc_analysis():
    """Example: Basic ROC analysis workflow."""
    # Initialize quantify engine
    engine = QuantifyEngine(
        truth_vcf="path/to/truth.vcf",
        query_vcf="path/to/query.vcf",
        enable_roc_analysis=True,
        quality_stratification=True
    )

    # Run complete analysis
    results = engine.quantify()

    # Access ROC curve data
    snp_roc = engine.roc_data["snp"]
    print(f"SNP ROC curve has {len(snp_roc['thresholds'])} points")

    # Access quality stratification
    quality_bins = engine.quality_metrics["bin_metrics"]
    for bin_name, metrics in quality_bins.items():
        print(f"{bin_name}: Precision={metrics['PRECISION']:.3f}, "
              f"Recall={metrics['RECALL']:.3f}")

    # Write results to files
    engine._write_roc_results("output/roc_analysis")
```

### Advanced Configuration Example

```python
def example_advanced_roc_configuration():
    """Example: Advanced ROC analysis configuration."""
    # Configure for high-precision analysis
    engine = QuantifyEngine(
        truth_vcf="path/to/truth.vcf",
        query_vcf="path/to/query.vcf",
        enable_roc_analysis=True,
        roc_bootstrap_samples=5000,  # More samples for better CI estimates
        quality_stratification=True,
        quantify_method="ga4gh"  # Use GA4GH standardized method
    )

    # Run analysis
    results = engine.quantify()

    # Generate detailed output with confidence intervals
    engine._write_roc_results("detailed_analysis")

    # Access confidence intervals
    ci_data = engine.bootstrap_confidence_intervals
    for variant_type, intervals in ci_data.items():
        if "precision_ci" in intervals:
            avg_precision_width = sum(
                ci["upper"] - ci["lower"]
                for ci in intervals["precision_ci"]
            ) / len(intervals["precision_ci"])
            print(f"{variant_type} average precision CI width: {avg_precision_width:.4f}")
```

This comprehensive integration test documentation provides the framework for validating the Phase 2 ROC analysis functionality across different scenarios and use cases.
