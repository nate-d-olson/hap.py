# Quantify Module Implementation Prompt

## Objective

Complete the implementation of the quantify module in the modernized hap.py codebase to restore full functionality for variant calling benchmarking and ROC analysis. This implementation should provide equivalent functionality to the original C++ quantify tool while leveraging modern Python practices and libraries.

## Current Implementation Status

### ✅ Completed Foundation (60% complete)
- **QuantifyEngine class framework** in `python_quantify.py`
- **VCF loading and basic parsing** with pysam integration
- **Basic stratification framework** by variant type, size, zygosity
- **Command-line interface** through `qfy.py`
- **Test infrastructure** with unit and integration test suites
- **Data models** and type definitions
- **Integration points** with vcfeval and metrics calculator

### ❌ Critical Missing Components (40% remaining)

#### 1. Core Variant Matching Algorithms (HIGH PRIORITY)
**File:** `src/hap_py/haplo/python_quantify.py`
**Method:** `QuantifyEngine._match_variants()`
**Current State:** Placeholder implementation only

**Required Implementation:**
```python
def _match_variants(self, truth_variants: List[Variant],
                   query_variants: List[Variant]) -> MatchResult:
    """
    Implement sophisticated variant matching logic equivalent to C++ XCMPQuantify.

    Requirements:
    - Complex variant matching with position and allele comparison
    - Superlocus analysis for complex variants
    - Multi-allelic variant support
    - Benchmarking decision tracking (BD, BVT, QQ fields)
    - Normalization and left-alignment handling
    """
    # IMPLEMENTATION NEEDED
```

**Algorithm Requirements:**
- **Position-based matching:** Variants at same genomic position
- **Allele sequence matching:** Exact or normalized allele comparison
- **Complex variant handling:** MNPs, complex indels, structural variants
- **Superlocus analysis:** Grouping nearby variants for comparison
- **Multi-allelic support:** Handling variants with multiple alternate alleles
- **Decision tracking:** Recording match/mismatch reasons (BD, BVT, QQ)

#### 2. ROC Analysis Implementation (HIGH PRIORITY)
**File:** `src/hap_py/haplo/quantify.py`
**Function:** `_generate_roc_curves()`
**Current State:** Basic structure, no algorithm implementation

**Required Implementation:**
```python
def _generate_roc_curves(results: MatchResult, quality_field: str) -> pd.DataFrame:
    """
    Generate ROC curves from variant matching results.

    Requirements:
    - Quality score threshold analysis
    - TP/FP/FN calculations at each threshold
    - Sensitivity and precision computation
    - Confidence interval calculation
    - Multiple quality field support
    """
    # IMPLEMENTATION NEEDED
```

**Algorithm Requirements:**
- **Quality score processing:** Extract and sort by quality values
- **Threshold analysis:** Calculate metrics at each quality threshold
- **ROC data generation:** Sensitivity vs (1-Precision) curves
- **Confidence intervals:** Statistical confidence bounds
- **Multiple metrics:** Support for QUAL, GQ, and custom quality fields

#### 3. Advanced Stratification (MEDIUM PRIORITY)
**File:** `src/hap_py/haplo/python_quantify.py`
**Method:** `QuantifyEngine._stratify_results()`
**Current State:** Basic type/size stratification only

**Enhancement Requirements:**
- **BED region stratification:** Analyze results by genomic regions
- **Custom stratification:** User-defined stratification categories
- **Feature-based filtering:** Filter by variant annotations
- **Population analysis:** Stratify by allele frequency or population

## Implementation Approach

### Phase 1: Core Algorithm Implementation

#### Step 1: Variant Matching Algorithm
**Priority:** CRITICAL
**Estimated Effort:** 2-3 weeks

**Implementation Strategy:**
1. **Study original C++ algorithm** in `src/c++/lib/quantify/XCMPQuantify.cpp`
2. **Implement basic position/allele matching**
3. **Add complex variant support** (MNPs, complex indels)
4. **Implement superlocus analysis**
5. **Add benchmarking decision tracking**
6. **Validate against C++ output** using test datasets

**Key Algorithms to Implement:**
```python
def _match_single_position(self, truth_vars: List[Variant],
                          query_vars: List[Variant]) -> List[Match]:
    """Match variants at a single genomic position."""
    pass

def _analyze_superlocus(self, variants: List[Variant]) -> SuperlocusAnalysis:
    """Analyze complex variants that span multiple positions."""
    pass

def _normalize_variants(self, variants: List[Variant]) -> List[Variant]:
    """Left-align and normalize variants for comparison."""
    pass
```

#### Step 2: ROC Analysis Implementation
**Priority:** CRITICAL
**Estimated Effort:** 1-2 weeks

**Implementation Strategy:**
1. **Extract quality scores** from VCF INFO/FORMAT fields
2. **Generate quality thresholds** (percentiles or fixed intervals)
3. **Calculate TP/FP/FN** at each threshold
4. **Compute sensitivity/precision** metrics
5. **Generate ROC curves** with confidence intervals
6. **Validate against original outputs**

**Key Functions to Implement:**
```python
def _extract_quality_scores(self, variants: List[Variant],
                           quality_field: str) -> List[float]:
    """Extract quality scores from variants."""
    pass

def _calculate_metrics_at_threshold(self, matches: MatchResult,
                                   threshold: float) -> ROCPoint:
    """Calculate TP/FP/FN at given quality threshold."""
    pass

def _compute_confidence_intervals(self, roc_data: pd.DataFrame) -> pd.DataFrame:
    """Add confidence intervals to ROC data."""
    pass
```

#### Step 3: Enhanced Stratification
**Priority:** MEDIUM
**Estimated Effort:** 1-2 weeks

**Implementation Requirements:**
- BED file parsing for region-based stratification
- Custom stratification category support
- Feature-based filtering capabilities
- Performance optimization for large stratifications

### Phase 2: Performance Optimization

#### Memory Efficiency
- **Streaming VCF processing** for large files
- **Memory-mapped file access** for reference data
- **Efficient data structures** for variant storage
- **Garbage collection optimization**

#### Processing Speed
- **Vectorized operations** using NumPy
- **Parallel processing** with multiprocessing
- **Optimized algorithms** for hotspot functions
- **Cython implementation** for critical paths (if needed)

**Performance Targets:**
- 100K variants: <1 minute
- 1M variants: <10 minutes
- Memory usage: <2GB for 1M variants

### Phase 3: Testing and Validation

#### Unit Testing
- **Algorithm correctness** with known input/output pairs
- **Edge case handling** (empty files, malformed data)
- **Performance regression** testing
- **Memory usage** validation

#### Integration Testing
- **End-to-end workflows** with real datasets
- **Output format compatibility** with downstream tools
- **Accuracy validation** against original C++ implementation

## Technical Specifications

### Input Requirements
- **Truth VCF:** Reference standard variant calls
- **Query VCF:** Test variant calls to evaluate
- **Reference FASTA:** Genome reference for normalization
- **BED files:** Optional regions for stratification
- **Configuration:** Quality fields, thresholds, output options

### Output Requirements
- **Summary TSV:** Overall metrics and stratified results
- **ROC curves:** Sensitivity vs precision data
- **Detailed TSV:** Per-variant match information
- **JSON output:** Machine-readable results (optional)

### Algorithm Specifications

#### Variant Matching Logic
```
For each genomic position:
  1. Group truth and query variants
  2. Normalize all variants (left-align, trim)
  3. Compare alleles for exact matches
  4. Handle multi-allelic cases
  5. Record match decisions (TP/FP/FN)
  6. Track benchmarking fields (BD, BVT, QQ)
```

#### ROC Analysis Logic
```
For each quality threshold:
  1. Filter variants by quality score
  2. Count TP/FP/FN from match results
  3. Calculate sensitivity = TP/(TP+FN)
  4. Calculate precision = TP/(TP+FP)
  5. Store ROC point (sensitivity, 1-precision)
  6. Compute confidence intervals
```

### Data Models

#### Core Classes
```python
@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    quality: Optional[float]
    info: Dict[str, Any]
    samples: List[Dict[str, Any]]

@dataclass
class Match:
    truth_variant: Optional[Variant]
    query_variant: Optional[Variant]
    match_type: MatchType  # TP, FP, FN
    decision: str  # BD, BVT, QQ
    distance: int  # For complex matches

@dataclass
class ROCPoint:
    threshold: float
    sensitivity: float
    precision: float
    tp_count: int
    fp_count: int
    fn_count: int
```

## Success Criteria

### Functional Requirements
- [ ] Variant matching produces identical results to original C++ on test datasets
- [ ] ROC analysis generates accurate sensitivity/precision curves
- [ ] All stratification categories work correctly
- [ ] Output formats match original specifications
- [ ] Integration with vcfeval pipeline works seamlessly

### Performance Requirements
- [ ] Processes 100K variants in <1 minute
- [ ] Handles 1M+ variants without memory issues
- [ ] Performance within 2x of original C++ implementation
- [ ] No memory leaks or excessive garbage collection

### Quality Requirements
- [ ] 100% unit test coverage for core algorithms
- [ ] Integration tests pass with real datasets
- [ ] Code follows project style guidelines
- [ ] Comprehensive documentation and examples

## Implementation Resources

### Reference Materials
- **Original C++ source:** `src/c++/lib/quantify/` and `src/c++/main/quantify.cpp`
- **Original Python wrapper:** `src/python/Haplo/quantify.py`
- **Test datasets:** `example/` directory with known good outputs
- **Documentation:** Original hap.py documentation and papers

### Development Tools
- **Profiling:** cProfile, line_profiler for performance analysis
- **Testing:** pytest with coverage reporting
- **Debugging:** pdb, logging for algorithm development
- **Validation:** diff tools for output comparison

### External Dependencies
- **pysam:** VCF/BCF file handling
- **pandas/numpy:** Data analysis and computation
- **matplotlib:** ROC curve visualization
- **scipy:** Statistical functions for confidence intervals

## Next Steps

1. **Set up development environment** with all dependencies
2. **Study original C++ algorithm** implementation details
3. **Implement variant matching** starting with simple cases
4. **Add comprehensive unit tests** for each algorithm component
5. **Implement ROC analysis** with validation against known outputs
6. **Optimize performance** through profiling and optimization
7. **Complete integration testing** with real-world datasets

This implementation will restore the quantify functionality to the modernized
hap.py codebase while maintaining compatibility with existing workflows and
providing the performance required for production genomics pipelines.
