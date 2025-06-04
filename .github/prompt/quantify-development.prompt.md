The objective is to complete the next task in the development of the `quantify` module for hap.py.

Steps to follow:
1. Review the current implementation status and identify completed, partially implemented, and missing critical components.
2. Determine changes required for the next task in the development plan.
3. Implement the necessary changes to the codebase.
4. Test the changes to ensure they meet the requirements and do not introduce new issues.
5. As appropriate update documentation, specifically `docs/quantify*md`, `.github/instructions/quantify-implementation.md`.
6. Document changes using an informative commit message.
7. Update `.github/prompt/quantify-development.prompt.md` to reflect the current state and next steps in the development process for use in the next development session.

# Quantify Module Development Plan

## Executive Summary

This document outlines the comprehensive development plan for implementing the
quantify module in the modernized hap.py codebase. Based on extensive analysis
of both the current Python implementation and the original C++ codebase,
this plan provides a strategic approach to complete the quantify functionality
while maintaining performance and compatibility.

## Current Implementation Status

### ✅ Completed Components

**Core Infrastructure:**
- `QuantifyEngine` class with basic VCF processing framework (776 lines in `python_quantify.py`)
- Variant loading and filtering capabilities (`_load_variants`, `_process_variant` methods)
- Basic stratification by variant type, indel size, and zygosity (`_stratify_results` method)
- Command-line interface through `qfy.py`
- Integration with `MetricsCalculator` for basic statistics
- Test framework covering both unit and integration scenarios
- ROC curve generation structure in `quantify.py`

**Data Models and Integration:**
- `QuantifyModels` with proper data structures
- Integration with vcfeval for variant comparison
- Basic VCF header validation and processing
- Configuration management through argparse

### 🔄 Partially Implemented

**Variant Processing:**
- Variant matching framework exists but needs core algorithm completion
- Basic variant type classification (SNP, INDEL, COMPLEX)
- Placeholder implementations for benchmarking superlocus handling

**Analysis and Metrics:**
- ROC analysis structure exists but needs proper quality scoring implementation
- Basic counting metrics implemented
- Stratification framework needs enhancement for regions and custom BED files

### ❌ Missing Critical Components

**Core Algorithms:**
- Sophisticated variant matching algorithms (equivalent to C++ XCMPQuantify/GA4GHQuantify)
- Benchmarking decision tracking (BD, BVT, QQ fields)
- Superlocus analysis and complex variant handling
- Multi-allelic variant processing

**Performance and Scalability:**
- Optimization for large VCF processing (>1M variants)
- Memory-efficient streaming for large datasets
- Multi-threaded processing capabilities

**Advanced Features:**
- Custom stratification by BED regions
- Feature-based analysis and filtering
- Advanced ROC curve generation with confidence intervals

## Development Strategy

### Phase 1: Core Algorithm Implementation (Priority: HIGH)

**Objective:** Complete the fundamental variant matching and analysis algorithms

**Tasks:**
1. **Implement Robust Variant Matching**
   - Complete `_match_variants` method in `QuantifyEngine`
   - Add support for complex variant matching (MNPs, complex indels)
   - Implement benchmarking superlocus handling
   - Add multi-allelic variant support

2. **Enhanced ROC Analysis**
   - Implement proper quality score processing
   - Add benchmarking decision tracking (BD, BVT, QQ fields)
   - Complete confidence interval calculations
   - Add support for custom quality thresholds

3. **Stratification Enhancement**
   - Implement BED region-based stratification
   - Add custom stratification categories
   - Enhance variant type classification accuracy
   - Add support for feature-based filtering

**Files to Modify:**
- `src/hap_py/haplo/python_quantify.py` (lines 400-600: variant matching)
- `src/hap_py/haplo/quantify.py` (ROC generation functions)
- `src/hap_py/haplo/quantify_models.py` (data models)

**Estimated Duration:** 3-4 weeks

### Phase 2: Performance Optimization (Priority: MEDIUM)

**Objective:** Optimize performance to handle large genomic datasets efficiently

**Performance Target:** Process 1M+ variants in <10 minutes (matching original C++ performance)

**Implementation Options:**

#### Option A: Pure Python with NumPy/Pandas Optimization
**Pros:**
- Maintains code simplicity and maintainability
- Leverages mature Python ecosystem
- Easier debugging and testing

**Cons:**
- May not achieve C++ performance levels
- Memory usage potentially higher

**Implementation:**
- Use NumPy vectorized operations for variant processing
- Leverage Pandas for efficient data manipulation
- Implement memory-mapped file reading for large VCFs
- Use multiprocessing for parallel variant processing

#### Option B: Hybrid Python/Cython Implementation
**Pros:**
- Near C++ performance for critical paths
- Maintains Python interface
- Gradual optimization possible

**Cons:**
- Additional build complexity
- Debugging more challenging
- Platform-specific compilation

**Implementation:**
- Identify performance bottlenecks through profiling
- Implement critical algorithms in Cython
- Use Cython memoryviews for efficient data access
- Maintain Python fallbacks for compatibility

**Recommendation:** Start with Option A (Pure Python optimization) and migrate to Option B only if performance benchmarks indicate necessity.

**Tasks:**
1. **Performance Profiling**
   - Benchmark current implementation with large datasets
   - Identify bottlenecks using cProfile and line_profiler
   - Compare performance against original C++ implementation

2. **Memory Optimization**
   - Implement streaming VCF processing
   - Use memory-mapped files for large datasets
   - Optimize data structures for memory efficiency

3. **Parallel Processing**
   - Implement multiprocessing for variant analysis
   - Add thread-safe variant matching
   - Optimize I/O operations

**Files to Create/Modify:**
- `src/hap_py/haplo/performance/` (new module for optimized algorithms)
- `src/hap_py/haplo/python_quantify.py` (optimization updates)
- `setup.py` or `pyproject.toml` (Cython build configuration if needed)

**Estimated Duration:** 2-3 weeks

### Phase 3: Advanced Features and Integration (Priority: LOW)

**Objective:** Complete advanced features and ensure seamless integration

**Tasks:**
1. **Advanced Stratification**
   - Custom BED region processing
   - Feature-based variant analysis
   - Population-specific stratification

2. **Enhanced Output Formats**
   - Complete TSV output formatting
   - Add JSON output option
   - Implement summary statistics

3. **Integration Testing**
   - Comprehensive end-to-end testing
   - Performance regression testing
   - Compatibility testing with existing workflows

**Files to Modify:**
- `src/hap_py/haplo/quantify.py` (output formatting)
- `tests/integration/test_integration_quantify.py` (comprehensive tests)
- `example/` (usage examples and documentation)

**Estimated Duration:** 2-3 weeks

## Technical Architecture

### Core Classes and Modules

```
src/hap_py/haplo/
├── quantify.py              # Main quantify interface and orchestration
├── python_quantify.py       # QuantifyEngine core implementation
├── quantify_models.py       # Data models and structures
└── performance/             # Performance-optimized algorithms (Phase 2)
    ├── __init__.py
    ├── variant_matching.py  # Optimized variant matching
    ├── roc_analysis.py      # Optimized ROC calculations
    └── streaming.py         # Memory-efficient processing
```

### Data Flow Architecture

```
VCF Input Files
       ↓
QuantifyEngine._load_variants()
       ↓
QuantifyEngine._match_variants()
       ↓
QuantifyEngine._stratify_results()
       ↓
ROC Analysis & Metrics Calculation
       ↓
Output Generation (TSV/JSON)
```

### Integration Points

1. **vcfeval Integration:** Quantify processes vcfeval output for detailed analysis
2. **MetricsCalculator:** Handles basic statistical calculations
3. **VCF Processing:** Leverages existing VCF parsing infrastructure
4. **Configuration Management:** Uses shared configuration system

## Performance Requirements and Benchmarks

### Target Performance Metrics

| Metric | Target | Current | Original C++ |
|--------|--------|---------|--------------|
| 100K variants | <1 min | TBD | <30s |
| 1M variants | <10 min | TBD | <5 min |
| Memory usage | <2GB | TBD | <1GB |
| Accuracy | 100% match | TBD | Reference |

### Benchmarking Strategy

1. **Create standardized test datasets** with varying sizes (10K, 100K, 1M variants)
2. **Implement performance test suite** to track regression
3. **Compare outputs** with original C++ implementation for accuracy validation
4. **Profile memory usage** under different workloads

## Testing Strategy

### Unit Testing
- **Variant matching algorithms** with known input/output pairs
- **ROC calculation accuracy** with synthetic datasets
- **Stratification logic** with edge cases
- **Performance regression tests** with timing assertions

### Integration Testing
- **End-to-end workflows** matching original hap.py behavior
- **Large dataset processing** for performance validation
- **Output format compatibility** with downstream tools
- **Error handling** for malformed inputs

### Test Data Requirements
- **Small test datasets** for unit testing (included)
- **Medium datasets** for integration testing (10K-100K variants)
- **Large datasets** for performance testing (1M+ variants)
- **Edge case datasets** for robustness testing

## Risk Mitigation

### Technical Risks

1. **Performance Gap vs C++**
   - **Mitigation:** Incremental optimization with profiling
   - **Fallback:** Cython implementation for critical paths

2. **Memory Usage for Large Datasets**
   - **Mitigation:** Streaming processing and memory mapping
   - **Fallback:** Chunked processing with temporary files

3. **Algorithm Complexity**
   - **Mitigation:** Incremental implementation with extensive testing
   - **Fallback:** Reference C++ implementation for validation

### Project Risks

1. **Timeline Overrun**
   - **Mitigation:** Phased implementation with clear deliverables
   - **Fallback:** Reduce scope for non-critical features

2. **Compatibility Issues**
   - **Mitigation:** Extensive testing against original outputs
   - **Fallback:** Migration guide for breaking changes

## Success Criteria

### Phase 1 Success Criteria
- [ ] All unit tests pass with >95% code coverage
- [ ] Core variant matching produces identical results to original C++
- [ ] ROC analysis generates correct curves for test datasets
- [ ] Integration tests demonstrate end-to-end functionality

### Phase 2 Success Criteria
- [ ] Performance within 2x of original C++ implementation
- [ ] Memory usage reasonable for target datasets
- [ ] No performance regression in existing functionality
- [ ] Scalability to 1M+ variant datasets

### Phase 3 Success Criteria
- [ ] Complete feature parity with original quantify
- [ ] All advanced features working correctly
- [ ] Documentation and examples complete
- [ ] Ready for production use

## Implementation Timeline

### Month 1: Core Algorithm Implementation
- **Week 1-2:** Variant matching algorithms
- **Week 3:** ROC analysis implementation
- **Week 4:** Enhanced stratification

### Month 2: Performance Optimization
- **Week 1:** Performance profiling and bottleneck identification
- **Week 2-3:** Python optimization implementation
- **Week 4:** Cython implementation (if needed)

### Month 3: Advanced Features and Polish
- **Week 1-2:** Advanced features implementation
- **Week 3:** Comprehensive testing and validation
- **Week 4:** Documentation and examples

## Conclusion

This development plan provides a strategic approach to completing the quantify module implementation. By focusing on core algorithms first, then optimizing performance, and finally adding advanced features, we ensure a solid foundation while maintaining development momentum.

The phased approach allows for early validation of core functionality and provides flexibility to adjust the implementation strategy based on performance requirements and technical constraints discovered during development.

The success of this implementation will restore the quantify functionality to the modernized hap.py codebase while maintaining the performance and accuracy expectations of the bioinformatics community.

## Development Objective and Process

The objective is to complete the next task in the development of the `quantify` module for hap.py.

Steps to follow:
1. Review the current implementation status and identify completed, partially implemented, and missing critical components.
2. Determine changes required for the next task in the development plan.
3. Implement the necessary changes to the codebase.
4. Test the changes to ensure they meet the requirements and do not introduce new issues.
5. As appropriate update documentation, specifically `docs/quantify*md`, `.github/instructions/quantify-implementation.md`.
6. Document changes using an informative commit message.
7. Update `.github/prompt/quantify-development.prompt.md` to reflect the current state and next steps in the development process for use in the next development session.
