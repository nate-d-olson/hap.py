# Implementation Guidelines for hap.py Development

## Phase Implementation Process

Each implementation phase should follow these guidelines:

1. **Analysis & Planning**
   - Review existing code and identify gaps
   - Create detailed implementation roadmap
   - Document requirements and acceptance criteria

2. **Implementation**
   - Follow Python 3 best practices and type hints
   - Add comprehensive tests for new functionality
   - Document code with Google-style docstrings

3. **Validation**
   - Create validation scripts to verify implementation
   - Compare results with original implementation
   - Document any differences or improvements

4. **Documentation**
   - Update module docstrings
   - Add examples to demonstration scripts
   - Update relevant sections of README.md

## Phase-Specific Guidelines

### Phase 1: Core Variant Matching

- Implement `_match_variants` method to handle both pandas Series and dictionary inputs
- Ensure consistent handling of allele compatibility
- Add robust tests for edge cases in variant classification

### Phase 2: ROC Analysis

- Implement methods for confidence interval calculation
- Add quality score stratification capabilities
- Generate output files consistent with original formats

### Phase 3: Superlocus Analysis

- Maintain API compatibility with original MultiSampleQuantifier
- Use region-based quantification for performance
- Ensure thread safety for parallel processing

### Phase 4: Performance Optimization

- Implement profiling to identify bottlenecks
- Consider numpy vectorization for performance-critical sections
- Benchmark against original implementation

### Phase 5: GA4GH Compliance

- Implement GA4GH formatter classes
- Add stratification capabilities required by the standard
- Generate standard-compliant metrics output

## Acceptance Criteria

Each phase is considered complete when:
- All tests pass
- Documentation is updated
- Original functionality is maintained
- Performance is comparable to original implementation
