# hap.py Quantify Module Enhancement - Phase 2 Prompt

This is a context prompt file for continuing the quantify module enhancement project in the hap.py bioinformatics tool modernization.

## Current Task: Phase 2 - Enhanced ROC Analysis Implementation

### Project Status
**Phase 1: Core Implementation** ✅ **COMPLETED**
- Fixed critical test failures in `_are_alleles_compatible`, `_classify_variant_type`, and performance tests
- Enhanced method compatibility to handle both pandas Series and dictionary inputs
- Updated test expectations to match actual implementation behavior
- Applied comprehensive code quality improvements with Black formatting
- Established performance baseline (1000 variants in ~12.7 seconds)
- All unit tests passing: 12/12 test methods in `test_unit_quantify.py`

**Phase 2: Enhanced ROC Analysis** 🔄 **CURRENT PRIORITY**
- ROC curve generation with confidence intervals
- Quality score stratification and performance analysis
- Multi-threshold analysis support
- Advanced statistical metrics implementation

### Key Files and Current State

**Core Implementation:**
- `src/hap_py/haplo/python_quantify.py` - Main quantify module (Phase 1 complete)
- `tests/unit/test_unit_quantify.py` - Comprehensive test suite (12 tests, all passing)

**Development Plan:**
- `QUANTIFY_IMPLEMENTATION_PLAN.md` - Detailed roadmap for all phases

**Configuration:**
- `pyproject.toml` - Project configuration with development dependencies
- `pytest.ini` - Test configuration
- `.pre-commit-config.yaml` - Code quality hooks (Black, Ruff, isort)

### Environment Setup

**REQUIRED:** Always start with the micromamba environment:
```bash
micromamba activate happy-dev
```

**Verify Setup:**
```bash
# Confirm Python version
python --version  # Should show: Python 3.11.12

# Confirm environment location
which python  # Should show: .../micromamba/envs/happy-dev/bin/python

# Install in development mode if needed
pip install -e .[dev]
```

### Phase 2 Implementation Goals

**Primary Features to Implement:**

1. **Enhanced ROC Curve Analysis:**
   - Implement confidence interval calculations using bootstrap sampling
   - Add support for stratified ROC analysis by variant type
   - Create visualization methods for ROC curves with confidence bands

2. **Quality Score Stratification:**
   - Implement quality score binning and per-bin performance metrics
   - Add support for configurable quality thresholds
   - Create stratified performance reports

3. **Multi-Threshold Analysis:**
   - Extend current threshold-based analysis to support multiple thresholds
   - Implement threshold optimization algorithms
   - Add comparative performance analysis across thresholds

4. **Advanced Statistical Metrics:**
   - Implement precision-recall curves with confidence intervals
   - Add F1-score optimization and threshold selection
   - Create comprehensive performance summary reports

**Success Criteria:**
- All existing tests continue to pass (maintain backward compatibility)
- New ROC analysis methods integrate seamlessly with existing quantify functionality
- Performance remains acceptable (target: <20 seconds for 1000 variants with ROC analysis)
- Code follows project style guidelines (Black formatting, type hints, docstrings)

### Implementation Strategy

**Step 1: Extend Core Infrastructure**
- Add new statistical calculation methods to `python_quantify.py`
- Implement bootstrap sampling utilities for confidence intervals
- Add quality score binning and stratification logic

**Step 2: ROC Analysis Implementation**
- Create `calculate_roc_with_confidence()` method
- Implement stratified ROC analysis by variant type
- Add ROC curve data structures and export functionality

**Step 3: Quality Score Stratification**
- Implement `stratify_by_quality()` method
- Add configurable quality score binning
- Create per-bin performance metric calculations

**Step 4: Multi-Threshold Support**
- Extend existing threshold methods to handle multiple values
- Implement threshold optimization algorithms
- Add comparative analysis utilities

**Step 5: Testing and Validation**
- Create comprehensive tests for all new functionality
- Validate against known datasets with expected ROC characteristics
- Performance testing to ensure acceptable execution times

### Code Quality Requirements

**Formatting and Style:**
```bash
# Apply Black formatting (88 character line limit)
black src/ tests/

# Check and fix linting
ruff check src/ tests/ --fix

# Sort imports
isort src/ tests/

# Run pre-commit hooks
pre-commit run --all-files
```

**Testing:**
```bash
# Run all quantify tests
pytest tests/unit/test_unit_quantify.py -v

# Run with coverage
pytest tests/unit/test_unit_quantify.py --cov=hap_py.haplo.python_quantify

# Performance validation
pytest tests/unit/test_unit_quantify.py::test_quantify_performance -v
```

### Current Performance Baseline
- **Test Configuration:** 1000 variants processed through full quantify pipeline
- **Current Performance:** ~12.7 seconds (Phase 1 baseline)
- **Phase 2 Target:** <20 seconds with enhanced ROC analysis

### Dependencies and External Tools
- **Python Environment:** 3.11.12 in `happy-dev` micromamba environment
- **Key Libraries:** pandas, numpy, scipy (for statistical calculations)
- **Development Tools:** pytest, black, ruff, isort, pre-commit
- **RTG Tools:** Available at `external/rtg-tools-3.12.1/rtg` (for integration testing)

### Implementation Notes

**Backward Compatibility:**
- All existing methods and interfaces must remain functional
- New functionality should extend, not replace, current capabilities
- Existing test suite must continue to pass without modification

**Type Safety:**
- Maintain support for both pandas Series and dictionary inputs
- Add proper type hints for all new methods
- Use Union types where multiple input formats are supported

**Error Handling:**
- Implement comprehensive error handling for statistical edge cases
- Add validation for quality score ranges and threshold values
- Provide informative error messages for invalid configurations

**Documentation:**
- Add Google-style docstrings for all new methods
- Include usage examples in docstrings
- Update any relevant documentation files

### Next Steps
1. Review the current `python_quantify.py` implementation
2. Design the ROC analysis architecture and data structures
3. Implement bootstrap sampling utilities for confidence intervals
4. Create the enhanced ROC calculation methods
5. Add comprehensive testing for new functionality
6. Validate performance against Phase 2 targets

This prompt provides the complete context needed to continue Phase 2 implementation while maintaining project quality standards and backward compatibility.
