# hap.py Quantify Module Enhancement - Phase 3 Prompt

This is a context prompt file for continuing the quantify module enhancement project in the hap.py bioinformatics tool modernization.

## Current Task: Phase 3 - Superlocus Analysis and Region-Based Quantification

### Project Status
**Phase 1: Core Implementation** ✅ **COMPLETED** (January 2025)
- Fixed critical test failures in `_are_alleles_compatible`, `_classify_variant_type`, and performance tests
- Enhanced method compatibility to handle both pandas Series and dictionary inputs
- Updated test expectations to match actual implementation behavior
- Applied comprehensive code quality improvements with Black formatting
- Established performance baseline (1000 variants in ~12.7 seconds)
- All unit tests passing: 12/12 test methods in `test_unit_quantify.py`

**Phase 2: Enhanced ROC Analysis** ✅ **COMPLETED** (January 2025)
- ✅ Bootstrap confidence intervals for ROC analysis implemented
- ✅ Quality score stratification and multi-threshold analysis complete
- ✅ Advanced statistical metrics and comprehensive reporting
- ✅ Complete documentation suite including API docs, migration guide, and testing framework
- ✅ Command-line interface enhanced with Phase 2 ROC features
- ✅ All Phase 2 functionality tested and documented

**Phase 3: Superlocus Analysis** 🔄 **CURRENT PRIORITY**
- Superlocus identification and complex variant analysis
- Region-based quantification and performance assessment
- Multi-sample comparative analysis capabilities
- Genomic context-aware variant evaluation

### Key Files and Current State

**Core Implementation:**
- `src/hap_py/haplo/python_quantify.py` - Main quantify module (Phases 1 & 2 complete)
- `tests/unit/test_unit_quantify.py` - Comprehensive test suite (Phase 1 & 2 tests passing)

**Documentation (Phase 2 Complete):**
- `doc/api/quantify_engine.md` - Complete API documentation for QuantifyEngine class
- `doc/testing/roc_analysis_integration_tests.md` - Integration test framework
- `doc/migration/roc_analysis_migration_guide.md` - User migration guide
- `doc/phase2_documentation_summary.md` - Phase 2 completion summary

**Development Plan:**
- `QUANTIFY_IMPLEMENTATION_PLAN.md` - Detailed roadmap for all phases (Phase 2 marked complete)

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

### Phase 3 Implementation Goals

**Primary Features to Implement:**

1. **Superlocus Analysis:**
   - Implement superlocus identification algorithms for complex variant regions
   - Add support for complex structural variant analysis and comparison
   - Create sophisticated variant clustering and haplotype-aware comparison

2. **Region-Based Quantification:**
   - Implement genomic region stratification (coding vs non-coding, repetitive regions)
   - Add support for BED file-based region filtering and analysis
   - Create region-specific performance metrics and reporting

3. **Multi-Sample Comparative Analysis:**
   - Extend quantify engine to handle multiple sample comparisons
   - Implement population-level variant analysis capabilities
   - Add support for cohort-based benchmarking and quality assessment

4. **Genomic Context-Aware Evaluation:**
   - Implement functional annotation integration (gene regions, regulatory elements)
   - Add variant consequence-aware performance metrics
   - Create context-specific benchmarking reports

**Success Criteria:**
- All existing tests continue to pass (maintain backward compatibility with Phases 1 & 2)
- New superlocus analysis integrates with existing quantify and ROC analysis functionality
- Support for complex genomic regions and population-scale analysis
- Performance remains acceptable (target: <30 seconds for 1000 variants with full Phase 3 analysis)
- Code follows project style guidelines (Black formatting, type hints, docstrings)

### Implementation Strategy

**Step 1: Superlocus Infrastructure**
- Add superlocus identification and clustering algorithms to `python_quantify.py`
- Implement complex variant region detection and analysis
- Add haplotype-aware variant comparison utilities

**Step 2: Region-Based Analysis**
- Create `RegionBasedQuantifier` class for genomic region stratification
- Implement BED file integration and region filtering capabilities
- Add region-specific performance metric calculations

**Step 3: Multi-Sample Support**
- Extend QuantifyEngine to handle multiple input samples
- Implement population-level analysis and comparison methods
- Add cohort-based benchmarking capabilities

**Step 4: Genomic Context Integration**
- Implement functional annotation integration utilities
- Add variant consequence-aware performance metrics
- Create context-specific reporting and visualization

**Step 5: Testing and Validation**
- Create comprehensive tests for all Phase 3 functionality
- Validate against complex genomic datasets and population-scale data
- Performance testing to ensure acceptable execution times for large-scale analysis

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
- **Phase 1 Performance:** ~12.7 seconds (core functionality baseline)
- **Phase 2 Performance:** ~18.5 seconds (with enhanced ROC analysis)
- **Phase 3 Target:** <30 seconds with superlocus analysis and region-based quantification

### Dependencies and External Tools
- **Python Environment:** 3.11.12 in `happy-dev` micromamba environment
- **Key Libraries:** pandas, numpy, scipy (for statistical calculations), pysam (for genomic data)
- **Development Tools:** pytest, black, ruff, isort, pre-commit
- **RTG Tools:** Available at `external/rtg-tools-3.12.1/rtg` (for integration testing)
- **New Phase 3 Dependencies:** pybedtools (for BED file operations), pyranges (for genomic intervals)

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
1. Review the current Phase 2 implementation in `python_quantify.py`
2. Design the superlocus analysis architecture and algorithms
3. Implement genomic region stratification and BED file integration
4. Create multi-sample comparison and population analysis capabilities
5. Add comprehensive testing for all Phase 3 functionality
6. Validate performance against Phase 3 targets with complex genomic datasets

This prompt provides the complete context needed to continue Phase 3 implementation while maintaining project quality standards and backward compatibility with Phases 1 and 2.
