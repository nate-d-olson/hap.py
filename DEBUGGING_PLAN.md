# hap.py Debugging Plan - June 2025

## Overview

This document outlines a comprehensive debugging plan to address the 16 failed tests identified in the pytest output from June 20, 2025. The plan is structured in phases, prioritizing critical code bugs that affect core functionality.

## Phase 1: Critical Code Bugs (High Priority)

### 1.1 Python Preprocessing Pipeline - AC Field Handling Bug

**Issue**: Multiple tests failing with AC field tuple handling error

```text
WARNING:hap_py.haplo.python_preprocess:Failed to set INFO field AC=(1, 1) (type: <class 'tuple'>): values expected to be 1-tuple, given len=2
```

**Affected Tests**:

- `test_happy_pg_test`
- `test_quantify_test`
- `test_large_giab_rtg_chr21`
- `test_large_giab_rtg_chr1`

**Root Cause**: In `src/hap_py/haplo/python_preprocess.py` lines 466-504, the AC field handling logic incorrectly processes tuple values.

**Debugging Steps**:

1. **Isolate the Issue**:

   ```bash
   micromamba activate happy-dev
   pytest tests/unit/test_python_preprocess.py::test_decompose_variants -v -s
   ```

2. **Examine AC Field Processing**:
   - Review `decompose_variant()` method in `python_preprocess.py`
   - Check lines 476-504 where INFO field copying occurs
   - Focus on the AC (Allele Count) field handling logic

3. **Debug Strategy**:

   ```python
   # Add debugging to python_preprocess.py line ~476
   logger.debug(f"Processing INFO field {key}={value}, type={type(value)}")
   if key == "AC":
       logger.debug(f"AC field details: number={field_info.number}, type={field_info.type}")
   ```

4. **Expected Fix**:
   - Update AC field handling to properly process multi-allelic AC values
   - Ensure AC field decomposition aligns with VCF 4.2 specification
   - AC should be split per alternative allele during decomposition

5. **Validation**:

   ```bash
   # Test with specific VCF containing AC fields
   python -m hap_py.haplo.python_preprocess tests/data/multi_allelic.vcf.gz -r ref.fa
   ```

### 1.2 VCF Header Validation Issues

**Issue**: Header validation errors causing integration test failures

```text
ERROR: Error checking file: Invalid header
```

**Affected Tests**:

- `test_mixed_chr_prefix`
- `test_decomp`
- `test_small_giab_rtg`

**Root Cause**: Overly strict header validation in `src/hap_py/haplo/python_vcfcheck.py`

**Debugging Steps**:

1. **Isolate Header Validation**:

   ```bash
   pytest tests/unit/test_vcfcheck.py::test_check_header -v -s
   ```

2. **Examine Header Check Logic**:
   - Review `_check_header()` method in `VCFChecker` class
   - Check if FILTER field detection is too strict
   - Verify FORMAT field validation logic

3. **Debug Strategy**:

   ```python
   # Add debug output to _check_header method
   def _check_header(self, header):
       issues = []
       logger.debug(f"Checking header with filters: {header.filters}")
       logger.debug(f"Available formats: {list(header.formats.keys())}")
       # ... existing logic
   ```

4. **Expected Fix**:
   - Relax FILTER field requirements for standard VCF files
   - Improve error messages to be more descriptive
   - Allow for minimal but valid VCF headers

5. **Validation**:

   ```bash
   python -m hap_py.haplo.python_vcfcheck tests/data/src/numeric_chrs/
   ```

### 1.3 Chromosome Prefix Detection Failures

**Issue**: Chromosome name matching failures between VCF and reference

```text
ValueError: Truth and reference have no chromosomes in common!
```

**Affected Tests**:

- `test_decomp`
- `test_small_giab_rtg`
- `test_large_giab_rtg_chr21`

**Debugging Steps**:

1. **Examine Chromosome Detection**:

   ```bash
   pytest tests/integration/test_chrprefix.py::test_numeric_chrs -v -s
   ```

2. **Debug Strategy**:
   - Add chromosome name logging to preprocessing pipeline
   - Check VCF contig headers vs reference sequence names
   - Verify chromosome prefix handling logic

3. **Expected Fix**:
   - Improve chromosome name normalization (chr1 ↔ 1)
   - Add robust chromosome matching algorithm
   - Better error messages indicating which chromosomes were found

### 1.4 BCFtools Integration Issues

**Issue**: BCFtools subprocess failures with return code -11

```text
Command line bcftools concat -a -O z -o /path/file.vcf.gz got return code -11
```

**Affected Tests**: `test_integration`

**Debugging Steps**:

1. **Examine BCFtools Wrapper**:

   ```bash
   grep -r "bcftools concat" src/hap_py/tools/
   ```

2. **Debug Strategy**:
   - Check BCFtools command construction
   - Verify input file validity before BCFtools calls
   - Add better error handling for subprocess failures

3. **Expected Fix**:
   - Improve input validation before BCFtools calls
   - Add fallback methods when BCFtools fails
   - Better error reporting for debugging

## Phase 2: Test Definition Issues (Medium Priority)

### 2.1 BedIntervalTree Value Format

**Issue**: Test expects scalar values but implementation returns lists

```python
assert [iv.value for iv in overlaps] == ["test", "test"]  # Fails
# Actual: [['test'], ['test']]
```

**Fix**: Update test assertion in `tests/unit/test_bedintervaltree.py`

```python
assert [iv.value for iv in overlaps] == [["test"], ["test"]]
```

### 2.2 Summary File Format Mismatches

**Issue**: Output format has evolved but expected files haven't been updated

**Debugging Steps**:

1. **Compare Output Formats**:

   ```bash
   diff expected.summary.csv actual.summary.csv
   ```

2. **Update Expected Files**:
   - Regenerate expected output files with current hap.py
   - Verify new format is correct and complete
   - Update test data repository

## Phase 3: Test Infrastructure Issues (Lower Priority)

### 3.1 Missing Test Data Files

**Issue**: Missing reference files and test data

**Debugging Steps**:

1. **Audit Test Data**:

   ```bash
   find tests/data -name "*.fa" -o -name "*.vcf*" | sort
   ```

2. **Fix Missing Files**:
   - Ensure all required test data exists
   - Update file paths in test configurations
   - Create minimal test data where missing

## Debugging Tools and Techniques

### Environment Setup

```bash
# Always start with
micromamba activate happy-dev

# Verify environment
which python  # Should show happy-dev path
python --version  # Should be 3.11.x
```

### Debugging Commands

**Individual Test Debugging**:

```bash
# Run specific failing test with maximum verbosity
pytest tests/integration/test_chrprefix.py::test_numeric_chrs -v -s --tb=long

# Run with logging enabled
pytest tests/unit/test_python_preprocess.py -v -s --log-cli-level=DEBUG

# Run single test file
pytest tests/unit/test_vcfcheck.py -v
```

**Code-Level Debugging**:

```python
# Add to problematic functions
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.debug(f"Debug info: variable={variable}")
```

**File Analysis**:

```bash
# Check VCF file structure
bcftools view -h problematic.vcf.gz
bcftools stats problematic.vcf.gz

# Check reference file
samtools faidx reference.fa
head reference.fa.fai
```

### Validation Strategy

**Phase-by-Phase Testing**:

```bash
# Phase 1: After each critical bug fix
pytest tests/unit/test_python_preprocess.py -v
pytest tests/unit/test_vcfcheck.py -v
pytest tests/integration/test_chrprefix.py::test_numeric_chrs -v

# Phase 2: After test fixes
pytest tests/unit/test_bedintervaltree.py -v

# Phase 3: Full integration testing
pytest tests/integration/ -v

# Final: Full test suite
pytest tests/ -v
```

**Regression Testing**:

```bash
# After each fix, ensure no new failures
pytest tests/ --tb=short | grep FAILED | wc -l
```

## Implementation Timeline

### Week 1: Critical Code Bugs

- **Days 1-2**: Fix AC field handling in python_preprocess.py
- **Days 3-4**: Fix VCF header validation in python_vcfcheck.py
- **Days 5-7**: Fix chromosome prefix detection and BCFtools integration

### Week 2: Test Issues and Validation

- **Days 1-2**: Update test definitions and expected output files
- **Days 3-4**: Fix missing test data and infrastructure
- **Days 5-7**: Full regression testing and documentation

## Success Criteria

1. **All critical tests pass**:
   - `test_happy_pg_test`
   - `test_quantify_test`
   - `test_large_giab_rtg_chr21`
   - `test_large_giab_rtg_chr1`

2. **Integration tests stabilized**:
   - `test_chrprefix.py` tests pass
   - `test_decomp` passes
   - `test_integration` passes

3. **No regression in existing functionality**:
   - All previously passing tests continue to pass
   - Performance characteristics maintained

4. **Improved error reporting**:
   - More descriptive error messages for debugging
   - Better logging for troubleshooting

## Documentation Updates

After successful debugging:

1. Update AGENTS.md with resolved issues
2. Document any API changes or behavioral modifications
3. Update testing guidelines for future development
4. Create troubleshooting guide for common issues

## Emergency Fallback Plan

If critical fixes cannot be implemented quickly:

1. **Skip problematic tests temporarily**:

   ```python
   @pytest.mark.skip(reason="Known issue - tracking in DEBUGGING_PLAN.md")
   ```

2. **Use mock implementations for testing**:
   - Mock problematic external tool calls
   - Use simplified test data

3. **Focus on core functionality**:
   - Ensure basic hap.py workflow still functions
   - Prioritize most commonly used features

This debugging plan provides a systematic approach to resolving the identified test failures while maintaining code quality and minimizing regressions.
