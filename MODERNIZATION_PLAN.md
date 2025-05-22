# hap.py Modernization Plan

## Overview

This document outlines a systematic plan to continue the Python 3 modernization of hap.py, with emphasis on replacing custom C++ code with established Python packages like pysam, maintaining clean code standards, and tracking progress effectively.

## Current Status Assessment

### ✅ Completed Modernization Steps

- [x] Basic Python 3.8+ compatibility established
- [x] Modern packaging with pyproject.toml (PEP 517/518)
- [x] Pre-commit hooks configured (black, ruff, isort, pyupgrade)
- [x] pysam dependency added for VCF/BCF handling
- [x] pytest framework for testing
- [x] Mock implementations for C++ fallbacks
- [x] Type hints infrastructure started
- [x] String handling utilities for Python 3 compatibility

### 🔄 In Progress

- [x] C++ dependency reduction
  - [x] blocksplit → Python implementation using pysam
  - [x] vcfcheck → Python implementation using pysam
  - [x] sequence utilities → Python implementation using BioPython
  - [x] quantify → Python implementation using pandas + pysam
  - [x] preprocess → Python implementation using pysam
  - [x] hapcmp → Python implementation using pysam (initial version complete)
- [x] Progress tracking and reporting infrastructure
- [ ] Shell script to pytest migration
- [ ] Complete type hint coverage
- [ ] Documentation modernization

### ❌ Not Started

- [ ] Systematic C++ component evaluation (for xcmp, scmp, multimerge)
- [ ] Performance benchmarking of Python vs C++ components (for hapcmp, xcmp, scmp, multimerge)
- [ ] CI/CD pipeline optimization
- [ ] Package distribution setup

## Modernization Progress Update

### Completed Components and Tools

The following components have been migrated from C++ to Python:

1. **blocksplit**
   - Created `python_blocksplit.py` using pysam
   - Added comprehensive tests in `test_python_blocksplit.py`
   - Features full compatibility with the original C++ version

2. **vcfcheck**
   - Created `python_vcfcheck.py` using pysam
   - Added tests in `test_vcfcheck.py`
   - Supports the same validation features as the C++ version

3. **Sequence Utilities**
   - Created `sequence_utils.py` using BioPython
   - Added tests in `test_sequence_utils.py`
   - Includes complement, reverse complement, and FASTA handling

4. **quantify**
   - Created `python_quantify.py` using pandas and pysam
   - Added tests in `test_quantify.py`
   - Supports variant comparison and stratification metrics

5. **hapcmp**
   - Created `python_hapcmp.py` using pysam
   - Added initial unit tests in `tests/unit/test_python_hapcmp.py`
   - Basic haplotype comparison functionality implemented.

### Progress Tracking Infrastructure

The following tools have been developed to track modernization progress:

1. **track_progress.py**
   - Monitors C++ to Python migration progress
   - Tracks Python 2 artifacts in the codebase
   - Reports on code quality metrics

2. **generate_report.py** and **generate_dashboard.py**
   - Creates detailed progress reports in Markdown
   - Generates visual dashboards of modernization status
   - Helps identify priority components for migration

3. **fix_python2_artifacts.py**
   - Automatically fixes common Python 2 artifacts
   - Handles print statements, xrange, string handling, etc.
   - Integrates with pre-commit hooks

4. **validate_replacements.py**
   - Validates Python replacements against C++ originals
   - Compares outputs for functional equivalence
   - Reports performance differences

5. **run_modernization.py**
   - Orchestrates the entire modernization process
   - Runs all tools in the correct sequence
   - Provides a simple interface for the modernization effort

### Next Steps

- Implement Python versions of `xcmp` and `scmp`
- Migrate remaining shell script tests to pytest
- Complete type hint coverage for all Python code
- Update documentation for Python 3 compatibility
- Create CI/CD pipeline for automated testing
- Refine and complete unit tests for `python_hapcmp.py`

## Phase 1: Code Quality and Standards (Weeks 1-2)

### 1.1 Automated Code Formatting and Linting

**Goal**: Establish consistent code quality across the entire codebase.

**Tasks**:

- [ ] Run pre-commit on all Python files and fix issues
- [ ] Configure ruff with project-specific rules
- [ ] Set up automated type checking with mypy
- [ ] Create code quality metrics baseline

**Commands to Execute**:

```bash
# Install and setup pre-commit
pip install pre-commit
pre-commit install

# Run all hooks on existing codebase
pre-commit run --all-files

# Fix any issues found
black src/python/
ruff check src/python/ --fix
isort src/python/
pyupgrade --py38-plus src/python/**/*.py
```

**Success Criteria**:

- All pre-commit hooks pass without errors
- Code coverage report generated
- No Python 2 syntax remaining

### 1.2 Dependency Modernization

**Goal**: Replace custom implementations with established packages.

**Priority Areas for Package Replacement**:

| Component | Current State | Replacement Package | Priority | Effort |
|-----------|---------------|-------------------|----------|---------|
| VCF/BCF parsing | Custom C++ + htslib | pysam | High | Medium |
| Sequence manipulation | Custom C++ | BioPython | High | Low |
| Statistical calculations | Custom C++ | numpy/scipy | Medium | Medium |
| File compression | Custom C++ | gzip/lzma modules | Medium | Low |
| Threading/parallel | Custom C++ | concurrent.futures | Low | Medium |

**Tasks**:

- [ ] Audit all C++ components in `src/c++/`
- [ ] Identify components that can be replaced with Python packages
- [ ] Create compatibility layer for gradual migration
- [ ] Performance benchmark critical paths

## Phase 2: C++ Dependency Reduction (Weeks 3-6)

### 2.1 VCF/BCF Processing Migration

**Goal**: Replace custom C++ VCF handling with pysam.

**Current C++ Tools to Replace**:

- `blocksplit` → Python implementation using pysam
- `quantify` → Python implementation using pandas + pysam
- `vcfcheck` → Python implementation using pysam validation
- `preprocess` → Python implementation using pysam + preprocessing logic
- `hapcmp` → Python implementation using pysam (initial version complete)
- `xcmp` → Python implementation (planned)
- `scmp` → Python implementation (planned)

**Implementation Strategy**:

```python
# Example: Replace blocksplit with Python
import pysam
import pandas as pd
from typing import Iterator, Tuple

def python_blocksplit(vcf_file: str, block_size: int = 1000) -> Iterator[Tuple[str, int, int]]:
    """
    Python replacement for C++ blocksplit using pysam.
    """
    with pysam.VariantFile(vcf_file) as vcf:
        current_block = []
        for record in vcf:
            current_block.append(record)
            if len(current_block) >= block_size:
                yield process_block(current_block)
                current_block = []
```

**Tasks**:

- [ ] Create Python equivalents for each C++ tool
- [ ] Implement compatibility wrappers
- [ ] Add comprehensive unit tests
- [ ] Performance comparison with C++ versions
- [ ] Gradual migration with feature flags

### 2.2 Sequence Processing Migration

**Goal**: Replace custom sequence manipulation with BioPython.

**Components to Replace**:

- Complement/reverse complement functions
- FASTA file handling
- Reference sequence processing

**Implementation**:

```python
from Bio.Seq import Seq
from Bio import SeqIO
from typing import Union

def complement_sequence(seq: Union[str, bytes]) -> str:
    """Replace C++ complement with BioPython implementation."""
    if isinstance(seq, bytes):
        seq = seq.decode('utf-8')
    return str(Seq(seq).complement())

def reverse_complement(seq: Union[str, bytes]) -> str:
    """Replace C++ reverse complement with BioPython implementation."""
    if isinstance(seq, bytes):
        seq = seq.decode('utf-8')
    return str(Seq(seq).reverse_complement())
```

## Phase 3: Testing and Validation (Weeks 7-8)

### 3.1 Comprehensive Test Suite

**Goal**: Ensure all functionality works correctly after modernization.

**Testing Strategy**:

- [ ] Unit tests for all new Python components
- [ ] Integration tests comparing Python vs C++ outputs
- [ ] Performance regression tests
- [ ] End-to-end validation with real datasets

**Test Infrastructure**:

```python
# tests/test_modernization.py
import pytest
from unittest.mock import patch
import tempfile
import os

class TestModernization:
    """Test suite for modernization components."""

    def test_python_vs_cpp_blocksplit(self):
        """Compare Python and C++ blocksplit outputs."""
        # Implementation here
        pass

    def test_pysam_integration(self):
        """Test pysam integration works correctly."""
        # Implementation here
        pass

    def test_performance_regression(self):
        """Ensure performance is acceptable."""
        # Implementation here
        pass
```

### 3.2 Migration Validation

**Goal**: Verify that migrated components produce identical results.

**Validation Process**:

- [ ] Output comparison tests
- [ ] Binary compatibility verification
- [ ] Performance benchmarking
- [ ] Memory usage analysis

## Phase 4: Documentation and Deployment (Weeks 9-10)

### 4.1 Documentation Updates

**Goal**: Update all documentation for Python 3 and new dependencies.

**Tasks**:

- [ ] Update README.md with new installation instructions
- [ ] Create migration guide for users
- [ ] Update API documentation
- [ ] Add troubleshooting guide

### 4.2 Distribution Setup

**Goal**: Enable easy installation via pip.

**Tasks**:

- [ ] Test PyPI package building
- [ ] Set up CI/CD for automated releases
- [ ] Create conda package
- [ ] Set up Docker containers

## Progress Tracking System

### 4.3 Tracking Mechanism

Create a progress tracking system to monitor modernization:

```python
# scripts/track_progress.py
import json
import subprocess
from pathlib import Path
from typing import Dict, List

class ModernizationTracker:
    """Track progress of modernization efforts."""

    def __init__(self, config_file: str = "modernization_progress.json"):
        self.config_file = Path(config_file)
        self.load_progress()

    def check_code_quality(self) -> Dict[str, bool]:
        """Check code quality metrics."""
        results = {}

        # Check pre-commit status
        try:
            subprocess.run(["pre-commit", "run", "--all-files"],
                         check=True, capture_output=True)
            results["pre_commit"] = True
        except subprocess.CalledProcessError:
            results["pre_commit"] = False

        # Check type coverage
        results["type_coverage"] = self.check_type_coverage()

        # Check test coverage
        results["test_coverage"] = self.check_test_coverage()

        return results

    def check_c_plus_plus_dependencies(self) -> Dict[str, int]:
        """Count remaining C++ dependencies."""
        cpp_files = list(Path("src/c++").rglob("*.cpp"))
        cpp_usage = {}

        for component in ["blocksplit", "quantify", "vcfcheck", "preprocess"]:
            cpp_usage[component] = len([f for f in cpp_files if component in f.name])

        return cpp_usage

    def generate_report(self) -> str:
        """Generate modernization progress report."""
        quality = self.check_code_quality()
        cpp_deps = self.check_c_plus_plus_dependencies()

        report = f"""
# Modernization Progress Report

## Code Quality Status
- Pre-commit hooks: {'✅ PASS' if quality['pre_commit'] else '❌ FAIL'}
- Type coverage: {quality.get('type_coverage', 0)}%
- Test coverage: {quality.get('test_coverage', 0)}%

## C++ Dependency Status
"""
        for component, count in cpp_deps.items():
            status = "✅ Migrated" if count == 0 else f"🔄 {count} files remaining"
            report += f"- {component}: {status}\n"

        return report
```

### 4.4 Automated Progress Monitoring

Create GitHub Actions workflow for continuous monitoring:

```yaml
# .github/workflows/modernization-progress.yml
name: Modernization Progress

on:
  push:
    branches: [ main, python3-migration ]
  pull_request:
    branches: [ main ]

jobs:
  track-progress:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3

    - name: Set up Python
      uses: actions/setup-python@v3
      with:
        python-version: '3.8'

    - name: Install dependencies
      run: |
        pip install -e .[dev]
        pip install pre-commit

    - name: Generate progress report
      run: |
        python scripts/track_progress.py > PROGRESS_REPORT.md

    - name: Update progress badge
      run: |
        # Update README badge with current progress
        python scripts/update_progress_badge.py

    - name: Comment PR with progress
      if: github.event_name == 'pull_request'
      uses: actions/github-script@v6
      with:
        script: |
          const fs = require('fs');
          const report = fs.readFileSync('PROGRESS_REPORT.md', 'utf8');
          github.rest.issues.createComment({
            issue_number: context.issue.number,
            owner: context.repo.owner,
            repo: context.repo.repo,
            body: report
          });
```

## Critical Review Points

### Weekly Review Checkpoints

1. **Week 2**: Code quality metrics established, no Python 2 artifacts remaining
2. **Week 4**: At least 50% of C++ components have Python equivalents
3. **Week 6**: All critical path components migrated and tested
4. **Week 8**: Performance benchmarks meet acceptance criteria
5. **Week 10**: Full migration complete, documentation updated

### Success Criteria for Each Phase

- **Phase 1**: 100% pre-commit compliance, zero Python 2 syntax
- **Phase 2**: Functional Python replacements for all target C++ components
- **Phase 3**: <5% performance regression on critical paths
- **Phase 4**: Successful pip installation from PyPI

### Risk Mitigation

- **Performance Risk**: Maintain C++ fallbacks during transition
- **Compatibility Risk**: Extensive regression testing
- **User Experience Risk**: Clear migration documentation and support

## Implementation Commands

### Daily Development Workflow

```bash
# Start work session
git checkout python3-migration
git pull origin python3-migration

# Run automated checks
pre-commit run --all-files
python scripts/track_progress.py

# Make changes and validate
python -m pytest tests/
python scripts/validate_migration.py

# Commit with progress update
git add .
git commit -m "feat: migrate blocksplit to pysam

Progress: C++ components remaining: 3/7"
git push origin python3-migration
```

### Weekly Progress Review

```bash
# Generate comprehensive report
python scripts/track_progress.py --detailed > weekly_report.md

# Run performance benchmarks
python scripts/benchmark_migration.py --compare-all

# Update documentation
python scripts/update_docs.py --auto-generate
```

## Next Immediate Steps (Week 1)

1. **Set up progress tracking**:

   ```bash
   pip install pre-commit
   pre-commit install
   python scripts/create_modernization_tracker.py
   ```

2. **Establish baseline metrics**:

   ```bash
   pre-commit run --all-files 2>&1 | tee baseline_issues.log
   python scripts/count_cpp_dependencies.py > baseline_cpp.json
   ```

3. **Create migration branch**:

   ```bash
   git checkout -b python3-migration-systematic
   git push -u origin python3-migration-systematic
   ```

4. **Start with highest-impact, lowest-risk migration**:
   - Begin with sequence processing functions (BioPython)
   - Create Python equivalents alongside C++ versions
   - Add feature flags for gradual rollout

This plan provides a systematic, trackable approach to continuing the modernization while prioritizing established Python packages over custom C++ implementations.
