#!/usr/bin/env python3
"""
GitHub Issues for hap.py Testing Problems

This script documents the issues that should be submitted to the repository
based on our test analysis.
"""


def main():
    print("GITHUB ISSUES TO SUBMIT - hap.py Testing Problems")
    print("=" * 60)

    issues = [
        {
            "title": "Integration Tests Hanging/Unresponsive During Execution",
            "priority": "HIGH",
            "labels": ["bug", "testing", "integration"],
            "description": """
## Problem Description

The integration test suite appears to hang or become unresponsive when running with `pytest tests/integration/`. This makes it difficult to complete the full test validation process.

## Observed Behavior

- Running `pytest tests/integration/ -v` causes the terminal to become unresponsive
- Individual tests may work but the full suite execution hangs
- Commands like `python -m pytest tests/integration/test_fastasize.py -v` don't return results in reasonable time

## Environment Details

- Python 3.11.12 in micromamba environment `happy-dev`
- pytest 8.3.5
- macOS (darwin platform)

## Impact

This prevents:
- Automated testing workflows
- Validation of the modernized codebase
- CI/CD pipeline implementation

## Potential Causes

1. **Test fixtures or setup issues**: Some tests may be waiting for external resources
2. **RTG tool path configuration**: Tests may be hanging while looking for RTG tools
3. **Temporary file/directory cleanup**: Tests may be stuck in filesystem operations
4. **Resource locks**: Tests may be competing for shared resources

## Suggested Solutions

1. **Add test timeouts**: Configure pytest with reasonable timeouts for integration tests
2. **Implement test isolation**: Ensure each test properly cleans up and doesn't affect others
3. **Fix external tool dependencies**: Ensure RTG and other tools are properly configured
4. **Add logging**: Increase verbosity to identify where tests are hanging

## Files Affected

- `tests/integration/` (all integration test files)
- `conftest.py` (pytest configuration)
- Individual test files that may have blocking operations
            """.strip(),
        },
        {
            "title": "GA4GH Integration Tests Failing with Import/Setup Errors",
            "priority": "MEDIUM",
            "labels": ["bug", "testing", "ga4gh"],
            "description": """
## Problem Description

Multiple GA4GH integration tests are showing ERROR status during execution, preventing validation of GA4GH compliance functionality.

## Observed Behavior

From partial test output, multiple tests in `test_ga4gh_integration.py` are failing with ERROR status:
- `test_init`
- `test_prepare_vcf_header`
- `test_transform_match_to_ga4gh`
- `test_annotate_vcf_record`
- `test_create_ga4gh_metrics`
- `test_write_ga4gh_metrics_file`

## Environment Details

- Python 3.11.12 in micromamba environment `happy-dev`
- All GA4GH unit tests are passing

## Impact

- Prevents validation of GA4GH compliance functionality
- Blocks testing of GA4GH integration with QuantifyEngine
- May indicate issues with modernized GA4GH implementation

## Potential Causes

1. **Import issues**: GA4GH modules may not be importing correctly in integration context
2. **Mock configuration**: Test fixtures and mock objects may be incorrectly configured
3. **Dependency issues**: Required libraries or modules may be missing
4. **Setup/teardown problems**: Test initialization may be failing

## Suggested Solutions

1. **Verify imports**: Check that all GA4GH modules import correctly in test environment
2. **Fix test fixtures**: Ensure mock objects and test fixtures are properly configured
3. **Add error handling**: Improve error reporting to identify specific failure points
4. **Validate dependencies**: Ensure all required dependencies are available

## Files Affected

- `tests/integration/test_ga4gh_integration.py`
- `src/hap_py/haplo/ga4gh_integration.py`
- `src/hap_py/haplo/ga4gh_compliance.py`
            """.strip(),
        },
        {
            "title": "VCF Processing Integration Tests Failing",
            "priority": "MEDIUM",
            "labels": ["bug", "testing", "vcf-processing"],
            "description": """
## Problem Description

Several integration tests related to core VCF processing functionality are failing, including chromosome prefix handling and variant decomposition.

## Observed Behavior

From partial test output, the following tests are failing:
- `test_chrprefix.py::test_numeric_chrs`
- `test_chrprefix.py::test_chr_prefixed`
- `test_chrprefix.py::test_mixed_chr_prefix`
- `test_decomp.py::test_decomp`
- `test_faulty_variants.py::test_faulty_variant_handling`

## Environment Details

- Python 3.11.12 in micromamba environment `happy-dev`
- RTG tools available at expected path
- Test data files appear to be present

## Impact

- Core VCF processing functionality may have regressions
- Chromosome handling logic may not be working correctly
- Variant decomposition may be broken

## Potential Causes

1. **RTG tool integration**: Path configuration or tool execution issues
2. **Test data**: Reference files or test data may be incorrect/missing
3. **VCF processing pipeline**: Modernization may have introduced bugs
4. **Path handling**: File path resolution issues in modernized code

## Suggested Solutions

1. **Verify RTG integration**: Ensure RTG tools are properly configured and accessible
2. **Check test data**: Validate that all required test data files exist and are correct
3. **Test individual components**: Run VCF processing components separately to isolate issues
4. **Review modernization**: Check for regressions introduced during Python 3 conversion

## Files Affected

- `tests/integration/test_chrprefix.py`
- `tests/integration/test_decomp.py`
- `tests/integration/test_faulty_variants.py`
- VCF processing modules in `src/hap_py/haplo/`
            """.strip(),
        },
    ]

    for i, issue in enumerate(issues, 1):
        print(f"\n{'='*60}")
        print(f"ISSUE #{i}: {issue['title']}")
        print(f"Priority: {issue['priority']}")
        print(f"Labels: {', '.join(issue['labels'])}")
        print(f"{'='*60}")
        print(issue["description"])
        print()

    print("\n" + "=" * 60)
    print("SUMMARY:")
    print(f"Total Issues to Submit: {len(issues)}")
    print("High Priority: 1")
    print("Medium Priority: 2")
    print("=" * 60)


if __name__ == "__main__":
    main()
