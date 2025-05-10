---
applyTo: "**/*{test,tests}*/**"
---
# Testing Standards for hap.py

## Test Data
- Use small representative test datasets when possible
- Include edge cases in test data (e.g., complex variants, edge of chromosomes)
- Document the source and meaning of test data files
- Prefer synthetic data over real patient data for tests

## Test Structure
- Migrate to pytest framework when updating existing tests
- Use fixtures for common setup and teardown
- Implement appropriate parameterization for testing multiple inputs
- Include both unit tests and integration tests

## Test Coverage
- Focus on critical variant comparison logic
- Test normalization and preprocessing functions thoroughly
- Include tests for edge cases in genomic data (e.g., chromosome ends, overlapping variants)
- Test performance with scaled-down but representative data

## Test Validation
- Include known-correct outputs for comparison
- Verify consistency with previous versions where appropriate
- Document expected behavior clearly
- Add regression tests for previously identified bugs
