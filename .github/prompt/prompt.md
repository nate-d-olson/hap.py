## Next Development Session: hap.py MVP Implementation

This session will focus on implementing Phase 1 of the hap.py MVP plan, which involves foundational stability and automation.

### Phase 1: Foundational Stability & Automation

**Objective**: Establish a reliable CI/CD pipeline and ensure core tests are robust.

*   **1.1 Configure and Enable CI/CD**:
    *   Review and activate existing GitHub Actions workflows (`.github/workflows/ci.yml`, `python-ci.yml`).
    *   Ensure all code quality checks (Black, Ruff, isort, mypy) are integrated and enforced in CI.
    *   Configure CI to run all unit and integration tests on every push/pull request.
    *   Set up automated reporting for test results and code coverage.
*   **1.2 Improve Unit Test Coverage**:
    *   Identify critical modules/functions in `src/hap_py/haplo/` (e.g., `vcfeval.py`, `python_preprocess.py`, `quantify.py`) that lack sufficient unit test coverage.
    *   Prioritize porting complex logic from existing integration tests into new, focused unit tests.
    *   Utilize `pytest` fixtures (`tests/conftest.py`) and mocking (e.g., `unittest.mock`) to isolate unit tests from external dependencies and large data.
    *   Aim for a target unit test coverage percentage (e.g., 80-90%) for core modules.
*   **1.3 Stabilize Existing Tests**:
    *   Address any flaky or consistently failing tests (unit or integration).
    *   Ensure integration tests correctly handle external tool availability (e.g., RTG, bcftools) by skipping if not present, rather than failing.
    *   Refine test data management to avoid issues with missing or outdated example data.

**Please begin by reviewing the existing CI/CD workflows and identifying any gaps or areas for improvement based on the objectives above.**
