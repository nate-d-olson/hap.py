---
mode: 'agent'
description: "Multi-Phase Modernization Roadmap: Python 3 Migration for hap.py"
tools: ["file_search", "semantic_search","read_file", "insert_edit_into_file","replace_string_in_file", "terminal", "terminalLastCommand", "python_style_guide", "store_memory", "planning","reasoning","task_management"]
---

# Multi-Phase Modernization Roadmap: Python 3 Migration for hap.py (CONTINUATION)

This document tracks the ongoing modernization of hap.py for Python 3, focusing on PEP 517/518 compliance, streamlined installation, reproducibility, and making releases consumable via standard Python tooling. This is a living document and should be updated as each phase is completed.
Core funcationality to run hap.py as described in [doc/vcfeval_best_practices.md](doc/vcfeval_best_practices.md) should be maintained. Future validation test will compare output generated using this command for the old and mondernized versions. This document focuses on the modernization of the installation and testing process.

## Recent Changes (as of May 19, 2025)

- **PEP 517/518 Build System:**
  - `pyproject.toml` now uses `setuptools.build_meta` as the build backend.
  - Added compatibility `setup.py` for backward compatibility.
  - Fixed critical install errors by removing `backend-path` and unused dependencies.
- **Dependency Consolidation:**
  - All runtime and optional dependencies are now declared in `pyproject.toml`.
  - Confirmed and removed all references to `pybedtools` (no longer required or present in code).
- **C/C++/Cython Build Integration:**
  - Updated `CMakeLists.txt` for Python 3 and scikit-build compatibility.
- **Continuous Integration:**
  - `.github/workflows/ci.yml` builds, lints, and tests across OSes and Python versions.
- **Documentation:**
  - `README.md` updated for new install/build instructions and markdown linting.
  - Created `MIGRATION_GUIDE.md` to help users transition to the new installation method.
- **Install Testing:**
  - Verified that `pip install .` now works after dependency cleanup.
- **Console Script Entry Points:**
  - Added entry points for all command-line tools in `pyproject.toml`.
- **Test Framework:**
  - Created `pytest.ini` and `conftest.py` for pytest configuration.
  - Migrated 12 shell tests to pytest format including leftshift_test.sh and quantify_stratification_test.sh.
  - Created test utilities in `tests/utils.py`.
  - Added migration script to help convert shell tests to pytest.
- **Type Annotations:**
  - Added type hints to core modules including Haplo/cython_compat.py.
  - Started adding type hints to Tools/bcftools.py and Haplo/quantify.py.
- **Deprecation:**
  - Added deprecation warnings to `install.py` for legacy installation method.

## Next Steps

1. **Complete Migration of Shell Tests to Pytest (Priority)**
   - Continue migrating the remaining 12 shell tests to pytest format.
   - Ensure all tests have appropriate markers (integration, cpp, etc.).
   - Add more test fixtures as needed for common test scenarios.

2. **Improve Type Checking**
   - Continue adding type hints to remaining core modules.
   - Run mypy checks to ensure type consistency.
   - Update docstrings to Google-style format.
   - Use black and ruff to format code.

3. **Enhance Documentation**
   - Complete API docs using Sphinx autodoc.
   - Consider setting up ReadTheDocs integration.
   - Create man pages for command-line tools.

4. **Refine Build System**
   - Test wheel building on Windows.
   - Configure for PyPI publication.
   - Document build requirements in README.

5. **Prepare for 1.0.0 Release**
   - Complete test coverage of critical components.
   - Remove deprecated functionality and installation methods.
   - Create release checklist and publish to PyPI.

3. **Improve Type Checking**
   - Add more type hints to core Python modules.
   - Configure mypy for strict checking where appropriate.
   - Document type-related design decisions.

4. **Improve Build System Integration**
   - Further refine scikit-build integration.
   - Add more helpful error messages for common build issues.
   - Test on additional platforms (Windows).

5. **Prepare for Release**
   - Complete all test migrations and ensure CI passes.
   - Verify build artifacts (sdist, wheel) work correctly.
   - Plan for v1.0.0 release with final removal of legacy code.

## Considerations and TODOs

- **Code Quality Tools:** Use pre-commit hooks to enforce code quality standards.
- **Type Annotation Progress:** Track the progress of adding type hints to key modules.
- **Automated Testing:** Improve CI/CD workflow to run tests on all platforms.
- **Windows Support:** Test and document Windows-specific installation steps.
- **Migration Guide:** Keep the migration guide updated with the latest best practices.
- **Python 3.7+ Compatibility:** Ensure all code is compatible with Python 3.7 and newer.

---
**This document is a continuation and should be updated as each phase of the modernization is completed.**

1. Replace ad-hoc install scripts with a PEP 517/518 build

What’s there today

install.py tall_py3.py wrappers

setup.cfg + an implicit setup.py buried in CMake integration

a bunch of shell scripts (configure.sh, diagnose_build.sh, etc.)

Where to go

Adopt a pyproject.toml that names setuptools (or Poetry/Flit) as the build backend.

Declare all build-time requirements (Cython, pytest-runner if needed) under [build-system].

Migrate the metadata (name, version, dependencies, extras) from setup.cfg into pyproject.toml or keep a minimal setup.cfg for declarative fields.

Benefits:

Users can now do

pip install .

(or even pip install git+https://…) without having to run custom scripts.

IDEs and editors immediately recognize your project as PEP-compliant.

2. Clean up and consolidate dependency declarations

Move everything from happy.requirements.txt and happy.requirements.py3.txt into the install_requires and extras_require sections in your build config.

If there are optional features (e.g. Cython/C++ integration, plotting), expose them via extras:

[project.optional-dependencies]
cython = ["cython>=0.x"]
docs   = ["sphinx", "sphinx-rtd-theme"]

Pin minimal versions only where strictly needed (e.g. ABI abstractions, binary wheels).

3. Leverage modern C/C++/Cython build integration

Currently, CMake drives the C++ bits and you have tests like test_cpp_integration.sh. To align with Python packaging:

Use setuptools’ CMake support (via scikit-build or cmake-python-distributions) so that pip wheel . will invoke your CMake step automatically.

Consolidate all Cython extensions under a single ext_modules list in your setup.py/pyproject.toml config.

Clearly document the toolchain requirements for users who need to build from source (e.g. “CMake ≥ 3.15, a C++14 compiler, Python headers”).

4. Continuous Integration & automated release artifacts

GitHub Actions workflow matrix for Ubuntu/macOS/Windows × Python 3.8–3.11:

Build sdist + wheel

Run your full test suite (unit + C++ integration + doctests)

Lint/check type hints (black, flake8, mypy)

Upon a successful tagged release, automatically publish to PyPI (using pypa/gh-action-pypi-publish) and/or GitHub Packages.

5. Improve the user-facing CLI and documentation

Declare console_scripts entry points in your build config, e.g.

[project.scripts]
hap = "happy.main:cli"

Update README to show:

Installation

pip install hap.py

Quick-start example using the hap command

Advanced build instructions only if users need to work from source.

Include an automatically generated man page or Markdown API guide via Sphinx’s autodoc and host it on ReadTheDocs.

6. Gradual deprecation of legacy artifacts

Once the PEP 517/518 flow is solid, mark install.py, install_py3.py, and related migration shell scripts as deprecated, with warnings directing users to the new workflow.

Remove Python 2–only branches and code paths. Ensure your CI matrix no longer tests on 2.x.

In a major version bump (e.g. v1.0.0), fully remove the old machinery.

Putting it all together: an implementation sketch

Week 1–2:

Draft pyproject.toml with minimal metadata.

Migrate and validate dependencies.

Smoke-test pip install . in a clean venv.

Week 3:

Integrate CMake via scikit-build or setuptools-cmake.

Adjust tests to run under the new build_ext; verify wheels build.

Week 4:

Stand up GitHub Actions: build + test + lint on each PR.

Add type checking & code formatting enforcement.

Week 5:

Publish test releases to TestPyPI.

Revise README/docs; update installation instructions.

Week 6:

Deprecate old scripts with warnings.

Tag a v1.0.0 “py3-modernized” release.

## Next Steps and Considerations

- MyPy Configuration: The MyPy step is a placeholder. You'll need to configure MyPy (e.g., in pyproject.toml or mypy.ini) and specify the correct paths for type checking (e.g., mypy src/python).
- Pytest Integration: The Pytest step is also a placeholder. You need to ensure your tests are structured in a way that Pytest can discover and run them. This might involve migrating existing tests from shell scripts to Pytest. The command pytest will run tests if they follow standard Pytest discovery rules.
- Refining Dependency Installation: The step "Install Python build dependencies" currently installs build tools like setuptools, scikit-build, etc., explicitly. Then, pip install .[dev] installs the project. This is generally fine. The commented-out line pip install -r pyproject.toml is not standard; pip install . or pip install -e . (for editable installs) along with extras like .[dev] is the correct way to install from a pyproject.toml-based project.
- C++ Integration Tests: The current workflow doesn't explicitly run the C++ integration tests (e.g., test_cpp_integration.sh). These will need to be either converted to Pytest or run as a separate step in the CI.
- Publishing to PyPI: To enable publishing to PyPI, you'll need to:
- Uncomment the "Publish to PyPI" step.
- Create an API token on PyPI and add it as a secret (e.g., PYPI_API_TOKEN) to your GitHub repository settings.
- Ensure your pyproject.toml has all the necessary metadata for a PyPI release.

This document tracks the ongoing modernization of hap.py for Python 3, focusing on PEP 517/518 compliance, streamlined installation, reproducibility, and making releases consumable via standard Python tooling. This is a living document and should be updated as each phase is completed. Use planning, reasoning, and task management along with websearch to find latest documentation for dependencies. Making sure the package passes all pre-commit hooks.
Deep thinking and planning while keeping track of progress is required to ensure that the modernization process does not break existing functionality. The goal is to make the package more maintainable, easier to install, and compliant with modern Python packaging standards.
