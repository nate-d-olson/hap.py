# Agents & Onboarding Guide

This document provides guidance for new contributors and automation agents working on the `happy` benchmarking tool.

## Workspace Setup

1. Clone the repository:
   ```bash
   git clone https://<repo-url>.git
   cd hap.py
   ```
2. Create and activate the dev environment:
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip setuptools wheel
   pip install -r requirements-dev.txt
   pre-commit install
   ```
   This installs pre-commit hooks, pytest, and other dev tooling.

## Pre-commit Checks

Run all hooks before pushing:
```bash
pre-commit run --all-files
```
This includes formatting (Black, isort), linting (Ruff), syntax upgrades (pyupgrade), and type checks (mypy scoped to core modules).

## Build & Installation

Install the Python package:
```bash
pip install -e .
# or for a regular install:
pip install .
```


Use the provided `nox` sessions, which install the package before testing or linting.

## Configuration Files

Key project configuration files:
- pyproject.toml: build-system requirements and project metadata
- setup.py / setup.cfg: package definitions and flake8 settings
- requirements-dev.txt: development dependencies for formatting, linting, and testing
- .pre-commit-config.yaml: pre-commit hooks configuration
- noxfile.py: reproducible sessions for linting, formatting, type checking, and tests
- Jenkinsfile: CI pipeline for automated builds and tests
 - Dockerfile and .dockerignore: Docker image setup for development or CI
- .codex/plan_*.md: Dated plan files for tracking progress

## CLI Entry Points

The following console scripts are available after installation:
- hap.py: Haplotype benchmarking driver
- qfy: Quantification driver
- pre: VCF preprocessing tool

Use `python -m happy.<cmd> --help` or the script name with `--help` to view usage. Smoke tests in `tests/test_cli.py` validate help output.

## Source Code Layout

The repository follows a multi-language layout under the `src` directory:
- src/python: core Python modules
  - Haplo: benchmarking and variant comparison logic
  - Tools: helper utilities (VCF parsing, metrics, etc.)
  - happy: CLI entry points (hap.py, qfy, and pre)
- src/c++: native C++ libraries and tools
- src/sh: shell scripts for integration tests and wrappers (run_tests.sh, rtg-wrapper.sh, etc.)
- src/R: R scripts for additional analyses and reports
- src/data: reference and test datasets used by examples and tests
Examples and integration test scenarios are provided in the `example/` directory.

## Testing

Unit tests and smoke tests are in `tests/`. Run:
```bash
pytest -q
```
(requires pytest installation in the env)

## Reproducible Environments

In addition to the standard setup, we provide `nox` sessions for reproducible linting, formatting, type checking, and testing environments. After activating the dev environment, you can use:

```bash
nox -s lint
nox -s format
nox -s type_check
nox -s tests
```

## Development Plan

While working on the project, track your milestones in dated plan files under `.codex/` (e.g., `.codex/plan_YYYY-MM-DD.md`) and update them as you make progress. For significant completed work, update `README.md` with a dated changelog entry referencing relevant documentation.

Key upcoming milestones:
- Continue dropping remaining Cython/C++ extension dependencies.
- Advance packaging to a PEP 517-only build system.
- Add unit and integration tests for pure-Python modules.
- Refine CLI scripts and expand test coverage.

## Contact & Resources
- Issues and PRs: use GitHub issues for bugs and feature requests.
- Documentation: `README.md` and `doc/` folder.

_Happy coding!_
