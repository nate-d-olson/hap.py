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

## CLI Entry Points

The following console scripts are available after installation:
- hap.py: Haplotype benchmarking driver
- qfy: Quantification (was `som.py`)
- pre: Preprocessing VCF files
- ftx: Somatic VCF feature extraction
- cnx: Advanced chromosome-based count\ns
- ovc: Overlap-based VCF comparisons

Use `python -m happy.<cmd> --help` or the script name with `--help` to view usage. Smoke tests in `tests/test_cli.py` validate help output.

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

Refer to `.codex/plan_2025-05-19.md` for the high-level modernization plan and Task breakdown. Key upcoming milestones:
- Type-hint CLI modules and remove skip directives
- Add end-to-end example-based tests
- Enforce full mypy and Ruff coverage on all code
- Finalize packaging (move metadata to `pyproject.toml`)
- Release v0.5.0 with Python 3-only support

## Contact & Resources
- Code authors: Peter Krusche `<pkrusche@illumina.com>` and team.
- Issues and PRs: use GitHub issues for bugs and feature requests.
- Documentation: `README.md` and `doc/` folder.

_Happy coding!_
