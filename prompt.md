This session will refine the RTG vcfeval integration flags, indexing, and JSON metrics support,
improve CLI help and documentation, optimize performance, and expand test coverage for comparison engine features.

Refer to `.codex/plan_2025-06-02.md` under "Next: Comparison Engine Integration" for the detailed plan status
and outline of remaining tasks.

Key next steps:
- Expose vcfeval options (e.g., `--roc`, `--threads`, `--Xloose-match-distance`) through the CLI.
- Generate and index the annotated VCF output (`.vcf.gz` with tabix) automatically.
- Add support for writing JSON metrics (`--write-json`) alongside CSV summaries.
- Update CLI help text, README.md, and AGENTS.md to document new flags and workflow.
- Optimize performance: manage temp directories, improve file I/O, and profile hotspots.
- Expand integration tests for vcfeval-specific flags, ROC output parsing, and JSON metrics.

## Repository Overview

- Transition to pure-Python implementation and a PEP 517 build.
- CLI entry points: `hap.py`, `qfy`, `pre`.
- Shell scripts in `src/sh` for integration tests.
- Documentation in `doc/`, `README.md`, and `AGENTS.md`.
- Prepared plan structure in the `.codex/` directory.
- Development Environment
   - Python 3.7+ with `requirements-dev.txt`, `pre-commit`, and `nox` for testing, linting, and formatting.
   - Build is now pure-Python; CMake-based builds and scripts are deprecated.

## Next Steps for Plan Implementation

- Improve test coverage:
  - Add and enforce unit, integration, and coverage thresholds (e.g., 90%).
- Enforce type safety:
  - Add type annotations, enable mypy in CI.
- Optimize performance:
  - Implement microbenchmarks, profile hotspots, and address bottlenecks.
- Enhance documentation and packaging:
  - Separate core logic from CLI, add docstrings, generate Sphinx docs, migrate metadata to pyproject.toml.
- CI/CD and developer experience:
  - Add GitHub Actions workflows for tests, linting, formatting, benchmarks, and coverage, and automate releases.
- Migration of native extensions:
  - Draft and prioritize refactoring remaining C++/Cython modules to pure-Python or C bindings.

After completing the next step of the development plan, or hit a blocker that requires revising the plan,
update the README.md, AGENTS.md, and the plan document as appropriate.
Next then commit changes to the repository using `git commit` with an appropriately detailed commit message.
Finally report development status and what, if any,
additional development work is required to test the modernized version with a
whole genome example variant callset comparison.
