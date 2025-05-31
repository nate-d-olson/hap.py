Please proceed with the next task in the development plan.
Refer to `.codex/plan_2025-05-30.md` for detailed development plan.
Before starting the next step review the codebase, `AGENTS.md`, and the plan file. To they accurately describe the next logical step to complete the project objective for modernization the hap.py bioinformatics tool codebase while maintaining key functionality.

## Repository Overview

- Transition to pure-Python implementation and a PEP 517 build.
- CLI entry points: `hap.py`, `qfy`, `pre`.
- Shell scripts in `src/sh` for integration tests.
- Documentation in `doc/`, `README.md`, and `AGENTS.md`.
- Prepared plan structure in the `.codex/` directory.
- Development Environment
   - Python 3.7+ with `requirements-dev.txt`, `pre-commit`, and `nox` for testing, linting, and formatting.
   - Build is now pure-Python; CMake-based builds and scripts are deprecated.

## Recent Changes (2025-05-29)

- Updated `AGENTS.md` with refreshed upcoming milestones.
- Created `.codex/plan_2025-05-30.md` outlining objectives.

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

When you are done, update the README.md, AGENTS.md, and the plan document as appropriate then commit the changes to the repository.
