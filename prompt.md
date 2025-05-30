Please proceed with the next task in the development plan. Review the codebase, `AGENTS.md`, and the plan file to align with project objectives.

## Repository Overview
- `happy`: Haplotype benchmarking tool for germline variant comparison.
- Transition to pure-Python implementation and a PEP 517 build.
- CLI entry points: `hap.py`, `qfy`, `pre`.
- Shell scripts in `src/sh` for integration tests.
- Documentation in `doc/`, `README.md`, and `AGENTS.md`.

## Recent Changes (2025-05-29)
- Updated `AGENTS.md` with refreshed upcoming milestones.
- Created `.codex/plan_2025-05-30.md` outlining objectives.
- Prepared plan structure in the `.codex/` directory.

## Plan for Next Session
Refer to `.codex/plan_2025-05-30.md` for detailed objectives for today's session.

## Development Environment
- Python 3.7+ with `requirements-dev.txt`, `pre-commit`, and `nox` for testing, linting, and formatting.
- Optional C++ components: use `./configure.sh` and CMake.

## Next Steps
1. Migrate the `sequence_utils` Cython extension to pure-Python.
2. Add unit tests for `sequence_utils` and validate functionality.
3. Run the full test suite (`pytest`) to ensure no regressions.
4. Update documentation (`AGENTS.md`, `README.md`) to reflect completed migrations.
5. Review remaining Cython/C++ modules and plan their migration.
