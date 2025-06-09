# hap.py Test Suite

This directory contains tests for the hap.py project. The tests are organized into two main directories:

- `unit/`: Unit tests for individual components
- `integration/`: Integration tests that verify the interaction between components

## Running Tests

Before running the test suite it is recommended to activate the conda environment
defined in `environment-dev.yml`:

```bash
conda env create -f environment-dev.yml    # one-time setup
conda activate happy-dev
```

Example:

```bash
conda activate happy-dev && pytest
```

Tests can then be run using pytest:

```bash
# Run all tests
pytest

# Run only unit tests
pytest tests/unit

# Run only integration tests
pytest tests/integration

# Run tests with specific markers
pytest -m "not integration"  # Skip integration tests
pytest -m "not cpp"          # Skip tests that require C++ components
pytest -m "not slow"         # Skip slow tests
```

## Test Markers

The tests use the following markers:

- `integration`: Marks tests as integration tests
- `cpp`: Marks tests that require C++ components
- `slow`: Marks tests that take more than a few seconds
- `external`: Marks tests that require external dependencies

## Test Utilities

Common utilities for tests are available in `tests/utils.py`.

## Adding New Tests

When adding new tests:

1. Follow the existing test organization
2. Use appropriate markers
3. Add docstrings explaining what each test does
4. Use the utilities in `tests/utils.py`

## Migrating Shell Tests

Shell tests are being migrated from `src/sh` to pytest tests. A helper script is available:

```bash
python scripts/migrate_test.py src/sh/run_test_name.sh tests/integration/test_name.py
```

This will generate a template pytest file that you can then complete.

## Test Environment

Tests assume that:

1. The hap.py package is installed or available in the Python path
2. C++ components have been built (for tests with the `cpp` marker)
3. A reference genome is available (either via `HGREF` environment variable or in the example directory)
4. `bgzip` and `tabix` executables are available in `build/bin` or on the `PATH`. If not, the helper functions in `tests/utils.py` fall back to `pysam` for compression and indexing.
5. The `rtg` executable from RTG Tools is available. Set the `RTG` or
   `RTGTOOLS_PATH` environment variable to its location if it is not on
   your `PATH`.

### Building the Project Before Running Tests

Before running tests, especially integration tests, you must build the C++ components. From the project root run:

```bash
cmake -B build -S .
cmake --build build
```

This places the compiled binaries in `build/bin`. Certain tests rely on these executables, so ensure they are built before running:

```bash
pytest
```

If the binaries are missing you may see errors like:

```
AssertionError: hap.py failed with error: No such file or directory
```

## GitHub Actions CI

The CI pipeline builds the project and runs the tests in the GitHub Actions environment. See the `.github/workflows` directory for the configuration.
