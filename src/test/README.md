# Test Data Overview

This file documents the organization of test data used by the hap.py test suite.
All data required for tests lives under `tests/data`.

```
tests/data/
├── example/      # files originally from the repository 'example/' directory
├── src/          # files previously in 'src/data'
└── common/       # shared reference files used by multiple tests
```

Each subdirectory mirrors the previous layout so that existing tests can
reference data using helpers from `tests.utils`. The `example` directory is
used by tests like `test_blocksplit` and `test_giab`. Data under `src` is
used by many integration tests including `test_faulty_variants`,
`test_fp_accuracy`, `test_leftshift`, `test_multimerge`, and
`test_pathtraversal`.

The `common` directory contains files shared between multiple datasets. For
instance `test.fa` and its index are referenced by the `open_indel` and
`pathtraversal` tests.

When adding new tests, place any required data in an appropriate
subdirectory under `tests/data` and document its use here.
