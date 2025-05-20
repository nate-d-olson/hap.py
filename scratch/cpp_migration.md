# C++ Module Migration Plan

This document tracks the migration of C++/Cython modules in `happy.Haplo` to pure Python
implementations or to rely on performant Python libraries (e.g., NumPy, pysam).

## Modules to Migrate

- **happyroc**: ROC curve processing (_happyroc.pyx_ / C++ `Roc.hh`).
- **variant_processor**: Variant-level processing (_variant_processor.pyx_ / C++ `variant/*.hh`).
- **sequence_utils**: Sequence handling (_sequence_utils.pyx_ / C++ helpers).

## Migration Phases

1. **happyroc Migration**
   - Replace Cython/C++ ROC extension with the existing PythonRocCurve class.
   - Remove `.pyx` and extension build for `happy.Haplo.happyroc`.
   - Drop C++ `Roc.hh` and corresponding implementation after verifying nothing else depends on it.
   - Ensure `create_roc_curve()` yields the pure-Python implementation by default.

2. **variant_processor Migration** (future)
   - Analyze functionality in C++ `variant_processor.pyx` and `variant/*.hh`.
   - Replace VCF-level operations with pysam-based Python code if possible.
   - Port remaining logic to Python, using NumPy/Pandas as needed.

3. **sequence_utils Migration** (future)
   - Review sequence operations in C++ extension.
   - Migrate to Python using Biopython/pysam or NumPy vectorization where appropriate.

4. **Cleanup**
   - Remove C++ include files, adjust setup.py to drop Cython extensions.
   - Update documentation and tests to use pure-Python implementations.

---

**Next Steps:**
- Begin Phase 1: migrate `happyroc` to pure Python. Document progress below.

### Phase 1 Progress

- [x] Add pure-Python `happyroc.py` module in `happy.Haplo`.
- [x] Update setup.py to remove `happy.Haplo.happyroc` Cython extension.
- [x] Delete `src/python/Haplo/happyroc.pyx` and Cython-generated `.cpp` file.
- [x] Remove `src/c++/include/helpers/Roc.hh` and `src/c++/lib/tools/Roc.cpp` if unused.
- [x] Run tests to ensure `create_roc_curve()` and fallback PythonRocCurve work as expected.

### Phase 2 Progress

- [x] Analyze functionality in C++ `variant_processor.pyx` and `variant/*.hh`.
- [x] Implement pure-Python `src/python/Haplo/variant_processor.py` using pysam and NumPy/Pandas.
- [x] Update setup.py to remove `happy.Haplo.variant_processor` Cython extension.
- [x] Delete `src/python/Haplo/variant_processor.pyx` and generated `variant_processor.cpp`.
- [ ] Update tests to cover new VariantProcessor functionality.
- [ ] Ensure tests for Cython modernization pass using the new implementation.
- [ ] Benchmark performance against Cython version and document results.

### Phase 3 Progress

- [x] Analyze functionality in C++ `sequence_utils.pyx` and its helper headers.
- [x] Implement pure-Python `src/python/Haplo/sequence_utils.py` using NumPy or pure Python.
- [x] Update setup.py to remove `happy.Haplo.sequence_utils` Cython extension.
- [x] Delete `src/python/Haplo/sequence_utils.pyx` and `src/python/Haplo/sequence_utils.cpp`.
- [ ] Update tests to cover new sequence_utils functionality.
- [ ] Benchmark performance of the Python implementation and document results.
