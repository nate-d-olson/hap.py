---
applyTo: "**"
---
# Cython Modernization Guidelines for Python 3

## Cython Extension Modules

### Build System Updates

1. **Modernized setup.py/pyproject.toml**
   ```python
   # Modern setup.py for Cython
   from setuptools import setup, Extension
   from Cython.Build import cythonize
   import numpy

   extensions = [
       Extension(
           "happy.cython_module",
           ["src/python/happy/cython_module.pyx"],
           include_dirs=[numpy.get_include(), "src/c++/include"],
           libraries=["haplotypes"],
           library_dirs=["build/lib"],
           extra_compile_args=["-std=c++11"],
       )
   ]

   setup(
       # ...other setup parameters...
       ext_modules=cythonize(extensions),
   )
   ```

2. **CMake Integration**
   ```cmake
   # CythonSupport.cmake usage
   add_cython_module(
       cython_module
       src/python/happy/cython_module.pyx
       LIBRARIES
           haplotypes
           ${HAPLOTYPES_ALL_LIBS}
       INCLUDES
           ${NUMPY_INCLUDE_DIRS}
           ${PYTHON_INCLUDE_DIRS}
   )
   ```

### Language Level Updates

1. **Set Appropriate Language Level**
   ```cython
   # cython: language_level=3
   # distutils: language=c++

   # Python 3 imports
   from __future__ import annotations  # Python 3.7+ forward annotation reference
   ```

2. **Unicode vs Bytes**
   ```cython
   # Python 2 style
   cdef char* get_seq(str python_str):
       return python_str

   # Python 3 style
   cdef char* get_seq(str python_str):
       py_bytes = python_str.encode('utf8')
       return py_bytes
   ```

## Memory Management

### Reference Counting

1. **Ownership Management**
   ```cython
   # Safety: explicitly keep references to Python objects when passing to C++
   cdef class VariantWrapper:
       cdef Variant* c_variant
       cdef object _owner  # Keeps reference to owner object

       def __cinit__(self, owner):
           self._owner = owner
           self.c_variant = get_variant_from_owner(owner)
   ```

2. **Temporary References**
   ```cython
   # Bad: temporary reference
   cdef Variant* var = convert_to_variant(some_function())

   # Good: keep reference
   temp_obj = some_function()
   cdef Variant* var = convert_to_variant(temp_obj)
   # Use var
   # temp_obj reference kept until no longer needed
   ```

### Memory Allocation

1. **Safe Memory Allocation**
   ```cython
   from libc.stdlib cimport malloc, free
   from libc.string cimport memcpy

   cdef char* safe_alloc(size_t size) except NULL:
       cdef char* buffer = <char*>malloc(size)
       if buffer == NULL:
           raise MemoryError("Failed to allocate memory")
       return buffer
   ```

2. **Exception-Safe Resource Management**
   ```cython
   cdef class Buffer:
       cdef char* data
       cdef size_t size

       def __cinit__(self, size_t size):
           self.size = size
           self.data = <char*>malloc(size)
           if self.data == NULL:
               raise MemoryError()

       def __dealloc__(self):
           if self.data != NULL:
               free(self.data)
               self.data = NULL
   ```

## String Handling

### UTF-8 Encoding

1. **String Conversion**
   ```cython
   # String to bytes
   def py_to_c_string(s):
       if isinstance(s, str):
           return s.encode('utf8')
       return s  # Assume already bytes

   # C string to Python string
   cdef str c_to_py_string(const char* c_str):
       return c_str.decode('utf8') if c_str != NULL else ""
   ```

2. **File Path Handling**
   ```cython
   # Python 3 style path handling
   cdef class FileHandler:
       cdef FILE* c_file

       def __cinit__(self, path):
           py_bytes = str(path).encode('utf8')
           self.c_file = fopen(py_bytes, "rb")
           if self.c_file == NULL:
               raise IOError(f"Could not open file {path}")

       def __dealloc__(self):
           if self.c_file != NULL:
               fclose(self.c_file)
   ```

### String Literals

1. **Python 3 String Literals**
   ```cython
   # Python 2
   py_str = 'string'
   py_unicode = u'unicode'

   # Python 3
   py_str = 'string'  # Unicode by default
   py_bytes = b'bytes'
   ```

## Type Handling

### Type Annotations

1. **Type Declarations**
   ```cython
   # Basic types
   cdef int i = 0
   cdef double d = 0.0
   cdef bint flag = True

   # Complex types
   from libcpp.vector cimport vector
   from libcpp.string cimport string

   cdef vector[int] vec
   cdef string cpp_str
   ```

2. **Using NumPy Arrays**
   ```cython
   # cython: language_level=3
   import numpy as np
   cimport numpy as np

   np.import_array()  # Initialize NumPy C API

   def process_array(np.ndarray[np.float64_t, ndim=2] arr not None):
       cdef Py_ssize_t i, j
       cdef double value
       cdef np.ndarray[np.float64_t, ndim=2] result = np.zeros_like(arr)

       for i in range(arr.shape[0]):
           for j in range(arr.shape[1]):
               value = arr[i, j]
               result[i, j] = process_value(value)

       return result
   ```

### Python Objects

1. **Object Handling**
   ```cython
   cdef process_variant(object variant):
       # Check for proper type
       if not isinstance(variant, Variant):
           raise TypeError("Expected Variant object")

       # Extract C++ pointer from Python object
       cdef Variant* c_var = get_c_variant(variant)
   ```

2. **Property Access**
   ```cython
   # Python 2 style
   cdef class GenomicInterval:
       cdef public str chrom
       cdef public int start, end

   # Python 3 style - same but with type annotations
   cdef class GenomicInterval:
       cdef public str chrom
       cdef public int start, end

       def __init__(self, str chrom, int start, int end):
           self.chrom = chrom
           self.start = start
           self.end = end
   ```

## Bioinformatics-Specific Patterns

### VCF Record Processing

1. **Variant Record Handling**
   ```cython
   # Efficient variant processing
   cdef class VariantProcessor:
       cdef vector[Variant*] c_variants

       def add_variant(self, variant):
           # Convert Python variant to C++ variant
           cdef Variant* c_variant = new Variant()

           # Copy data from Python to C++
           py_bytes = variant.chrom.encode('utf8')
           c_variant.chrom = string(py_bytes)
           c_variant.position = variant.position

           # Store in vector
           self.c_variants.push_back(c_variant)

       def __dealloc__(self):
           # Clean up C++ objects
           for i in range(self.c_variants.size()):
               if self.c_variants[i] != NULL:
                   del self.c_variants[i]
   ```

2. **FASTA/Sequence Processing**
   ```cython
   # Fast sequence handling
   cdef bytes complement_sequence(bytes seq):
       cdef:
           size_t i, n = len(seq)
           char* c_seq = seq
           char* result = <char*>malloc(n + 1)

       if result == NULL:
           raise MemoryError()

       try:
           for i in range(n):
               if c_seq[i] == 'A':
                   result[i] = 'T'
               elif c_seq[i] == 'T':
                   result[i] = 'A'
               elif c_seq[i] == 'G':
                   result[i] = 'C'
               elif c_seq[i] == 'C':
                   result[i] = 'G'
               else:
                   result[i] = c_seq[i]

           result[n] = 0  # Null terminator
           return result[:n]  # Convert back to Python bytes
       finally:
           free(result)
   ```

### Thread Safety

1. **GIL Management**
   ```cython
   def process_variants(variants, threads=1):
       # Prepare for parallel processing
       cdef:
           size_t n_variants = len(variants)
           vector[Variant*] c_variants = prepare_variants(variants)
           vector[Result*] c_results

       with nogil:  # Release GIL for true parallelism
           c_results = process_variants_parallel(c_variants, threads)

       # Convert results back to Python objects
       return [create_py_result(c_results[i]) for i in range(c_results.size())]
   ```

2. **Thread Pool**
   ```cython
   from cython.parallel import prange

   def parallel_process(np.ndarray[np.float64_t, ndim=2] data):
       cdef:
           Py_ssize_t i, j, rows, cols
           np.ndarray[np.float64_t, ndim=2] result

       rows = data.shape[0]
       cols = data.shape[1]
       result = np.zeros((rows, cols), dtype=np.float64)

       # Process in parallel with OpenMP
       for i in prange(rows, nogil=True, schedule='static'):
           for j in range(cols):
               result[i, j] = process_value_nogil(data[i, j])

       return result
   ```

## Testing Cython Modules

### Unit Testing

1. **Test Python Interface**
   ```python
   def test_variant_processor():
       # Create test data
       variants = [
           Variant("chr1", 100, "A", "T"),
           Variant("chr1", 200, "G", "C")
       ]

       # Test Cython module
       processor = VariantProcessor()
       for v in variants:
           processor.add_variant(v)

       results = processor.process_all()
       assert len(results) == 2
       assert results[0].chrom == "chr1"
       assert results[1].position == 200
   ```

2. **Memory Leak Testing**
   ```python
   def test_no_memory_leaks():
       import gc
       import psutil
       import os

       process = psutil.Process(os.getpid())

       # Get initial memory usage
       gc.collect()
       initial = process.memory_info().rss

       # Run computation many times
       for _ in range(1000):
           large_data = create_large_dataset()
           result = cython_process_function(large_data)
           del large_data
           del result

       # Check final memory usage
       gc.collect()
       final = process.memory_info().rss

       # Assert memory growth is minimal
       assert (final - initial) < 1024 * 1024  # Less than 1MB growth
   ```

### Integration Testing

1. **End-to-End Testing**
   ```python
   def test_variant_comparison_pipeline():
       # Set up test data
       truth_vcf = "test_data/truth.vcf"
       query_vcf = "test_data/query.vcf"

       # Run full pipeline
       result = run_comparison(truth_vcf, query_vcf)

       # Validate outputs
       assert os.path.exists(result.summary_csv)
       assert result.precision > 0.9
       assert result.recall > 0.9
   ```

## Migration Checklist

1. **Build System Updates**
   - [ ] Update pyproject.toml for Cython modules
   - [ ] Configure Language level to Python 3 in all .pyx files
   - [ ] Update CMake integration for Cython

2. **Code Updates**
   - [ ] Fix string handling (Unicode vs bytes)
   - [ ] Update memory management patterns
   - [ ] Fix object reference management
   - [ ] Update NumPy array handling

3. **Testing**
   - [ ] Add unit tests for all Cython modules
   - [ ] Test with various Python 3 versions
   - [ ] Verify no memory leaks
   - [ ] Benchmark performance against Python 2 version
