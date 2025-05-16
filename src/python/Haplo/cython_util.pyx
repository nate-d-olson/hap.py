# cython: language_level=3
# distutils: language=c++
# Helper functions for string encoding/decoding and safe C++ interactions
from libc.stdlib cimport malloc, free
from libc.stdio cimport FILE, fopen, fclose

cdef bytes _ensure_bytes(s):
    """Convert a Python string to bytes."""
    if isinstance(s, str):
        return s.encode('utf8')
    return s  # Assume already bytes

cdef str _ensure_str(const char* c_str):
    """Convert a C string to Python string."""
    if c_str == NULL:
        return ""
    return c_str.decode('utf8')

# String handling for file paths
cdef class FileHandler:
    cdef FILE* c_file

    def __cinit__(self, path):
        py_bytes = _ensure_bytes(str(path))
        self.c_file = fopen(py_bytes, "rb")
        if self.c_file == NULL:
            raise IOError(f"Could not open file {path}")

    def __dealloc__(self):
        if self.c_file != NULL:
            fclose(self.c_file)

# Safe memory allocation
cdef char* safe_alloc(size_t size) except NULL:
    cdef char* buffer = <char*>malloc(size)
    if buffer == NULL:
        raise MemoryError("Failed to allocate memory")
    return buffer

# Exception-safe resource management
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
            
# Example of proper reference handling
cdef class VariantWrapper:
    cdef Variant* c_variant
    cdef object _owner  # Keeps reference to owner object

    def __cinit__(self, owner):
        self._owner = owner  # Store reference to prevent garbage collection
        self.c_variant = get_variant_from_owner(owner)
        
    def __dealloc__(self):
        # No need to free c_variant as it's owned by _owner
        pass