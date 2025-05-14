---
applyTo: "**"
---

# Python 3 Migration Guidelines with C++ Integration

## Current Status (Updated May 14, 2025)

- Total Python files: 48
- Fully migrated files: 21 (43.8%)
- Partially migrated files: 27 (56.2%)
- Total remaining issues: 90

## Priority Areas

1. Haplo module (32 issues) - Highest priority
2. Tools module (25 issues)
3. Somatic module (12 issues)
4. Core files (hap.py, som.py, qfy.py) - 13 issues total

## Migration Workflow

3. **Parallel Development Strategy**
   - While addressing C++ build issues:
     - Work on standalone Python modules that don't depend on C++ components
     - Update unit tests to support Python 3
     - Create Python 3 versions of utility functions
   - After C++ builds are fixed:
     - Update Python/C++ integration points
     - Test full system with Python 3

4. **Testing Approach**
   - Create mocks for C++ components during Python migration
   - Test Python components independently when possible
   - Verify Python/C++ integration with small test cases
   - Run full test suite once both systems are working

## Python-C++ Integration Points

### Handling Extension Modules
- Update Python C extension code carefully:

```python
# In Python 2, this might work:
from happy_ext import SomeExtension

# In Python 3, ensure extensions are properly built and imported:
try:
  from happy_ext import SomeExtension
except ImportError:
  # Fallback to pure Python implementation or show clear error
  raise ImportError("C++ extensions failed to build. Please check build logs.")
```

### Mock Objects for Testing

Create Python mock objects for C++ components during development:

```python
# Mock of C++ extension class for testing
class MockVariantProcessor:
    def __init__(self):
        self.variants = []
        
    def add_variant(self, variant):
        self.variants.append(variant)
        
    def process(self):
        # Simple Python implementation of C++ algorithm for testing
        return [self._process_one(v) for v in self.variants]
        
    def _process_one(self, variant):
        # Simplified processing logic
        return {"chrom": variant.chrom, "processed": True}
```

## Python 3 Migration Steps for Core Components

1. Tools Module Updates
  * Fix src/python/Tools/__init__.py first
  * Update src/python/Tools/bcftools.py for Python 3 string handling
  * Create tests for each utility function

2. Data Processing Components
  * Update `src/python/pre.py` to handle file paths and encoding properly
  * Fix string vs bytes handling in preprocessing functions
  * Ensure compatibility with both VCF and BCF formats

3. Main Scripts
  * Update `src/python/qfy.py` after Tools and pre.py
  * Update `src/python/hap.py` last after all dependencies are fixed
  * Handle command-line arguments and file path encoding

## Genomic Data-Specific Considerations

### VCF/BCF File Handling
- Replace direct file reads with proper encoding-aware operations
  ```python
  # Python 2
  with open(filename) as f:
      header = f.readline()
      
  # Python 3
  with open(filename, 'rt', encoding='utf-8') as f:
      header = f.readline()
  ```

- Update binary file handling for BCF files
  ```python
  # Python 2
  with open(filename, 'rb') as f:
      magic = f.read(3)
  
  # Python 3
  with open(filename, 'rb') as f:
      magic = f.read(3)  # Returns bytes in Python 3
      if magic == b'BCF':  # Use bytes literal
          # Process BCF file
  ```

### Sequence Data Processing

- Fix string vs bytes handling in sequence data
  ```python
  # Python 2
  def reverse_complement(seq):
      comp_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
      return ''.join([comp_dict[base] for base in reversed(seq)])
  
  # Python 3
  def reverse_complement(seq):
      # Handle both string and bytes
      is_bytes = isinstance(seq, bytes)
      if is_bytes:
          seq = seq.decode('ascii')
      
      comp_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
      result = ''.join([comp_dict[base] for base in reversed(seq)])
      
      if is_bytes:
          return result.encode('ascii')
      return result
  ```

### Sorting and Comparison
- Update chromosome sorting logic for Python 3
  ```python
  # Python 2
  chrom_list.sort(lambda a, b: cmp_chromosomes(a, b))
  
  # Python 3
  from functools import cmp_to_key
  chrom_list.sort(key=cmp_to_key(cmp_chromosomes))
  ```

### Binary Data Structures
- Update struct packing/unpacking for genomic coordinate data
  ```python
  # Python 2
  import struct
  def pack_position(chrom_id, position):
      return struct.pack("iI", chrom_id, position)
  
  # Python 3
  import struct
  def pack_position(chrom_id, position):
      return struct.pack("iI", chrom_id, position)  # Unchanged, but returns bytes
  ```

### Hashing and Equality
- Update hash-based collections with proper Python 3 behavior
  ```python
  # Python 2
  variant_dict = {variant: count for variant, count in variants}
  
  # Python 3
  # Ensure __hash__ and __eq__ methods are properly implemented
  class Variant:
      def __hash__(self):
          return hash((self.chrom, self.pos, self.ref, self.alt))
      
      def __eq__(self, other):
          if not isinstance(other, Variant):
              return False
          return (self.chrom == other.chrom and 
                  self.pos == other.pos and
                  self.ref == other.ref and
                  self.alt == other.alt)
  ```

## pysam-Specific Updates

### Updates for pysam ≥0.15
- Replace deprecated alignment methods
  ```python
  # Python 2 with older pysam
  for read in samfile.fetch():
      if read.is_proper_pair and not read.is_unmapped:
          tid = read.tid
  
  # Python 3 with newer pysam
  for read in samfile.fetch():
      if read.is_proper_pair and not read.is_unmapped:
          tid = read.reference_id  # tid is deprecated
  ```

- Update VCF handling for pysam API changes
  ```python
  # Python 2 with older pysam
  bcf_in = pysam.VariantFile(filename)
  for rec in bcf_in:
      sample = rec.samples[sample_name]
      gt = sample['GT']
  
  # Python 3 with newer pysam
  bcf_in = pysam.VariantFile(filename)
  for rec in bcf_in:
      sample = rec.samples[sample_name]
      gt = sample['GT']  # Same API, but returns NumPy array-like in newer versions
      gt_list = list(gt)  # Convert to regular list if needed
  ```

### Handling File Objects
- Update file context management with pysam
  ```python
  # Python 2
  samfile = pysam.AlignmentFile(filename, "rb")
  try:
      # process samfile
  finally:
      samfile.close()
  
  # Python 3
  with pysam.AlignmentFile(filename, "rb") as samfile:
      # process samfile
  ```

## Cython Integration

### Memory Management
- Update memory allocation patterns in Cython
  ```python
  # Python 2/Cython
  cdef char* c_string = <char*>malloc(length + 1)
  # ...use c_string...
  free(c_string)
  
  # Python 3/Cython
  from libc.stdlib cimport malloc, free
  cdef char* c_string = <char*>malloc(length + 1)
  # ...use c_string...
  free(c_string)
  ```

### Object Lifetime
- Fix object lifetime issues in Python 3
  ```cython
  # Python 2/Cython
  cdef object py_obj
  py_obj = get_python_object()
  c_obj = convert_to_c(py_obj)
  # py_obj might be garbage collected here, invalidating c_obj
  
  # Python 3/Cython
  cdef object py_obj
  py_obj = get_python_object()
  c_obj = convert_to_c(py_obj)
  # Keep reference until no longer needed
  some_function(c_obj)
  # Now py_obj can be garbage collected
  ```

### String Handling in Cython
- Update string handling in Cython for Python 3
  ```cython
  # Python 2/Cython
  cdef char* c_str = python_str
  
  # Python 3/Cython
  # For Python str -> C char*
  py_bytes = python_str.encode('utf8')
  cdef char* c_str = py_bytes
  
  # For C char* -> Python str
  cdef char* c_str = get_c_string()
  py_str = c_str.decode('utf8')
  ```

## Common Issues & Solutions

### Iterator Changes
- Update iteration patterns for Python 3
  ```python
  # Python 2
  vcf_records = list(vcf_file)  # Materializes all records
  vcf_it = iter(vcf_file)
  while True:
      try:
          rec = vcf_it.next()  # .next() method
      except StopIteration:
          break
  
  # Python 3
  vcf_records = list(vcf_file)  # Still works, but consider generator approach
  vcf_it = iter(vcf_file)
  while True:
      try:
          rec = next(vcf_it)  # Built-in next() function
      except StopIteration:
          break
  ```

### Metadata Processing
- Update metadata handling for Python 3 compatibility
  ```python
  # Python 2
  vcf_meta = {str(k): str(v) for k, v in metadata.iteritems()}
  
  # Python 3
  vcf_meta = {str(k): str(v) for k, v in metadata.items()}
  ```

### Integer Division
- Fix integer division in coordinate calculations
  ```python
  # Python 2
  mid_point = (start + end) / 2  # Integer division
  
  # Python 3
  mid_point = (start + end) // 2  # Explicit integer division
  ```

## Testing Strategies for Bioinformatics Data

1. **Sequence Equivalence Tests**
   - Compare sequence strings byte-by-byte
   - Verify correctness of complement/reverse operations

2. **Coordinate Validation**
   - Test genomic coordinate calculations
   - Verify correct sorting of chromosomes

3. **VCF/BCF Consistency**
   - Compare VCF records before and after migration
   - Check genotype encoding/decoding

4. **Statistical Output Validation**
   - Compare numerical results with tolerance for floating-point
   - Verify statistical calculations match between versions

5. **Memory Usage Monitoring**
   - Test with incrementally larger datasets
   - Verify no memory leaks in Python 3 version
