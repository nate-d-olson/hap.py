# Modern setup for Cython extension modules
import os

import numpy
from setuptools import Extension, setup

# Try to import Cython, fall back to no extensions if not available
try:
    from Cython.Build import cythonize

    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False
    print("WARNING: Cython not available, installing without C extensions")


# Find include directories
def get_include_dirs():
    include_dirs = [
        numpy.get_include(),
        os.path.join(os.path.dirname(__file__), "src/c++/include"),
    ]
    return include_dirs


# Find library directories
def get_library_dirs():
    # Get build directory from command line arguments or use default
    build_dir = os.environ.get(
        "HAPPY_BUILD_DIR", os.path.join(os.path.dirname(__file__), "build")
    )
    library_dirs = [
        os.path.join(build_dir, "lib"),
    ]
    return library_dirs


# Configure extensions with proper Python 3 support
ext_modules = []
if HAVE_CYTHON:
    extensions = [
        Extension(
            "happy.Haplo.happyroc",
            ["src/python/Haplo/happyroc.pyx"],
            include_dirs=get_include_dirs(),
            libraries=["haplotypes"],
            library_dirs=get_library_dirs(),
            extra_compile_args=["-std=c++11"],
            language="c++",
        ),
        Extension(
            "happy.Haplo.variant_processor",
            ["src/python/Haplo/variant_processor.pyx"],
            include_dirs=get_include_dirs(),
            libraries=["haplotypes"],
            library_dirs=get_library_dirs(),
            extra_compile_args=["-std=c++11"],
            language="c++",
        ),
        Extension(
            "happy.Haplo.sequence_utils",
            ["src/python/Haplo/sequence_utils.pyx"],
            include_dirs=get_include_dirs(),
            libraries=["haplotypes"],
            library_dirs=get_library_dirs(),
            extra_compile_args=["-std=c++11"],
            language="c++",
        ),
    ]

    # Use Cython compiler directives for Python 3 compatibility
    cython_directives = {
        "language_level": "3",  # Python 3 mode
        "boundscheck": False,  # Turn off bounds-checking
        "wraparound": False,  # Turn off negative indexing
    }

    ext_modules = cythonize(extensions, compiler_directives=cython_directives)

setup(
    name="happy",
    version="0.4.0",
    description="Haplotype Comparison Tools",
    packages=["happy", "happy.Haplo", "happy.Tools"],
    package_dir={"happy": "src/python"},
    ext_modules=ext_modules,
    install_requires=[
        "numpy>=1.15.0",
        "pysam>=0.15.0",
        "scipy>=1.0.0",
        "pandas>=0.23.0",
        "cython>=0.29.0",
    ],
    python_requires=">=3.6",
)
