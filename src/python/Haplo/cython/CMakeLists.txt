# CMakeLists.txt for building Cython modules in Python 3
# This file should be placed in the src/python/Haplo/cython directory

cmake_minimum_required(VERSION 3.10)

# Find Python3 and its components
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
include(${CMAKE_SOURCE_DIR}/src/cmake/CythonSupport.cmake)

# Include directories for C++ headers
include_directories(
    ${CMAKE_SOURCE_DIR}/src/c++/include
    ${CMAKE_BINARY_DIR}/include
)

# Add the Cython module
add_cython_module(
    internal
    _internal.pyx
    LIBRARIES
        haplotypes
        ${HAPLOTYPES_ALL_LIBS}
    INCLUDES
        ${CMAKE_SOURCE_DIR}/src/c++/include 
        ${CMAKE_BINARY_DIR}/include
)

# Install the Python module files
install(FILES 
    __init__.py 
    DESTINATION ${CMAKE_BINARY_DIR}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/Haplo/cython
)
