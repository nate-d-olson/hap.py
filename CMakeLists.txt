# Top-level Cmake file for Haplotype Tools (Python 3 version)
#
# Author: Peter Krusche <pkrusche@illumina.com>
# Python 3 Migration: [Your name]
#
cmake_minimum_required(VERSION 3.10)
project(HAPLOTYPES VERSION 0.3.15 LANGUAGES CXX C)

# Set scikit-build specific variables
set(SKBUILD_MODULE_NAME hap_py)
set(SKBUILD_PACKAGE_DIR src)

# Set CMAKE policies
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# Output directories
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include custom CMake modules
include(cxx)
include(configureFiles)
include(CythonSupport)

# Python version detection
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Found Python: ${Python3_EXECUTABLE} (version: ${Python3_VERSION})")
if(Python3_VERSION VERSION_LESS 3.7)
  message(FATAL_ERROR "Python 3.7 or later is required")
endif()

# NumPy detection and inclusion for Cython
find_package(NumPy REQUIRED)
include_directories(${NumPy_INCLUDE_DIRS})
message(STATUS "Found NumPy: ${NumPy_VERSION} at ${NumPy_INCLUDE_DIRS}")

# RTG tools / vcfeval configuration
if(BUILD_VCFEVAL)
    message(STATUS "Will build rtg-tools.")
    set(ENV{BUILD_VCFEVAL} "1")
    set(ENV{VCFEVAL_WRAPPER} ${VCFEVAL_WRAPPER})
    set(VCFEVAL_AVAILABLE 1)
else()
    set(VCFEVAL_AVAILABLE 0)
endif()

# External dependencies are now managed through modern CMake
# Check for pre-built dependencies or build them if needed
if(NOT EXISTS "${CMAKE_BINARY_DIR}/include/boost/version.hpp")
    message(STATUS "External dependencies not found. Please run: cmake --build build --target external_deps")
endif()

# Boost configuration
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

# Find system Boost or use local built version
if(DEFINED ENV{BOOST_ROOT})
    set(BOOST_ROOT $ENV{BOOST_ROOT})
    message(STATUS "Using system Boost from BOOST_ROOT environment: ${BOOST_ROOT}")
else()
    set(BOOST_ROOT ${CMAKE_BINARY_DIR})
    message(STATUS "Using built Boost at ${CMAKE_BINARY_DIR}")
endif()

# Find Boost components
find_package(Boost 1.55.0 COMPONENTS thread iostreams regex filesystem system program_options REQUIRED)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find Boost libraries. Please install or set BOOST_ROOT.")
endif()

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries: ${Boost_LIBRARIES}")

# Include directories
include_directories(${CMAKE_BINARY_DIR}/include)
link_directories(${CMAKE_BINARY_DIR}/lib)

# Find other dependencies
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

# Include subdirectories
add_subdirectory(src/hap_py/haplo/cython)

message(STATUS "Configuration completed successfully. Use 'pip install .' to build and install.")
