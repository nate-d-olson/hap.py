# Top-level Cmake file for Haplotype Tools (Python 3 version)
#
# Author: Peter Krusche <pkrusche@illumina.com>
# Python 3 Migration: [Your name]
#
cmake_minimum_required(VERSION 3.10)
project(HAPLOTYPES VERSION 0.3.15 LANGUAGES CXX C)

# Set scikit-build specific variables
# This helps scikit-build find the Python package directory
set(SKBUILD_MODULE_NAME hap_py) # Changed from Haplo
set(SKBUILD_PACKAGE_DIR src)    # Changed from src/python, as pyproject.toml now specifies package_dir = {"": "src"} essentially via packages = [{include = "hap_py", from = "src"}]

# Set CMAKE policies
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)  # Use <PackageName>_ROOT variables
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
if(Python3_VERSION VERSION_LESS 3.7) # Updated to 3.7 as per pyproject.toml
  message(FATAL_ERROR "Python 3.7 or later is required")
endif()

# NumPy detection and inclusion for Cython
# scikit-build handles this more seamlessly, but explicit find_package can be useful
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

# Build external dependencies
# This step might need to be adapted or removed if scikit-build handles dependencies differently
# For now, assume it's still needed for non-Python dependencies
execute_process(
    COMMAND ${CMAKE_SOURCE_DIR}/external/make_dependencies.sh # Assuming this script is Python 3 compatible or updated
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    RESULT_VARIABLE EXTERNAL_SUCCESS)

if(NOT "${EXTERNAL_SUCCESS}" STREQUAL "0")
    message(FATAL_ERROR "Building external dependencies has failed. See error messages above.")
endif()

# Boost configuration
set(Boost_USE_STATIC_LIBS ON)  # only find static libs
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)  # Changed to OFF for better compatibility with Python

# Find system Boost or use local built version
if(DEFINED ENV{BOOST_ROOT})
    set(BOOST_ROOT $ENV{BOOST_ROOT})
    message(STATUS "Using system Boost from BOOST_ROOT environment: ${BOOST_ROOT}")
    if(EXISTS "${BOOST_ROOT}/include/boost/version.hpp")
        message(STATUS "Boost include directory found at ${BOOST_ROOT}/include/boost")
    else()
        message(WARNING "Boost headers not found at ${BOOST_ROOT}/include/boost")
    endif()
else()
    set(BOOST_ROOT ${CMAKE_BINARY_DIR})
    message(STATUS "Using our own Boost, which was built at ${CMAKE_SOURCE_DIR}/external/boost_subset_1_58_0")
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/include/boost/version.hpp")
        message(WARNING "Built Boost headers not found at ${CMAKE_BINARY_DIR}/include/boost")
    endif()
endif()

# Find Boost components
find_package(Boost 1.55.0 COMPONENTS thread iostreams regex filesystem system program_options REQUIRED)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find Boost libraries. Please install or set BOOST_ROOT, or ensure it is built by external/make_dependencies.sh.") # Clarified error
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
add_subdirectory(src/c++/lib)
add_subdirectory(src/hap_py/haplo/cython) # Updated path

# Install rules for Python package files (IMPORTANT for scikit-build)
# This tells CMake to install the Python parts of the package.
# scikit-build uses these install rules to package your project correctly.

# Example: Install Python modules
# install(DIRECTORY ${SKBUILD_PACKAGE_DIR}/${SKBUILD_MODULE_NAME}/
#         DESTINATION ${SKBUILD_MODULE_NAME}
#         FILES_MATCHING PATTERN "*.py"
# )

# Example: Install package data (if any, e.g. non-code files within your package)
# install(DIRECTORY ${SKBUILD_PACKAGE_DIR}/${SKBUILD_MODULE_NAME}/data/
#         DESTINATION ${SKBUILD_MODULE_NAME}/data
#         FILES_MATCHING PATTERN "*"
# )

# The Cython modules built via add_cython_target from CythonSupport.cmake
# should be automatically handled by scikit-build if they are linked to a target
# that scikit-build is aware of, or if their output location is standard.
# Ensure your Cython modules are correctly defined as libraries and linked.

message(STATUS "Configuration completed successfully. Use 'pip install .' to build and install.")
