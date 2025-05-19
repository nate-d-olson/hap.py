# Top-level Cmake file for Haplotype Tools (Python 3 version)
#
# Author: Peter Krusche <pkrusche@illumina.com>
# Python 3 Migration: [Your name]
#
cmake_minimum_required(VERSION 3.10)
project(HAPLOTYPES VERSION 0.3.15 LANGUAGES CXX C)

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
if(Python3_VERSION VERSION_LESS 3.6)
  message(FATAL_ERROR "Python 3.6 or later is required")
endif()

# NumPy detection is handled in CythonSupport.cmake for Cython modules

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
execute_process(
    COMMAND ${CMAKE_SOURCE_DIR}/external/make_dependencies_py3.sh
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
add_subdirectory(src/c++/lib) # Ensure correct path to library CMakeLists
# Add the directory containing Cython modules
add_subdirectory(src/python/Haplo/cython) # Ensure this path is correct

# Status message at the end
message(STATUS "Configuration completed successfully. Run 'make' to build.")
