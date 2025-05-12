# Common properties for hap.py
# This module defines version numbers and other global constants

# Project information
set(HAPPY_PROJECT_NAME "hap.py")
set(HAPPY_PROJECT_DESCRIPTION "Haplotype-based variant calling assessment tools")
set(HAPPY_PROJECT_HOMEPAGE "https://github.com/illumina/hap.py")

# Version numbers (should be updated for releases)
set(HAPPY_VERSION_MAJOR 1)
set(HAPPY_VERSION_MINOR 0)
set(HAPPY_VERSION_PATCH 0)
set(HAPPY_VERSION_SUFFIX "")  # e.g., "-beta1"

# Derived version strings
set(HAPPY_VERSION "${HAPPY_VERSION_MAJOR}.${HAPPY_VERSION_MINOR}.${HAPPY_VERSION_PATCH}${HAPPY_VERSION_SUFFIX}")
set(HAPPY_VERSION_SHORT "${HAPPY_VERSION_MAJOR}.${HAPPY_VERSION_MINOR}")
set(HAPPY_SOVERSION "${HAPPY_VERSION_MAJOR}")

# Use Git version if available (for development builds)
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
    RESULT_VARIABLE GIT_DESCRIBE_RESULT
    ERROR_VARIABLE GIT_DESCRIBE_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  
  if(GIT_DESCRIBE_RESULT EQUAL 0)
    # Use Git version for development builds
    set(HAPLOTYPES_VERSION "${GIT_DESCRIBE_VERSION}")
  else()
    # Use CMake version as fallback
    set(HAPLOTYPES_VERSION "${HAPPY_VERSION}")
  endif()
else()
  # No Git, use CMake version
  set(HAPLOTYPES_VERSION "${HAPPY_VERSION}")
endif()

# Dependency requirements (consistent minimum versions)
set(BOOST_MIN_VERSION 1.74.0)
set(HTSLIB_MIN_VERSION 1.15.1)
set(BCFTOOLS_MIN_VERSION 1.15.1)
set(ZLIB_MIN_VERSION 1.2.12)
set(JSONCPP_MIN_VERSION 1.9.5)

# Set up Python version requirements
set(PYTHON_MIN_VERSION 3.8)
set(CYTHON_MIN_VERSION 0.29.24)

# Build options (can be overridden)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(BUILD_TESTS "Build test programs" ON)
option(BUILD_DOCS "Build documentation" OFF)
option(BUILD_PYTHON "Build Python bindings" ON)
option(BUILD_VCFEVAL "Build with vcfeval support" OFF)
option(USE_SGE "Enable SGE support" OFF)
option(USE_SYSTEM_LIBS "Try to use system libraries instead of downloading" ON)
option(FORCE_PACKAGED_LIBS "Always use packaged libraries even if system ones exist" OFF)

# Platform-specific settings
if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif()

# Output settings for debugging
message(STATUS "hap.py version: ${HAPPY_VERSION} (${HAPLOTYPES_VERSION})")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Python integration: ${BUILD_PYTHON}")
message(STATUS "vcfeval integration: ${BUILD_VCFEVAL}")
