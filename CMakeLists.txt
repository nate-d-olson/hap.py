# Top-level Cmake file for Haplotype Tools
#
# Author: Peter Krusche <pkrusche@illumina.com>
#
cmake_minimum_required(VERSION 3.14...3.27)
project(HAPLOTYPES
    VERSION 1.0.0
    DESCRIPTION "Benchmarking tool for small variant calls"
    LANGUAGES CXX C)

# Add our custom CMake modules
list(APPEND CMAKE_MODULE_PATH 
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake"
)

# Include our modern CMake modules
include(VersionCheck)
include(PackageHelpers)
include(ModernDependencies)

# Output directories
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Legacy includes
include(cxx)
include(configureFiles)

if(BUILD_VCFEVAL)
    message("Will build rtg-tools.")
    set( ENV{BUILD_VCFEVAL} "1" )
    set( ENV{VCFEVAL_WRAPPER} ${VCFEVAL_WRAPPER} )
    set (VCFEVAL_AVAILABLE 1)
else()
    set (VCFEVAL_AVAILABLE 0)
endif()

# Dependencies configuration
include(FetchContent)

# Modern Boost setup
set(BOOST_MIN_VERSION 1.74.0)
find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS thread iostreams regex unit_test_framework filesystem system program_options)
if(NOT Boost_FOUND)
    message(STATUS "Boost ${BOOST_MIN_VERSION}+ not found in system, fetching from source")
    FetchContent_Declare(
        Boost
        URL https://boostorg.jfrog.io/artifactory/main/release/1.74.0/source/boost_1_74_0.tar.bz2
        URL_HASH SHA256=83bfc1507731a0906e387fc28b7ef5417d591429e51e788417fe9ff025e116b1
    )
    FetchContent_MakeAvailable(Boost)
endif()

# ZLIB configuration
find_package(ZLIB)
if(NOT ZLIB_FOUND)
    message(STATUS "ZLIB not found in system, fetching from source")
    FetchContent_Declare(
        ZLIB
        URL https://zlib.net/zlib-1.2.13.tar.gz
        URL_HASH SHA256=b3a24de97a8fdbc835b9833169501030b8977031bcb54b3b3ac13740f846ab30
    )
    FetchContent_MakeAvailable(ZLIB)
endif()

# HTSlib configuration
find_package(HTSlib)
if(NOT HTSlib_FOUND)
    message(STATUS "HTSlib not found in system, fetching from source")
    FetchContent_Declare(
        htslib
        GIT_REPOSITORY https://github.com/samtools/htslib.git
        GIT_TAG 1.17  # Latest stable release
    )
    FetchContent_MakeAvailable(htslib)
endif()

# BCFtools configuration (if needed)
if(BUILD_VCFEVAL)
    find_package(BCFtools)
    if(NOT BCFtools_FOUND)
        message(STATUS "BCFtools not found in system, fetching from source")
        FetchContent_Declare(
            bcftools
            GIT_REPOSITORY https://github.com/samtools/bcftools.git
            GIT_TAG 1.17  # Latest stable release
        )
        FetchContent_MakeAvailable(bcftools)
    endif()
endif()

# Include directories
include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${Boost_INCLUDE_DIRS})

# Additional external library configuration
FetchContent_Declare(
    jsoncpp
    GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
    GIT_TAG 1.9.5  # Latest stable version
)
FetchContent_MakeAvailable(jsoncpp)

# Include local headers
include_directories(${HAPLOTYPES_SOURCE_DIR}/external/klib)
include_directories(${HAPLOTYPES_SOURCE_DIR}/external/intervaltree)
include_directories(${HAPLOTYPES_SOURCE_DIR}/src/c++/include)

# SGE support configuration
if(USE_SGE)
    message(STATUS "SGE mode selected, you will need to use --force-interactive to run hap.py")
    set(SGE_AVAILABLE 1)
else()
    set(SGE_AVAILABLE 0)
endif()

# Required system packages
find_package(Threads REQUIRED)

# Build our libraries
add_subdirectory(src/c++/lib)

# Configure library dependencies
set(HAPLOTYPES_ALL_LIBS 
    ${HAPLOTYPES_LIBRARY}
    ${GSSW_LIBRARY}
    ${KLIB_LIBRARY}
    jsoncpp_lib  # Modern jsoncpp target
    ${Boost_LIBRARIES}
    ${HTSLIB_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT})


execute_process(COMMAND git describe --tags --always
    OUTPUT_VARIABLE HAPLOTYPES_VERSION
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/c++/include/Version.hh.in"
               "${CMAKE_BINARY_DIR}/include/Version.hh")

add_subdirectory (src/c++/main)
add_subdirectory (src/data)
add_subdirectory (src/c++/test)
add_subdirectory (src/python)

# install
install (DIRECTORY "${CMAKE_BINARY_DIR}/bin/" DESTINATION bin FILES_MATCHING PATTERN "*" PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_READ WORLD_EXECUTE)
install (DIRECTORY "${CMAKE_BINARY_DIR}/lib/" DESTINATION lib)
if(BUILD_VCFEVAL)
    install (DIRECTORY "${CMAKE_BINARY_DIR}/libexec/" DESTINATION libexec USE_SOURCE_PERMISSIONS)
endif()
