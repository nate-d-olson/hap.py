include(FetchContent)
include(ExternalProject)

# Handle macOS Homebrew paths
if(APPLE)
    execute_process(
        COMMAND brew --prefix
        OUTPUT_VARIABLE HOMEBREW_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(HOMEBREW_PREFIX)
        list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}")
        list(APPEND CMAKE_LIBRARY_PATH "${HOMEBREW_PREFIX}/lib")
        list(APPEND CMAKE_INCLUDE_PATH "${HOMEBREW_PREFIX}/include")
    endif()
endif()

# Set minimum versions for dependencies
set(BOOST_MIN_VERSION 1.74.0)  # Modern version with good C++11/14 support
set(HTSLIB_MIN_VERSION 1.15.1) # Modern HTSlib version
set(BCFTOOLS_MIN_VERSION 1.15.1)
set(ZLIB_MIN_VERSION 1.2.12)   # macOS system zlib version

# Add option to force using packaged dependencies
option(USE_SYSTEM_LIBS "Try to use system libraries instead of downloading" ON)
option(FORCE_PACKAGED_LIBS "Always use packaged libraries even if system ones exist" OFF)

# Configure Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_NO_BOOST_CMAKE ON)
set(Boost_NO_SYSTEM_PATHS OFF)

# First try to find system installations
find_package(ZLIB)
find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS thread iostreams regex unit_test_framework filesystem system program_options)

if(NOT ZLIB_FOUND)
    message(STATUS "System zlib not found, will fetch and build from source")
    FetchContent_Declare(
        zlib
        URL https://github.com/madler/zlib/releases/download/v${ZLIB_MIN_VERSION}/zlib-${ZLIB_MIN_VERSION}.tar.gz
        URL_HASH SHA256=b3a24de97a8fdbc835b9833169501030b8977031bcb54b3b3ac13740f846ab30
    )
    FetchContent_MakeAvailable(zlib)
endif()

if(NOT Boost_FOUND)
    message(STATUS "System Boost not found, will fetch and build from source")
    set(BOOST_VERSION 1.82.0)
    string(REPLACE "." "_" BOOST_VERSION_UNDERSCORE ${BOOST_VERSION})
    FetchContent_Declare(
        Boost
        URL https://boostorg.jfrog.io/artifactory/main/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION_UNDERSCORE}.tar.gz
        URL_HASH SHA256=66a469b6e608a154af0b0bae8dd6de9641e43c46c4615bd1c3295251bd6e3534
    )
    FetchContent_MakeAvailable(Boost)
endif()

# HTSlib
FetchContent_Declare(
    htslib
    GIT_REPOSITORY https://github.com/samtools/htslib.git
    GIT_TAG ${HTSLIB_MIN_VERSION}
)
FetchContent_MakeAvailable(htslib)

# BCFtools
FetchContent_Declare(
    bcftools
    GIT_REPOSITORY https://github.com/samtools/bcftools.git
    GIT_TAG ${BCFTOOLS_MIN_VERSION}
)
FetchContent_MakeAvailable(bcftools)

# Configure include and link directories
include_directories(
    ${ZLIB_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${htslib_SOURCE_DIR}
    ${bcftools_SOURCE_DIR}
)

link_directories(
    ${ZLIB_LIBRARY_DIRS}
    ${Boost_LIBRARY_DIRS}
    ${htslib_BINARY_DIR}
    ${bcftools_BINARY_DIR}
)
