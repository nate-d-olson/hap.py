# Top-level CMake file for Haplotype Tools (hap.py)
cmake_minimum_required(VERSION 3.14...3.27)

# Set CMake policies
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW) # use PROFILE directly for Xcode
endif()
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW) # allow options to be set without overriding
endif()
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # use <PackageName>_ROOT variables
endif()

# Include GNUInstallDirs early (needs to be before project command)
include(GNUInstallDirs)

# Add our custom CMake modules
list(APPEND CMAKE_MODULE_PATH 
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake"
)

# Include core properties module before project() to set version info
include(HappyProperties)

# Project declaration with metadata
project(HAPLOTYPES
    VERSION ${HAPPY_VERSION}
    DESCRIPTION "${HAPPY_PROJECT_DESCRIPTION}"
    HOMEPAGE_URL "${HAPPY_PROJECT_HOMEPAGE}"
    LANGUAGES CXX C)

# Include our modern CMake modules
include(FetchContent)
include(CMakePackageConfigHelpers)
include(VersionCheck)
include(PackageHelpers)
include(MacOSPaths)
include(CppSettings)      # C++ specific settings
include(CythonSupport)    # Cython and Python integration
include(AppleCodeSign)    # macOS code signing support

# Set output directories for consistency across platforms
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Legacy includes - to be modernized in future PRs
include(cxx)
include(configureFiles)

# Configure vcfeval support
if(BUILD_VCFEVAL)
    message(STATUS "Building with vcfeval support")
    set(ENV{BUILD_VCFEVAL} "1")
    set(ENV{VCFEVAL_WRAPPER} ${VCFEVAL_WRAPPER})
    set(VCFEVAL_AVAILABLE 1)
else()
    set(VCFEVAL_AVAILABLE 0)
endif()

# Include our modern dependency management
include(ModernDependencies)

# Required system packages
find_package(Threads REQUIRED)

# Configure JsonCPP (using modern FetchContent)
FetchContent_Declare(
    jsoncpp
    GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
    GIT_TAG 1.9.5
    CMAKE_CACHE_ARGS
        -DJSONCPP_WITH_TESTS:BOOL=OFF
        -DJSONCPP_WITH_POST_BUILD_UNITTEST:BOOL=OFF
        -DBUILD_SHARED_LIBS:BOOL=ON
)
FetchContent_MakeAvailable(jsoncpp)

# Configure include directories
include_directories(
    ${CMAKE_BINARY_DIR}/include
    ${HAPLOTYPES_SOURCE_DIR}/external/klib
    ${HAPLOTYPES_SOURCE_DIR}/external/intervaltree
    ${HAPLOTYPES_SOURCE_DIR}/src/c++/include
)

# Configure SGE support
if(USE_SGE)
    message(STATUS "SGE mode selected, you will need to use --force-interactive to run hap.py")
    set(SGE_AVAILABLE 1)
else()
    set(SGE_AVAILABLE 0)
endif()

# Extract version from Git
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always
        OUTPUT_VARIABLE HAPLOTYPES_VERSION
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()

# If Git isn't available or fails, use a default version
if(NOT DEFINED HAPLOTYPES_VERSION OR HAPLOTYPES_VERSION STREQUAL "")
    set(HAPLOTYPES_VERSION "${CMAKE_PROJECT_VERSION}")
endif()

# Configure the version header
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/c++/include/Version.hh.in"
               "${CMAKE_BINARY_DIR}/include/Version.hh"
               @ONLY)

# Build our libraries
add_subdirectory(src/c++/lib)

# Configure library dependencies
set(HAPLOTYPES_ALL_LIBS 
    ${HAPLOTYPES_LIBRARY}
    ${GSSW_LIBRARY}
    ${KLIB_LIBRARY}
    jsoncpp_lib
    ${Boost_LIBRARIES}
    ${HTSLIB_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT})


# Add subdirectories for components
add_subdirectory(src/c++/main)
add_subdirectory(src/data)

# Build tests conditionally
if(BUILD_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(src/c++/test)
endif()

# Build Python bindings conditionally
if(BUILD_PYTHON)
    add_subdirectory(src/python)
endif()

# Documentation build
if(BUILD_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        add_subdirectory(doc)
    else()
        message(STATUS "Doxygen not found, documentation will not be built")
    endif()
endif()

# Installation rules
install(DIRECTORY "${CMAKE_BINARY_DIR}/bin/"
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        USE_SOURCE_PERMISSIONS
        FILES_MATCHING PATTERN "*")

install(DIRECTORY "${CMAKE_BINARY_DIR}/lib/"
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(BUILD_VCFEVAL)
    install(DIRECTORY "${CMAKE_BINARY_DIR}/libexec/"
            DESTINATION ${CMAKE_INSTALL_LIBEXECDIR}
            USE_SOURCE_PERMISSIONS)
endif()

# Package configuration files
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/packages/hap.py-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/hap.py-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hap.py
)

# Version config file
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/hap.py-config-version.cmake
    VERSION ${CMAKE_PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Install package config files
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/hap.py-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/hap.py-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hap.py
)

# Configure CPack for packaging
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPackConfig.cmake)
