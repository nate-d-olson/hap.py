# FindBCFtools.cmake
# Finds the BCFtools library and executables
#
# This will define:
# BCFTOOLS_FOUND - True if BCFtools was found
# BCFTOOLS_INCLUDE_DIRS - BCFtools include directories
# BCFTOOLS_LIBRARIES - BCFtools libraries
# BCFTOOLS_EXECUTABLE - BCFtools executable

# Try to find BCFtools using pkg-config first
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_BCFTOOLS QUIET bcftools)
endif()

# Find bcftools executable
find_program(BCFTOOLS_EXECUTABLE
    NAMES bcftools
    PATHS
        ${PC_BCFTOOLS_PREFIX}/bin
        /usr/local/bin
        /usr/bin
        /opt/homebrew/bin
)

# Find the include directory
find_path(BCFTOOLS_INCLUDE_DIR
    NAMES bcftools/bcf.h
    PATHS
        ${PC_BCFTOOLS_INCLUDEDIR}
        ${PC_BCFTOOLS_INCLUDE_DIRS}
        /usr/local/include
        /usr/include
        /opt/homebrew/include
)

# Find the library
find_library(BCFTOOLS_LIBRARY
    NAMES bcftools
    PATHS
        ${PC_BCFTOOLS_LIBDIR}
        ${PC_BCFTOOLS_LIBRARY_DIRS}
        /usr/local/lib
        /usr/lib
        /opt/homebrew/lib
)

# Set version string
if(BCFTOOLS_EXECUTABLE)
    execute_process(
        COMMAND ${BCFTOOLS_EXECUTABLE} --version
        OUTPUT_VARIABLE BCFTOOLS_VERSION_OUTPUT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(BCFTOOLS_VERSION_OUTPUT MATCHES "bcftools [0-9\\.]+")
        string(REGEX REPLACE ".*bcftools ([0-9\\.]+).*" "\\1" BCFTOOLS_VERSION_STRING "${BCFTOOLS_VERSION_OUTPUT}")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BCFtools
    REQUIRED_VARS
        BCFTOOLS_EXECUTABLE
        BCFTOOLS_INCLUDE_DIR
    VERSION_VAR BCFTOOLS_VERSION_STRING
)

if(BCFTOOLS_FOUND)
    set(BCFTOOLS_LIBRARIES ${BCFTOOLS_LIBRARY})
    set(BCFTOOLS_INCLUDE_DIRS ${BCFTOOLS_INCLUDE_DIR})
    
    # Add imported target
    if(NOT TARGET BCFtools::BCFtools)
        add_library(BCFtools::BCFtools UNKNOWN IMPORTED)
        set_target_properties(BCFtools::BCFtools PROPERTIES
            IMPORTED_LOCATION "${BCFTOOLS_EXECUTABLE}"
            INTERFACE_INCLUDE_DIRECTORIES "${BCFTOOLS_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(
    BCFTOOLS_INCLUDE_DIR
    BCFTOOLS_LIBRARY
    BCFTOOLS_EXECUTABLE
)
