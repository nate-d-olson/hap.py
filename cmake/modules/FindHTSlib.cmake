# FindHTSlib.cmake
# Finds the HTSlib library
#
# This will define:
# HTSLIB_FOUND - True if HTSlib was found
# HTSLIB_INCLUDE_DIRS - HTSlib include directories
# HTSLIB_LIBRARIES - HTSlib libraries

# Try to find HTSlib using pkg-config first
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_HTSLIB QUIET htslib)
endif()

# Find the include directory
find_path(HTSLIB_INCLUDE_DIR
    NAMES htslib/sam.h
    PATHS
        ${PC_HTSLIB_INCLUDEDIR}
        ${PC_HTSLIB_INCLUDE_DIRS}
        /opt/homebrew/include
        /opt/homebrew/opt/htslib/include
        /usr/local/include
        /usr/include
    PATH_SUFFIXES htslib
)

# Find the library
find_library(HTSLIB_LIBRARY
    NAMES hts libhts.1
    PATHS
        ${PC_HTSLIB_LIBDIR}
        ${PC_HTSLIB_LIBRARY_DIRS}
        /opt/homebrew/lib
        /opt/homebrew/opt/htslib/lib
        /usr/local/lib
        /usr/lib
)

# Set version string
set(HTSLIB_VERSION_STRING ${PC_HTSLIB_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HTSlib
    REQUIRED_VARS
        HTSLIB_LIBRARY
        HTSLIB_INCLUDE_DIR
    VERSION_VAR HTSLIB_VERSION_STRING
)

if(HTSLIB_FOUND)
    set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY})
    set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR})
    
    # Add imported target
    if(NOT TARGET HTSlib::HTSlib)
        add_library(HTSlib::HTSlib UNKNOWN IMPORTED)
        set_target_properties(HTSlib::HTSlib PROPERTIES
            IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(HTSLIB_INCLUDE_DIR HTSLIB_LIBRARY)
