# MacOSPaths.cmake
# Handle macOS specific paths and dependencies

if(APPLE)
    # Try to find Homebrew prefix
    execute_process(
        COMMAND brew --prefix
        OUTPUT_VARIABLE HOMEBREW_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(HOMEBREW_PREFIX)
        # Add Homebrew paths
        list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}")
        list(APPEND CMAKE_LIBRARY_PATH "${HOMEBREW_PREFIX}/lib")
        list(APPEND CMAKE_INCLUDE_PATH "${HOMEBREW_PREFIX}/include")
        
        # Add package-specific paths
        list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}/opt/zlib")
        list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}/opt/boost")
        list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}/opt/htslib")
        list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}/opt/bcftools")
        list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}/opt/jsoncpp")
    endif()
endif()
