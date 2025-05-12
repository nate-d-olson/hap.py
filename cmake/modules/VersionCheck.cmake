# Check for minimum required versions of build tools and dependencies

# CMake version check
if(CMAKE_VERSION VERSION_LESS 3.14)
    message(FATAL_ERROR "CMake >= 3.14 required. Current version: ${CMAKE_VERSION}")
endif()

# C++ compiler version check
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0)
        message(FATAL_ERROR "GCC >= 7.0 required. Current version: ${CMAKE_CXX_COMPILER_VERSION}")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message(FATAL_ERROR "Clang >= 6.0 required. Current version: ${CMAKE_CXX_COMPILER_VERSION}")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
        message(FATAL_ERROR "AppleClang >= 10.0 required. Current version: ${CMAKE_CXX_COMPILER_VERSION}")
    endif()
endif()

# Python version check (needed for build scripts)
find_package(Python3 3.6 REQUIRED COMPONENTS Interpreter Development)
