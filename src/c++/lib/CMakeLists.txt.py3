# CMakeLists.txt for Haplo library (Python 3 version)
# This file builds the C++ components of hap.py and makes them available to Python

cmake_minimum_required(VERSION 3.10)

# Include necessary CMake modules
include(${CMAKE_SOURCE_DIR}/src/cmake/CythonSupport.cmake)

# Collect source files
file(GLOB_RECURSE HAPLO_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/align/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/diploidgraphs/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/quantify/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scmp/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tools/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/variant/*.cpp)

# Create the static library
add_library(haplotypes STATIC ${HAPLO_SRCS})

# Include directories
target_include_directories(haplotypes PUBLIC
    ${CMAKE_SOURCE_DIR}/src/c++/include
    ${CMAKE_BINARY_DIR}/include
    ${Boost_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(haplotypes PUBLIC
    ${Boost_LIBRARIES}
    ${HTSLIB_LIBRARIES}
)

# Set C++ standard
target_compile_features(haplotypes PUBLIC cxx_std_11)

# Add Python extension module
add_cython_module(
    _internal
    ${CMAKE_SOURCE_DIR}/src/python/Haplo/cython/_internal.pyx
    LIBRARIES
        haplotypes
        ${HAPLOTYPES_ALL_LIBS}
    INCLUDES
        ${CMAKE_SOURCE_DIR}/src/c++/include
        ${CMAKE_BINARY_DIR}/include
)

# Install library
install(TARGETS haplotypes
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
