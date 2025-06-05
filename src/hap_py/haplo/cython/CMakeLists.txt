# CMakeLists.txt for building Cython modules in Python 3
# This file should be placed in the src/python/Haplo/cython directory

cmake_minimum_required(VERSION 3.10)

# Find Python3 and its components
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(NumPy REQUIRED)

# Include directories
include_directories(
    ${CMAKE_BINARY_DIR}/include
    ${Python3_INCLUDE_DIRS}
    ${NumPy_INCLUDE_DIRS}
)

# Set module output directory
set(CYTHON_MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Define Cython modules - update these paths based on actual .pyx files
# Note: These should match the actual Cython source files in your project

# Check if the .pyx files exist and add them conditionally
set(CYTHON_SOURCES
    happyroc.pyx
    variant_processor.pyx
    sequence_utils.pyx
)

foreach(CYTHON_SOURCE ${CYTHON_SOURCES})
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${CYTHON_SOURCE}")
        get_filename_component(MODULE_NAME ${CYTHON_SOURCE} NAME_WE)

        add_cython_module(
            ${MODULE_NAME} ${CYTHON_SOURCE}
            INCLUDES
                ${CMAKE_BINARY_DIR}/include
        )

        message(STATUS "Added Cython module: ${MODULE_NAME}")
    else()
        message(WARNING "Cython source file not found: ${CYTHON_SOURCE}")
    endif()
endforeach()
