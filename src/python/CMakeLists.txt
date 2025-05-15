# Python module build configuration for hap.py with Python 3 support
# Updated version of the original CMakeLists.txt

# Get Python information through CMake's FindPython3
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Include Cython support
include(${CMAKE_SOURCE_DIR}/src/cmake/CythonSupport.cmake)

# Define Python destination directories
set(PYTHON3_LIB_DIR "${CMAKE_BINARY_DIR}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")
message(STATUS "Will install Python modules to: ${PYTHON3_LIB_DIR}")

# Copy Python scripts to bin directory
file(GLOB COPY_FILES *.py *.py.py3)
add_custom_target(python_bin ALL COMMENT "Copying Python scripts...")

foreach(FILENAME ${COPY_FILES})
    get_filename_component(TARGETNAME ${FILENAME} "NAME")
    # # If it's a .py3 file, we'll copy it as a .py file
    # if(TARGETNAME MATCHES "\.py3$")
    #     string(REGEX REPLACE "\.py3$" ".py" TARGETNAME ${TARGETNAME})
    # endif()
    set(SRC "${FILENAME}")
    set(DST "${CMAKE_BINARY_DIR}/bin/${TARGETNAME}")

    add_custom_command(
        TARGET python_bin
        COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
    )
endforeach(FILENAME)

# Copy Python package modules
SET(PYTHON_SRC "${CMAKE_SOURCE_DIR}/src/python/Haplo")
SET(PYTHON_DST "${PYTHON3_LIB_DIR}/Haplo")
add_custom_target(python_src1 ALL COMMENT "Copying ${PYTHON_SRC} -> ${PYTHON_DST}")

add_custom_command(TARGET python_src1
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PYTHON_SRC} ${PYTHON_DST}
)

# SET(PYTHON_SRC "${CMAKE_SOURCE_DIR}/src/python/Somatic")
# SET(PYTHON_DST "${PYTHON3_LIB_DIR}/Somatic")
# add_custom_target(python_src2 ALL COMMENT "Copying ${PYTHON_SRC} -> ${PYTHON_DST}")

# add_custom_command(TARGET python_src2
#                    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PYTHON_SRC} ${PYTHON_DST}
# )

SET(PYTHON_SRC "${CMAKE_SOURCE_DIR}/src/python/Tools")
SET(PYTHON_DST "${PYTHON3_LIB_DIR}/Tools")
add_custom_target(python_src3 ALL COMMENT "Copying ${PYTHON_SRC} -> ${PYTHON_DST}")

add_custom_command(TARGET python_src3

# Add Cython module build
add_subdirectory(Haplo/cython)
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${PYTHON_SRC} ${PYTHON_DST}
)

# Generate version file
configure_file("Haplo/version.py.in"
               "${PYTHON3_LIB_DIR}/Haplo/version.py")

# Cython integration
# We need to find where the C++ extension module is built and ensure it's properly
# linked with Python 3

# Create cython directory if it doesn't exist
file(MAKE_DIRECTORY "${PYTHON3_LIB_DIR}/Haplo/cython")

# Create a simple __init__.py file in the cython directory
file(WRITE "${PYTHON3_LIB_DIR}/Haplo/cython/__init__.py" "# Cython package for Haplo\n")
