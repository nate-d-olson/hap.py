# CythonSupport.py3.cmake - Python 3 compatible Cython support for CMake
# Author: [Your Name]
# Based on the original CythonSupport.cmake with updates for Python 3

include(CMakeParseArguments)

find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# Find Cython executable for Python 3
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import Cython.Build; print(Cython.Build.cython_executable)"
    OUTPUT_VARIABLE CYTHON_EXECUTABLE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

if(NOT CYTHON_EXECUTABLE)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import sys; sys.stdout.write(sys.exec_prefix)"
        OUTPUT_VARIABLE PYTHON_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    find_program(CYTHON_EXECUTABLE 
        NAMES cython cython.py cython3
        HINTS ${PYTHON_PREFIX}/bin
    )
endif()

if(NOT CYTHON_EXECUTABLE)
    message(FATAL_ERROR "Cython not found. Please install Cython for Python 3.")
endif()

message(STATUS "Found Cython: ${CYTHON_EXECUTABLE}")

# Detect NumPy include directory
if(Python3_NumPy_FOUND)
    set(NUMPY_INCLUDE_DIR "${Python3_NumPy_INCLUDE_DIRS}")
    message(STATUS "NumPy include directory: ${NUMPY_INCLUDE_DIR}")
else()
    message(WARNING "NumPy headers not found. Some Cython modules may not build correctly.")
endif()

# Function to add a Cython module
function(add_cython_module _name _source)
    cmake_parse_arguments(CY "" "" "INCLUDES;LIBRARIES;EMBED" ${ARGN})
    
    message(STATUS "Adding Cython module ${_name}")
    
    get_filename_component(_source_ext ${_source} EXT)
    get_filename_component(_source_path ${_source} PATH)
    get_filename_component(_source_name ${_source} NAME_WE)

    # Make directories in the binary tree
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_source_path})

    # C++ output file
    set(SOURCE_CPP ${CMAKE_CURRENT_BINARY_DIR}/${_source_path}/${_source_name}.cpp)

    # Cython language level and flags
    set(CYTHON_FLAGS "--cplus" "--directive" "language_level=3")
    
    # Add Cython compilation command
    add_custom_command(
        OUTPUT ${SOURCE_CPP}
        COMMAND ${CYTHON_EXECUTABLE}
        ARGS ${CYTHON_FLAGS} -o ${SOURCE_CPP} ${CMAKE_CURRENT_SOURCE_DIR}/${_source} 
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_source}
        COMMENT "Cythonizing ${_source}"
    )

    set_source_files_properties(
        ${SOURCE_CPP}
        PROPERTIES GENERATED TRUE
    )

    # Set up include directories
    set(INCLUDE_DIRS
        ${Python3_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/${_source_path} # For .pxd files
        ${CY_INCLUDES}
    )
    
    # Add NumPy include directory if available
    if(NUMPY_INCLUDE_DIR)
        list(APPEND INCLUDE_DIRS ${NUMPY_INCLUDE_DIR})
    endif()

    # Add module depending on module type
    if(CY_EMBED)
        # Static library for embedding
        add_library(${_name} STATIC ${SOURCE_CPP} ${CY_EMBED})
        target_include_directories(${_name} PRIVATE ${INCLUDE_DIRS})
        target_link_libraries(${_name} PRIVATE ${CY_LIBRARIES})
    else()
        # Python extension module
        add_library(${_name} MODULE ${SOURCE_CPP})
        target_include_directories(${_name} PRIVATE ${INCLUDE_DIRS})
        target_link_libraries(${_name} PRIVATE ${CY_LIBRARIES})
        
        # Don't add a lib prefix - Python expects modules without it
        set_target_properties(${_name} PROPERTIES
            PREFIX ""
            OUTPUT_NAME "${_name}"
        )
        
        # Platform-specific extensions and settings
        if(WIN32)
            set_target_properties(${_name} PROPERTIES SUFFIX ".pyd")
        elseif(APPLE)
            set_target_properties(${_name} PROPERTIES 
                SUFFIX ".so"
                LINK_FLAGS "-undefined dynamic_lookup"
            )
        else()
            set_target_properties(${_name} PROPERTIES SUFFIX ".so")
        endif()
    endif()

    # Install destination for the Python module
    if(NOT CY_EMBED)
        set(PYTHON_MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${_source_path}")
        string(REPLACE "/" "." PYTHON_MODULE_NAME "${_source_path}/${_source_name}")
        string(REGEX REPLACE "^src\\.python\\." "" PYTHON_MODULE_NAME "${PYTHON_MODULE_NAME}")
        
        message(STATUS "Will install Python module ${PYTHON_MODULE_NAME} from ${PYTHON_MODULE_DIR}")
        
        # Figure out the installation path
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import sys; import distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=True, standard_lib=False))"
            OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        
        install(TARGETS ${_name} 
            LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/${_source_path}"
            RUNTIME DESTINATION "${PYTHON_SITE_PACKAGES}/${_source_path}"
        )
    endif()
endfunction()
