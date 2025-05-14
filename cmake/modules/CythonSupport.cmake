# Cython and Python integration module for hap.py
# This module provides utilities for working with Cython components

# Find Python interpreter
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Find Cython if not already found
if(NOT DEFINED CYTHON_EXECUTABLE)
  # Try to find cython using Python
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import cython; print(cython.__path__[0])"
    RESULT_VARIABLE CYTHON_RESULT
    OUTPUT_VARIABLE CYTHON_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  
  if(CYTHON_RESULT EQUAL 0)
    find_program(CYTHON_EXECUTABLE "cython" 
      HINTS "${CYTHON_PATH}/../bin"
            "${CYTHON_PATH}/bin"
      PATH_SUFFIXES bin
    )
  endif()
  
  if(NOT CYTHON_EXECUTABLE)
    # Fallback to standard find_program
    find_program(CYTHON_EXECUTABLE NAMES cython cython3)
  endif()
  
  if(NOT CYTHON_EXECUTABLE)
    message(STATUS "Cython not found, Python bindings will not be built")
    set(BUILD_PYTHON OFF CACHE BOOL "Build Python bindings" FORCE)
  else()
    message(STATUS "Found Cython: ${CYTHON_EXECUTABLE}")
  endif()
endif()

# Function to add a Cython module
function(add_cython_module module_name sources)
  if(NOT CYTHON_EXECUTABLE)
    message(WARNING "Cython not found, skipping module ${module_name}")
    return()
  endif()
  
  # Parse additional arguments
  set(options )
  set(oneValueArgs OUTPUT_DIR INCLUDE_DIRECTORIES)
  set(multiValueArgs LIBRARIES DEPENDS)
  cmake_parse_arguments(ACM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  
  # Default output directory
  if(NOT ACM_OUTPUT_DIR)
    set(ACM_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
  
  # Create output directory if it doesn't exist
  file(MAKE_DIRECTORY ${ACM_OUTPUT_DIR})
  
  # For each source file
  foreach(cython_source ${sources})
    get_filename_component(cython_module_name "${cython_source}" NAME_WE)
    
    # Generate C++ source
    set(generated_file "${ACM_OUTPUT_DIR}/${cython_module_name}.cpp")
    
    add_custom_command(
      OUTPUT ${generated_file}
      COMMAND ${CYTHON_EXECUTABLE} --cplus -3 --fast-fail
              -o ${generated_file} ${cython_source}
      DEPENDS ${cython_source} ${ACM_DEPENDS}
      COMMENT "Cythonizing ${cython_source}"
    )
    
    # Build module
    Python3_add_library(${module_name} MODULE ${generated_file})
    
    # Set include directories
    if(ACM_INCLUDE_DIRECTORIES)
      target_include_directories(${module_name} PRIVATE ${ACM_INCLUDE_DIRECTORIES})
    endif()
    
    # Add Python include directories
    target_include_directories(${module_name} PRIVATE ${Python3_INCLUDE_DIRS})
    
    # Set module properties
    set_target_properties(${module_name} PROPERTIES 
      PREFIX ""
      CXX_STANDARD 11
      CXX_STANDARD_REQUIRED YES
      CXX_EXTENSIONS NO
    )
    
    if(APPLE)
      set_target_properties(${module_name} PROPERTIES
        SUFFIX ".so"  # Use .so even on macOS for Python compatibility
      )
    endif()
    
    # Link libraries
    if(ACM_LIBRARIES)
      target_link_libraries(${module_name} PRIVATE ${ACM_LIBRARIES})
    endif()
  endforeach()
endfunction()
