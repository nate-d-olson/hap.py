# C++ Settings Module for hap.py
# This module provides consistent C++ settings across the project

# Set C++ standard and compiler options
macro(happy_set_cpp_standard target)
  set_target_properties(${target} PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE YES
  )
  
  # Platform-specific compiler flags
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE 
      -Wall 
      -Wextra
      -Wpedantic
      -Wno-unused-parameter  # Common in bioinformatics API code
    )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(${target} PRIVATE 
      /W4
      /MP  # Multi-processor compilation
    )
  endif()
endmacro()

# Add genomics-specific compiler definitions
macro(happy_add_genomics_definitions target)
  target_compile_definitions(${target} PRIVATE
    HAPLOTYPE_COMPARISON_TOOLS  # Main application identifier
    _FILE_OFFSET_BITS=64        # Support large files
    $<$<CONFIG:DEBUG>:DEBUG>    # Debug mode flag
  )
endmacro()

# Set optimization levels appropriate for genomic applications
macro(happy_set_optimization_flags target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE
      $<$<CONFIG:RELEASE>:-O3>          # High optimization for Release
      $<$<CONFIG:RELEASE>:-ffast-math>  # Fast math for genomics calculations
      $<$<CONFIG:DEBUG>:-O0>            # No optimization for Debug
      $<$<CONFIG:DEBUG>:-g3>            # Max debug info
    )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(${target} PRIVATE
      $<$<CONFIG:RELEASE>:/O2>  # High optimization
      $<$<CONFIG:DEBUG>:/Od>    # No optimization
      $<$<CONFIG:DEBUG>:/Zi>    # Debug info
    )
  endif()
endmacro()

# Configure a hap.py C++ executable with standard settings
macro(happy_add_executable target sources)
  add_executable(${target} ${sources})
  happy_set_cpp_standard(${target})
  happy_add_genomics_definitions(${target})
  happy_set_optimization_flags(${target})
  
  # Apply code signing on Apple platforms if module is included
  if(COMMAND apple_codesign)
    apple_codesign(${target})
  endif()
endmacro()

# Configure a hap.py C++ library with standard settings
macro(happy_add_library target type sources)
  add_library(${target} ${type} ${sources})
  happy_set_cpp_standard(${target})
  happy_add_genomics_definitions(${target})
  happy_set_optimization_flags(${target})
endmacro()

# Helper function to setup include directories for targets
function(happy_setup_includes target)
  target_include_directories(${target}
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/c++/include>
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:include>
    PRIVATE
      ${CMAKE_SOURCE_DIR}/external/klib
      ${CMAKE_SOURCE_DIR}/external/intervaltree
  )
endfunction()

# Helper function to link common dependencies
function(happy_link_common_libs target)
  target_link_libraries(${target}
    PUBLIC
      ${HAPLOTYPES_LIBRARY}
      jsoncpp_lib
      ${Boost_LIBRARIES}
      ${HTSLIB_LIBRARIES}
      ${ZLIB_LIBRARIES}
    PRIVATE
      ${CMAKE_THREAD_LIBS_INIT}
  )
endfunction()
