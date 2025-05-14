# Apple code signing support for hap.py
# This module provides functions for code signing on macOS

# Check if we're on macOS
if(NOT APPLE)
  return()
endif()

# Find the codesign utility
find_program(CODESIGN_EXECUTABLE NAMES codesign)

if(NOT CODESIGN_EXECUTABLE)
  message(STATUS "codesign utility not found, binaries will not be signed")
  return()
endif()

# Define code signing identity option
set(APPLE_CODE_SIGN_IDENTITY "" CACHE STRING 
    "Code signing identity for macOS binaries (empty to disable signing)")

# Function to sign executables
function(apple_codesign target)
  if(NOT APPLE OR NOT CODESIGN_EXECUTABLE OR NOT APPLE_CODE_SIGN_IDENTITY)
    return()
  endif()
  
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND ${CODESIGN_EXECUTABLE} --force --timestamp --options runtime
            --sign "${APPLE_CODE_SIGN_IDENTITY}"
            $<TARGET_FILE:${target}>
    COMMENT "Code signing ${target} with identity '${APPLE_CODE_SIGN_IDENTITY}'"
    VERBATIM
  )
endfunction()

# Output status message if signing is enabled
if(APPLE_CODE_SIGN_IDENTITY)
  message(STATUS "macOS code signing enabled with identity: ${APPLE_CODE_SIGN_IDENTITY}")
endif()
