
# Modern CMake configuration for haplotypes library
# Define library components with their sources
set(LIBRARY_COMPONENTS
    align
    diploidgraphs
    quantify
    scmp
    tools
    variant
)

# Create a list for library sources
set(LIBRARY_SOURCES "")

# Add sources from each component directory
foreach(component ${LIBRARY_COMPONENTS})
    file(GLOB_RECURSE COMPONENT_SOURCES "${component}/*.cpp" "${component}/*.c")
    list(APPEND LIBRARY_SOURCES ${COMPONENT_SOURCES})
endforeach()

# Define the library using our standardized macro
happy_add_library(haplotypes STATIC "${LIBRARY_SOURCES}")

# Define include directories for this target
target_include_directories(haplotypes
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link dependencies with our standardized function
target_link_libraries(haplotypes
    PUBLIC
        jsoncpp_lib
        ${Boost_LIBRARIES}
        ${HTSLIB_LIBRARIES}
        ${ZLIB_LIBRARIES}
    PRIVATE
        ${CMAKE_THREAD_LIBS_INIT}
)

# Export the target name to parent scope
set(HAPLOTYPES_LIBRARY haplotypes PARENT_SCOPE)

# Install the library and header files
install(TARGETS haplotypes
    EXPORT hap_py_targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)