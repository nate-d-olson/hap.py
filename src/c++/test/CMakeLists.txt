# Find test sources
file(GLOB TEST_SOURCES test_*.cpp test_*.c)

# Create test executable with our standardized macro
happy_add_executable(test_haplotypes "${TEST_SOURCES}")

# Setup standard includes
happy_setup_includes(test_haplotypes)

# Link against libraries
target_link_libraries(test_haplotypes 
    PRIVATE
        ${HAPLOTYPES_ALL_LIBS}
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

# For each test file, add a CTest test
file(GLOB TEST_FILES test_*.cpp)
foreach(TEST_FILE ${TEST_FILES})
    # Extract the test name from the filename (without extension)
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    
    # Register the test with CTest
    add_test(
        NAME ${TEST_NAME}
        COMMAND test_haplotypes --run_test=${TEST_NAME} --log_level=test_suite
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    
    # Set test properties
    set_tests_properties(${TEST_NAME} PROPERTIES
        ENVIRONMENT "BOOST_TEST_LOG_LEVEL=test_suite"
    )
endforeach()

# Add a convenience target to run all tests
add_custom_target(check 
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS test_haplotypes
    COMMENT "Running all tests"
)

# Install test executable if requested
option(INSTALL_TESTS "Install test executables" OFF)
if(INSTALL_TESTS)
    install(TARGETS test_haplotypes
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/tests
    )
endif()
