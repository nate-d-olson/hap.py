# Helper functions for finding system packages

include(FindPackageHandleStandardArgs)

# Helper function to find system libraries
function(find_system_package NAME VERSION)
    if(USE_SYSTEM_LIBS AND NOT FORCE_PACKAGED_LIBS)
        find_package(${NAME} ${VERSION} QUIET)
        if(${NAME}_FOUND)
            message(STATUS "Using system ${NAME} version ${${NAME}_VERSION}")
            return()
        endif()
    endif()
    
    if(FORCE_PACKAGED_LIBS)
        message(STATUS "Forced to use packaged ${NAME}")
    else()
        message(STATUS "System ${NAME} >= ${VERSION} not found, will use packaged version")
    endif()
endfunction()

# Helper function to handle git repositories
function(fetch_git_repository NAME URL TAG)
    FetchContent_Declare(
        ${NAME}
        GIT_REPOSITORY ${URL}
        GIT_TAG ${TAG}
        GIT_SHALLOW TRUE
    )
    
    FetchContent_GetProperties(${NAME})
    if(NOT ${NAME}_POPULATED)
        message(STATUS "Fetching ${NAME} from ${URL}")
        FetchContent_Populate(${NAME})
    endif()
endfunction()
