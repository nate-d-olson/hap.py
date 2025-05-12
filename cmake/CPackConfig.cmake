# CPack configuration for hap.py
include(InstallRequiredSystemLibraries)

# Version information from CMake
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})

# Package information
set(CPACK_PACKAGE_NAME "hap.py")
set(CPACK_PACKAGE_VENDOR "Illumina")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Haplotype comparison tools for benchmarking variant calls")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "hap.py")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

# Source package settings
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /build
    /*.build
    /*.cbp
    /*.user
    .DS_Store
)

# Platform-specific settings
if(UNIX)
    if(CMAKE_SYSTEM_NAME MATCHES Linux)
        # Linux specific settings
        set(CPACK_GENERATOR "TGZ;DEB;RPM")
        set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Illumina")
        set(CPACK_DEBIAN_PACKAGE_DEPENDS "python3, zlib1g, libboost-all-dev, htslib-dev")
        set(CPACK_RPM_PACKAGE_REQUIRES "python3, zlib, boost-devel, htslib-devel")
    elseif(APPLE)
        # macOS specific settings
        set(CPACK_GENERATOR "TGZ;productbuild")
        set(CPACK_PRODUCTBUILD_IDENTIFIER "com.illumina.hap.py")
    endif()
else()
    # Windows specific settings
    set(CPACK_GENERATOR "ZIP;NSIS")
endif()

include(CPack)
