# hap.py Modern Build System Documentation

## Overview

The hap.py build system has been modernized using CMake best practices and structured in a modular way to improve maintainability, portability, and extensibility. This document explains the new structure and how to use it effectively.

## Build System Components

### Core CMake Files

- **CMakeLists.txt** - Main project configuration
- **CMakePresets.json** - Standard build configurations for different platforms and build types
- **cmake/modules/** - Specialized CMake modules

### CMake Modules

The build system uses several custom modules to maintain consistency and simplify configuration:

#### 1. HappyProperties.cmake

Centralizes version information and global settings:

- Project metadata and version numbers
- Dependency version requirements
- Build options with defaults

#### 2. CppSettings.cmake

Standardizes C++ compilation across all components:

- Compiler flags and warning settings
- Optimization levels
- Standard C++ settings
- Helper macros for adding libraries and executables

#### 3. CythonSupport.cmake

Manages Python/C++ integration:

- Python and Cython detection
- Utilities for building Cython modules
- Platform-specific settings for Python extensions

#### 4. Find Modules

Custom modules to locate external dependencies:

- FindHTSlib.cmake - Locates the HTSlib bioinformatics library
- FindBCFtools.cmake - Locates BCFtools for variant manipulation

#### 5. Support Modules

Additional utility modules:

- MacOSPaths.cmake - macOS-specific path handling
- AppleCodeSign.cmake - Code signing support for macOS
- PackageHelpers.cmake - Package configuration utilities
- VersionCheck.cmake - Version compatibility checking

## Using the Build System

### Basic Build

To build the project using the default configuration:

```bash
mkdir build
cd build
cmake ..
make -j8
```

### Using CMake Presets

For standard configurations, use CMake presets:

```bash
# Debug build
cmake --preset=debug
cmake --build --preset=debug

# Release build
cmake --preset=release
cmake --build --preset=release
```

### Available Options

Set these options with `-D<OPTION>=<VALUE>` when running CMake:

| Option | Default | Description |
|--------|---------|-------------|
| BUILD_TESTS | ON | Build test programs |
| BUILD_DOCS | OFF | Build documentation |  
| BUILD_PYTHON | ON | Build Python bindings |
| BUILD_VCFEVAL | OFF | Build with vcfeval support |
| USE_SGE | OFF | Enable SGE support |
| USE_SYSTEM_LIBS | ON | Use system libraries if available |
| FORCE_PACKAGED_LIBS | OFF | Always use packaged libraries |

Example:

```bash
cmake -DBUILD_PYTHON=OFF -DBUILD_TESTS=ON ..
```

### Dependencies

The build system automatically handles dependencies through:

1. System libraries (if USE_SYSTEM_LIBS=ON)
2. Fetching and building dependencies (if system libraries aren't found)
3. Using packaged dependencies (if FORCE_PACKAGED_LIBS=ON)

### Cross-Platform Support

The build system supports:

- Linux (GCC, Clang)
- macOS (Apple Clang)
- Windows (experimental, MSVC)

## Advanced Usage

### Creating Custom Build Configurations

Add custom configurations to CMakePresets.json:

```json
{
  "name": "custom-config",
  "displayName": "Custom Configuration",
  "generator": "Ninja",
  "binaryDir": "${sourceDir}/build/custom",
  "cacheVariables": {
    "CMAKE_BUILD_TYPE": "RelWithDebInfo",
    "BUILD_PYTHON": "ON",
    "USE_SGE": "ON"
  }
}
```

### Extending the Build System

To add new C++ components:

1. Create a new directory in src/c++/
2. Add CMakeLists.txt using the provided macros:

```cmake
# Source files
set(SOURCES
    component1.cpp
    component2.cpp
)

# Create library
happy_add_library(my_component STATIC "${SOURCES}")
happy_setup_includes(my_component)
target_link_libraries(my_component ${HAPLOTYPES_ALL_LIBS})
```

### Adding Python Extensions

For new Cython modules:

```cmake
add_cython_module(
    my_module
    my_module.pyx
    LIBRARIES 
        haplotypes
        ${HAPLOTYPES_ALL_LIBS}
)
```

## Contributing

When modifying the build system:

1. Keep component-specific settings in their own CMakeLists.txt files
2. Add reusable functionality as modules in cmake/modules/
3. Maintain backward compatibility with existing builds
4. Document significant changes
5. Test on multiple platforms before submitting PRs
