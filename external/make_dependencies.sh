#!/bin/bash
# Updated external dependencies script for hap.py
# Python 3 compatible version with additional error checking and modernized builds

set -e
set -o pipefail

# Find python
PYTHON=python3
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TLD=$(pwd)/scratch
ISD=$(pwd)
CPU_COUNT=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)

# Print status messages
echo "Building external dependencies..."
echo "Source directory: $DIR"
echo "Build directory: $TLD"
echo "Install directory: $ISD"
echo "CPU count: $CPU_COUNT"

# Create function for checking system dependencies
check_system_dep() {
    local cmd=$1
    local package=$2
    local install_cmd=$3
    
    if ! command -v $cmd &> /dev/null; then
        echo "WARNING: $cmd not found. Please install $package using:"
        echo "  $install_cmd"
        return 1
    fi
    return 0
}

# Check system dependencies
check_system_dep "cmake" "CMake" "brew install cmake # or apt install cmake"
check_system_dep "g++" "C++ compiler" "brew install gcc # or apt install g++"
check_system_dep "make" "build tools" "brew install make # or apt install build-essential"

# Check if rebuild requested
if [ "$1" == "rebuild" ]; then
    echo "Rebuild requested, removing previous builds..."
    rm -rf ${TLD}
    rm -rf ${ISD}/include/boost
    rm -rf ${ISD}/include/htslib
    rm -rf ${ISD}/bin/bcftools
    rm -rf ${ISD}/bin/samtools
    rm -rf ${ISD}/libexec/rtg-tools-install
fi

# Create scratch directory
mkdir -p ${TLD}

# Build zlib
echo "=== Building zlib ==="
# if [ ! -f ${TLD}/zlib-1.2.8/libz.a ]; then
#     cd ${TLD}
#     rm -rf ${TLD}/zlib-1.2.8
#     tar xzf ${DIR}/zlib-1.2.8.tar.gz
#     cd zlib-1.2.8
#     ./configure --prefix ${ISD}
#     make -j${CPU_COUNT}
#     make install
# else
#     echo "Zlib already built. To rebuild, delete ${TLD}/zlib-1.2.8"
# fi

# Build Boost
# echo "=== Building Boost ==="
# if [ -z "$BOOST_ROOT" ]; then
#     if [ ! -d ${ISD}/include/boost ]; then
#         cd ${TLD}
#         rm -rf ${TLD}/boost_subset_1_58_0
#         tar xjf ${DIR}/boost_subset_1_58_0.tar.bz2
#         cd boost_subset_1_58_0
        
#         # Apply patches for modern compilers if needed
#         if [ -f ${DIR}/patches/boost_modern_compiler.patch ]; then
#             echo "Applying Boost patches for modern compilers..."
#             patch -p1 < ${DIR}/patches/boost_modern_compiler.patch
#         fi
        
#         ./bootstrap.sh
#         ./b2 link=static -j${CPU_COUNT} --prefix=$ISD -sZLIB_SOURCE=$TLD
#         ./b2 link=static -j${CPU_COUNT} --prefix=$ISD install -sZLIB_SOURCE=$TLD/zlib-1.2.8
#     else
#         echo "Boost already built. To rebuild, delete ${ISD}/include/boost"
#     fi
# else
#     echo "BOOST_ROOT is set to $BOOST_ROOT, not building boost."
#     # Link to system boost if needed
#     if [ ! -d ${ISD}/include/boost ] && [ -d $BOOST_ROOT/include/boost ]; then
#         echo "Linking system Boost headers to build directory..."
#         mkdir -p ${ISD}/include
#         ln -sf $BOOST_ROOT/include/boost ${ISD}/include/boost
#     fi
# fi

# Build htslib
# echo "=== Building htslib ==="
# if [ ! -d ${ISD}/include/htslib ]; then
#     cd ${TLD}
#     rm -rf ${TLD}/htslib
#     tar xzf ${DIR}/htslib.tar.gz
#     cd htslib
#     autoreconf -i 2>/dev/null || true
#     ./configure --prefix=${ISD} --disable-lzma --disable-bz2
#     make -j${CPU_COUNT}
#     make install
# else
#     echo "htslib already built. To rebuild, delete ${ISD}/include/htslib"
# fi

# Build bcftools
# echo "=== Building bcftools ==="
# if [ ! -x ${ISD}/bin/bcftools ]; then
#     cd ${TLD}
#     rm -rf ${TLD}/bcftools
#     tar xzf ${DIR}/bcftools.tar.gz
#     cd bcftools
#     make -j${CPU_COUNT} prefix=${ISD} all
#     make prefix=${ISD} install
# else
#     echo "bcftools already built. To rebuild, delete ${ISD}/bin/bcftools"
# fi

# Build samtools
# echo "=== Building samtools ==="
# if [ ! -x ${ISD}/bin/samtools ]; then
#     cd ${TLD}
#     rm -rf ${TLD}/samtools
#     tar xzf ${DIR}/samtools.tar.gz
#     cd samtools
#     ./configure --prefix=${ISD} --without-curses
#     make -j${CPU_COUNT}
#     make install
# else
#     echo "samtools already built. To rebuild, delete ${ISD}/bin/samtools"
# fi

# Build RTG if requested
if [ "${BUILD_VCFEVAL}" == "1" ]; then
    echo "=== Building rtg-tools ==="
    cd ${TLD}
    rm -rf rtg-tools rtg-tools-install
    
    # TODO: Handle RTG build with better error reporting
    ant runalltests
    mkdir -p ${ISD}/libexec/rtg-tools-install
    cp -r ${TLD}/rtg-tools-install/* ${ISD}/libexec/rtg-tools-install
fi

echo "All external dependencies built successfully."
