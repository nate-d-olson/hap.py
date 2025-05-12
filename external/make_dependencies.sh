#!/bin/bash

set -e

echo "This script is deprecated. Dependencies are now handled by CMake using FetchContent."
echo "Please build the project using standard CMake commands:"
echo "  mkdir build && cd build"
echo "  cmake .."
echo "  cmake --build ."

exit 0

if [ ! -d ${ISD}/include/htslib ] ;
then
    cd ${TLD}
    rm -rf ${TLD}/htslib
    tar xzf ${DIR}/htslib.tar.gz
    cd htslib
    ./configure --prefix=${ISD} \
        CFLAGS=-I${ISD}/include\ -g \
        CXXFLAGS=-I${ISD}/include\ -g \
        LDFLAGS=-L${ISD}/lib \
        --disable-plugins \
        --disable-libcurl \
        --disable-lzma \
        --disable-bz2
    make -j4
    make -j4 install

    # windows shared folder workaround
    if [ -e "${ISD}/lib/libhts.so" ] && [ ! -e "${ISD}/lib/libhts.so.1" ] ;
    then
        cp ${ISD}/lib/libhts.so ${ISD}/lib/libhts.so.1
    fi
else
    echo "HTSLIB already built. To rebuild, delete ${ISD}/include/htslib"
fi

if [ ! -f ${ISD}/bin/bcftools ];
then
    cd ${TLD}
    rm -rf ${TLD}/bcftools
    tar xzf ${DIR}/bcftools.tar.gz
    cd bcftools
    make -j4 prefix=${ISD}
    make -j4 prefix=${ISD} install
else
    echo "bcftools already built. To rebuild, delete ${ISD}/bin/bcftools"
fi

if [ ! -f ${ISD}/bin/samtools ];
then
    cd ${TLD}
    rm -rf ${TLD}/samtools
    tar xzf ${DIR}/samtools.tar.gz
    cd samtools
    autoconf -Wno-syntax || autoconf -Wno-syntax
    ./configure --prefix=${ISD} \
        --with-htslib=${TLD}/htslib \
        --without-curses \
        CFLAGS=-I${ISD}/include \
        CPPFLAGS=-I${ISD}/include \
        LDFLAGS=-L${ISD}/lib
    make -j4
    make -j4 install
else
    echo "samtools already built. To rebuild, delete ${ISD}/bin/samtools"
fi

# get vcfeval
# https://github.com/RealTimeGenomics/rtg-tools/archive/ga4gh-test.zip
if [[ ! -z $BUILD_VCFEVAL ]]; then
    if [[ ! -d ${ISD}/libexec/rtg-tools-install ]]; then
        echo "Building rtg-tools"
        cd ${TLD}
        mkdir -p ${TLD}/rtg-tools
        cd rtg-tools
        wget http://github.com/RealTimeGenomics/rtg-tools/archive/3.12.1.tar.gz -O ${TLD}/rtg-tools/rtg-tools.tar.gz
        tar xvf rtg-tools.tar.gz
        cd rtg-tools-3.12.1

        if [[ ! -z ${ANT_HOME} ]]; then
            $ANT_HOME/bin/ant zip-nojre
        else
            ant zip-nojre
        fi
        cd ..

        RTG_ZIPFILE=$(ls rtg-tools-3.12.1/dist/*-nojre.zip | head -1)
        RTG_BASE=$(basename $RTG_ZIPFILE -nojre.zip)
        jar xvf $RTG_ZIPFILE
        mv $RTG_BASE ${ISD}/libexec/rtg-tools-install
        cp ${DIR}/rtg.cfg ${ISD}/libexec/rtg-tools-install
        chmod +x ${ISD}/libexec/rtg-tools-install/rtg
    else
        echo "rtg-tools is already built. To rebuild, delete ${ISD}/libexec/rtg-tools-install"
    fi
    if [[ -f ${VCFEVAL_WRAPPER} ]]; then
        echo "using wrapper for rtg-tools: ${VCFEVAL_WRAPPER}"
        cp ${VCFEVAL_WRAPPER} ${ISD}/libexec/rtg-tools-install/rtg-wrapper.sh
        chmod +x ${ISD}/libexec/rtg-tools-install/rtg-wrapper.sh
    fi
fi

