#!/bin/bash
#
# Detect common shell variables for testing
#
# Author: Peter Krusche <pkrusche@illumina.com>

if [[ -d "/illumina/thirdparty/graphviz/latest/bin" ]]; then
	export PATH=$PATH:/illumina/thirdparty/graphviz/latest/bin
fi

# fallback HG19 locations
if [[ ! -f "$HG19" ]]; then
	# HGREF variable can point to reference (hap.py asks for this, see src/python/Tools/__init__.py)
	export HG19="${HGREF}"
fi

if [[ ! -f "$HG19" ]]; then
	export HG19=/illumina/development/iSAAC/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
fi

if [[ ! -f "$HG19" ]]; then
	export HG19=~/workspace/human_genome/hg19.fa
fi

# fallback for finding hc binaries
if [[ ! -f "${HCDIR}/hap.py" ]]; then
	export HCDIR="$(pwd)/bin"
fi

if [[ ! -f "${HCDIR}/hap.py" ]]; then
	export HCDIR="$(pwd)/build/bin"
fi

if [[ ! -f "${HCDIR}/hap.py" ]]; then
	export HCDIR="$(dirname $(pwd))/build/bin"
fi

if [[ ! -f "${HCDIR}/hap.py" ]]; then
	# Try to find it relative to this script
	export HCDIR="$(cd ${DIR}/../../build/bin && pwd)"
fi

if [[ ! -f "${HCDIR}/hap.py" ]]; then
	echo "Cannot find HC binaries. Set HCDIR to the bin directory."
	exit 1
fi

export PATH="$PATH:${HCDIR}"

# htslib is a shared object
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${HCDIR}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HCDIR}/lib

# detect or use python
if [[ -z ${PYTHON} ]]; then
	DEFAULT_PYTHON=1
	if [ -f "/illumina/sync/software/groups/hap.py/latest/python-ve/bin/python-wrapper.sh" ]; then
	    export PYTHON=/illumina/sync/software/groups/hap.py/latest/python-ve/bin/python-wrapper.sh
	fi
else
	DEFAULT_PYTHON=0
fi

export PYTHON=${PYTHON:-python}

PYVERSION=$(${PYTHON} --version 2>&1)
if [[ "$PYVERSION" != "Python 2.7."* ]] && [[ "$PYVERSION" != "Python 3."* ]] && [[ $DEFAULT_PYTHON == 1 ]]; then
	PYTHON=python3
	if ! command -v $PYTHON &> /dev/null; then
		PYTHON=python2.7
	fi
fi


export HCVERSION=`${PYTHON} ${HCDIR}/hap.py --version 2>/dev/null || echo "unknown"`
if [[ "$HCVERSION" == "" ]]; then
    echo "Warning: Cannot run hap.py to extract version information. Continuing anyway."
    export HCVERSION="unknown"
fi
