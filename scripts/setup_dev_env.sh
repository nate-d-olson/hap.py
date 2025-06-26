#!/bin/bash
# This script sets up the development environment for hap.py.
# It uses micromamba for external tools and a Python venv for libraries.

set -e

# --- Helper Functions ---
print_info() {
    echo "INFO: $1"
}

print_error() {
    echo "ERROR: $1" >&2
    exit 1
}

# --- Micromamba Setup ---
print_info "Setting up micromamba for external dependencies..."

if ! command -v micromamba &> /dev/null; then
    print_error "Micromamba is not installed. Please install it first."
fi

# Detect OS architecture
ARCH=$(uname -m)
OS=$(uname -s)

SUBDIR=""
if [ "$OS" = "Darwin" ]; then
    if [ "$ARCH" = "arm64" ]; then
        SUBDIR="osx-arm64"
    else
        SUBDIR="osx-64"
    fi
elif [ "$OS" = "Linux" ]; then
    if [ "$ARCH" = "aarch64" ]; then
        SUBDIR="linux-aarch64"
    else
        SUBDIR="linux-64"
    fi
fi

print_info "Cleaning micromamba cache..."
micromamba clean -a -y

print_info "Uninstalling existing rtg-tools and bcftools..."
micromamba uninstall -y -n base rtg-tools bcftools openjdk || true

print_info "Installing rtg-tools and bcftools with micromamba for platform $SUBDIR..."
if [ -n "$SUBDIR" ]; then
    micromamba install -y -n base -c bioconda -c conda-forge --platform $SUBDIR rtg-tools bcftools
else
    micromamba install -y -n base -c bioconda -c conda-forge rtg-tools bcftools
fi


# --- Python Virtual Environment Setup ---
print_info "Setting up Python virtual environment..."

if [ -d ".venv" ]; then
    print_info "Virtual environment .venv already exists. Skipping creation."
else
    python3 -m venv .venv
fi

print_info "Activating the virtual environment..."
source .venv/bin/activate

# --- Python Dependencies Installation ---
print_info "Installing Python dependencies from pyproject.toml..."
pip install -e ".[dev,bio,viz]"

print_info "Development environment setup complete."
