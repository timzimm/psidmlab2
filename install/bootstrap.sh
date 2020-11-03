#!/bin/bash
# Installs (i.e. downloads and potantially compiles) all 3rd party dependencies.

HOME=$(pwd)

command -v vcpkg > /dev/null
if [[ $? == 1 ]]; then
    echo "Please add vcpkg to your PATH"
    exit 1
fi
VCPKG_REPO=$(dirname $(which vcpkg))
cd $VCPKG_REPO
# up-to-date dependencies
./vcpkg install boost-math
./vcpkg install boost-variant 
./vcpkg install boost-units 
./vcpkg install hdf5
./vcpkg install gsl
./vcpkg install gtest

# nlohmann json version 3.7.3
git checkout ef90276f8d8032e932d94090240d16a6801470d6 -- ports/nlohmann-json
./vcpkg install nlohmann-json

# Blas Library
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cat /proc/cpuinfo | grep Intel > /dev/null
    if [[ $? == 0 ]]; then
        # Anything newer than MKL 2018 will do
        git checkout ea2e5d7b006fc5bf6adc8789c5929d638bbf3e97 -- ports/intel-mkl
        ./vcpkg install intel-mkl
    else
        ./vcpkg install openblas
    fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
    # We live in a pre-Apple-Silicon world, so Intel is assumed implictely
    # Anything newer than MKL 2018 will do
    git checkout ea2e5d7b006fc5bf6adc8789c5929d638bbf3e97 -- ports/intel-mkl
    ./vcpkg install intel-mkl
fi

# blaze 3.8
git checkout e78db0b6d7efb3e2072417274ccd1903072de28f -- ports/blaze
./vcpkg install blaze

# MKL ships with FFTW3 interface
./vcpkg list | grep intel-mkl > /dev/null
if [[ $? == 1 ]]; then
    ./vcpkg install fftw3[avx2]
fi

cd $HOME
