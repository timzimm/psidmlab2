#!/bin/bash
# Installs (i.e. downloads and potantially compiles) all 3rd party dependencies.
# We assume vcpkg is already installed and integrated. Make sure <VCKPG_ROOT> is
# in $PATH.

vcpkg install boost-math
vcpkg install boost-variant 
vcpkg install boost-core 
vcpkg install fftw3
vcpkg install hdf5
vcpkg install --head blaze
vcpkg install nlohmann-json
vcpkg install intel-mkl
