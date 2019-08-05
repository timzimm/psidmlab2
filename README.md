# psiDMLab

## What is psiDMLab

psiDMLab is meant as a numerical playground for the Schroedinger-Poisson equations, 
a model that utilises quantum mechanical techniques to approximate the numerically 
intractable Vlasov-Poisson system for the evolution of dark matter. 
Currently, our focus lies on modularity not utmost performance. The idea is to draft 
a numerically stable and accurate scheme that can later on be optimized to run 
on large scale systems and hopefully is capable of simulating 3+1 dimensional 
cosmological setups comparable to what is known for N-Body simulations.

## Theory behind psiDMLab

psiDM aims to model the dynamics of (cold) dark matter by employing a time evolution
described by the Schroedinger-Poisson (SP) PDE system:
```math
    i \partial_t \psi(x,t) = \left[-\frac{\hbar^2}{2m} \nabla^2 +
    mV(x,t)\right]\psi(x,t)
```
```math
    \nabla^2 V = \frac{4\pi G\rho_{m0}}{a}\left(\frac{|\psi|^2}{\langle|\psi|^2\rangle_V } - 1\right)
```
From a cosmological point of view, two distinct applications of wavelike dark
matter can be identified.
### Ultralight Axion Dark Matter
Axions are considered a viable dark matter particle candidate. In this
framework, axions (spin 0 bosons) form a cosmic Bose-Einstein condensate that
obeys SP, i.e. evolves quantum mechanically under its self gravity. For
sufficiently small axion masses, structure formation on the kpc-scale is controlled
by quantum mechanical effects such as Heissenberg's uncertainty princlple acting
as a "regularizer" of the emerging density structures. This solves long standing
difficulites of the standard $`\Lambda`$CDM model, such as the core cusp problem
or the predicition of diverging densities in cosmic halos. On scales vastly
larger than the axion de-Broglie wavelength $`\psi`$DM recovers the standard CDM
description.
### Approximation to the Vlasov-Poisson system (VP)
Standard CDM obeys a collisionless Boltzmann equation (also known as Vlasov
equation in plasma physics) coupled to Poisson's
equation. Assuming a collisionless setting, this can be shown to be true even
for the one-particle distribution function $`f(x,p,t)`$. Unfortunately, solving
VP is hard both numerically and analytically. Therefore, people (almost always)
rely on N-Body simulations yielding a noisy sample from the true phase space
distributionn. Partial analytical results, however, show a good correspondence
between SP and VP which is not fully understood. Numerical results support those
analytical findings. For instance, the associated potential of the wavefunction
was shown to converge to the classical potential $`\propto{\left(\frac{\hbar}{m}\right)^2}`$.
Moreover, choosing Husimi's distribution to compute the associated phase
space distribution of $`\psi`$ to obtain a "smoothed eveolution" produces
distribution moments (such as matter densities or veloctiy dispersions) that
resemble closely what we expect in the classical case.
### 

## What Dependenicies does the Simulation have
Currently psiDMLab depends on:
* **[BLAZE (HEAD)](https://bitbucket.org/blaze-lib/blaze/src/master/)**: 
    A (smart) expression template based, high performance linear algebra library written in C++14. 
* **[Boost 1.70](http://www.boost.org)**: For root finding, boost::variant etc.
* **[FFTW 3.3.8](http://www.fftw.org)**: For, well, DFTs.
* **[nlohmann's json parser](https://github.com/nlohmann/json)**: 
    For simulation parameter management.
* **[HDF5 1.10.5](https://www.hdfgroup.org/solutions/hdf5/)**: 
    A high-performance data management and storage suite for I/O
* **[Intel MKL](https://software.intel.com/en-us/mkl)**: 
    For BLAS and LAPACK(E).
* **C++17 compiler**: like **[clang++](https://llvm.org),
    [ic++](https://software.intel.com/en-us/c-compilers),
    [g++](https://gcc.gnu.org)**
* **cmake 3.14**: build file generator

## How to Install
Installing all of the above manually is a mess. You don't want this. Trust me.
The most common case is that you develop locally, but run large scale
simulations on a cluster. Currently both "modes of operation" are dealt
differently when it comes to installation. This will change in the future when
[singularity](https://www.sylabs.io/singularity/) for macOS is less buggy.

### Development
We provide a bootstrap script that works together with [vcpkg](https://github.com/microsoft/vcpkg)
to install all third party dependencies. vcpkg works on Linux, macOS and
Windows. Note that you still need to provide a C++17 compiler and cmake. Moreover, make
sure to install [Intel MKL](https://software.intel.com/en-us/mkl) before you
procede with the steps below.

That said, start by installing vcpkg:
```bash
~$ git clone https://github.com/Microsoft/vcpkg.git
~$ cd vcpkg
~$ ./bootstrap-vcpkg.sh
~$ ./vcpkg integrate install
~$ echo 'export PATH=~/vcpkg:$PATH' >> ~/.bashrc
```

To get psiDMLab, clone this repo and run the bootstrap script:
```bash
~$ cd
~$ git clone git@gitlab.com:ttz/psidm2.git
~$ cd psidm2/install
~/psidm2/install$ source /opt/intel/mkl/bin/mklvars.sh intel64
~/psidm2/install$ ./bootstrap.sh
# Coffee break!
```

That's it. You can now generate the makefile with cmake by executing
```bash
~/psidm2$ mkdir build
~/psidm2$ cmake -DCMAKE_TOOLCHAIN_FILE=~/vcpkg/scripts/buildsystems/vcpkg.cmake \
                -DCMAKE_BUILD_TYPE=Release . -B build
~/psidm2$ cd build
~/psidm2/build$ make
```

### Cluster
TODO

## How to Run
TODO

## A Worked Example
TODO

## psiDM in the literature
TODO
