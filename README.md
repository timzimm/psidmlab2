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
    Provided as submodule.
* **[Boost 1.65](http://www.boost.org)**: For root finding, boost::variant etc.
* **[FFTW 3.3.7](http://www.fftw.org)**: For, well, DFTs.
* **[nlohmann's json parser](https://github.com/nlohmann/json)**: 
    For simulation parameter management. Provided as submodule.
* **[HDF5 1.10.0](https://www.hdfgroup.org/solutions/hdf5/)**: 
    A high-performance data management and storage suite for I/O
* **[LAPACK(E)](https://software.intel.com/en-us/mkl)**: 
    Any implementation of LAPACK(E) will do, but if you work on
    Intel chips **MKL** is _highly_ recommended.
* **[BLAS](https://software.intel.com/en-us/mkl)**: 
    For optimized matrix-vector operations that blaze defers to in its
    backend. Again, any implementation is fine
    (**[ATLAS](http://math-atlas.sourceforge.net),
    [openBLAS](https://www.openblas.net), ...**) but
    **MKL** is the way to go, if you are working with Intel chips.
* **C++17 compiler**: like **[clang++](https://llvm.org),
    [ic++](https://software.intel.com/en-us/c-compilers),
    [g++](https://gcc.gnu.org)**
* **cmake 3.10**: build file generator

## How to Install
Installing all of the above is a mess. You don't want this. Trust me.

For the sake of consistency, and to make it easier for you to get started, we
provide a [singularity](https://www.sylabs.io/singularity/) container defintion file
as well as a prebuild container for x86_64 hosted on Sylab's Cloud Library. 
The container includes all the dependencies mentioned above and can directly be used 
for all development and execution purposes. MKL is preinstalled
which makes it quite heavy in terms of disk space.

**Their is no need to install any dependencies on your own and performance will not 
be harmed by running psiDMLab from inside the container**

To get psiDMLab start by cloning this repo and its submodules:
```bash
# ssh ...
$ git clone --recurse-submodules git@gitlab.com:ttz/psidm2.git
# ...or https
$ git clone --recurse-submodules https://gitlab.com/ttz/psidm2.git
```
Next, install singularity by following 
[these](https://www.sylabs.io/guides/3.2/user-guide/installation.html#) instructions.
### Linux
Pull down the development container from Sylabs Cloud:
```bash
$ singularity pull <your container name> library://ttz/psi_dm_lab/dev
```

### MacOS
Singularity is still in its alpha stage on macOS. Currently some sort of
file permission issues stop me from compiling files using the software packaged
in the container. I opened an [issue](https://github.com/sylabs/singularity/issues/3636) 
for that purpose.

### Windows
Don't know. Install Linux :)

That's it! From now on compiling and running the code will happen by invoking
commands on the container.

## How to Compile
We assume the container resides in the projects root and is called dev, i.e.
```bash
$ ls
CMakeLists.txt          dev                    run                  tags
README.md               inc                    script
install                 src                    modules                     
```
To compile the code create a out-of-source build folder:
```bash
$ mkdir build; cd build
```
Now invoke cmake inside the container to create a makefile and start building
```bash
$ singularity exec dev cmake -DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=Intel10_64lp ..
$ singularity exec dev make
```

## How to Run
TODO

## A Worked Example
TODO

## psiDM in the literature
TODO
