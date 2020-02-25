cmake_minimum_required(VERSION 3.13)
project(psiDM2 VERSION 0.1
    DESCRIPTION "DM simulation with QM techniques"
    LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules/")

# Command line options
option(PSIDMLAB_SMP "Enable OpenMP for Blaze and FFTW3")

find_package(Boost REQUIRED)
find_package(hdf5 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(blaze CONFIG REQUIRED)

# MKL ships with FFTW3 wrappers. In this case we can omit the external FFTW
# dependency. Note that MKL does not provide all the functionality of FFTW.
# See https://software.intel.com/en-us/mkl-developer-reference-c-using-fftw3-wrappers
if(NOT DEFINED ENV{MKLROOT})
    find_package(FFTW3 CONFIG REQUIRED)
endif()

# Get rid of any predefined OpenMP dependencies. Default blaze is now serial.
set_property(TARGET blaze::blaze 
    PROPERTY INTERFACE_LINK_LIBRARIES "")
set_property(TARGET blaze::blaze 
    PROPERTY INTERFACE_COMPILE_OPTIONS "")

# Blaze configuration. 
# This adds compile definitions to blaze::blaze in order to override the 
# defaults set in blaze/config/*.h. 
# It is a modifaction of Blaze_Import shipped with blaze.
if(NOT DEFINED ENV{MKLROOT})
    set(CBLAS_HEADER <cblas.h>)
else()
    set(CBLAS_HEADER <mkl_cblas.h>)
endif()
if(${CMAKE_BUILD_TYPE} STREQUAL "Release" AND ${PSIDMLAB_SMP})
    set(SMP_BACKEND OpenMP)
else()
    set(SMP_BACKEND off)
endif()

include(blaze_config)
Blaze_Config(
    BLAS                    ON
    BLAS_PARALLEL           ON
    BLAS_MV                 ON
    BLAS_MM                 ON
    BLAS_INCLUDE            ${CBLAS_HEADER}
    LAPACK                  ON
    THREADING               ${SMP_BACKEND}
    CACHE_SIZE              auto
    VECTORIZATION           ON
    TRANSPOSE_FLAG          columnVector
    STORAGE_ORDER           rowMajor
)

add_executable(psidm 
    src/main.cc 
    src/state.cc 
    src/fftw.cc 
    src/cosmology.cc
    src/ic.cc
    src/observables.cc
    src/convolution_functions.cc
    src/interaction/poisson_fft.cc
    src/evolution/kinetic.cc
    src/evolution/interaction_potential_trapezodial.cc
    src/evolution/cap_poisson_potential_delta.cc
    src/evolution/concrete_splittings.cc
)

target_compile_features(psidm PRIVATE cxx_std_17)

# don't use -isystem, just -I
set_target_properties(psidm PROPERTIES 
    NO_SYSTEM_FROM_IMPORTED TRUE
)

target_compile_options(psidm
    PRIVATE $<$<CONFIG:Release>:
            -fno-omit-frame-pointer
            -ffast-math
            -march=native>
            $<$<CONFIG:Debug>:
            -fno-omit-frame-pointer
            -fsanitize=undefined
            -fsanitize=address>
)

target_link_options(psidm
    PRIVATE $<$<CONFIG:Debug>:
            -fsanitize=undefined
            -fsanitize=address>
)

if(${PSIDMLAB_SMP})
    target_compile_definitions(psidm PRIVATE PSIDMLAB_SMP) 
endif()

target_include_directories(psidm PRIVATE inc)

target_link_libraries(psidm 
            Boost::boost
            blaze::blaze
            hdf5::hdf5-static 
            hdf5::hdf5_hl-static
            nlohmann_json::nlohmann_json
)
