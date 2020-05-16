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
find_package(GSL REQUIRED)

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
    FUNCTION_TRACES         OFF
)

include_directories(
    ${CMAKE_SOURCE_DIR}/inc
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)
add_subdirectory(test)
