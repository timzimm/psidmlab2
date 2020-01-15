cmake_minimum_required(VERSION 3.13)
project(psiDM2 VERSION 0.1
    DESCRIPTION "DM simulation with QM techniques"
    LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules/")

# External dependencies
find_package(hdf5 CONFIG REQUIRED)
find_package(Boost REQUIRED)
find_package(FFTW3 CONFIG REQUIRED)
find_package(blaze CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(CBLAS)
find_package(OpenMP)

# Require C++17 for everything
set (CMAKE_CXX_STANDARD 17)

# don't use -isystem, just -I (macOS seems to get confused otherwise)
set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)

# If we find cblas on the target system, we configure blaze to use it
if(CBLAS_FOUND)
    if ("${CBLAS_INCLUDE_DIRS}" STREQUAL "")
        message(FATAL_ERROR "CBLAS was found but its include dir not. Pass it in via
        -DCBLAS_INCLUDE_DIRS")
    endif()
    # Each vendor names cblas differently. This queries the header file name
    # TODO Check SMT
    file(GLOB HEADERNAME RELATIVE "${CBLAS_INCLUDE_DIRS}"
        ${CBLAS_INCLUDE_DIRS}/*cblas.h)

    message(STATUS "CBLAS header found: ${HEADERNAME}")
endif()

set(BLAZE_DEFINITIONS BLAZE_BLAS_MODE
                      BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
                      BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION 
                      BLAZE_BLAS_IS_PARALLEL
                      BLAZE_BLAS_INCLUDE_FILE="${HEADERNAME}")

# TODO: Generator expression for build type selection
add_compile_definitions(${BLAZE_DEFINITIONS})

add_subdirectory(test)
add_subdirectory(src)
