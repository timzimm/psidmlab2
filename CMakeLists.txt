cmake_minimum_required(VERSION 3.13)
project(psiDM2 VERSION 0.1
    DESCRIPTION "DM simulation with QM techniques"
    LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules/")

find_package(hdf5 CONFIG REQUIRED)
find_package(Boost REQUIRED)
find_package(FFTW3 CONFIG REQUIRED)
find_package(blaze CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(CBLAS)
find_package(OpenMP)

add_executable(psidm 
    src/main.cc 
    src/state.cc 
    src/cosmology.cc
    src/ic.cc
    src/observables.cc
    src/convolution_functions.cc
    src/interaction/poisson_fft.cc
    src/interaction/gp.cc
    src/evolution/kinetic.cc
    src/evolution/interaction_potential_trapezodial.cc
    src/evolution/cap_poisson_potential_delta.cc
    src/evolution/concrete_splittings.cc
    src/evolution/uso_kdk_itp.cc
)

target_compile_features(psidm PRIVATE cxx_std_17)

# don't use -isystem, just -I
set_target_properties(psidm PROPERTIES NO_SYSTEM_FROM_IMPORTED TRUE)

#Link time optimization for release build - if supported by the compiler
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
    set_property(TARGET psidm PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
else()
    message(WARNING "IPO is not supported: ${output}")
endif()

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

target_compile_definitions(psidm
    PRIVATE $<$<CONFIG:Release>:
            -DBLAZE_BLAS_MODE
            -DBLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
            -DBLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION 
            -DBLAZE_BLAS_IS_PARALLEL
            -DBLAZE_BLAS_INCLUDE_FILE="${HEADERNAME}">
)

target_compile_options(psidm
    PRIVATE $<$<CONFIG:Release>:
            -g
            -fno-omit-frame-pointer
            -ffast-math
            -march=native>
    PRIVATE $<$<CONFIG:Debug>:
            -fno-omit-frame-pointer
            -fsanitize=undefined
            -fsanitize=address>
)

target_link_options(psidm
    PRIVATE $<$<CONFIG:Debug>:
            -fsanitize=undefined
            -fsanitize=address>
)

target_include_directories(psidm 
    PUBLIC  inc
            $<$<CONFIG:Release>:
            ${CBLAS_INCLUDE_DIRS}>
)

target_link_libraries(psidm 
    PRIVATE ${BLAS_LIBRARIES}
            Boost::boost
            blaze::blaze
            $<$<BOOL:OpenMP_CXX_FOUND>:OpenMP::OpenMP_CXX>
            FFTW3::fftw3
            hdf5::hdf5-static 
            hdf5::hdf5_hl-static
            nlohmann_json::nlohmann_json
)
