cmake_minimum_required(VERSION 3.9)
project(psiDM2 VERSION 0.1
    DESCRIPTION "DM simulation with QM techniques"
    LANGUAGES C CXX)

# Add modules/ folder for Find* scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")

# TODO ADD LAPACKE. Right now we implicitly assume that BLAS comes with LAPACK
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Boost REQUIRED)
find_package(FFTW REQUIRED)
find_package(CBLAS)
find_package(nlohmann_json 3.2.0 REQUIRED)

add_executable(psidm 
    src/main.cc 
    src/state.cc 
    src/cosmology.cc
    src/ic.cc
    src/observables.cc
    src/convolution_functions.cc
    src/poisson/fft.cc
    src/poisson/fd.cc
    src/poisson/identity.cc
    src/schroedinger/uso_dkd.cc
    src/schroedinger/uso_kdk.cc
    src/schroedinger/pccn.cc
    src/schroedinger/pc_cayley.cc
)

# If we find cblas on the target system, we configure blaze to use it
if(CBLAS_FOUND)
    # Each vendor names cblas differently. This queries the header file name
    # TODO Check SMT
    file(GLOB HEADERNAME RELATIVE ${CBLAS_INCLUDE_DIRS} "${CBLAS_INCLUDE_DIRS}/*cblas.h")

endif()

# We use if constexpr
target_compile_features(psidm PRIVATE cxx_std_17)

# BLAS backend on release
target_compile_definitions(psidm
    PRIVATE $<$<CONFIG:Release>:
            -DBLAZE_BLAS_MODE
            -DBLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
            -DBLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION 
            -DBLAZE_BLAS_IS_PARALLEL
            -DBLAZE_BLAS_INCLUDE_FILE="${HEADERNAME}">
)

#Link time optimization for release build - if supported by the compiler
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
    set_property(TARGET psidm PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
else()
  message(WARNING "IPO is not supported: ${output}")
endif()

target_compile_options(psidm
    PRIVATE $<$<CONFIG:Release>:
            -Ofast
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
            inc/blaze
            ${HDF5_C_INCLUDE_DIRS}
            ${Boost_INCLUDE_DIRS} 
            ${FFTW_INCLUDE_DIRS} 
)

target_link_libraries(psidm 
    PRIVATE ${HDF5_C_LIBRARIES} 
            ${FFTW_LIBRARIES}
            ${BLAS_LIBRARIES}
            nlohmann_json::nlohmann_json
)


