cmake_minimum_required(VERSION 3.1)
project(psiDM2 VERSION 0.1
    DESCRIPTION "DM simulation with QM techniques"
    LANGUAGES C CXX)

# Add modules/ folder for Find* scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")

find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Boost REQUIRED)
find_package(FFTW REQUIRED)
find_package(CBLAS)

add_executable(psidm 
    src/main.cc 
    src/state.cc 
    src/cosmology.cc
    src/ic.cc
    src/poisson/fft.cc
    src/poisson/fd.cc
    src/schroedinger/uso_dkd.cc
    src/schroedinger/uso_kdk.cc
    src/schroedinger/pccn.cc
)

# If we find cblas on the target system, we configure blaze to use it
if(CBLAS_FOUND)
    # Each vendor names cblas differently. This queries the header file name
    # TODO Check SMT
    file(GLOB HEADERNAME RELATIVE ${CBLAS_INCLUDE_DIRS} "${CBLAS_INCLUDE_DIRS}/*cblas.h")

    add_definitions(
        -DBLAZE_BLAS_MODE
        -DBLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
        -DBLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION 
        -DBLAZE_BLAS_IS_PARALLEL
        -DBLAZE_BLAS_INCLUDE_FILE="${HEADERNAME}")
endif()

target_compile_features(psidm PUBLIC cxx_std_14)

target_include_directories(psidm 
    PUBLIC  inc
            inc/blaze
            ${HDF5_C_INCLUDE_DIRS}
            ${Boost_INCLUDE_DIRS} 
            ${FFTW_INCLUDE_DIRS} 
)
target_link_libraries(psidm ${HDF5_C_LIBRARIES} ${FFTW_LIBRARIES}
    ${BLAS_LIBRARIES})


