cmake_minimum_required(VERSION 3.1)
project(psiDM2 VERSION 0.1
    DESCRIPTION "DM simulation with QM techniques"
    LANGUAGES C CXX)

# Add modules/ folder for Find* scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")

# Work with MKL for the moment
set(BLA_VENDOR Intel10_64lp)

find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Boost REQUIRED)
find_package(FFTW REQUIRED)
find_package(BLAS REQUIRED)

# message(STATUS ${BLAS_LINKER_FLAGS})
# message(STATUS ${BLAS_LIBRARIES})

# # Create target for FFTW
# message(STATUS "${FFTW_LIBRARIES}")
# message(STATUS "${FFTW_INCLUDE_DIRS}")
# add_library(fftw3 SHARED IMPORTED)
# set_target_properties(fftw3 
#     PROPERTIES IMPORTED_LOCATION "${/usr/local/lib/libfftw3.dylib}"
#     INTERFACE_INCLUDE_DIRECTORIES ${FFTW_INCLUDE_DIRS}
#     INTERFACE_COMPILE_OPTIONS ${FFTW_C_FLAGS}
#     INTERFACE_LINK_OPTIONS ${FFTW_LINKER_FLAGS}
# )

add_executable(psidm 
    src/main.cc 
    src/io.cc 
    src/common.cc 
    src/cosmology.cc
    src/ic.cc
    src/poisson/fft.cc
    src/poisson/fd.cc
    src/schroedinger/uso_dkd.cc
    src/schroedinger/uso_kdk.cc
)

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


