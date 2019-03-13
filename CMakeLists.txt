cmake_minimum_required(VERSION 3.1)
project(psiDM2 VERSION 0.1
    DESCRIPTION "DM simulation with QM techniques"
    LANGUAGES C CXX)

find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Boost REQUIRED)

add_executable(psidm src/main.cc src/io.cc src/common.cc)
target_compile_features(psidm PUBLIC cxx_std_11)

target_include_directories(psidm 
    PUBLIC  ${HDF5_C_INCLUDE_DIRS}
            ${Boost_INCLUDE_DIRS} 
            inc
)

target_link_libraries(psidm ${HDF5_C_LIBRARIES})

