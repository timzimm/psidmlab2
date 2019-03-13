cmake_minimum_required(VERSION 3.1)
project(psiDM2 VERSION 0.1
    DESCRIPTION "DM simulation with QM techniques"
    LANGUAGES C CXX)

find_package(HDF5 REQUIRED COMPONENTS C HL)

add_executable(iotest src/main.cc src/io.cc)
target_compile_features(iotest PUBLIC cxx_std_11)

target_include_directories(iotest 
    PRIVATE ${HDF5_C_INCLUDE_DIRS} inc)

target_link_libraries(iotest PRIVATE ${HDF5_C_LIBRARIES})

