enable_testing()
find_package(GTest MODULE REQUIRED)
 
include_directories(
    ${CMAKE_SOURCE_DIR}/inc
    $<$<CONFIG:Release>:${CBLAS_INCLUDE_DIRS}>
)

link_libraries(GTest::GTest GTest::Main)

add_executable(numbers numbers.cc)
add_executable(strings strings.cc)
 
gtest_discover_tests(strings)
gtest_discover_tests(numbers)
