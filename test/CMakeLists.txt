enable_testing()
find_package(GTest MODULE REQUIRED)
 
link_libraries(GTest::GTest GTest::Main)

add_executable(numbers numbers.cc)
add_executable(strings strings.cc)
 
gtest_discover_tests(strings)
gtest_discover_tests(numbers)
