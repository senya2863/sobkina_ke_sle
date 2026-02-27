cmake_minimum_required(VERSION 3.16)
project(TridiagonalMatrix)

set(CMAKE_CXX_STANDARD 17)

find_package(GTest REQUIRED)
enable_testing()

add_library(DiagThree src/tridiagonal.cpp)
target_include_directories(DiagThree PUBLIC src)

add_executable(DiagThreeTest tests/tridiagonal_test.cpp)
target_link_libraries(DiagThreeTest PRIVATE DiagThree GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(DiagThreeTest)