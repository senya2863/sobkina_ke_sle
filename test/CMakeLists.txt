find_package(GTest REQUIRED)

add_executable(tridiagonal_test tridiagonal_test.cpp)
target_link_libraries(tridiagonal_test PRIVATE MatrixLib GTest::gtest GTest::gtest_main)
target_include_directories(tridiagonal_test PRIVATE ${CMAKE_SOURCE_DIR}/include)

include(GoogleTest)
gtest_discover_tests(tridiagonal_test)
