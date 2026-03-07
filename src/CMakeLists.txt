cmake_minimum_required(VERSION 3.16)
project(MatrixProject)

set(CMAKE_CXX_STANDARD 17)

add_executable(MatrixTest
        src/dense.cpp
        src/csr.cpp
        src/matrix_test.cpp
)