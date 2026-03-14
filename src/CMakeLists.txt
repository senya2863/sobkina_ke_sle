cmake_minimum_required(VERSION 3.16)
project(MatrixProject)

set(CMAKE_CXX_STANDARD 17)


include_directories(include)

add_executable(matrix_test
        src/main.cpp
        src/dense.cpp
        src/csr.cpp
        src/vector_ops.cpp
)