# Общая библиотека: матрицы + векторные операции + триди-решатель
add_library(MatrixLib
    csr.cpp
    dense.cpp
    vector_ops.cpp
    tridiagonal.cpp
)
target_include_directories(MatrixLib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Бенчмарк: dense vs CSR умножение
add_executable(matrix_bench matrix_bench.cpp)
target_link_libraries(matrix_bench PRIVATE MatrixLib)

# QR-разложение (Хаусхолдер)
add_executable(qr qr.cpp)
target_link_libraries(qr PRIVATE MatrixLib)

# Итерационные методы: Jacobi / Gauss-Seidel / простая итерация
add_executable(gauss_jacobi gauss_jacobi.cpp)
target_link_libraries(gauss_jacobi PRIVATE MatrixLib)

#ДЗ 5: Чебышев и сравнение методов
add_executable(compare_methods compare_methods.cpp chebyshev_mpi.cpp)
target_link_libraries(compare_methods PRIVATE MatrixLib)

#ДЗ 6
# сравнение с симметризованными методами
add_executable(compare_methods_2 compare_methods_2.cpp
        symmetr_methods.cpp chebyshev_mpi.cpp)
target_link_libraries(compare_methods_2 PRIVATE MatrixLib)