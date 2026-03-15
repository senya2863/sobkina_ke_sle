#!/bin/bash
# reorganize.sh — запускать из корня репозитория
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "=== Reorganizing repo ==="

# ── 1. Временная папка — собираем всё что надо сохранить ─────────
mkdir -p _tmp/include _tmp/src _tmp/test _tmp/measures

# Заголовки
cp src/src/csr.h            _tmp/include/csr.h
cp src/src/dense.h          _tmp/include/dense.h
cp src/src/vector_ops.h     _tmp/include/vector_ops.h
cp tridiagonal/src/tridiagonal.h  _tmp/include/tridiagonal.h

# Исходники
cp src/src/csr.cpp              _tmp/src/csr.cpp
cp src/src/dense.cpp            _tmp/src/dense.cpp
cp src/src/vector_ops.cpp       _tmp/src/vector_ops.cpp
cp tridiagonal/src/tridiagonal.cpp  _tmp/src/tridiagonal.cpp
cp src/src/main.cpp             _tmp/src/matrix_bench.cpp
cp mpi_gauss_jacobi/main.cpp   _tmp/src/gauss_jacobi.cpp
cp qr/main.cpp                  _tmp/src/qr.cpp

# Тесты
cp tridiagonal/tests/tridiagonal_test.cpp  _tmp/test/tridiagonal_test.cpp

# Данные и графики
cp mpi_gauss_jacobi/plot.py  _tmp/measures/gauss_jacobi_plot.py
cp src/src/plot.py           _tmp/measures/matrix_plot.py
cp src/src/times.txt         _tmp/measures/times.txt
cp src/src/graph.png         _tmp/measures/graph.png
cp src/src/graph.pdf         _tmp/measures/graph.pdf

# ── 2. Убрать старые папки из git-индекса (файлы остаются на диске) ─
echo ""
echo "=== Removing build/idea folders from git index ==="
for d in \
  mpi_gauss_jacobi/build \
  mpi_gauss_jacobi/cmake-build-debug \
  mpi_gauss_jacobi/.idea \
  qr/cmake-build-debug \
  qr/.idea \
  src/build \
  src/cmake-build-debug \
  src/.idea \
  tridiagonal/build \
  tridiagonal/cmake-build-debug \
  tridiagonal/.idea; do
  git rm -r --cached "$d" 2>/dev/null && echo "  removed from index: $d" || true
done

# ── 3. Удалить старые папки с диска ─────────────────────────────
rm -rf mpi_gauss_jacobi qr src tridiagonal

# ── 4. Поставить новую структуру из temp ─────────────────────────
mv _tmp/include .
mv _tmp/src     .
mv _tmp/test    .
mv _tmp/measures .
rm -rf _tmp

# ── 5. Записать CMakeLists.txt ───────────────────────────────────
cat > CMakeLists.txt << 'EOF'
cmake_minimum_required(VERSION 3.16)
project(sobkina_ke_sle CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_SOURCE_DIR}/include)

enable_testing()

add_subdirectory(src)
add_subdirectory(test)
EOF

cat > src/CMakeLists.txt << 'EOF'
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
EOF

cat > test/CMakeLists.txt << 'EOF'
find_package(GTest REQUIRED)

add_executable(tridiagonal_test tridiagonal_test.cpp)
target_link_libraries(tridiagonal_test PRIVATE MatrixLib GTest::gtest GTest::gtest_main)
target_include_directories(tridiagonal_test PRIVATE ${CMAKE_SOURCE_DIR}/include)

include(GoogleTest)
gtest_discover_tests(tridiagonal_test)
EOF

# ── 6. .gitignore ────────────────────────────────────────────────
cat > .gitignore << 'EOF'
# Build artifacts
build/
cmake-build-debug/
cmake-build-release/

# CLion
.idea/

# Compiled
*.o
*.a
*.out

# CMake generated
CMakeCache.txt
CMakeFiles/
cmake_install.cmake
Makefile
build.ninja
*.ninja
CTestTestfile.cmake
Testing/
EOF

# ── 7. git add ───────────────────────────────────────────────────
git add .

echo ""
echo "======================================"
echo "  Готово! Структура:"
echo ""
echo "  include/   — заголовки (.h)"
echo "  src/       — исходники (.cpp) + CMakeLists.txt"
echo "  test/      — тесты + CMakeLists.txt"
echo "  measures/  — графики, CSV"
echo "  CMakeLists.txt"
echo "  .gitignore"
echo ""
echo "  Сборка:"
echo "    mkdir build"
echo "    cmake -S . -B build"
echo "    cmake --build build"
echo ""
echo "  Тесты:"
echo "    cd build && ctest"
echo ""
echo "  Запуск:"
echo "    ./build/src/gauss_jacobi"
echo "    ./build/src/qr"
echo "    ./build/src/matrix_bench"
echo ""
echo "  Закоммить:"
echo "    git commit -m 'refactor: restructure into single project'"
echo "    git push"
echo "======================================"