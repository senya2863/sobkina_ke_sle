#pragma once
#include <vector>
#include "iterative_methods.h"  // <-- включаем SparseMatrix отсюда

class EllipticGenerator {
public:
    static SparseMatrix generate2D(int nx, int ny);
    static std::vector<double> generateRHS(int nx, int ny);
};