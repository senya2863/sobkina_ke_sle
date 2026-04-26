#pragma once
#include <vector>
#include "iterative_methods.h"

class IncompleteCholesky {
public:
    static SparseMatrix decompose(const SparseMatrix& A);
    static std::vector<double> solveLower(const SparseMatrix& L, const std::vector<double>& b);
    static std::vector<double> solveUpper(const SparseMatrix& L, const std::vector<double>& y);
};