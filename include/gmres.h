#pragma once
#include <vector>
#include <cmath>
#include "iterative_methods.h"

class GMRES {
private:
    SparseMatrix& A;
    std::vector<double> b;
    int restart;

    std::vector<double> multiply(const std::vector<double>& x);
    double dot(const std::vector<double>& a, const std::vector<double>& b);
    double norm(const std::vector<double>& x);
    std::vector<double> subtract(const std::vector<double>& a, const std::vector<double>& b);

public:
    GMRES(SparseMatrix& matrix, const std::vector<double>& rhs, int restart = 30);

    int solve(std::vector<double>& x, int max_iter, double eps, 
              std::vector<double>& errors);
};