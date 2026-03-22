#pragma once
#include <vector>
#include <cmath>

struct SparseMatrix {
    int n;
    std::vector<std::vector<std::pair<int, double>>> rows;
};

class ChebyshevMPI {
private:
    SparseMatrix& A;
    std::vector<double> b;
    double lambda_min, lambda_max;

    void estimateSpectrum();

public:
    ChebyshevMPI(SparseMatrix& matrix, const std::vector<double>& rhs);
    //мпи с ускорением
    int solve(std::vector<double>& x, int max_iter, double eps,
              std::vector<double>& errors);
    //обычный мпи
    int solveSimple(std::vector<double>& x, int max_iter, double eps,
                    std::vector<double>& errors, double tau);
};