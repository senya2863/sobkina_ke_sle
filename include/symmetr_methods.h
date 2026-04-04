#pragma once
#include <vector>
#include <cmath>
#include <utility>

struct sparsematrix {
    int n;
    std::vector<std::vector<std::pair<int, double>>> rows;
};

class symmetrizedmethods {
private:
    sparsematrix& a;
    std::vector<double> b;
    double lambda_min, lambda_max;
    void estimatespectrum();

public:
    symmetrizedmethods(sparsematrix& matrix, const std::vector<double>& rhs);

    int solvesymmetricgs(std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors);
    int solvechebyshevsymmetrized(std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors);
};