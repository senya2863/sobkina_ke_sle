#pragma once
#include <vector>
#include <cmath>

struct SparseMatrix {
    int n;
    std::vector<std::vector<std::pair<int, double>>> rows;
};

class IterativeMethods {
private:
    SparseMatrix& A;
    std::vector<double> b;

    std::vector<double> multiply(const std::vector<double>& x);
    double dot(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> subtract(const std::vector<double>& a, const std::vector<double>& b);
    
public:
    IterativeMethods(SparseMatrix& matrix, const std::vector<double>& rhs);

    //гаусс
    int gaussSeidel(std::vector<double>& x, int max_iter, double eps, 
                    std::vector<double>& errors);

    //SOR
    int sor(std::vector<double>& x, int max_iter, double eps, 
            std::vector<double>& errors, double omega);

    //наискорейший спуск
    int steepestDescent(std::vector<double>& x, int max_iter, double eps,
                        std::vector<double>& errors);

    //CG
    int conjugateGradient(std::vector<double>& x, int max_iter, double eps,
                          std::vector<double>& errors);
};