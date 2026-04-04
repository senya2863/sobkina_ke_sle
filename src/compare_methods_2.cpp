#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <algorithm>
#include "chebyshev_mpi.h"
#include "symmetr_methods.h"

//разреженной матрицы
struct element {
    int col;
    double value;
};

struct sparse {
    int n;
    std::vector<std::vector<element>> rows;
};

double norm(const std::vector<double>& a, const std::vector<double>& b) {
    double max_err = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        max_err = std::max(max_err, std::abs(a[i] - b[i]));
    }
    return max_err;
}

// якоби
int jacobi(const sparse& a, const std::vector<double>& b, std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors) {
    int n = a.n;
    std::vector<double> x_new(n);
    
    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& e : a.rows[i]) {
                if (e.col >= n) continue;
                if (e.col == i) diag = e.value;
                else sum += e.value * x[e.col];
            }
            x_new[i] = (b[i] - sum) / diag;
        }
        
        double err = norm(x_new, x);
        errors.push_back(err);
        x = x_new;
        
        if (err < eps) return iter + 1;
    }
    return max_iter;
}

// гаусс-зейдель
int gaussseidel(const sparse& a, const std::vector<double>& b, std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors) {
    int n = a.n;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> old = x;
        
        for (int i = 0; i < n; ++i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& e : a.rows[i]) {
                if (e.col >= n) continue;
                if (e.col == i) diag = e.value;
                else sum += e.value * x[e.col];
            }
            x[i] = (b[i] - sum) / diag;
        }
        
        double err = norm(x, old);
        errors.push_back(err);
        
        if (err < eps) return iter + 1;
    }
    return max_iter;
}

// симметричная трёхдиагональная матрица
sparse createsymmetricmatrix(int n) {
    sparse a;
    a.n = n;
    a.rows.resize(n);
    
    // симметричная матрица: a_ii = 4, a_ij = -1 для соседних
    double alpha = 4.0;
    double beta = -1.0;
    
    for (int i = 0; i < n; ++i) {
        a.rows[i].push_back({i, alpha});
        if (i > 0) {
            a.rows[i].push_back({i-1, beta});
        }
        if (i < n-1) {
            a.rows[i].push_back({i+1, beta});
        }
    }
    
    return a;
}

// правая часть: b = a * x_true
std::vector<double> createrhs(const sparse& a) {
    int n = a.n;
    std::vector<double> b(n, 0.0);
    std::vector<double> x_true(n);
    
    for (int i = 0; i < n; ++i) {
        x_true[i] = (i + 1) * 5.0;
    }
    
    for (int i = 0; i < n; ++i) {
        b[i] = 0.0;
        for (const auto& e : a.rows[i]) {
            if (e.col < n) {
                b[i] += e.value * x_true[e.col];
            }
        }
    }
    
    return b;
}

int main() {
    int n = 500;
    int max_iter = 1000;
    double eps = 1e-6;

    sparse a = createsymmetricmatrix(n);
    std::vector<double> b = createrhs(a);

    std::vector<double> err_jacobi, err_gauss, err_sym_gauss, err_cheb_sym;
    
    //начальные приближения
    std::vector<double> x_jacobi(n, 0.0);
    std::vector<double> x_gauss(n, 0.0);
    std::vector<double> x_sym_gauss(n, 0.0);
    std::vector<double> x_cheb_sym(n, 0.0);
    
    //jacobi
    auto t0 = std::chrono::high_resolution_clock::now();
    int iter_jacobi = jacobi(a, b, x_jacobi, max_iter, eps, err_jacobi);
    auto t1 = std::chrono::high_resolution_clock::now();
    double time_jacobi = std::chrono::duration<double, std::milli>(t1 - t0).count();
    
    //gauss-seidel
    auto t2 = std::chrono::high_resolution_clock::now();
    int iter_gauss = gaussseidel(a, b, x_gauss, max_iter, eps, err_gauss);
    auto t3 = std::chrono::high_resolution_clock::now();
    double time_gauss = std::chrono::duration<double, std::milli>(t3 - t2).count();

    sparsematrix a_mpi;
    a_mpi.n = a.n;
    a_mpi.rows.resize(a.n);
    for (int i = 0; i < a.n; ++i) {
        for (const auto& e : a.rows[i]) {
            a_mpi.rows[i].push_back({e.col, e.value});
        }
    }
    
    //symmetric gauss-seidel
    symmetrizedmethods sym(a_mpi, b);
    auto t4 = std::chrono::high_resolution_clock::now();
    int iter_sym_gauss = sym.solvesymmetricgs(x_sym_gauss, max_iter, eps, err_sym_gauss);
    auto t5 = std::chrono::high_resolution_clock::now();
    double time_sym_gauss = std::chrono::duration<double, std::milli>(t5 - t4).count();
    
    //chebyshev symmetrized
    symmetrizedmethods cheb(a_mpi, b);
    auto t6 = std::chrono::high_resolution_clock::now();
    int iter_cheb_sym = cheb.solvechebyshevsymmetrized(x_cheb_sym, max_iter, eps, err_cheb_sym);
    auto t7 = std::chrono::high_resolution_clock::now();
    double time_cheb_sym = std::chrono::duration<double, std::milli>(t7 - t6).count();

    //csv
    std::ofstream fout("convergence2.csv");
    fout << "iter,jacobi_error,gaussseidel_error,symmetricgs_error,chebyshevsym_error" << std::endl;
    
    size_t max_iters = std::max({err_jacobi.size(), err_gauss.size(), err_sym_gauss.size(), err_cheb_sym.size()});
    
    for (size_t i = 0; i < max_iters; ++i) {
        fout << i + 1 << ",";
        fout << (i < err_jacobi.size() ? std::to_string(err_jacobi[i]) : "") << ",";
        fout << (i < err_gauss.size() ? std::to_string(err_gauss[i]) : "") << ",";
        fout << (i < err_sym_gauss.size() ? std::to_string(err_sym_gauss[i]) : "") << ",";
        fout << (i < err_cheb_sym.size() ? std::to_string(err_cheb_sym[i]) : "") << std::endl;
    }
    fout.close();

    std::ofstream time_out("timing2.csv");
    time_out << "method,iterations,time_ms" << std::endl;
    time_out << "jacobi," << iter_jacobi << "," << time_jacobi << std::endl;
    time_out << "gaussseidel," << iter_gauss << "," << time_gauss << std::endl;
    time_out << "symmetricgs," << iter_sym_gauss << "," << time_sym_gauss << std::endl;
    time_out << "chebyshevsym," << iter_cheb_sym << "," << time_cheb_sym << std::endl;
    time_out.close();
    
    std::cout << "\ndata saved to convergence2.csv and timing2.csv" << std::endl;
    
    return 0;
}