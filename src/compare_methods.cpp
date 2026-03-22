#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <algorithm>
#include "chebyshev_mpi.h"

// Структуры для разреженной матрицы (как в gauss_jacobi.cpp)
struct Element {
    int col;
    double value;
};

struct Sparse {
    int n;
    std::vector<std::vector<Element>> rows;
};

// Норма (как в gauss_jacobi.cpp)
double norm(const std::vector<double>& a, const std::vector<double>& b) {
    double max_err = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        max_err = std::max(max_err, std::abs(a[i] - b[i]));
    }
    return max_err;
}

// Якоби
int jacobi(const Sparse& A, const std::vector<double>& b, std::vector<double>& x,
           int max_iter, double eps, std::vector<double>& errors) {
    int n = A.n;
    std::vector<double> x_new(n);

    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& e : A.rows[i]) {
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

// Гаусс-Зейдель
int gaussSeidel(const Sparse& A, const std::vector<double>& b, std::vector<double>& x,
                int max_iter, double eps, std::vector<double>& errors) {
    int n = A.n;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> old = x;

        for (int i = 0; i < n; ++i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& e : A.rows[i]) {
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

// Трёхдиагональная матрица с плохой обусловленностью
Sparse createTridiagonalMatrix(int n) {
    Sparse A;
    A.n = n;
    A.rows.resize(n);

    double alpha = 2.0;   // главная диагональ
    double beta = -0.99;  // побочные (близко к -1 для медленной сходимости)

    for (int i = 0; i < n; ++i) {
        A.rows[i].push_back({i, alpha});
        if (i > 0) {
            A.rows[i].push_back({i-1, beta});
        }
        if (i < n-1) {
            A.rows[i].push_back({i+1, beta});
        }
    }
    return A;
}

// Правая часть: b = A * x_true, где x_true далеко от нуля
std::vector<double> createRHS(const Sparse& A) {
    int n = A.n;
    std::vector<double> b(n, 0.0);
    std::vector<double> x_true(n);

    for (int i = 0; i < n; ++i) {
        x_true[i] = (i + 1) * 10.0;  // 10, 20, 30, ..., n*10
    }

    for (int i = 0; i < n; ++i) {
        b[i] = 0.0;
        for (const auto& e : A.rows[i]) {
            if (e.col < n) {
                b[i] += e.value * x_true[e.col];
            }
        }
    }
    return b;
}

int main() {
    int n = 500;          // Размер матрицы (увеличен для наглядности)
    int max_iter = 1000;
    double eps = 1e-6;

    // Создаём задачу
    Sparse A = createTridiagonalMatrix(n);
    std::vector<double> b = createRHS(A);

    // Векторы ошибок
    std::vector<double> err_jacobi, err_gauss, err_mpi, err_cheb;

    // Начальные приближения (нули — далеко от решения)
    std::vector<double> x_jacobi(n, 0.0);
    std::vector<double> x_gauss(n, 0.0);
    std::vector<double> x_mpi(n, 0.0);
    std::vector<double> x_cheb(n, 0.0);

    // === Jacobi ===
    std::cout << "Running Jacobi..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    int iter_jacobi = jacobi(A, b, x_jacobi, max_iter, eps, err_jacobi);
    auto t1 = std::chrono::high_resolution_clock::now();
    double time_jacobi = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // === Gauss-Seidel ===
    std::cout << "Running Gauss-Seidel..." << std::endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    int iter_gauss = gaussSeidel(A, b, x_gauss, max_iter, eps, err_gauss);
    auto t3 = std::chrono::high_resolution_clock::now();
    double time_gauss = std::chrono::duration<double, std::milli>(t3 - t2).count();

    // Подготовка матрицы для MPI-классов
    SparseMatrix A_mpi;
    A_mpi.n = A.n;
    A_mpi.rows.resize(A.n);
    for (int i = 0; i < A.n; ++i) {
        for (const auto& e : A.rows[i]) {
            A_mpi.rows[i].push_back({e.col, e.value});
        }
    }

    // === Simple MPI ===
    std::cout << "Running Simple MPI..." << std::endl;
    ChebyshevMPI mpi(A_mpi, b);
    auto t4 = std::chrono::high_resolution_clock::now();
    // tau=0.2 — малый шаг для медленной сходимости (чтобы график был плавным)
    int iter_mpi = mpi.solveSimple(x_mpi, max_iter, eps, err_mpi, 0.2);
    auto t5 = std::chrono::high_resolution_clock::now();
    double time_mpi = std::chrono::duration<double, std::milli>(t5 - t4).count();

    // === Chebyshev MPI ===
    std::cout << "Running Chebyshev MPI..." << std::endl;
    ChebyshevMPI cheb(A_mpi, b);
    auto t6 = std::chrono::high_resolution_clock::now();
    int iter_cheb = cheb.solve(x_cheb, max_iter, eps, err_cheb);
    auto t7 = std::chrono::high_resolution_clock::now();
    double time_cheb = std::chrono::duration<double, std::milli>(t7 - t6).count();

    // === Вывод результатов ===
    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "Jacobi:        " << iter_jacobi << " iterations, "
              << time_jacobi << " ms" << std::endl;
    std::cout << "Gauss-Seidel:  " << iter_gauss << " iterations, "
              << time_gauss << " ms" << std::endl;
    std::cout << "Simple MPI:    " << iter_mpi << " iterations, "
              << time_mpi << " ms" << std::endl;
    std::cout << "Chebyshev MPI: " << iter_cheb << " iterations, "
              << time_cheb << " ms" << std::endl;

    std::ofstream fout("build/convergence.csv");
    fout << "Iter,MPI_error,Jacobi_error,GaussSeidel_error,Chebyshev_error" << std::endl;

    size_t max_iters = err_jacobi.size();
    max_iters = std::max(max_iters, err_gauss.size());
    max_iters = std::max(max_iters, err_mpi.size());
    max_iters = std::max(max_iters, err_cheb.size());

    for (size_t i = 0; i < max_iters; ++i) {
        fout << i + 1 << ",";
        fout << (i < err_mpi.size() ? std::to_string(err_mpi[i]) : "") << ",";
        fout << (i < err_jacobi.size() ? std::to_string(err_jacobi[i]) : "") << ",";
        fout << (i < err_gauss.size() ? std::to_string(err_gauss[i]) : "") << ",";
        fout << (i < err_cheb.size() ? std::to_string(err_cheb[i]) : "") << std::endl;
    }
    fout.close();

    // Данные по времени
    std::ofstream time_out("build/timing.csv");
    time_out << "Method,Iterations,Time_ms" << std::endl;
    time_out << "MPI," << iter_mpi << "," << time_mpi << std::endl;
    time_out << "Jacobi," << iter_jacobi << "," << time_jacobi << std::endl;
    time_out << "GaussSeidel," << iter_gauss << "," << time_gauss << std::endl;
    time_out << "Chebyshev," << iter_cheb << "," << time_cheb << std::endl;
    time_out.close();

    std::cout << "\nData saved to build/convergence.csv and build/timing.csv" << std::endl;
    
    return 0;
}