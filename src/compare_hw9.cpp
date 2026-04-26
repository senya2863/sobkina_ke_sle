#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <algorithm>
#include "level_matrix.h"
#include "icholesky.h"
#include "iterative_methods.h"
#include "elliptic_generator.h"

std::vector<double> multiplySparse(const SparseMatrix& A, const std::vector<double>& x) {
    int n = A.n;
    std::vector<double> result(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (const auto& [col, val] : A.rows[i]) {
            if (col < n) {
                result[i] += val * x[col];
            }
        }
    }
    return result;
}

double norm(const std::vector<double>& v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

//CG без предобуславливателя
int conjugateGradientPlain(
    const SparseMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    int max_iter,
    double eps,
    std::vector<double>& errors
) {
    int n = A.n;

    std::vector<double> r(n);
    std::vector<double> Ax(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (const auto& [col, val] : A.rows[i]) {
            if (col < n) Ax[i] += val * x[col];
        }
    }
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - Ax[i];
    }

    std::vector<double> p = r;

    double rs_old = 0.0;
    for (int i = 0; i < n; ++i) {
        rs_old += r[i] * r[i];
    }

    double err = std::sqrt(rs_old);
    errors.push_back(err);
    if (err < eps) {
        return 0;
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> Ap(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (const auto& [col, val] : A.rows[i]) {
                if (col < n) Ap[i] += val * p[col];
            }
        }
        // pAp = (p, Ap)
        double pAp = 0.0;
        for (int i = 0; i < n; ++i) {
            pAp += p[i] * Ap[i];
        }
        if (std::abs(pAp) < 1e-15) {
            break;
        }
        double alpha = rs_old / pAp;
        // x = x + alpha * p
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
        }
        // r = r - alpha * Ap
        for (int i = 0; i < n; ++i) {
            r[i] -= alpha * Ap[i];
        }
        // rs_new = (r, r)
        double rs_new = 0.0;
        for (int i = 0; i < n; ++i) {
            rs_new += r[i] * r[i];
        }
        err = std::sqrt(rs_new);
        errors.push_back(err);

        if (err < eps) {
            return iter + 1;
        }

        double beta = rs_new / rs_old;
        rs_old = rs_new;
        // p = r + beta * p
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }
    }
    return max_iter;
}

// Метод сопряжённых градиентов с предобуславливателем
int conjugateGradientPrecond(
    const SparseMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    int max_iter,
    double eps,
    std::vector<double>& errors,
    const SparseMatrix& L  // предобуславливатель: A ≈ L × Lᵀ
) {
    int n = A.n;

    // r = b - Ax
    std::vector<double> r = b;
    std::vector<double> Ax = multiplySparse(A, x);
    for (int i = 0; i < n; ++i) {
        r[i] -= Ax[i];
    }
    std::vector<double> y = IncompleteCholesky::solveLower(L, r);
    std::vector<double> z = IncompleteCholesky::solveUpper(L, y);
    std::vector<double> p = z;
    double rs_old = 0.0;
    for (int i = 0; i < n; ++i) {
        rs_old += r[i] * z[i];
    }
    double err = norm(r);
    errors.push_back(err);

    if (err < eps) return 0;
    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> Ap = multiplySparse(A, p);
        double alpha = rs_old;
        double pAp = 0.0;
        for (int i = 0; i < n; ++i) {
            pAp += p[i] * Ap[i];
        }
        if (std::abs(pAp) < 1e-15) break;
        alpha /= pAp;

        // x = x + alpha × p
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
        }

        // r = r - alpha × Ap
        for (int i = 0; i < n; ++i) {
            r[i] -= alpha * Ap[i];
        }

        y = IncompleteCholesky::solveLower(L, r);
        z = IncompleteCholesky::solveUpper(L, y);

        double rs_new = 0.0;
        for (int i = 0; i < n; ++i) {
            rs_new += r[i] * z[i];
        }
        double beta = rs_new / rs_old;
        rs_old = rs_new;
        // p = z + beta × p
        for (int i = 0; i < n; ++i) {
            p[i] = z[i] + beta * p[i];
        }
        err = norm(r);
        errors.push_back(err);

        if (err < eps) return iter + 1;
    }

    return max_iter;
}

int main() {
    int nx = 20, ny = 20;
    int max_iter = 500;
    double eps = 1e-6;

    SparseMatrix A = EllipticGenerator::generate2D(nx, ny);
    std::vector<double> b = EllipticGenerator::generateRHS(nx, ny);

    int n = A.n;
    std::cout << "Matrix size: " << n << "x" << n << std::endl;

    // === Пункт 1: Матрица уровней ===
    std::vector<int> levels = computeLevelMatrix(A);
    std::ofstream level_file("levels_hw9.csv");
    level_file << "row,level" << std::endl;
    for (int i = 0; i < n; ++i) {
        level_file << i << "," << levels[i] << std::endl;
    }
    level_file.close();
    SparseMatrix L = IncompleteCholesky::decompose(A);

    std::vector<double> err_cg_plain, err_cg_prec;

    // CG без предобуславливателя
    std::vector<double> x_plain(n, 0.0);
    auto t0 = std::chrono::high_resolution_clock::now();
    int iter_plain = conjugateGradientPlain(A, b, x_plain, max_iter, eps, err_cg_plain);
    auto t1 = std::chrono::high_resolution_clock::now();
    double time_plain = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // CG С предобуславливателем
    std::vector<double> x_prec(n, 0.0);
    auto t2 = std::chrono::high_resolution_clock::now();
    int iter_prec = conjugateGradientPrecond(A, b, x_prec, max_iter, eps, err_cg_prec, L);
    auto t3 = std::chrono::high_resolution_clock::now();
    double time_prec = std::chrono::duration<double, std::milli>(t3 - t2).count();

    std::cout << "CG (plain):           " << iter_plain << " iter, " << time_plain << " ms" << std::endl;
    std::cout << "CG (IC(0) precond):   " << iter_prec << " iter, " << time_prec << " ms" << std::endl;

    std::ofstream fout("convergence_hw9.csv");
    fout << "Iter,CG_plain_error,CG_prec_error" << std::endl;

    size_t max_iters = std::max(err_cg_plain.size(), err_cg_prec.size());
    for (size_t i = 0; i < max_iters; ++i) {
        fout << i + 1 << ",";
        fout << (i < err_cg_plain.size() ? std::to_string(err_cg_plain[i]) : "") << ",";
        fout << (i < err_cg_prec.size() ? std::to_string(err_cg_prec[i]) : "") << std::endl;
    }
    fout.close();

    std::ofstream time_out("timing_hw9.csv");
    time_out << "method,iterations,time_ms" << std::endl;
    time_out << "CG_plain," << iter_plain << "," << time_plain << std::endl;
    time_out << "CG_prec," << iter_prec << "," << time_prec << std::endl;
    time_out.close();

    std::cout << "\nData saved to convergence_hw9.csv and timing_hw9.csv" << std::endl;

    return 0;
}