#define _USE_MATH_DEFINES
#include "iterative_methods.h"
#include <algorithm>

IterativeMethods::IterativeMethods(SparseMatrix& matrix, const std::vector<double>& rhs)
    : A(matrix), b(rhs) {}

std::vector<double> IterativeMethods::multiply(const std::vector<double>& x) {
    std::vector<double> result(A.n, 0.0);
    for (int i = 0; i < A.n; ++i) {
        for (const auto& [col, val] : A.rows[i]) {
            if (col < A.n) {
                result[i] += val * x[col];
            }
        }
    }
    return result;
}

double IterativeMethods::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

std::vector<double> IterativeMethods::subtract(const std::vector<double>& a,
                                                const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

int IterativeMethods::gaussSeidel(std::vector<double>& x, int max_iter, double eps,
                                   std::vector<double>& errors) {
    int n = A.n;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> old = x;

        for (int i = 0; i < n; ++i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& [col, val] : A.rows[i]) {
                if (col >= n) continue;
                if (col == i) diag = val;
                else sum += val * x[col];
            }
            x[i] = (b[i] - sum) / diag;
        }

        double err = 0.0;
        for (int i = 0; i < n; ++i) {
            err = std::max(err, std::abs(x[i] - old[i]));
        }
        errors.push_back(err);

        if (err < eps) return iter + 1;
    }
    return max_iter;
}

int IterativeMethods::sor(std::vector<double>& x, int max_iter, double eps,
                          std::vector<double>& errors, double omega) {
    int n = A.n;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> old = x;

        for (int i = 0; i < n; ++i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& [col, val] : A.rows[i]) {
                if (col >= n) continue;
                if (col == i) diag = val;
                else sum += val * x[col];
            }
            double x_gs = (b[i] - sum) / diag;
            x[i] = (1.0 - omega) * old[i] + omega * x_gs;
        }

        double err = 0.0;
        for (int i = 0; i < n; ++i) {
            err = std::max(err, std::abs(x[i] - old[i]));
        }
        errors.push_back(err);

        if (err < eps) return iter + 1;
    }
    return max_iter;
}

int IterativeMethods::steepestDescent(std::vector<double>& x, int max_iter, double eps,
                                       std::vector<double>& errors) {
    int n = x.size();

    // r = b - Ax (начальная невязка)
    std::vector<double> r = subtract(b, multiply(x));

    // Записываем начальную ошибку
    double initial_err = 0.0;
    for (int i = 0; i < n; ++i) {
        initial_err += r[i] * r[i];
    }
    initial_err = std::sqrt(initial_err);
    errors.push_back(initial_err);

    if (initial_err < eps) return 0;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> old_x = x;
        std::vector<double> Ar = multiply(r);

        // alpha = (r, r) / (r, Ar)
        double r_dot_r = dot(r, r);
        double alpha = r_dot_r / dot(r, Ar);

        // x = x + alpha * r
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * r[i];
        }

        // Новая невязка r = b - Ax
        r = subtract(b, multiply(x));

        // === ОШИБКА = норма невязки ===
        double err = 0.0;
        for (int i = 0; i < n; ++i) {
            err += r[i] * r[i];
        }
        err = std::sqrt(err);
        errors.push_back(err);

        if (err < eps) return iter + 1;
    }
    return max_iter;
}

int IterativeMethods::conjugateGradient(std::vector<double>& x, int max_iter, double eps,
                                         std::vector<double>& errors) {
    int n = x.size();

    // r = b - Ax
    std::vector<double> r = subtract(b, multiply(x));
    std::vector<double> p = r;
    double rs_old = dot(r, r);

    // Записываем начальную ошибку = норма невязки
    double err = std::sqrt(rs_old);
    errors.push_back(err);

    // Если уже сошлись
    if (err < eps) return 0;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> Ap = multiply(p);
        double pAp = dot(p, Ap);

        if (std::abs(pAp) < 1e-15) {
            break;  // защита от деления на 0
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

        double rs_new = dot(r, r);

        // === ОШИБКА = норма невязки ===
        err = std::sqrt(rs_new);
        errors.push_back(err);

        if (err < eps) return iter + 1;

        double beta = rs_new / rs_old;
        rs_old = rs_new;

        // p = r + beta * p
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }
    }
    return max_iter;
}