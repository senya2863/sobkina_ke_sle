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

std::vector<double> IterativeMethods::subtract(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

int IterativeMethods::gaussSeidel(std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors) {
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

int IterativeMethods::sor(std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors, double omega) {
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

int IterativeMethods::steepestDescent(std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors) {
    // r = b - Ax
    std::vector<double> r = subtract(b, multiply(x));
    int n = x.size();

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> old_x = x;
        std::vector<double> Ar = multiply(r);

        // alpha = (r, r) / (r, Ar)
        double alpha = dot(r, r) / dot(r, Ar);

        // x = x + alpha * r
        for (int i = 0; i < n; ++i) {
            x[i] = x[i] + alpha * r[i];
        }

        r = subtract(b, multiply(x));

        double err = 0.0;
        for (int i = 0; i < n; ++i) {
            err = std::max(err, std::abs(x[i] - old_x[i]));
        }
        errors.push_back(err);

        if (err < eps) return iter + 1;
    }
    return max_iter;
}

int IterativeMethods::conjugateGradient(std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors) {
    // r = b - Ax
    std::vector<double> r = subtract(b, multiply(x));
    std::vector<double> p = r;  // d₀ = r₀
    double rs_old = dot(r, r);
    int n = x.size();

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> old_x = x;
        std::vector<double> Ap = multiply(p);

        // alpha = (r, r) / (p, Ap)
        double alpha = rs_old / dot(p, Ap);

        // x = x + alpha * p
        for (int i = 0; i < n; ++i) {
            x[i] = x[i] + alpha * p[i];
        }

        // r = r - alpha * Ap
        for (int i = 0; i < n; ++i) {
            r[i] = r[i] - alpha * Ap[i];
        }

        double rs_new = dot(r, r);

        // beta = (r_new, r_new) / (r_old, r_old)
        double beta = rs_new / rs_old;
        rs_old = rs_new;

        // p = r + beta * p
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }

        double err = std::sqrt(rs_new);
        errors.push_back(err);

        if (err < eps) return iter + 1;
    }
    return max_iter;
}
