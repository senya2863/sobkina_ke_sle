#define _USE_MATH_DEFINES
#include "gmres.h"
#include <algorithm>
#include <cmath>

GMRES::GMRES(SparseMatrix& matrix, const std::vector<double>& rhs, int restart)
    : A(matrix), b(rhs), restart(restart) {}

std::vector<double> GMRES::multiply(const std::vector<double>& x) {
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

double GMRES::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

double GMRES::norm(const std::vector<double>& x) {
    return std::sqrt(dot(x, x));
}

std::vector<double> GMRES::subtract(const std::vector<double>& a, 
                                     const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

int GMRES::solve(std::vector<double>& x, int max_iter, double eps, 
                 std::vector<double>& errors) {
    int n = A.n;

    std::vector<double> r = subtract(b, multiply(x));
    double beta = norm(r);

    errors.push_back(beta);
    if (beta < eps) return 0;
    
    int total_iter = 0;
    while (total_iter < max_iter) {
        // Базис Крылова
        std::vector<std::vector<double>> V(restart + 1, std::vector<double>(n, 0.0));
        // Матрица Хессенберга
        std::vector<std::vector<double>> H(restart + 1, std::vector<double>(restart, 0.0));
        // Вращения Гивенса
        std::vector<double> cs(restart, 0.0);
        std::vector<double> sn(restart, 0.0);
        std::vector<double> g(restart + 1, 0.0);
        g[0] = beta;

        double norm_r = norm(r);
        if (norm_r < 1e-12) break;
        
        for (int i = 0; i < n; ++i) {
            V[0][i] = r[i] / norm_r;
        }

        for (int j = 0; j < restart && total_iter < max_iter; ++j) {
            total_iter++;

            std::vector<double> w = multiply(V[j]);

            for (int i = 0; i <= j; ++i) {
                H[i][j] = dot(w, V[i]);
                for (int k = 0; k < n; ++k) {
                    w[k] -= H[i][j] * V[i][k];
                }
            }
            H[j + 1][j] = norm(w);

            if (H[j + 1][j] > 1e-10) {
                for (int k = 0; k < n; ++k) {
                    V[j + 1][k] = w[k] / H[j + 1][j];
                }
            }

            for (int i = 0; i < j; ++i) {
                double temp = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                H[i][j] = temp;
            }

            double rho = std::sqrt(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]);
            if (rho < 1e-12) rho = 1e-12;
            cs[j] = H[j][j] / rho;
            sn[j] = H[j + 1][j] / rho;

            H[j][j] = cs[j] * H[j][j] + sn[j] * H[j + 1][j];
            H[j + 1][j] = 0.0;
            
            double temp_g = cs[j] * g[j] + sn[j] * g[j + 1];
            g[j + 1] = -sn[j] * g[j] + cs[j] * g[j + 1];
            g[j] = temp_g;
            
            // === ОШИБКА: реальная норма невязки ===
            std::vector<double> Ax = multiply(x);
            std::vector<double> residual(n);
            for (int i = 0; i < n; ++i) {
                residual[i] = b[i] - Ax[i];
            }
            double err = norm(residual);
            errors.push_back(err);
            
            if (err < eps) {
                std::vector<double> y(j + 1, 0.0);
                for (int i = j; i >= 0; --i) {
                    y[i] = g[i];
                    for (int k = i + 1; k <= j; ++k) {
                        y[i] -= H[i][k] * y[k];
                    }
                    if (std::abs(H[i][i]) > 1e-12) {
                        y[i] /= H[i][i];
                    }
                }

                for (int i = 0; i <= j; ++i) {
                    for (int k = 0; k < n; ++k) {
                        x[k] += y[i] * V[i][k];
                    }
                }
                
                return total_iter;
            }
        }

        std::vector<double> y(restart, 0.0);
        for (int i = restart - 1; i >= 0; --i) {
            y[i] = g[i];
            for (int k = i + 1; k < restart; ++k) {
                y[i] -= H[i][k] * y[k];
            }
            if (std::abs(H[i][i]) > 1e-12) {
                y[i] /= H[i][i];
            }
        }
        
        for (int i = 0; i < restart; ++i) {
            for (int k = 0; k < n; ++k) {
                x[k] += y[i] * V[i][k];
            }
        }

        r = subtract(b, multiply(x));
        beta = norm(r);
        
        if (beta < eps) {
            errors.push_back(beta);
            return total_iter;
        }
    }
    
    return total_iter;
}