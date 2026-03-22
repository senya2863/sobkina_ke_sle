#include "chebyshev_mpi.h"
#include <algorithm>
#include <numeric>

ChebyshevMPI::ChebyshevMPI(SparseMatrix& matrix, const std::vector<double>& rhs)
    : A(matrix), b(rhs) {
    estimateSpectrum();
}

void ChebyshevMPI::estimateSpectrum() {
  //смотрим через диагонали для разреж или трехдиагональной матрицы
    double row_sum_max = 0.0;
    
    for (int i = 0; i < A.n; ++i) {
        double diag = 0.0, off_diag_sum = 0.0;
        
        for (const auto& [col, val] : A.rows[i]) {
            if (col == i) {
                diag = std::abs(val);
            } else {
                off_diag_sum += std::abs(val);
            }
        }

        double radius = off_diag_sum;
        row_sum_max = std::max(row_sum_max, diag + radius);
    }
    
    //приближ границы
    lambda_min = 0.1 * row_sum_max;
    lambda_max = row_sum_max;
}

int ChebyshevMPI::solveSimple(std::vector<double>& x, int max_iter, double eps,
                               std::vector<double>& errors, double tau) {
    std::vector<double> x_new(A.n);
    std::vector<double> Ax(A.n);
    
    for (int iter = 0; iter < max_iter; ++iter) {
        // Вычисляем Ax
        for (int i = 0; i < A.n; ++i) {
            Ax[i] = 0.0;
            for (const auto& [col, val] : A.rows[i]) {
                if (col < A.n) {
                    Ax[i] += val * x[col];
                }
            }
        }
        
        // x_new = x + tau * (b - Ax)
        double max_err = 0.0;
        for (int i = 0; i < A.n; ++i) {
            x_new[i] = x[i] + tau * (b[i] - Ax[i]);
            max_err = std::max(max_err, std::abs(x_new[i] - x[i]));
        }
        
        errors.push_back(max_err);
        x = x_new;
        
        if (max_err < eps) {
            return iter + 1;
        }
    }
    
    return max_iter;
}

int ChebyshevMPI::solve(std::vector<double>& x, int max_iter, double eps,
                        std::vector<double>& errors) {
    // Чебышёвские параметры
    double sigma = (lambda_max - lambda_min) / (lambda_max + lambda_min);
    double rho = (1.0 - std::sqrt(1.0 - sigma * sigma)) / sigma;
    
    std::vector<double> x_new(A.n);
    std::vector<double> Ax(A.n);
    std::vector<double> r(A.n); //невязка
    
    //обычный мпи начало
    double tau_0 = 2.0 / (lambda_max + lambda_min);
    
    for (int i = 0; i < A.n; ++i) {
        Ax[i] = 0.0;
        for (const auto& [col, val] : A.rows[i]) {
            if (col < A.n) {
                Ax[i] += val * x[col];
            }
        }
        r[i] = b[i] - Ax[i];
        x_new[i] = x[i] + tau_0 * r[i];
    }
    
    double max_err = 0.0;
    for (int i = 0; i < A.n; ++i) {
        max_err = std::max(max_err, std::abs(x_new[i] - x[i]));
    }
    errors.push_back(max_err);
    x = x_new;
    
    if (max_err < eps) return 1;

     //iterations
    for (int k = 1; k < max_iter; ++k) {
        double theta_k = 2.0 / (lambda_max + lambda_min + 
                     (lambda_max - lambda_min) * std::cos(M_PI * (2*k + 1) / (2*k + 2)));
        
        //Ax
        for (int i = 0; i < A.n; ++i) {
            Ax[i] = 0.0;
            for (const auto& [col, val] : A.rows[i]) {
                if (col < A.n) {
                    Ax[i] += val * x[col];
                }
            }
            r[i] = b[i] - Ax[i];
        }
        
        // x_new = x + theta_k * r
        max_err = 0.0;
        for (int i = 0; i < A.n; ++i) {
            x_new[i] = x[i] + theta_k * r[i];
            max_err = std::max(max_err, std::abs(x_new[i] - x[i]));
        }
        
        errors.push_back(max_err);
        x = x_new;
        
        if (max_err < eps) {
            return k + 1;
        }
    }
    
    return max_iter;
}