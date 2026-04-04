#include "symmetr_methods.h"
#include <algorithm>
#include <numeric>

symmetrizedmethods::symmetrizedmethods(sparsematrix& matrix, const std::vector<double>& rhs) : a(matrix), b(rhs) {
    estimatespectrum();
}

void symmetrizedmethods::estimatespectrum() { //через суммы по строкам
    double row_sum_max = 0.0;
    
    for (int i = 0; i < a.n; ++i) {
        double diag = 0.0, off_diag_sum = 0.0;
        
        for (const auto& [col, val] : a.rows[i]) {
            if (col == i) {
                diag = std::abs(val);
            } else {
                off_diag_sum += std::abs(val);
            }
        }
        
        row_sum_max = std::max(row_sum_max, diag + off_diag_sum);
    }
    
    lambda_min = 0.1 * row_sum_max;
    lambda_max = row_sum_max;
}

int symmetrizedmethods::solvesymmetricgs(std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors) {
    int n = a.n;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> x_old = x;
        
        //прямой проход (как обычный гаусс)
        for (int i = 0; i < n; ++i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& [col, val] : a.rows[i]) {
                if (col >= n) continue;
                if (col == i) diag = val;
                else sum += val * x[col];
            }
            x[i] = (b[i] - sum) / diag;
        }
        
        //обратный проход (от n-1 до 0)
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& [col, val] : a.rows[i]) {
                if (col >= n) continue;
                if (col == i) diag = val;
                else sum += val * x[col];
            }
            x[i] = (b[i] - sum) / diag;
        }
        
        //ошибки
        double max_err = 0.0;
        for (int i = 0; i < n; ++i) {
            max_err = std::max(max_err, std::abs(x[i] - x_old[i]));
        }
        errors.push_back(max_err);
        
        if (max_err < eps) {
            return (iter + 1) * 2; //без этого не работает
        }
    }
    
    return max_iter * 2;
}

int symmetrizedmethods::solvechebyshevsymmetrized(std::vector<double>& x, int max_iter, double eps, std::vector<double>& errors) {
    // чебышёвские параметры
    double sigma = (lambda_max - lambda_min) / (lambda_max + lambda_min);
    
    int n = a.n;
    std::vector<double> x_old(n);
    std::vector<double> r(n);

    double tau_0 = 2.0 / (lambda_max + lambda_min);
    
    // один шаг симметричного гаусс-зейделя для получения невязки
    for (int i = 0; i < n; ++i) {
        double sum = 0.0, diag = 0.0;
        for (const auto& [col, val] : a.rows[i]) {
            if (col >= n) continue;
            if (col == i) diag = val;
            else sum += val * x[col];
        }
        x_old[i] = (b[i] - sum) / diag;
    }
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0, diag = 0.0;
        for (const auto& [col, val] : a.rows[i]) {
            if (col >= n) continue;
            if (col == i) diag = val;
            else sum += val * x_old[col];
        }
        x_old[i] = (b[i] - sum) / diag;
    }
    
    //невязка r = b - ax
    for (int i = 0; i < n; ++i) {
        double ax_i = 0.0;
        for (const auto& [col, val] : a.rows[i]) {
            if (col < n) ax_i += val * x_old[col];
        }
        r[i] = b[i] - ax_i;
    }
    
    // x = x + tau_0 * r
    double max_err = 0.0;
    for (int i = 0; i < n; ++i) {
        x[i] = x_old[i] + tau_0 * r[i];
        max_err = std::max(max_err, std::abs(x[i] - x_old[i]));
    }
    errors.push_back(max_err);
    
    if (max_err < eps) return 1;
    
    // чебыевские итерации
    for (int k = 1; k < max_iter; ++k) {
        x_old = x;
        
        // чебышевский параметр
        double theta_k = 2.0 / (lambda_max + lambda_min + (lambda_max - lambda_min) * std::cos(M_PI * (2*k + 1) / (2*k + 2)));
        
        // шаг симметричного гаусс-зейделя
        for (int i = 0; i < n; ++i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& [col, val] : a.rows[i]) {
                if (col >= n) continue;
                if (col == i) diag = val;
                else sum += val * x[col];
            }
            x[i] = (b[i] - sum) / diag;
        }
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0.0, diag = 0.0;
            for (const auto& [col, val] : a.rows[i]) {
                if (col >= n) continue;
                if (col == i) diag = val;
                else sum += val * x[col];
            }
            x[i] = (b[i] - sum) / diag;
        }
        
        // невязка
        for (int i = 0; i < n; ++i) {
            double ax_i = 0.0;
            for (const auto& [col, val] : a.rows[i]) {
                if (col < n) ax_i += val * x[col];
            }
            r[i] = b[i] - ax_i;
        }
        
        // чебышёвская коррекция
        max_err = 0.0;
        for (int i = 0; i < n; ++i) {
            double x_new = x[i] + theta_k * r[i];
            max_err = std::max(max_err, std::abs(x_new - x[i]));
            x[i] = x_new;
        }
        errors.push_back(max_err);
        
        if (max_err < eps) {
            return k + 1;
        }
    }
    
    return max_iter;
}