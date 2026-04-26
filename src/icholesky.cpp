#include "icholesky.h"
#include <cmath>
#include <algorithm>

SparseMatrix IncompleteCholesky::decompose(const SparseMatrix& A) {
    int n = A.n;
    SparseMatrix L;
    L.n = n;
    L.rows.resize(n);

    // Копируем структуру из A (нижняя треугольная часть)
    for (int i = 0; i < n; ++i) {
        for (const auto& [col, val] : A.rows[i]) {
            if (col <= i) {
                L.rows[i].push_back({col, 0.0});
            }
        }
    }

    // Заполняем L
    for (int i = 0; i < n; ++i) {
        // Диагональ
        double sum_diag = 0.0;
        for (const auto& [k, l_ik] : L.rows[i]) {
            if (k < i) {
                sum_diag += l_ik * l_ik;
            }
        }

        double a_ii = 0.0;
        for (const auto& [col, val] : A.rows[i]) {
            if (col == i) {
                a_ii = val;
                break;
            }
        }

        double l_ii = std::sqrt(std::max(a_ii - sum_diag, 1e-12));

        for (auto& [col, val] : L.rows[i]) {
            if (col == i) {
                val = l_ii;
                break;
            }
        }

        // Недиагональные элементы
        for (auto& [j, val] : L.rows[i]) {
            if (j < i) {
                double sum_off = 0.0;
                for (const auto& [k, l_ik] : L.rows[i]) {
                    if (k >= j) continue;
                    for (const auto& [kj, l_jk] : L.rows[j]) {
                        if (kj == k) {
                            sum_off += l_ik * l_jk;
                            break;
                        }
                    }
                }

                double a_ij = 0.0;
                for (const auto& [col, a_val] : A.rows[i]) {
                    if (col == j) {
                        a_ij = a_val;
                        break;
                    }
                }

                double l_jj = 0.0;
                for (const auto& [col, l_val] : L.rows[j]) {
                    if (col == j) {
                        l_jj = l_val;
                        break;
                    }
                }

                val = (a_ij - sum_off) / std::max(l_jj, 1e-12);
            }
        }
    }

    return L;
}

std::vector<double> IncompleteCholesky::solveLower(const SparseMatrix& L, const std::vector<double>& b) {
    int n = L.n;
    std::vector<double> y(n);
    
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (const auto& [col, val] : L.rows[i]) {
            if (col < i) {
                sum += val * y[col];
            }
        }
        
        double l_ii = 0.0;
        for (const auto& [col, val] : L.rows[i]) {
            if (col == i) {
                l_ii = val;
                break;
            }
        }
        
        y[i] = (b[i] - sum) / l_ii;
    }
    
    return y;
}

std::vector<double> IncompleteCholesky::solveUpper(const SparseMatrix& L, const std::vector<double>& y) {
    int n = L.n;
    std::vector<double> x(n);
    
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;

        for (int j = i + 1; j < n; ++j) {
            for (const auto& [col, val] : L.rows[j]) {
                if (col == i) {
                    sum += val * x[j];
                    break;
                }
            }
        }
        
        double l_ii = 0.0;
        for (const auto& [col, val] : L.rows[i]) {
            if (col == i) {
                l_ii = val;
                break;
            }
        }
        
        x[i] = (y[i] - sum) / l_ii;
    }
    
    return x;
}