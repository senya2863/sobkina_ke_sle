#include "elliptic_generator.h"
#include <cmath>

SparseMatrix EllipticGenerator::generate2D(int nx, int ny) {
    int n = nx * ny;
    SparseMatrix A;
    A.n = n;
    A.rows.resize(n);
    
    double h = 1.0 / (nx + 1);
    double coeff = 1.0 / (h * h);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = j * nx + i;

            A.rows[idx].push_back({idx, 4.0 * coeff});

            if (i > 0) {
                A.rows[idx].push_back({idx - 1, -coeff});
            }

            if (i < nx - 1) {
                A.rows[idx].push_back({idx + 1, -coeff});
            }

            if (j > 0) {
                A.rows[idx].push_back({idx - nx, -coeff});
            }

            if (j < ny - 1) {
                A.rows[idx].push_back({idx + nx, -coeff});
            }
        }
    }
    
    return A;
}

std::vector<double> EllipticGenerator::generateRHS(int nx, int ny) {
    int n = nx * ny;
    std::vector<double> b(n, 0.0);
    
    double h = 1.0 / (nx + 1);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = j * nx + i;
            double x = (i + 1) * h;
            double y = (j + 1) * h;

            b[idx] = 2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
        }
    }
    
    return b;
}