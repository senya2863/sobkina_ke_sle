#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <algorithm>
#include "gmres.h"
#include "iterative_methods.h"
#include "elliptic_generator.h"

int main() {
    int nx = 20, ny = 20;
    int max_iter = 500;
    double eps = 1e-6;

    SparseMatrix A = EllipticGenerator::generate2D(nx, ny);
    std::vector<double> b = EllipticGenerator::generateRHS(nx, ny);

    int n = A.n;
    std::cout << "Matrix size: " << n << "x" << n << std::endl;

    std::vector<double> x_test(n, 0.0);
    std::vector<double> Ax_test(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (const auto& [col, val] : A.rows[i]) {
            if (col < n) Ax_test[i] += val * x_test[col];
        }
    }
    double init_residual = 0.0;
    for (int i = 0; i < n; ++i) {
        double r_i = b[i] - Ax_test[i];
        init_residual += r_i * r_i;
    }
    init_residual = std::sqrt(init_residual);
    std::cout << "DEBUG: Initial residual norm = " << init_residual << std::endl;
    std::cout << "DEBUG: eps = " << eps << std::endl;
    std::cout << "DEBUG: b[0] = " << b[0] << ", b[100] = " << b[100] << std::endl;
    std::cout << "DEBUG: b[200] = " << b[200] << ", b[300] = " << b[300] << std::endl;

    std::vector<double> err_gs, err_sor, err_cg, err_gmres;

    std::vector<double> x_gs(n, 0.0);
    std::vector<double> x_sor(n, 0.0);
    std::vector<double> x_cg(n, 0.0);
    std::vector<double> x_gmres(n, 0.0);

    IterativeMethods methods(A, b);
    GMRES gmres(A, b, 30);

    // Gauss-Seidel
    auto t0 = std::chrono::high_resolution_clock::now();
    int iter_gs = methods.gaussSeidel(x_gs, max_iter, eps, err_gs);
    auto t1 = std::chrono::high_resolution_clock::now();
    double time_gs = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // SOR
    auto t2 = std::chrono::high_resolution_clock::now();
    int iter_sor = methods.sor(x_sor, max_iter, eps, err_sor, 1.7);
    auto t3 = std::chrono::high_resolution_clock::now();
    double time_sor = std::chrono::duration<double, std::milli>(t3 - t2).count();

    // CG
    auto t4 = std::chrono::high_resolution_clock::now();
    int iter_cg = methods.conjugateGradient(x_cg, max_iter, eps, err_cg);
    auto t5 = std::chrono::high_resolution_clock::now();
    double time_cg = std::chrono::duration<double, std::milli>(t5 - t4).count();

    // GMRES
    auto t6 = std::chrono::high_resolution_clock::now();
    int iter_gmres = gmres.solve(x_gmres, max_iter, eps, err_gmres);
    auto t7 = std::chrono::high_resolution_clock::now();
    double time_gmres = std::chrono::duration<double, std::milli>(t7 - t6).count();

    std::cout << "Gauss-Seidel:  " << iter_gs << " iter, " << time_gs << " ms" << std::endl;
    std::cout << "SOR:   " << iter_sor << " iter, " << time_sor << " ms" << std::endl;
    std::cout << "CG:            " << iter_cg << " iter, " << time_cg << " ms" << std::endl;
    std::cout << "GMRES:  " << iter_gmres << " iter, " << time_gmres << " ms" << std::endl;

    //CSV
    std::ofstream fout("convergence_hw8.csv");
    fout << "Iter,GS_error,SOR_error,CG_error,GMRES_error" << std::endl;

    size_t max_iters = std::max({err_gs.size(), err_sor.size(),
                                  err_cg.size(), err_gmres.size()});

    for (size_t i = 0; i < max_iters; ++i) {
        fout << i + 1 << ",";
        fout << (i < err_gs.size() ? std::to_string(err_gs[i]) : "") << ",";
        fout << (i < err_sor.size() ? std::to_string(err_sor[i]) : "") << ",";
        fout << (i < err_cg.size() ? std::to_string(err_cg[i]) : "") << ",";
        fout << (i < err_gmres.size() ? std::to_string(err_gmres[i]) : "") << std::endl;
    }
    fout.close();

    std::ofstream time_out("timing_hw8.csv");
    time_out << "method,iterations,time_ms" << std::endl;
    time_out << "GS," << iter_gs << "," << time_gs << std::endl;
    time_out << "SOR," << iter_sor << "," << time_sor << std::endl;
    time_out << "CG," << iter_cg << "," << time_cg << std::endl;
    time_out << "GMRES," << iter_gmres << "," << time_gmres << std::endl;
    time_out.close();
    
    std::cout << "\nData saved to convergence_hw8.csv and timing_hw8.csv" << std::endl;
    
    return 0;
}