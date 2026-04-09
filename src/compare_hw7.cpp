#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <algorithm>
#include "iterative_methods.h"
#include "elliptic_generator.h"

int main() {
    int nx = 20, ny = 20;
    int max_iter = 1000;
    double eps = 1e-8;

    SparseMatrix A = EllipticGenerator::generate2D(nx, ny);
    std::vector<double> b = EllipticGenerator::generateRHS(nx, ny);

    int n = A.n;
    std::cout << "Matrix size: " << n << "x" << n << std::endl;

    std::vector<double> err_gs, err_sor, err_sor_cheb, err_sd, err_cg;

    std::vector<double> x_gs(n, 0.0);
    std::vector<double> x_sor(n, 0.0);
    std::vector<double> x_sor_cheb(n, 0.0);
    std::vector<double> x_sd(n, 0.0);
    std::vector<double> x_cg(n, 0.0);

    IterativeMethods methods(A, b);

    std::cout << "Running Gauss-Seidel..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    int iter_gs = methods.gaussSeidel(x_gs, max_iter, eps, err_gs);
    auto t1 = std::chrono::high_resolution_clock::now();
    double time_gs = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "Running SOR..." << std::endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    int iter_sor = methods.sor(x_sor, max_iter, eps, err_sor, 1.5);
    auto t3 = std::chrono::high_resolution_clock::now();
    double time_sor = std::chrono::duration<double, std::milli>(t3 - t2).count();

    std::cout << "Running SOR + Chebyshev..." << std::endl;
    auto t4 = std::chrono::high_resolution_clock::now();
    int iter_sor_cheb = methods.sor(x_sor_cheb, max_iter, eps, err_sor_cheb, 1.8);
    auto t5 = std::chrono::high_resolution_clock::now();
    double time_sor_cheb = std::chrono::duration<double, std::milli>(t5 - t4).count();

    std::cout << "Running Steepest Descent..." << std::endl;
    auto t6 = std::chrono::high_resolution_clock::now();
    int iter_sd = methods.steepestDescent(x_sd, max_iter, eps, err_sd);
    auto t7 = std::chrono::high_resolution_clock::now();
    double time_sd = std::chrono::duration<double, std::milli>(t7 - t6).count();

    std::cout << "Running Conjugate Gradient..." << std::endl;
    auto t8 = std::chrono::high_resolution_clock::now();
    int iter_cg = methods.conjugateGradient(x_cg, max_iter, eps, err_cg);
    auto t9 = std::chrono::high_resolution_clock::now();
    double time_cg = std::chrono::duration<double, std::milli>(t9 - t8).count();

    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "Gauss-Seidel:     " << iter_gs << " iter, " << time_gs << " ms" << std::endl;
    std::cout << "SOR (w=1.5):      " << iter_sor << " iter, " << time_sor << " ms" << std::endl;
    std::cout << "SOR+Cheb (w=1.8): " << iter_sor_cheb << " iter, " << time_sor_cheb << " ms" << std::endl;
    std::cout << "Steepest Descent: " << iter_sd << " iter, " << time_sd << " ms" << std::endl;
    std::cout << "CG:               " << iter_cg << " iter, " << time_cg << " ms" << std::endl;

    std::ofstream fout("convergence_hw7.csv");
    fout << "Iter,GS_error,SOR_error,SORCheb_error,SD_error,CG_error" << std::endl;

    size_t max_iters = std::max({err_gs.size(), err_sor.size(), err_sor_cheb.size(),
                                  err_sd.size(), err_cg.size()});

    for (size_t i = 0; i < max_iters; ++i) {
        fout << i + 1 << ",";
        fout << (i < err_gs.size() ? std::to_string(err_gs[i]) : "") << ",";
        fout << (i < err_sor.size() ? std::to_string(err_sor[i]) : "") << ",";
        fout << (i < err_sor_cheb.size() ? std::to_string(err_sor_cheb[i]) : "") << ",";
        fout << (i < err_sd.size() ? std::to_string(err_sd[i]) : "") << ",";
        fout << (i < err_cg.size() ? std::to_string(err_cg[i]) : "") << std::endl;
    }
    fout.close();

    std::ofstream time_out("timing_hw7.csv");
    time_out << "method,iterations,time_ms" << std::endl;
    time_out << "GS," << iter_gs << "," << time_gs << std::endl;
    time_out << "SOR," << iter_sor << "," << time_sor << std::endl;
    time_out << "SORCheb," << iter_sor_cheb << "," << time_sor_cheb << std::endl;
    time_out << "SD," << iter_sd << "," << time_sd << std::endl;
    time_out << "CG," << iter_cg << "," << time_cg << std::endl;
    time_out.close();

    std::cout << "\nData saved to convergence_hw7.csv and timing_hw7.csv" << std::endl;
    
    return 0;
}