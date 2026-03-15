#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <algorithm>

struct el
{
    int col;
    double value;
};

struct sparse
{
    int n;
    std::vector<std::vector<el>> rows;
};

double norm(const std::vector<double>& a, const std::vector<double>& b)
{
    double max = 0.0;
    for (size_t i = 0; i < a.size(); i++)
    {
        double diff = std::abs(a[i] - b[i]);
        if (diff > max)
            max = diff;
    }
    return max;
}

// Простая итерация: x_new = x + tau*(b - A*x)
int simple_iter(const sparse& A, const std::vector<double>& b, std::vector<double>& x,
                int max_iter, double eps, double tau, std::vector<double>& errors)
{
    int n = A.n;
    std::vector<double> x_new(n, 0.0);

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            double Ax = 0.0;
            for (const el& e : A.rows[i])
                if (e.col >= 0 && e.col < n)
                    Ax += e.value * x[e.col];
            x_new[i] = x[i] + tau * (b[i] - Ax);
        }
        double err = norm(x_new, x);
        errors.push_back(err);
        x = x_new;
        if (err < eps) return iter + 1;
    }
    return max_iter;
}

int jacobi(const sparse& A, const std::vector<double>& b, std::vector<double>& x,
           int max_iter, double eps, std::vector<double>& errors)
{
    int n = A.n;
    std::vector<double> x_new(n);

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0, diag = 0.0;
            for (const el& e : A.rows[i]) {
                if (e.col >= n) continue;
                if (e.col == i) diag = e.value;
                else            sum += e.value * x[e.col];
            }
            x_new[i] = (b[i] - sum) / diag;
        }
        double err = norm(x_new, x);  // считаем один раз
        errors.push_back(err);
        x = x_new;
        if (err < eps) return iter + 1;
    }
    return max_iter;
}

// Метод Гаусса-Зейделя
int gauss_seidel(const sparse& A, const std::vector<double>& b, std::vector<double>& x,
                 int max_iter, double eps, std::vector<double>& errors)
{
    int n = A.n;

    for (int iter = 0; iter < max_iter; iter++) {
        std::vector<double> old = x;
        for (int i = 0; i < n; i++) {
            double sum = 0.0, diag = 0.0;
            for (const el& e : A.rows[i]) {
                if (e.col >= n) continue;
                if (e.col == i) diag = e.value;
                else            sum += e.value * x[e.col];
            }
            x[i] = (b[i] - sum) / diag;
        }
        double err = norm(x, old);
        errors.push_back(err);
        if (err < eps) return iter + 1;
    }
    return max_iter;
}

int main()
{
    sparse A;
    A.n = 3;
    A.rows.resize(3);
    A.rows[0] = {{0, 4},  {1, -1}};
    A.rows[1] = {{0, -1}, {1, 4},  {2, -1}};
    A.rows[2] = {{1, -1}, {2, 3}};

    std::vector<double> b = {15, 10, 10};
    int    max_iter = 100;
    double eps      = 1e-6;
    double tau      = 0.1;

    std::vector<double> x_jacobi(A.n, 0.0);
    std::vector<double> x_gauss (A.n, 0.0);
    std::vector<double> x_iter  (A.n, 0.0);
    std::vector<double> err_jacobi, err_gauss, err_iter;

    auto t0 = std::chrono::high_resolution_clock::now();
    int n_jacobi = jacobi     (A, b, x_jacobi, max_iter, eps, err_jacobi);
    auto t1 = std::chrono::high_resolution_clock::now();
    int n_gauss  = gauss_seidel(A, b, x_gauss,  max_iter, eps, err_gauss);
    auto t2 = std::chrono::high_resolution_clock::now();
    int n_iter   = simple_iter (A, b, x_iter,   max_iter, eps, tau, err_iter);
    auto t3 = std::chrono::high_resolution_clock::now();

    std::cout << "Jacobi        iterations: " << n_jacobi
              << ", time: " << std::chrono::duration<double>(t1-t0).count() << " s\n";
    std::cout << "Gauss-Seidel  iterations: " << n_gauss
              << ", time: " << std::chrono::duration<double>(t2-t1).count() << " s\n";
    std::cout << "Simple iter   iterations: " << n_iter
              << ", time: " << std::chrono::duration<double>(t3-t2).count() << " s\n";


    std::ofstream fout("convergence.csv");
    fout << "Iter,Jacobi,GaussSeidel,SimpleIter\n";
    size_t rows = std::max({err_jacobi.size(), err_gauss.size(), err_iter.size()});
    for (size_t i = 0; i < rows; i++) {
        fout << i + 1 << ",";
        fout << (i < err_jacobi.size() ? std::to_string(err_jacobi[i]) : "") << ",";
        fout << (i < err_gauss.size()  ? std::to_string(err_gauss[i])  : "") << ",";
        fout << (i < err_iter.size()   ? std::to_string(err_iter[i])   : "") << "\n";
    }
    fout.close();

    return 0;
}