#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

//sparse
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
        {
            max = diff;
        }
    }
    return max;
}

//easy method
int mpi(const sparse &A, const std::vector<double> &b, std::vector<double> &x, int max_iter, double eps, double t) {
    int n= A.n;
    std:: vector<double> x_new(n);

    for (int iter =0; iter<max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            double Ax=0.0;
            for (const el &e: A.rows[i]) {
                Ax += e.value * x[e.col];

            }
            x_new[i]=x[i]+t*(b[i] - Ax);
        }
        if (norm(x_new, x) < eps)
            return iter + 1;
        x = x_new;
    }
    return max_iter;
}

//jacobi
int jacobi (const sparse&A, const std::vector<double>& b, std::vector<double>& x,  int max_iter, double eps) {
    int n= A.n;

    std:: vector<double> x_new(n);

    for (int iter =0; iter<max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            double diag = 0.0;
            for (const el &e: A.rows[i]) {
                if (e.col == i) {
                    diag = e.value;
                }
                else {
                    sum += e.value*x[e.col];
                }
            }
            x_new[i]=(b[i]-sum)/diag;
        }
        if (norm(x_new,x)<eps) {
            return iter+1;
        }
        x = x_new;
    }
    return max_iter;
}

//gauss
int gauss (const sparse&A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double eps) {
    int n= A.n;

    for (int iter =0; iter<max_iter; iter++) {
        std:: vector<double> old=x;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            double diag = 0.0;
            for (const el &e: A.rows[i]) {
                if (e.col == i) {
                    diag = e.value;
                }
                else {
                    sum += e.value*x[e.col];
                }
            }
            x[i]=(b[i]-sum)/diag;
        }
        if (norm(x,old)<eps) {
            return iter+1;
        }
    }
    return max_iter;
}

int main()
{
    sparse A;
    A.n = 3;
    A.rows.resize(3);

    A.rows[0] = {{0, 4}, {1, -1}};
    A.rows[1] = {{0, -1}, {1, 4}, {2, -1}};
    A.rows[2] = {{1, -1}, {2, 3}};

    std::vector<double> b = {15, 10, 10};

    double eps = 1e-6;
    int max_iter = 10000;

    std::vector<double> x_jacobi(3, 0.0);
    std::vector<double> x_gauss(3, 0.0);
    std::vector<double> x_mpi(3, 0.0);

    auto start1 = std::chrono::high_resolution_clock::now();
    int iter_jacobi = jacobi(A, b, x_jacobi, max_iter, eps);
    auto end1 = std::chrono::high_resolution_clock::now();

    auto start2 = std::chrono::high_resolution_clock::now();
    int iter_gauss = gauss(A, b, x_gauss, max_iter, eps);
    auto end2 = std::chrono::high_resolution_clock::now();

    auto start3 = std::chrono::high_resolution_clock::now();
    int iter_mpi = mpi(A, b, x_mpi, max_iter, eps, 0.1);
    auto end3 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> t_jacobi = end1 - start1;
    std::chrono::duration<double> t_gauss = end2 - start2;
    std::chrono::duration<double> t_mpi = end3 - start3;

    std::cout << "jacobi iterations: " << iter_jacobi << "\n";
    std::cout << "jacobi time: " << t_jacobi.count() << " s\n\n";

    std::cout << "gauss iterations: " << iter_gauss << "\n";
    std::cout << "gauss time: " << t_gauss.count() << " s\n\n";

    std::cout << "mpi iterations: " << iter_mpi << "\n";
    std::cout << "mpi time: " << t_mpi.count() << " s\n";

    return 0;
}
// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.