#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>

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
int mpi(const sparse&A, const std::vector<double>& b, std::vector<double>& x,int max_iter, double eps, double tau, std::vector<double>& errors)
{
    int n = A.n;
    std::vector<double> x_new(n, 0.0);

    for (int iter=0; iter<max_iter; iter++) {
        for (int i=0; i<n; i++) {
            double Ax=0.0;
            for (const el &e : A.rows[i]) {
                if (e.col >= 0 && e.col < n)
                    Ax += e.value * x[e.col];
            }
            x_new[i] = x[i] + tau*(b[i]-Ax);
        }
        double err = norm(x_new, x);
        errors.push_back(err);

        if (err < eps) return iter+1;
        x = x_new;
    }
    return max_iter;
}

//jacobi
int jacobi (const sparse&A, const std::vector<double>& b, std::vector<double>& x,  int max_iter, double eps, std::vector<double> &errors) {
    int n= A.n;

    std:: vector<double> x_new(n);

    for (int iter =0; iter<max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            double diag = 0.0;
            for (const el &e: A.rows[i]) {
                if (e.col >= n) continue;
                if (e.col == i) {
                    diag = e.value;
                }
                else {
                    sum += e.value*x[e.col];
                }
            }
            x_new[i]=(b[i]-sum)/diag;
        }
        double err = norm(x_new, x);
        errors.push_back(err);

        if (norm(x_new,x)<eps) {
            return iter+1;
        }
        x = x_new;
    }
    return max_iter;
}

//gauss
int gauss (const sparse&A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double eps, std::vector<double> &errors) {
    int n= A.n;

    for (int iter =0; iter<max_iter; iter++) {
        std:: vector<double> old=x;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            double diag = 0.0;
            for (const el &e: A.rows[i]) {
                if (e.col >= n) continue;
                if (e.col == i) {
                    diag = e.value;
                }
                else {
                    sum += e.value*x[e.col];
                }
            }
            x[i]=(b[i]-sum)/diag;
        }
        double err = norm(x, old);
        errors.push_back(err);

        if (norm(x,old)<eps) {
            return iter+1;
        }
    }
    return max_iter;
}

int main()
{
    sparse A;
    A.n=3;
    A.rows.resize(3);
    A.rows[0] = {{0,4},{1,-1}};
    A.rows[1] = {{0,-1},{1,4},{2,-1}};
    A.rows[2] = {{1,-1},{2,3}};

    std::vector<double> b = {15,10,10};
    int max_iter=100;
    double eps=1e-6;
    double tau=0.1;

    std::vector<double> x_jacobi(A.n,0.0);
    std::vector<double> x_gauss(A.n,0.0);
    std::vector<double> x_mpi(A.n,0.0);
    std::vector<double> err_jacobi, err_gauss, err_mpi;

    auto start = std::chrono::high_resolution_clock::now();
    int iter_jacobi = jacobi(A,b,x_jacobi,max_iter,eps,err_jacobi);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_jacobi = end-start;

    start = std::chrono::high_resolution_clock::now();
    int iter_gauss = gauss(A,b,x_gauss,max_iter,eps,err_gauss);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_gauss = end-start;

    start = std::chrono::high_resolution_clock::now();
    int iter_mpi = mpi(A,b,x_mpi,max_iter,eps,tau,err_mpi);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_mpi = end-start;

    std::cout << "Jacobi iterations: " << iter_jacobi << ", time: " << t_jacobi.count() << " s\n";
    std::cout << "Gauss-Seidel iterations: " << iter_gauss << ", time: " << t_gauss.count() << " s\n";
    std::cout << "MPI iterations: " << iter_mpi << ", time: " << t_mpi.count() << " s\n";

    std::ofstream fout("convergence.csv");
    fout << "Iter,Jacobi,GaussSeidel,MPI\n";
    for (size_t i=0; i<err_jacobi.size(); i++)
        fout << i+1 << "," << err_jacobi[i] << "," << err_gauss[i] << "," << err_mpi[i] << "\n";
    fout.close();

    return 0;
}

// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.