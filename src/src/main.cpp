#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include "dense.h"
#include "csr.h"
#include "vector_ops.h"

double rnd() {
    static std::mt19937 gen(42);
    static std::uniform_real_distribution<> dis(0.0,1.0);
    return dis(gen);
}

densematrix generatedense(size_t n,double sparsity) {
    densematrix A(n,n);

    for(size_t i=0;i<n;i++)
        for(size_t j=0;j<n;j++)
            if(rnd()<sparsity)
                A.at(i,j)=rnd();

    return A;
}

csr generatecsr(size_t n,double sparsity) {
    csr A(n,n);

    for(size_t i=0;i<n;i++)
        for(size_t j=0;j<n;j++)
            if(rnd()<sparsity)
                A.addvalue(i,j,rnd());

    A.final();
    return A;
}

std::vector<double> randomvector(size_t n) {
    std::vector<double> v(n);

    for(size_t i=0;i<n;i++)
        v[i]=rnd();

    return v;
}

int main() {

    std::vector<size_t> sizes = {100,500,1000,2000};
    double sparsity = 0.01;

    std::cout<<"size dense csr\n";

    for(size_t n : sizes) {

        densematrix A = generatedense(n,sparsity);
        csr B = generatecsr(n,sparsity);
        auto x = randomvector(n);

        auto start = std::chrono::high_resolution_clock::now();
        auto y1 = A.multiply(x);
        auto end = std::chrono::high_resolution_clock::now();

        double dtime =std::chrono::duration<double>(end-start).count();

        start = std::chrono::high_resolution_clock::now();
        auto y2 = B.multiply(x);
        end = std::chrono::high_resolution_clock::now();

        double ctime = std::chrono::duration<double>(end-start).count();

        std::cout<<n<<" "<<dtime<<" "<<ctime<<"\n";
    }
}