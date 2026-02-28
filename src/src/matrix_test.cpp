#include "dense.h"
#include "csr.h"
#include <iostream>
#include <vector>

int main() {
  //dense
    densematrix D(2, 2);
    D.at(0, 0) = 1.0;
    D.at(0, 1) = 2.0;
    D.at(1, 0) = 3.0;
    D.at(1, 1) = 4.0;

    std::vector<double> x = {1.0, 1.0};
    auto y_dense = D.multiply(x);
    if (y_dense[0] == 3.0 && y_dense[1] == 7.0)
        std::cout << "dense good\n";
    else
        std::cout << "dense error\n";
//csr
    csr A(2, 2);
    A.addvalue(0, 0, 1.0);
    A.addvalue(0, 1, 2.0);
    A.addvalue(1, 0, 3.0);
    A.addvalue(1, 1, 4.0);
    A.final();
    auto y_csr = A.multiply(x);
    if (y_csr[0] == 3.0 && y_csr[1] == 7.0)
        std::cout << "csr good\n";
    else
        std::cout << "csr error\n";

    return 0;
}