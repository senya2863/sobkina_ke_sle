#include "dense.h"

densematrix::densematrix(size_t r, size_t c) {
    rows = r;
    cols = c;
    data.resize(r * c);
}
double& densematrix::at(size_t i, size_t j) {
    return data[i*cols+ j];
}
double densematrix::at(size_t i, size_t j) const {
    return data[i *cols + j];
}

std::vector<double> densematrix::multiply(const std::vector<double>& vec) const {
    std::vector<double> result(rows);
    for (size_t i = 0; i < rows; i++) {
        result[i] = 0;
        for (size_t j = 0; j < cols; j++) {
            result[i] += at(i, j) * vec[j];
        }
    }
    return result;
}