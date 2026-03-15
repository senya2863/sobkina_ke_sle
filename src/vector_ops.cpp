#include "vector_ops.h"

std::vector<double> add(const std::vector<double>& a,const std::vector<double>& b) {
    std::vector<double> res(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        res[i] = a[i] + b[i];
    }
    return res;
}

double dot(const std::vector<double>& a,const std::vector<double>& b) {
    double sum = 0;

    for (size_t i = 0; i < a.size(); i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

std::vector<double> scalarmul(const std::vector<double>& v,double k) {
    std::vector<double> res(v.size());

    for (size_t i = 0; i < v.size(); i++) {
        res[i] = v[i] * k;
    }
    return res;
}