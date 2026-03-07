#include <vector>

class densematrix {
private:
    size_t rows, cols;
    std::vector<double> data;
public:
    densematrix(size_t r, size_t c);

    double& at(size_t i, size_t j);
    double at(size_t i, size_t j) const;

    std::vector<double> multiply(const std::vector<double>& vec) const;
};