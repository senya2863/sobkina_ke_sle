#include <vector>

class csr {
private:
    size_t rows, cols;
    std::vector<double> values;
    std::vector<size_t> colindex;
    std::vector<size_t> rowindex;

public:
    csr(size_t r, size_t c);
    void addvalue(size_t row, size_t col, double value);
    void final();
    std::vector<double> multiply(const std::vector<double>& vec) const;
};