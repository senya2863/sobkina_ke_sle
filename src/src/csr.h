#include <vector>

class csr{
  private:
    size_t rows;
    size_t cols;

    std::vector<double> values;
    std::vector<size_t> colindex;
    std::vector<size_t> rowindex;

    public:
      csr(size_t rows, size_t cols);
      void addvalue(size_t row, size_t col, double value);
      void final();
      std
}
