#include "csr.h"

csr::csr(size_t r, size_t c) : rows(r), cols(c) {
  rowindex.resize(rows + 1);
}

void csr::addvalue(size_t row, size_t col, double value) {
  values.push_back(value);
  colindex.push_back(col);
  rowindex[row + 1]++;
}

void csr::final() {
  for (size_t i = 1; i <= rows; ++i) {
    rowindex[i] += rowindex[i - 1];
  }
}

std::vector<double> csr::multiply(const std::vector<double>& vec) const {
  std::vector<double> result(rows);
  for (size_t i = 0; i < rows; i++) {
    result[i]=0;
    for (size_t j = rowindex[i]; j < rowindex[i + 1]; j++) {
      result[i] += values[j] * vec[colindex[j]];
    }
  }
  return result;
}