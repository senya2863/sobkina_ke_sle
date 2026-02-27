#include "dense.h"

densematrix:: densematrix(int r, int c){
  rows =r;
  cols =c;
  data.resize(r*c);
}

double& DenseMatrix::at(size_t i, size_t j)
{
  return data[i * cols + j];
}

double DenseMatrix::at(size_t i, size_t j) const
{
  return data[i * cols + j];
}

size_t DenseMatrix::getRows() const
{
  return rows;
}

size_t DenseMatrix::getCols() const
{
  return cols;
}

std::vector<double> DenseMatrix::multiply(const std::vector<double>& vec) const
{
  std::vector<double> result(rows);
  for (size_t i = 0; i < rows; i++)
  {
    result[i] = 0;
    for (size_t j = 0; j < cols; j++)
    {
      result[i] += at(i, j) * vec[j];
    }
  }
  return result;
}