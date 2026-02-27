#include <vector>

class densematrix{
  private:
    size_t rows, cols;
    std::vector<double> data;
   public:
     densematrix(int r, int c);

     double &at(int i, int j); //change
     double at(int i, int j) const;

     int getrows() const;
     int getcols() const;

     std::vector<double> multiply(const std::vector<double>& vec) const;
 };