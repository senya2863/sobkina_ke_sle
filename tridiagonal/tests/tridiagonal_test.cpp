#include <gtest/gtest.h>
#include "tridiagonal.h"

TEST(progonka, bebebe) {
    std::vector<double> a = {0,1,1};
    std::vector<double> b = {4,4,4};
    std::vector<double> c = {1,1,0};
    std::vector<double> d = {6,6,6};

    auto x = solve(a,b,c,d);

    // Просто проверяем, что вектор нужного размера
    EXPECT_EQ(x.size(), 3);
}