#include <gtest/gtest.h>

#include <vector>

#include "los.h"
#include "s21_matrix.h"

#define PRESICION_LOS 1e-3

std::tuple<s21::S21Matrix<long double>, s21::S21Matrix<long double>>
MatrixForTestLOS() {
  s21::S21Matrix<long double> test_matrix(4, 4);
  test_matrix(0, 0) = 4.5;
  test_matrix(0, 1) = 4.21;
  test_matrix(0, 2) = -5.21;
  test_matrix(0, 3) = 1.23;

  test_matrix(1, 0) = 2.0;
  test_matrix(1, 1) = 1.87;
  test_matrix(1, 2) = 0.15;
  test_matrix(1, 3) = -2.11;

  test_matrix(2, 0) = 0.12;
  test_matrix(2, 1) = 0.73;
  test_matrix(2, 2) = 0.31;
  test_matrix(2, 3) = 0.83;

  test_matrix(3, 0) = -1.42;
  test_matrix(3, 1) = 3.41;
  test_matrix(3, 2) = 1.73;
  test_matrix(3, 3) = 0.69;

  s21::S21Matrix<long double> test_free_matrix(4, 1);
  test_free_matrix(0, 0) = 4.7;
  test_free_matrix(1, 0) = 1.84;
  test_free_matrix(2, 0) = 2.98;
  test_free_matrix(3, 0) = 4.31;
  return std::make_tuple(test_matrix, test_free_matrix);
}

TEST(los, solve) {
  auto matrices = MatrixForTestLOS();
  auto koef_matrix = std::get<0>(matrices);
  auto free_matrix = std::get<1>(matrices);
  std::cout << "Koef matrix: \n" << koef_matrix << '\n';
  std::cout << "Free members matrix: \n" << free_matrix << '\n';

  for (double koef = 0.1; koef <= 1.1; koef += 0.1) {
    auto result = s21::LOS<long double>::Solve(koef_matrix, free_matrix, koef);
    auto result_m = std::get<0>(result);
    auto epoches = std::get<1>(result);
    std::cout << "Result with a = " << koef << ": " << epoches << " epoches\n";
    for (const auto &it : result_m) std::cout << it << " ";
    std::cout << "\n\n";
    long double dif_sum = 0;
    for (size_t row = 0; row < koef_matrix.GetRows(); ++row) {
      long double res = 0;
      for (size_t col = 0; col < koef_matrix.GetCols(); ++col)
        res += koef_matrix(row, col) * result_m[col];
      long double dif = fabs(res - free_matrix(row, 0));
      std::cout << "Difference for row " << row + 1 << " = "
                << fabs(res - free_matrix(row, 0)) << '\n';
      dif_sum += dif;
    }
    std::cout << "\nTotal difference = " << dif_sum << "\n\n";
    std::cout << "-----------------------------------\n";
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
