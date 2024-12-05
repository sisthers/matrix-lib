#pragma once

#include <tuple>
#include <vector>

#include "s21_matrix.h"

namespace s21 {

template <typename T>
class LOS {
 public:
  static std::tuple<std::vector<T>, uint32_t> Solve(
      const S21Matrix<T>& koef_matrix, const S21Matrix<T>& free_members_matrix,
      const T ac_koef);
};

template <typename T>
std::tuple<std::vector<T>, uint32_t> LOS<T>::Solve(
    const S21Matrix<T>& koef_matrix, const S21Matrix<T>& free_members_matrix,
    const T ac_koef) {
  T precision = 1e-3;
  auto C_matrix = koef_matrix.Transpose() * koef_matrix;

  T max_row_sum = 0;
  for (uint32_t row = 0; row < C_matrix.GetRows(); ++row) {
    T sum = 0;
    for (uint32_t col = 0; col < C_matrix.GetCols(); ++col)
      sum += C_matrix(row, col);
    if (sum >= max_row_sum) max_row_sum = sum;
  }

  T max_col_sum = 0;
  for (uint32_t col = 0; col < C_matrix.GetCols(); ++col) {
    T sum = 0;
    for (uint32_t row = 0; row < C_matrix.GetRows(); ++row)
      sum += C_matrix(row, col);
    if (sum >= max_col_sum) max_col_sum = sum;
  }

  T sum_sqrt = 0;
  for (uint32_t row = 0; row < C_matrix.GetRows(); ++row)
    for (uint32_t col = 0; col < C_matrix.GetCols(); ++col)
      sum_sqrt += C_matrix(row, col) * C_matrix(row, col);

  sum_sqrt = sqrt(sum_sqrt);

  T min_delta = max_row_sum;
  if (max_col_sum < min_delta) min_delta = max_col_sum;
  if (sum_sqrt < min_delta) min_delta = sum_sqrt;

  auto E_matrix = S21Matrix<T>(C_matrix.GetRows());
  for (uint32_t row = 0; row < C_matrix.GetRows(); ++row)
    E_matrix(row, row) = 1;

  auto G_matrix = E_matrix - ac_koef * C_matrix * (1.0 / min_delta);
  auto F_matrix = ac_koef * koef_matrix.Transpose() * free_members_matrix *
                  (1.0 / min_delta);

  auto X_matrix = F_matrix;
  uint32_t epoch = 1;
  while (true) {
    auto result = G_matrix * X_matrix + F_matrix;
    auto dif = koef_matrix * result - free_members_matrix;
    T dif_norm = 0;
    for (uint32_t row = 0; row < dif.GetRows(); ++row)
      dif_norm += fabs(dif(row, 0));
    if (dif_norm <= precision) {
      std::vector<T> result_vector(result.GetRows());
      for (size_t row = 0; row < result_vector.size(); ++row)
        result_vector[row] = result(row, 0);

      return std::make_tuple(result_vector, epoch);
    }
    X_matrix = result;
    ++epoch;
  }
}

}  // namespace s21