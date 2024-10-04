#pragma once

#include <vector>

#include "s21_matrix.h"

namespace s21 {

template <typename T>
class Gauss {
 public:
  static std::vector<T> SolveSimple(const S21Matrix<T>& koef_matrix,
                                    const S21Matrix<T>& free_members_matrix);
  static std::vector<T> SolveWithMainElementChoice(
      const S21Matrix<T>& koef_matrix, const S21Matrix<T>& free_members_matrix);

 private:
  static constexpr double kPrecision = 1e-6;
};

template <typename T>
std::vector<T> Gauss<T>::SolveSimple(const S21Matrix<T>& koef_matrix,
                                     const S21Matrix<T>& free_members_matrix) {
  if (koef_matrix.GetCols() != koef_matrix.GetRows())
    throw std::invalid_argument("Koef matrix must be quadratic");
  if (koef_matrix.GetRows() != free_members_matrix.GetRows() ||
      free_members_matrix.GetCols() != 1)
    throw std::invalid_argument("Invalid free members matrix");

  S21Matrix<T> temp_matrix = koef_matrix;
  temp_matrix.SetCols(temp_matrix.GetCols() + 1);
  for (size_t row = 0; row < temp_matrix.GetRows(); ++row)
    temp_matrix(row, temp_matrix.GetCols() - 1) = free_members_matrix(row, 0);

  for (size_t row = 0; row < temp_matrix.GetRows(); ++row) {
    T koef = temp_matrix(row, row);
    for (size_t col = row; col < temp_matrix.GetCols(); ++col)
      temp_matrix(row, col) /= koef;

    for (size_t next_row = row + 1; next_row != temp_matrix.GetRows();
         ++next_row) {
      T koef = temp_matrix(next_row, row) / temp_matrix(row, row);
      for (size_t col = row; col < temp_matrix.GetCols(); ++col)
        temp_matrix(next_row, col) -= temp_matrix(row, col) * koef;
    }
  }

  for (int64_t row = temp_matrix.GetRows() - 1; row >= 0; --row) {
    for (int64_t prev_row = row - 1; prev_row >= 0; --prev_row) {
      T koef = temp_matrix(prev_row, row);
      for (int64_t col = temp_matrix.GetCols() - 1; col >= row; --col)
        temp_matrix(prev_row, col) -= temp_matrix(row, col) * koef;
    }
  }
  std::vector<double> result(temp_matrix.GetRows());
  for (size_t i = 0; i < result.size(); ++i) {
    result[i] = temp_matrix(i, temp_matrix.GetCols() - 1);
    std::cout << result[i] << ' ';
  }
  std::cout << "\n\n";
  std::cout << temp_matrix;
  return std::vector<T>();
}

template <typename T>
std::vector<T> Gauss<T>::SolveWithMainElementChoice(
    const S21Matrix<T>& koef_matrix, const S21Matrix<T>& free_members_matrix) {
  if (koef_matrix.GetCols() != koef_matrix.GetRows())
    throw std::invalid_argument("Koef matrix must be quadratic");
  if (koef_matrix.GetRows() != free_members_matrix.GetRows() ||
      free_members_matrix.GetCols() != 1)
    throw std::invalid_argument("Invalid free members matrix");

  S21Matrix<T> temp_matrix = koef_matrix;
  temp_matrix.SetCols(temp_matrix.GetCols() + 1);

  std::vector<T> change_vector(koef_matrix.GetRows());
  for (size_t row = 0; row < temp_matrix.GetRows(); ++row)
    temp_matrix(row, temp_matrix.GetCols() - 1) = free_members_matrix(row, 0);

  for (size_t col = 0; col < koef_matrix.GetCols(); ++col)
    change_vector[col] = col;

  for (size_t row = 0; row < temp_matrix.GetRows(); ++row) {
    T max_el = temp_matrix(row, row);
    size_t res_row = row, res_col = row;
    for (size_t f_row = row; f_row < temp_matrix.GetRows(); ++f_row)
      for (size_t f_col = row; f_col < temp_matrix.GetRows(); ++f_col) {
        if (fabs(temp_matrix(f_row, f_col)) > max_el) {
          max_el = fabs(temp_matrix(f_row, f_col));
          res_row = f_row;
          res_col = f_col;
        }
      }

    for (size_t f_col = row; f_col < temp_matrix.GetCols(); ++f_col) {
      std::swap(temp_matrix(row, f_col), temp_matrix(res_row, f_col));
    }

    for (size_t f_row = row; f_row < temp_matrix.GetRows(); ++f_row) {
      std::swap(temp_matrix(f_row, row), temp_matrix(f_row, res_col));
      if (f_row == row) std::swap(change_vector[row], change_vector[res_col]);
    }

    T koef = temp_matrix(row, row);
    for (size_t col = row; col < temp_matrix.GetCols(); ++col)
      temp_matrix(row, col) /= koef;

    for (size_t next_row = row + 1; next_row != temp_matrix.GetRows();
         ++next_row) {
      T koef = temp_matrix(next_row, row) / temp_matrix(row, row);
      for (size_t col = row; col < temp_matrix.GetCols(); ++col)
        temp_matrix(next_row, col) -= temp_matrix(row, col) * koef;
    }
  }

  for (int64_t row = temp_matrix.GetRows() - 1; row >= 0; --row) {
    for (int64_t prev_row = row - 1; prev_row >= 0; --prev_row) {
      T koef = temp_matrix(prev_row, row);
      for (int64_t col = temp_matrix.GetCols() - 1; col >= row; --col)
        temp_matrix(prev_row, col) -= temp_matrix(row, col) * koef;
    }
  }

  std::vector<T> result(temp_matrix.GetRows());
  for (size_t i = 0; i < result.size(); ++i)
    result[change_vector[i]] = temp_matrix(i, temp_matrix.GetCols() - 1);

  for (const auto& it : result) {
    std::cout << it << " ";
  }

  std::cout << "\n\n";
  std::cout << temp_matrix;
  return result;
}

}  // namespace s21