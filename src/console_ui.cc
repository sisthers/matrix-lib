#include "gauss.h"
#include "s21_matrix.h"

int main() {
  s21::S21Matrix<long double> test_matrix(4, 4);
  test_matrix(0, 0) = 4.3;
  test_matrix(0, 1) = 2.8;
  test_matrix(0, 2) = 1.7;
  test_matrix(0, 3) = 5 * 1.73;

  test_matrix(1, 0) = 2.0;
  test_matrix(1, 1) = 1.3;
  test_matrix(1, 2) = 2.5 * 2.64;
  test_matrix(1, 3) = -3.4;

  test_matrix(2, 0) = -3.9;
  test_matrix(2, 1) = 5.4 * 2.23;
  test_matrix(2, 2) = -1.3;
  test_matrix(2, 3) = 0.6;

  test_matrix(3, 0) = 1.9 * 1.41;
  test_matrix(3, 1) = -2;
  test_matrix(3, 2) = 4.8;
  test_matrix(3, 3) = -1.6;

  s21::S21Matrix<long double> test_free_matrix(4, 1);
  test_free_matrix(0, 0) = -0.8;
  test_free_matrix(1, 0) = 1.7;
  test_free_matrix(2, 0) = -3.8;
  test_free_matrix(3, 0) = 2;
  s21::Gauss<long double>::SolveSimple(test_matrix, test_free_matrix);
}

/*(s21::S21Matrix<long double> test_matrix(4, 4);
  test_matrix(0, 0) = 4.3;
  test_matrix(0, 1) = 2.8;
  test_matrix(0, 2) = 1.7;
  test_matrix(0, 3) = 5 * SQRT_3;

  test_matrix(1, 0) = 2.0;
  test_matrix(1, 1) = 1.3;
  test_matrix(1, 2) = 2.5 * SQRT_7;
  test_matrix(1, 3) = -3.4;

  test_matrix(2, 0) = -3.9;
  test_matrix(2, 1) = 5.4 * SQRT_5;
  test_matrix(2, 2) = -1.3;
  test_matrix(2, 3) = 0.6;

  test_matrix(2, 0) = 1.9 * SQRT_2;
  test_matrix(2, 1) = -2;
  test_matrix(2, 2) = 4.8;
  test_matrix(2, 3) = -1.6;

  s21::S21Matrix<long double> test_free_matrix(4, 1);
  test_free_matrix(0, 0) = -0.8;
  test_free_matrix(1, 0) = 1.7;
  test_free_matrix(2, 0) = -3.8;
  test_free_matrix(3, 0) = 2;
  s21::Gauss<long double>::SolveSimple(test_matrix, test_free_matrix);)*/