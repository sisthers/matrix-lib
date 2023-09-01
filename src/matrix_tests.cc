#include <gtest/gtest.h>

#include <vector>

#include "s21_matrix.h"

TEST(funcs, sum) {
  s21::S21Matrix<double> m1(2, 5);
  int c = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) m1(i, j) = c++;
  s21::S21Matrix<double> m2(m1);
  m2.SumMatrix(m1);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) ASSERT_NEAR(m1(i, j) * 2, m2(i, j), 1e-6);
}
TEST(funcs, sum_fail) {
  s21::S21Matrix<double> m1(2, 5);
  s21::S21Matrix<double> m2(2, 4);
  bool catched = false;
  try {
    m1.SumMatrix(m2);
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(funcs, sub) {
  s21::S21Matrix<double> m1(2, 5);
  int c = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) m1(i, j) = c++;
  s21::S21Matrix<double> m2(m1);
  m2.SubMatrix(m1);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) ASSERT_NEAR(0.0, m2(i, j), 1e-6);
}

TEST(funcs, sub_fail) {
  s21::S21Matrix<double> m1(2, 5);
  s21::S21Matrix<double> m2(2, 4);
  bool catched = false;
  try {
    m1.SubMatrix(m2);
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(funcs, muln) {
  s21::S21Matrix<double> m1(2, 5);
  int c = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) m1(i, j) = c++;
  s21::S21Matrix<double> m2(m1);
  m1.MulNumber(12.5);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) ASSERT_NEAR(m1(i, j), m2(i, j) * 12.5, 1e-6);
}

TEST(funcs, mulm) {
  s21::S21Matrix<double> m1;
  s21::S21Matrix<double> m2;
  s21::S21Matrix<double> m3;
  int c = 1;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++, c++) {
      m1(i, j) = c;
      m2(i, j) = c * c;
    }
  m1.MulMatrix(m2);
  m3(0, 0) = 180;
  m3(0, 1) = 246;
  m3(0, 2) = 324;
  m3(1, 0) = 378;
  m3(1, 1) = 525;
  m3(1, 2) = 702;
  m3(2, 0) = 576;
  m3(2, 1) = 804;
  m3(2, 2) = 1080;
  ASSERT_EQ(m3.EqMatrix(m1), true);

  s21::S21Matrix<double> m4(1, 5);
  s21::S21Matrix<double> m5(5, 2);
  s21::S21Matrix<double> m6(1, 2);
  for (int i = 0, c = 1; i < 5; i++, c++) m4(0, i) = c;
  for (int i = 0, c = 1; i < 5; i++)
    for (int j = 0; j < 2; j++, c++) m5(i, j) = c;
  m6(0, 0) = 95;
  m6(0, 1) = 110;
  m4.MulMatrix(m5);
  ASSERT_EQ(m6.EqMatrix(m4), true);
}

TEST(funcs, mulm_fail) {
  s21::S21Matrix<double> m1(3, 5);
  s21::S21Matrix<double> m2(4, 3);
  bool catched = false;
  try {
    m1.MulMatrix(m2);
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(funcs, trans) {
  s21::S21Matrix<double> m1(2, 4);
  s21::S21Matrix<double> m2(4, 2);
  m1(0, 0) = 14;
  m1(0, 1) = -12;
  m1(0, 2) = 45;
  m1(0, 3) = 12;
  m1(1, 0) = 3;
  m1(1, 1) = 0;
  m1(1, 2) = 2;
  m1(1, 3) = 5;
  m2(0, 0) = 14;
  m2(0, 1) = 3;
  m2(1, 0) = -12;
  m2(1, 1) = 0;
  m2(2, 0) = 45;
  m2(2, 1) = 2;
  m2(3, 0) = 12;
  m2(3, 1) = 5;
  s21::S21Matrix<double> m3 = m1.Transpose();
  ASSERT_EQ(m3.EqMatrix(m2), true);
}

TEST(funcs, det) {
  s21::S21Matrix<double> m1(5, 5);
  m1(0, 0) = 15.16;
  m1(0, 1) = 11.51;
  m1(0, 2) = -4.37;
  m1(0, 3) = 6.70;
  m1(0, 4) = 2.66;
  m1(1, 0) = -0.51;
  m1(1, 1) = 0.94;
  m1(1, 2) = 9.83;
  m1(1, 3) = -9.44;
  m1(1, 4) = 7.91;
  m1(2, 0) = -5.64;
  m1(2, 1) = -5.80;
  m1(2, 2) = 15.42;
  m1(2, 3) = 7.41;
  m1(2, 4) = 5.54;
  m1(3, 0) = -9.79;
  m1(3, 1) = 19.11;
  m1(3, 2) = 14.98;
  m1(3, 3) = -4.88;
  m1(3, 4) = 4.76;
  m1(4, 0) = 7.07;
  m1(4, 1) = -2.91;
  m1(4, 2) = 9.48;
  m1(4, 3) = 5.90;
  m1(4, 4) = 17.21;
  ASSERT_NEAR(m1.Determinant(), 805689.261085, 1e-6);

  s21::S21Matrix<double> m2(2, 2);
  m2(0, 0) = -8;
  m2(0, 1) = -7;
  m2(1, 0) = 10;
  m2(1, 1) = 3;
  ASSERT_NEAR(m2.Determinant(), 46, 1e-6);
}

TEST(funcs, det_fail) {
  s21::S21Matrix<double> m1(2, 3);
  bool catched = false;
  try {
    m1.Determinant();
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(funcs, calcc) {
  s21::S21Matrix<double> m1(4, 4);
  s21::S21Matrix<double> m2(4, 4);
  int a = 1;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++, a++) m1(i, j) = a * a;
  m2(0, 0) = -512;
  m2(0, 1) = 1536;
  m2(0, 2) = -1536;
  m2(0, 3) = 512;
  m2(1, 0) = 1536;
  m2(1, 1) = -4608;
  m2(1, 2) = 4608;
  m2(1, 3) = -1536;
  m2(2, 0) = -1536;
  m2(2, 1) = 4608;
  m2(2, 2) = -4608;
  m2(2, 3) = 1536;
  m2(3, 0) = 512;
  m2(3, 1) = -1536;
  m2(3, 2) = 1536;
  m2(3, 3) = -512;
  s21::S21Matrix<double> m3 = m1.CalcComplements();
  ASSERT_EQ(m2.EqMatrix(m3), true);

  s21::S21Matrix<double> m4(3, 3);
  s21::S21Matrix<double> m5(3, 3);
  m4(0, 0) = 1;
  m4(0, 1) = 2;
  m4(0, 2) = 3;
  m4(1, 0) = 0;
  m4(1, 1) = 4;
  m4(1, 2) = 2;
  m4(2, 0) = 5;
  m4(2, 1) = 2;
  m4(2, 2) = 1;
  m5(0, 0) = 0;
  m5(0, 1) = 10;
  m5(0, 2) = -20;
  m5(1, 0) = 4;
  m5(1, 1) = -14;
  m5(1, 2) = 8;
  m5(2, 0) = -8;
  m5(2, 1) = -2;
  m5(2, 2) = 4;
  s21::S21Matrix<double> m6 = m4.CalcComplements();
  ASSERT_EQ(m5.EqMatrix(m6), true);

  s21::S21Matrix<double> m7(2, 2);
  s21::S21Matrix<double> m8(2, 2);
  m7(0, 0) = 1;
  m7(0, 1) = 2;
  m7(1, 0) = 3;
  m7(1, 1) = 4;
  m8(0, 0) = 4;
  m8(0, 1) = -3;
  m8(1, 0) = -2;
  m8(1, 1) = 1;
  s21::S21Matrix<double> m9 = m7.CalcComplements();
  ASSERT_EQ(m8.EqMatrix(m9), true);
}

TEST(funcs, calcc_fail) {
  s21::S21Matrix<double> m1(2, 3);
  bool catched = false;
  try {
    s21::S21Matrix<double> m2 = m1.CalcComplements();
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(funcs, eq) {
  s21::S21Matrix<double> m1;
  s21::S21Matrix<double> m2;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      m1(i, j) = i + j;
      m2(i, j) = i + j;
    }
  ASSERT_EQ(m1.EqMatrix(m2), true);

  s21::S21Matrix<double> m3 = m1;
  ASSERT_EQ(m3.EqMatrix(m1), true);

  s21::S21Matrix<double> m4(m3);
  ASSERT_EQ(m4.EqMatrix(m1), true);

  m4(0, 0) = 10;
  ASSERT_EQ(m4.EqMatrix(m1), false);

  s21::S21Matrix<double> m5(2, 4);
  s21::S21Matrix<double> m6(2, 3);
  ASSERT_EQ(m5.EqMatrix(m6), false);
}

TEST(funcs, inv) {
  s21::S21Matrix<double> m1(4, 4);
  s21::S21Matrix<double> m2(4, 4);
  int c = 1;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++, c++) {
      m1(i, j) = c;
      if (i == j)
        m2(i, j) = 1;
      else
        m2(i, j) = 0;
    }
  m1(0, 0) = 17;
  m1(3, 3) = 17;
  s21::S21Matrix<double> m3 = m1.InverseMatrix();
  m1.MulMatrix(m3);
  ASSERT_EQ(m2.EqMatrix(m1), true);

  s21::S21Matrix<double> m4;
  s21::S21Matrix<double> m5;
  m4(0, 0) = 1;
  m4(0, 1) = 2;
  m4(0, 2) = 3;
  m4(1, 0) = 0;
  m4(1, 1) = 4;
  m4(1, 2) = 2;
  m4(2, 0) = 5;
  m4(2, 1) = 2;
  m4(2, 2) = 1;
  m5(0, 0) = (double)0;
  m5(0, 1) = (double)-1 / 10;
  m5(0, 2) = (double)1 / 5;
  m5(1, 0) = (double)-1 / 4;
  m5(1, 1) = (double)7 / 20;
  m5(1, 2) = (double)1 / 20;
  m5(2, 0) = (double)1 / 2;
  m5(2, 1) = (double)-1 / 5;
  m5(2, 2) = (double)-1 / 10;
  s21::S21Matrix<double> m6 = m4.InverseMatrix();
  ASSERT_EQ(m5.EqMatrix(m6), true);

  s21::S21Matrix<double> m7(1, 1);
  s21::S21Matrix<double> m8(1, 1);
  m7(0, 0) = 10;
  m8(0, 0) = 1;
  s21::S21Matrix<double> m9 = m7.InverseMatrix();
  m7.MulMatrix(m9);
  ASSERT_EQ(m8.EqMatrix(m7), true);

  s21::S21Matrix<double> m10(2, 2);
  s21::S21Matrix<double> m11(2, 2);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
      if (i == j)
        m11(i, j) = 1;
      else
        m11(i, j) = 0;
    }
  m10(0, 0) = 15;
  m10(0, 1) = 2;
  m10(1, 0) = 21;
  m10(1, 1) = 3;
  s21::S21Matrix<double> m12 = m10.InverseMatrix();
  m10.MulMatrix(m12);
  ASSERT_EQ(m11.EqMatrix(m10), true);
}

TEST(funcs, inv_fail) {
  s21::S21Matrix<double> m1(2, 3);
  bool catched = false;
  try {
    s21::S21Matrix<double> m2 = m1.InverseMatrix();
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);

  s21::S21Matrix<double> m3(4, 4);
  for (int i = 0, c = 0; i < 4; i++)
    for (int j = 0; j < 4; j++, c++) m3(i, j) = c * c;
  catched = false;
  try {
    s21::S21Matrix<double> m4 = m3.InverseMatrix();
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(opers, bnplus) {
  s21::S21Matrix<double> m1;
  s21::S21Matrix<double> m2;
  s21::S21Matrix<double> m3;
  int c1 = 1, c2 = 2;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      m3(i, j) = c1 + c2;
      m1(i, j) = c1++;
      m2(i, j) = c2++;
    }
  s21::S21Matrix<double> m4 = m1 + m2;
  ASSERT_EQ(m3.EqMatrix(m4), true);
}
TEST(opers, unplus) {
  s21::S21Matrix<double> m1(2, 2);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) m1(i, j) = i + j + 1;
  s21::S21Matrix<double> m2 = +m1;
  ASSERT_EQ(m1.EqMatrix(m2), true);
}

TEST(opers, bnminus) {
  s21::S21Matrix<double> m1;
  s21::S21Matrix<double> m2;
  s21::S21Matrix<double> m3;
  int c1 = 1, c2 = 2;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      m3(i, j) = c1 - c2;
      m1(i, j) = c1++;
      m2(i, j) = c2++;
    }
  s21::S21Matrix<double> m4 = m1 - m2;
  ASSERT_EQ(m3.EqMatrix(m4), true);
}

TEST(opers, unminus) {
  s21::S21Matrix<double> m1(2, 2);
  s21::S21Matrix<double> m2(2, 2);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
      m1(i, j) = i + j + 1;
      m2(i, j) = -(i + j + 1);
    }
  s21::S21Matrix<double> m3 = -m1;
  ASSERT_EQ(m2.EqMatrix(m3), true);
}

TEST(opers, mulm) {
  s21::S21Matrix<double> m1(1, 5);
  s21::S21Matrix<double> m2(5, 2);
  s21::S21Matrix<double> m3(1, 2);
  for (int i = 0, c = 1; i < 5; i++, c++) m1(0, i) = c;
  for (int i = 0, c = 1; i < 5; i++)
    for (int j = 0; j < 2; j++, c++) m2(i, j) = c;
  m3(0, 0) = 95;
  m3(0, 1) = 110;
  s21::S21Matrix<double> m4 = m1 * m2;
  ASSERT_EQ(m4.EqMatrix(m3), true);
}

TEST(opers, muln) {
  s21::S21Matrix<double> m1(2, 5);
  int c = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) m1(i, j) = c++;
  s21::S21Matrix<double> m2 = m1 * 12.5;
  s21::S21Matrix<double> m3 = 12.5 * m1;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) {
      ASSERT_NEAR(m1(i, j) * 12.5, m2(i, j), 1e-6);
      ASSERT_NEAR(m1(i, j) * 12.5, m3(i, j), 1e-6);
    }
}

TEST(opers, plas) {
  s21::S21Matrix<double> m1;
  s21::S21Matrix<double> m2;
  s21::S21Matrix<double> m3;
  int c1 = 1, c2 = 2;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      m3(i, j) = c1 + c2;
      m1(i, j) = c1++;
      m2(i, j) = c2++;
    }
  m1 += m2;
  ASSERT_EQ(m1.EqMatrix(m3), true);
}

TEST(opers, mnas) {
  s21::S21Matrix<double> m1;
  s21::S21Matrix<double> m2;
  s21::S21Matrix<double> m3;
  int c1 = 1, c2 = 2;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      m3(i, j) = c1 - c2;
      m1(i, j) = c1++;
      m2(i, j) = c2++;
    }
  m1 -= m2;
  ASSERT_EQ(m1.EqMatrix(m3), true);
}

TEST(opers, mulnas) {
  s21::S21Matrix<double> m1(2, 5);
  int c = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) m1(i, j) = c++;
  m1 *= 12.5;
  c = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++) {
      ASSERT_NEAR(m1(i, j), c++ * 12.5, 1e-6);
    }
}

TEST(opers, mulmas) {
  s21::S21Matrix<double> m1(1, 5);
  s21::S21Matrix<double> m2(5, 2);
  s21::S21Matrix<double> m3(1, 2);
  for (int i = 0, c = 1; i < 5; i++, c++) m1(0, i) = c;
  for (int i = 0, c = 1; i < 5; i++)
    for (int j = 0; j < 2; j++, c++) m2(i, j) = c;
  m3(0, 0) = 95;
  m3(0, 1) = 110;
  m1 *= m2;
  ASSERT_EQ(m1.EqMatrix(m3), true);
}

TEST(opers, eq) {
  s21::S21Matrix<double> m1;
  s21::S21Matrix<double> m2;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      m1(i, j) = i + j;
      m2(i, j) = i + j;
    }
  ASSERT_EQ(m1 == m2, true);
  m2(0, 0) = 100;
  ASSERT_EQ(m1 == m2, false);
}

TEST(opers, par) {
  s21::S21Matrix<double> m1(1, 1);
  m1(0, 0) = 123.123;
  ASSERT_NEAR(m1(0, 0), 123.123, 1e-6);
}

TEST(opers, par_fail) {
  s21::S21Matrix<double> m1(1, 1);
  bool catched = false;
  try {
    m1(2, 0) = 123.123;
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(opers, parconst) {
  s21::S21Matrix<double> m1(1, 1);
  m1(0, 0) = 123.123;
  const s21::S21Matrix<double> m2(m1);
  ASSERT_NEAR(m2(0, 0), 123.123, 1e-6);
}

TEST(opers, parconst_fail) {
  s21::S21Matrix<double> m1(1, 1);
  m1(0, 0) = 123.123;
  const s21::S21Matrix<double> m2(m1);
  bool catched = false;
  try {
    m1(0, 0) = m2(-1, 0);
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(opers, as) {
  s21::S21Matrix<double> m1;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m1(i, j) = i + j;
  s21::S21Matrix<double> m2 = m1.Transpose();
  m1 = m1.Transpose();
  ASSERT_EQ(m1.EqMatrix(m2), true);
}

TEST(cons, defaultcon) {
  s21::S21Matrix<double> m1;
  ASSERT_EQ(m1.GetRows(), 3);
  ASSERT_EQ(m1.GetCols(), 3);
}

TEST(cons, paramscon) {
  s21::S21Matrix<double> m1(2, 9);
  ASSERT_EQ(m1.GetRows(), 2);
  ASSERT_EQ(m1.GetCols(), 9);
}

TEST(cons, paramscon_fail) {
  bool catched = false;
  try {
    s21::S21Matrix<double> m1(0, -5);
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(cons, copycon) {
  s21::S21Matrix<double> m1(2, 2);
  s21::S21Matrix<double> m2(m1);
  ASSERT_EQ(m1.EqMatrix(m2), true);
}

TEST(cons, movecon) {
  std::vector<s21::S21Matrix<double>> vec;
  vec.push_back(s21::S21Matrix<double>(4, 4));
  ASSERT_EQ(vec[0].EqMatrix(s21::S21Matrix<double>(4, 4)), true);
}

TEST(cons, rowsgetset) {
  s21::S21Matrix<double> m1(2, 2);
  s21::S21Matrix<double> m2(1, 2);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) m1(i, j) = i + j;
  m2(0, 0) = 0;
  m2(0, 1) = 1;
  m1.SetRows(1);
  ASSERT_EQ(m2.EqMatrix(m1), true);
  ASSERT_EQ(m1.GetRows(), 1);

  m1.SetRows(3);
  s21::S21Matrix<double> m3(3, 2);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 2; j++) m3(i, j) = 0;
  m3(0, 1) = 1;
  ASSERT_EQ(m3.EqMatrix(m1), true);
  ASSERT_EQ(m1.GetRows(), 3);
}

TEST(cons, rowsgetset_fail) {
  s21::S21Matrix<double> m1(2, 2);
  bool catched = false;
  try {
    m1.SetRows(0);
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

TEST(cons, colsgetset) {
  s21::S21Matrix<double> m1(2, 2);
  s21::S21Matrix<double> m2(2, 1);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) m1(i, j) = i + j;
  m2(0, 0) = 0;
  m2(1, 0) = 1;
  m1.SetCols(1);
  ASSERT_EQ(m2.EqMatrix(m1), true);
  ASSERT_EQ(m1.GetCols(), 1);

  m1.SetCols(3);
  s21::S21Matrix<double> m3(2, 3);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 3; j++) m3(i, j) = 0;
  m3(1, 0) = 1;
  ASSERT_EQ(m3.EqMatrix(m1), true);
  ASSERT_EQ(m1.GetCols(), 3);
}

TEST(cons, colsgetset_fail) {
  s21::S21Matrix<double> m1(2, 2);
  bool catched = false;
  try {
    m1.SetCols(0);
  } catch (std::out_of_range &exp) {
    catched = true;
  }
  ASSERT_EQ(catched, true);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
