#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "nr3.h"
#include <print>
#include <string>

namespace util {
void print(MatDoub mat, std::string symbol = "") {
  if (symbol.compare(""))
    std::println("{} , mat_size {}x{}", symbol, mat.nrows(), mat.ncols());

  for (int m = 0; m < mat.nrows(); m++) {
    for (int n = 0; n < mat.ncols(); n++) {
      std::print("{:^15}\t", mat[m][n]);
    }
    std::println("");
  }
  std::println("");
}

void printDiag(MatDoub mat, std::string symbol = "") {
  if (symbol.compare(""))
    std::println("{} , mat_diag_size {}x{}", symbol, mat.nrows(), mat.ncols());
  double nmax = mat.nrows() < mat.nrows() ? mat.nrows() : mat.nrows();
  for (int n = 0; n < nmax; n++) {
    std::print("{:^15}\t", mat[n][n]);
  }
  std::println("");
}

MatDoub diag(VecDoub &V) {
  double m = V.size();
  MatDoub M(m, m);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      M[i][j] = 0;

  for (int i = 0; i < m; i++)
    M[i][i] = V[i];

  return M;
}

void print(VecDoub vec, std::string symbol = "") {
  if (symbol.compare("")) {
    std::println("{} , vec_size {}", symbol, vec.size());
  }

  for (int m = 0; m < vec.size(); m++) {
    std::print("{:^15}\t", vec[m]);
  }
  std::println("");
}

MatDoub Transpose(const MatDoub &Mat) {
  MatDoub res(Mat.ncols(), Mat.nrows());
  for (int n = 0; n < res.nrows(); n++) {
    for (int m = 0; m < res.ncols(); m++) {
      res[n][m] = Mat[m][n];
    }
  }
  return res;
}

MatDoub T(const MatDoub &Mat) { return Transpose(Mat); }

double norm(VecDoub &v) {
  double res = 0;
  for (int i = 0; i < v.size(); i++) {
    res += v[i] * v[i];
  }
  return sqrt(res);
}
} // namespace util

MatDoub operator*(const MatDoub &A1, const MatDoub &A2) {
  if (A1.ncols() != A2.nrows()) {
    std::println(stderr, "in prod: the number of rows in A1 is not equal to "
                         "the number of cols in A2");
  }

  MatDoub res(A1.nrows(), A2.ncols());
  for (int n = 0; n < A1.nrows(); n++) {
    for (int m = 0; m < A2.ncols(); m++) {
      double temp = 0;
      for (int i = 0; i < A1.ncols(); i++) {
        temp += A1[n][i] * A2[i][m];
      }
      res[n][m] = temp;
    }
  }
  return res;
}

VecDoub operator*(const MatDoub &A, const VecDoub &b) {
  if (A.ncols() != b.size()) {
    std::println(stderr, "in prod: the number of rows in A is not equal to the "
                         "size of vector b");
  }

  VecDoub res(A.nrows());
  for (int n = 0; n < A.nrows(); n++) {
    double temp = 0;
    for (int m = 0; m < A.ncols(); m++) {
      temp += A[n][m] * b[m];
    }
    res[n] = temp;
  }
  return res;
}

VecDoub operator-(const VecDoub &v1, const VecDoub &v2) {
  if (v1.size() != v2.size()) {
    std::println(
        stderr,
        "in minus: the size of  cclv1 is not equal to the size of vector b");
  }
  VecDoub res(v1.size());
  for (int i = 0; i < v1.size(); i++) {
    res[i] = v1[i] - v2[i];
  }
  return res;
}

VecDoub operator+(const VecDoub &a, const VecDoub &b) {
  if (a.size() != b.size()) {
    std::println(stderr, "in prod: the number of rows in A is not equal to the "
                         "size of vector b");
  }
  VecDoub res(a.size());
  for (int i = 0; i < a.size(); i++) {
    res[i] = a[i] + b[i];
  }
  return res;
}

VecDoub operator/(const VecDoub &v, double s) {
  VecDoub res(v.size());
  for (int i = 0; i < v.size(); i++) {
    res[i] = v[i] / s;
  }
  return res;
}

VecDoub operator*(const VecDoub &v, double s) {
  VecDoub res(v.size());
  for (int i = 0; i < v.size(); i++) {
    res[i] = v[i] * s;
  }
  return res;
}

VecDoub operator*(double s, const VecDoub &v) { return v * s; }

#endif /* UTILITIES_H_ */
