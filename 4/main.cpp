#include "nr3.h"
#include "svd.h"
#include "utilities.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

std::pair<std::vector<double>, std::vector<double>>
readData(const std::string &filename) { // Læse data (y, x)
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Could not open file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<double> x_vals, y_vals;
  double x, y;
  while (file >> y >> x) {
    x_vals.push_back(x);
    y_vals.push_back(y);
  }

  return {x_vals, y_vals};
}

void residualError(MatDoub A, VecDoub b, VecDoub x) {
  VecDoub top = A * x - b;
  Doub error = util::norm(top) / util::norm(b);
  std::cout << "\nRelative residual error epsilon_residual = " << error
            << std::endl;
}

void randomFittingError(const MatDoub &A) {
  int m = A.nrows(), n = A.ncols();
  Doub result = std::sqrt((m - n) / static_cast<Doub>(m));
  std::cout << "\nRandom fit Error= " << result << std::endl;
}

void sigErrorEstimate(MatDoub A, VecDoub b, VecDoub x, const SVD &svd) {
  int M = x.size();

  VecDoub sigma(M);

  for (int j = 0; j < M; j++) {
    double sum = 0.0;
    for (int i = 0; i < M; i++) {
      if (svd.w[i] > svd.tsh) { // values above threshold
        double term = svd.v[j][i] / svd.w[i];
        sum += term * term;
      }
    }
    sigma[j] = sqrt(sum);
  }

  std::cout << "\nAnd now the error:" << std::endl;
  std::cout << "sigma\tVector " << M << "D:" << std::endl;
  util::print(sigma);
}

void solvePontius(const std::string &filename) { // y = a0 + a1*x + a2*x^2
  auto [x_vals, y_vals] = readData(filename);
  const int N = x_vals.size();
  const int M = 3;

  MatDoub A(N, M);
  VecDoub b(N);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      A[i][j] = std::pow(x_vals[i], j);
    }
    b[i] = y_vals[i];
  }

  VecDoub x(M);
  SVD svd(A);

  svd.solve(b, x, svd.eps);

  std::cout << "\nBest fit parameters for Pontius using SVD: " << std::endl;
  util::print(x);

  sigErrorEstimate(A, b, x, svd);

  residualError(A, b, x);

  randomFittingError(A);
}

void solveFilip(
    const std::string &filename) { // y = B0 + B1*x + B2*(x**2) + ... +
                                   // B9*(x**9) + B10*(x**10) + e
  auto [x_vals, y_vals] = readData(filename);
  const int N = x_vals.size();
  const int M = 11;

  MatDoub A(N, M);
  VecDoub b(N);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      A[i][j] = std::pow(x_vals[i], j);
    }
    b[i] = y_vals[i];
  }

  VecDoub x(M);
  SVD svd(A);

  svd.solve(b, x, svd.eps);

  std::cout << "\nBest fit parameters for Filip using SVD: " << std::endl;
  util::print(x);

  sigErrorEstimate(A, b, x, svd);

  residualError(A, b, x);

  randomFittingError(A);
}

int main() {
  const std::string pontius_path =
      "/Users/patrickandersen/Desktop/6 semester/Numeriske "
      "Metoder/Code/4/PontiusData.dat";
  solvePontius(pontius_path);

  const std::string filip_path =
      "/Users/patrickandersen/Desktop/6 semester/Numeriske "
      "Metoder/Code/4/FilipData.dat";
  solveFilip(filip_path);

  return 0;
}
