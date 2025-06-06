#include "nr3.h"
#include "svd.h"
#include "utilities.h"
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

  std::cout << "\nBest fit parameters for Prontius using SVD: " << std::endl;
  util::print(x);
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

  // util::print(svd.u);

  svd.solve(b, x, svd.eps);

  std::cout << "\nBest fit parameters for Filip using SVD: " << std::endl;
  util::print(x);
}

int main() {
  const std::string pontius_path =
      "/Users/patrickandersen/Desktop/6 semester/Numeriske "
      "Metoder/Code/3/PontiusData.dat";
  solvePontius(pontius_path);

  const std::string filip_path =
      "/Users/patrickandersen/Desktop/6 semester/Numeriske "
      "Metoder/Code/3/FilipData.dat";
  solveFilip(filip_path);

  return 0;
}
