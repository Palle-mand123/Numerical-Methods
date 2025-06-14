#include "nr3.h"
#include "svd.h"
#include "utilities.h"
#include <print>

MatDoub readMatrix(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  int M, N;
  file >> M >> N;

  MatDoub matrix(M, N);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      file >> matrix[i][j];
    }
  }

  file.close();
  return matrix;
}

VecDoub readVector(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  int M, N;
  file >> M >> N;
  if (N != 1) {
    std::cerr << "Expected a vector (N=1), but got N=" << N << std::endl;
    exit(1);
  }

  VecDoub vec(M);
  for (int i = 0; i < M; i++) {
    file >> vec[i];
  }

  file.close();
  return vec;
}

void problemI(MatDoub A) {
  SVD svd(A);
  std::cout << "Diagonal elements of W:\n" << std::endl;
  util::print(svd.w);
}

void problemII(MatDoub A) {
  SVD svd(A);

  util::print(svd.nullspace(1e-10));
}

VecDoub problemIII(VecDoub b, MatDoub A) {
  VecDoub x(A.ncols());

  SVD svd(A);

  svd.solve(b, x);

  std::cout << "Best fit solution parameters x:\n" << std::endl;
  util::print(x);
  return x;
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
  std::cout << "\nRandom fit Error = " << result << std::endl;
}

void sigErrorEstimate(MatDoub A, VecDoub b, VecDoub x) {
  SVD svd(A);
  VecDoub sigma = VecDoub(svd.w.size());

  for (int j = 0; j < svd.w.size(); j++) {
    double sum = 0.0;
    for (int i = 0; i < svd.w.size(); i++) {
      if (svd.w[i] > svd.tsh) {
        sum += pow(svd.v[j][i] / svd.w[i], 2);
      }
    }
    sigma[j] = sqrt(sum);
  }

  std::cout << "sigma\tVector " << svd.w.size() << "D:" << std::endl;
  util::print(sigma);
}

void problemIV(MatDoub A, VecDoub b, VecDoub x) {
  residualError(A, b, x);
  randomFittingError(A);
  std::cout << "\n" << std::endl;
  sigErrorEstimate(A, b, x);
  std::cout << "\n" << std::endl;
}

int main() {

  MatDoub A = readMatrix("/Users/patrickandersen/Desktop/6 semester/Numeriske "
                         "Metoder/Code/Eksamen/1/NUM S25_Ex1A.dat");
  VecDoub b = readVector("/Users/patrickandersen/Desktop/6 semester/Numeriske "
                         "Metoder/Code/Eksamen/1/NUM S25_Ex1b.dat");

  std::cout << "\n------------------ Problem i ------------------------\n"
            << std::endl;
  problemI(A);
  std::cout << "\n------------------ Problem ii ------------------------\n"
            << std::endl;
  problemII(A);
  std::cout << "\n------------------ Problem iii ------------------------\n"
            << std::endl;
  VecDoub x = problemIII(b, A);
  std::cout << "\n------------------ Problem iV ------------------------\n"
            << std::endl;
  problemIV(A, b, x);

  return 0;
}
