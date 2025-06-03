#include "nr3.h"
#include "svd.h"
#include "utilities.h"
#include <print>
#include <utility>

MatDoub readMatrix(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(1);
  }

  int M, N;
  file >> M >> N; // Read number of rows and columns

  MatDoub matrix(M, N); // Create matrix with M rows and N columns
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      file >> matrix[i][j]; // Read each element
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
  if (N != 1) { // Check if its a vector
    std::cerr << "Expected a vector (N=1), but got N=" << N << std::endl;
    exit(1);
  }

  VecDoub vec(M); // Create vector with M elements
  for (int i = 0; i < M; i++) {
    file >> vec[i]; // Read each element
  }

  file.close();
  return vec;
}

void problemI(MatDoub A) {
  SVD svd(A);
  std::cout << "Diagonal elements of W:" << std::endl;
  util::print(svd.w);
}

void problemII(VecDoub b, MatDoub A) {

  VecDoub x(A.ncols());

  SVD svd(A);

  svd.solve(b, x);

  std::cout << "Best fit solution parameters x:" << std::endl;
  util::print(x);
}

void problemIII(VecDoub b, MatDoub A) {
  int m = A.nrows(), n = A.ncols();

  // The least squares solution x
  VecDoub x(n);
  SVD svd(A);
  svd.solve(b, x);

  // residual vector r = Ax - b
  VecDoub Ax = A * x;
  VecDoub r(m);
  for (int i = 0; i < m; i++) {
    r[i] = Ax[i] - b[i];
  }

  // The norm of the residual ||Ax - b||
  Doub norm_r = 0.0;
  for (int i = 0; i < m; i++) {
    norm_r += r[i] * r[i];
  }
  norm_r = std::sqrt(norm_r);

  // The norm of b ||b||
  Doub norm_b = 0.0;
  for (int i = 0; i < m; i++) {
    norm_b += b[i] * b[i];
  }
  norm_b = std::sqrt(norm_b);

  // The relative residual error
  Doub epsilon_residual =
      (norm_b != 0.0) ? (norm_r / norm_b) : 0.0; // Avoid division by zero
  std::cout << "Relative residual error epsilon_residual = " << epsilon_residual
            << std::endl;

  // The expected error for a random model
  Doub random_model_error = std::sqrt((m - n) / static_cast<Doub>(m));
  std::cout << "Expected error for a random model = " << random_model_error
            << std::endl;

  // Error on parameter estimates (std)
  VecDoub delta_x(n);
  for (int j = 0; j < n; j++) {
    Doub sum = 0.0;
    for (int i = 0; i < n; i++) {
      if (svd.w[i] > 1e-10) { // Avoid division by zero
        Doub v_ji = svd.v[j][i];
        sum += (v_ji * v_ji) / (svd.w[i] * svd.w[i]);
      }
    }
    delta_x[j] = std::sqrt(sum);
  }

  std::cout << "Accuracy: Standard deviations of x: " << std::endl;
  util::print(delta_x);
}

void problemIV(VecDoub b, MatDoub A) {
  int m = A.nrows(), n = A.ncols();

  VecDoub x(n);
  SVD svd(A);
  svd.solve(b, x);

  // Residual vector r = Ax - b
  VecDoub Ax = A * x;
  VecDoub r(m);
  for (int i = 0; i < m; i++) {
    r[i] = Ax[i] - b[i];
  }

  std::cout << "Residual vector:" << std::endl;
  util::print(r);
}

void problemV(VecDoub b, MatDoub A) {
  int m = A.nrows(), n = A.ncols();

  VecDoub x(n);
  SVD svd(A);
  svd.solve(b, x);

  // Residual vector r = Ax - b
  VecDoub Ax = A * x;
  VecDoub r(m);
  for (int i = 0; i < m; i++) {
    r[i] = Ax[i] - b[i];
  }

  // sigma_i
  VecDoub sigma(m);
  for (int i = 0; i < m; i++) {
    sigma[i] = std::max(1.0, std::abs(r[i]));
  }

  // Prime A'
  MatDoub A_prime(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      A_prime[i][j] = A[i][j] / sigma[i];
    }
  }

  // Prime b'
  VecDoub b_prime(m);
  for (int i = 0; i < m; i++) {
    b_prime[i] = b[i] / sigma[i];
  }

  std::cout << "New A[0][0] = " << A_prime[0][0] << std::endl;
  std::cout << "New b[6] = " << b_prime[6] << std::endl;
}

void problemVI(VecDoub b, MatDoub A) {
  int m = A.nrows(), n = A.ncols();

  // Compute the initial solution x using SVD
  VecDoub x(n);
  SVD svd(A);
  svd.solve(b, x);

  // Residual vector r = Ax - b
  VecDoub Ax = A * x;
  VecDoub r(m);
  for (int i = 0; i < m; i++) {
    r[i] = Ax[i] - b[i];
  }

  // sigma_i
  VecDoub sigma(m);
  for (int i = 0; i < m; i++) {
    sigma[i] = std::max(1.0, std::abs(r[i]));
  }

  // Prime A'
  MatDoub A_prime(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      A_prime[i][j] = A[i][j] / sigma[i];
    }
  }
  // Prime b'
  VecDoub b_prime(m);
  for (int i = 0; i < m; i++) {
    b_prime[i] = b[i] / sigma[i];
  }

  // Solve the new system A' x = b' using SVD
  SVD svd_prime(A_prime);
  VecDoub x_new(n);
  svd_prime.solve(b_prime, x_new);

  // Output new solution
  std::cout << "New solution x:" << std::endl;
  util::print(x_new);
}

int main(int argc, char *argv[]) {

  MatDoub A = readMatrix("/Users/patrickandersen/Desktop/6 semester/Numeriske "
                         "Metoder/Code/man1/Ex1A.dat");
  VecDoub b = readVector("/Users/patrickandersen/Desktop/6 semester/Numeriske "
                         "Metoder/Code/man1/Ex1b.dat");

  // Verify dimensions
  std::cout << "Matrix A: " << A.nrows() << " x " << A.ncols() << std::endl;
  std::cout << "Vector b: " << b.size() << std::endl;

  std::cout << "\nMatrix A:" << std::endl;
  util::print(A);
  std::cout << "\nVector b:" << std::endl;
  util::print(b);

  // I use \n to make space inbetween

  std::cout << "\nProblem I" << std::endl;
  problemI(A);
  
  std::cout << "\nProblem II" << std::endl;
  problemII(b, A);

  std::cout << "\nProblem III" << std::endl;
  problemIII(b, A);

  std::cout << "\nProblem IV" << std::endl;
  problemIV(b, A);

  std::cout << "\nProblem V" << std::endl;
  problemV(b, A);

  std::cout << "\nProblem VI" << std::endl;
  problemVI(b, A);

  return 0;
}
