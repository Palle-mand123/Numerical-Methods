#include "cholesky.h"
#include "ludcmp.h"
#include "nr3.h"
#include "utilities.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

int main() {

  ifstream file("/Users/patrickandersen/Desktop/6 semester/Numeriske "
                "Metoder/Code/2/PontiusData.dat");
  if (!file.is_open()) {
    cerr << "Could not open data file." << endl;
    return 1;
  }

  vector<double> x_vals, y_vals;
  double x, y;
  while (file >> y >> x) {
    x_vals.push_back(x);
    y_vals.push_back(y);
  }
  file.close();

  int N = x_vals.size();
  const int M = 3; // a0, a1, a2

  MatDoub A(N, M);
  VecDoub b(N);

  for (int i = 0; i < N; ++i) {
    A[i][0] = 1.0;
    A[i][1] = x_vals[i];
    A[i][2] = x_vals[i] * x_vals[i];
    b[i] = y_vals[i];
  }

  cout << "A Matrix:" << endl;
  util::print(A); // A matrix

  cout << "\nB Vector:" << endl;
  util::print(b); // B vector

  //  At*A and At*b
  MatDoub AtA(M, M, 0.0);
  VecDoub Atb(M, 0.0);
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      for (int k = 0; k < N; ++k) {
        AtA[i][j] += A[k][i] * A[k][j];
      }
    }
    for (int k = 0; k < N; ++k) {
      Atb[i] += A[k][i] * b[k];
    }
  }

  cout << "\n A_T * b vector:" << endl;
  util::print(Atb); // A^T*b

  cout << "\n A_T * A Matrix:" << endl;
  util::print(AtA); // A^T*A Matrix

  VecDoub x_lu(M);
  LUdcmp lu_solver(AtA);
  lu_solver.solve(Atb, x_lu);

  cout << "LU decomposition:" << endl;
  util::print(x_lu);

  VecDoub x_chol(M);
  Cholesky chol_solver(AtA);
  chol_solver.solve(Atb, x_chol);

  cout << "\nCholesky decomposition:" << endl;
  util::print(x_chol);

  return 0;
}
