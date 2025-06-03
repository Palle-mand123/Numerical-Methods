#include "cholesky.h"
#include "ludcmp.h"
#include "nr3.h"
#include "utilities.h"
#include <fstream>
#include <print>

void pontSolution(VecDoub x, VecDoub y, int parameters) {

  int size = x.size();
  MatDoub A(size, parameters);
  VecDoub b(size);
  int k = 1;

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < parameters; j++) {
      A[i][j] = pow(x[i], j) / k;
    }
    b[i] = y[i] / k;
  }

  MatDoub A_T = util::Transpose(A);

  // Lektion 2 slide 10
  MatDoub C = A_T * A;
  VecDoub c = A_T * b;
  util::print(A);

  // use to find a: Solve C * a = c

  // LU decomposition for at finde best values: a0, a1, a2
  auto ludcmp_solver = LUdcmp(C);
  VecDoub a_ludcmp(parameters);
  ludcmp_solver.solve(c, a_ludcmp);
  util::print(a_ludcmp, "LU decomposition");

  std::cout << " " << std::endl;

  // Cholesky for at finde best values: a0, a1, a2
  VecDoub a_cholesky(parameters);
  auto cholesky_solver = Cholesky(C);
  cholesky_solver.solve(c, a_cholesky);
  util::print(a_cholesky, "Cholesky");

  std::cout << " " << std::endl;
}

int main() {

  // Load PontiusData.dat
  int pontDataLength = 40;
  VecDoub pontX(pontDataLength);
  VecDoub pontY(pontDataLength);

  ifstream readPont("/Users/patrickandersen/Desktop/6 semester/Numeriske "
                    "Metoder/Code/2/PontiusData.dat");

  if (!readPont) {
    cerr << "Cant Open PontiusData.dat" << endl;
    return 1;
  }

  for (int i = 0; i < pontDataLength; i++) {
    readPont >> pontX[i];
    readPont >> pontY[i];
  }

  pontSolution(pontX, pontY, 3);

  return 0;
}
