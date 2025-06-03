#include "ludcmp.h"
#include "nr3.h"
#include "utilities.h"
#include <print>


int main() {
  // Exercise 1:
  // Solve A x = b using LU decomposition, and print the result.

  MatDoub A(3, 3);
  A[0][0] = 1.0;
  A[0][1] = 2.0;
  A[0][2] = 3.0;
  A[1][0] = 2.0;
  A[1][1] = -4.0;
  A[1][2] = 6.0;
  A[2][0] = 3.0;
  A[2][1] = -9.0;
  A[2][2] = -3.0;

  // Another way to implement A:
  double vals[] = {1.0, 2.0, 3.0, 2.0, -4.0, 6.0, 3.0, -9.0, -3.0};
  MatDoub A_new(3, 3, vals);

  VecDoub b(3);
  b[0] = 5.0;
  b[1] = 18.0;
  b[2] = 6.0;

  VecDoub x(3);

  // evaluate x
  auto LUdcmp_solver = LUdcmp(A);
  LUdcmp_solver.solve(b, x);

  MatDoub L(3, 3);
  MatDoub U(3, 3);
  LUdcmp_solver.decompose(L, U);

  // print x
  println("A:");
  util::print(A);

  println("b:");
  util::print(b);

  println("x:");
  util::print(x);

  println("Lower:");
  util::print(L);

  println("Upper:");
  util::print(U);

  LUdcmp_solver.mprove(b, x);
  println("Verify x: ");
  util::print(x);

  return 0;
}
