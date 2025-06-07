#include "nr3.h"
#include "roots.h"
#include <cmath>
#include <print>

struct equation {
  double operator()(double x) const { return x - std::cos(x); }
};

// For Newton's method, you need both function and derivative
struct equationWithDerivative {
  double operator()(double x) const { return x - std::cos(x); }
  double df(double x) const { return 1.0 + std::sin(x); }
};

void solveBisection(double a, double b, double tol, equation &f) {
  double prev_mid = 0.0;
  double dx = b - a;
  int k = 1;

  std::print("Bisection\n");
  std::print("k  {:>10}  {:>10}  {:>10}\n", "xmin", "xmax", "dx");

  // Bisection loop
  while (true) {
    double mid = 0.5 * (a + b);
    dx = b - a;

    std::print("{:2}  {:10.6f}  {:10.6f}  {:10.6f}\n", k, a, b, dx);

    double fmid = f(mid);
    if (dx <= tol || std::fabs(fmid) <= tol) {
      prev_mid = mid;
      break;
    }

    if (f(a) * fmid < 0) {
      b = mid;
    } else {
      a = mid;
    }

    prev_mid = mid;
    k++;
  }

  std::print("Result: {:.6f}\n", prev_mid);
}

void numericalRecipesBisection(double a, double b, double tol, equation &f) {
  try {
    double root = rtbis(f, a, b, tol);
    std::print("Root found: {:.8f}\n", root);
    std::print("Verification: f({:.8f}) = {:.8f}\n", root, f(root));
  } catch (const char *error) {
    std::print("Error: {}\n", error);
  }
}

void solveNewton(double x0, double tol, equationWithDerivative &f,
                 int maxIter = 100) {
  double x = x0;

  std::print("Newton\n");
  std::print("{:>5} {:>15} {:>15}\n", "k", "x", "dx");

  for (int k = 1; k <= maxIter; k++) {
    double fx = f(x);
    double dfx = f.df(x);

    if (std::abs(dfx) < 1e-15) {
      std::print("Error: Derivative too small at iteration {}\n", k);
      return;
    }

    double dx = -fx / dfx; // Newton step: dx = -f(x)/f'(x)

    std::print("{:5} {:15.6f} {:15.6g}\n", k, x, dx);

    if (std::abs(dx) <= tol || std::abs(fx) <= tol) {
      x += dx;
      std::print("Result: {:.6f}\n", x);
      return;
    }

    x += dx;
  }

  std::print("Newton method failed to converge after {} iterations\n", maxIter);
}

void numericalRecipesNewton(double x1, double x2, double tol,
                            equationWithDerivative &f) {
  try {
    double root = rtnewt(f, x1, x2, tol);
    std::print("Newton Root found: {:.8f}\n", root);
    std::print("Newton Verification: f({:.8f}) = {:.8f}\n", root, f(root));
  } catch (const char *error) {
    std::print("Newton Error: {}\n", error);
  }
}

void solveSecant(double x0, double x1, double tol, equation &f,
                 int maxIter = 100) {
  double x_prev = x0;
  double x_curr = x1;

  std::print("Secant Method\n");
  std::print("{:>5} {:>15} {:>15}\n", "k", "x", "dx");

  for (int k = 1; k <= maxIter; k++) {
    double f_prev = f(x_prev);
    double f_curr = f(x_curr);

    if (std::abs(f_curr - f_prev) < 1e-15) {
      std::print("Error: Function values too close at iteration {}\n", k);
      return;
    }

    double x_new = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);
    double dx = x_new - x_curr;

    std::print("{:5} {:15.6f} {:15.6g}\n", k, x_curr, std::abs(dx));

    if (std::abs(dx) <= tol || std::abs(f_curr) <= tol) {
      std::print("Result: {:.6f}\n", x_new);
      return;
    }

    x_prev = x_curr;
    x_curr = x_new;
  }

  std::print("Secant method failed to converge after {} iterations\n", maxIter);
}

void numericalRecipesSecant(double x1, double x2, double tol, equation &f) {
  try {
    double root = rtsec(f, x1, x2, tol);
    std::print("Secant Root found: {:.8f}\n", root);
    std::print("Secant Verification: f({:.8f}) = {:.8f}\n", root, f(root));
  } catch (const char *error) {
    std::print("Secant Error: {}\n", error);
  }
}

int main() {
  equation f;
  equationWithDerivative fWithDerivative;

  double a = 0.0;
  double b = M_PI / 2;
  double tol = 1e-8;

  std::print("=== Bisection Method ===\n");
  solveBisection(a, b, tol, f);

  std::print("\n=== Numerical Recipes Bisection (library) ===\n");
  numericalRecipesBisection(a, b, tol, f);

  std::print("\n---------------------------------------------------------------"
             "----\n");

  std::print("\n=== Newton-Raphson Method (with table) ===\n");
  solveNewton(a, tol, fWithDerivative);

  // Your custom Newton is "pure" and only needs one guess, while the library
  // Newton is "safe" and needs bracketing for robustness. The library version
  // is more like a combination of Newton + Bisection methods.

  std::print("\n=== Newton-Raphson Method (library) ===\n");
  numericalRecipesNewton(a, b, tol, fWithDerivative);

  std::print("\n---------------------------------------------------------------"
             "----\n");

  std::print("\n=== Secant Method (with table) ===\n");
  solveSecant(a, b, tol, f);

  std::print("\n=== Secant Method (library) ===\n");
  numericalRecipesSecant(a, b, tol, f);

  return 0;
}
