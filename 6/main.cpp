#include "nr3.h"
#include "roots.h"
#include <cmath>
#include <cstdlib>
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

  std::print("\n k  {:>10}  {:>10}  {:>10}\n", "xmin", "xmax", "dx");

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

  std::print("\n{:>5} {:>15} {:>15}\n", "k", "x", "dx");

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

void solveSecant(double x0, double tol, equation &f, int maxIter = 100) {
  double x_prev = x0;
  double x_curr = 0.611015; // random initial guess

  std::print("\n{:>5} {:>15} {:>15}\n", "k", "x", "dx");

  for (int k = 1; k <= maxIter; k++) {
    double f_prev = f(x_prev);
    double f_curr = f(x_curr);

    if (std::abs(f_curr - f_prev) < 1e-15) {
      std::print("Error: Function values too close at iteration {}\n", k);
      return;
    }

    double x_new = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);
    double dx = x_curr - x_prev;

    std::print("{:5} {:15.6f} {:15.6g}\n", k, x_prev, std::abs(dx));

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

void solveFalsePosition(double a, double b, double tol, equation &f,
                        int maxIter = 100) {
  double fa = f(a);
  double fb = f(b);

  if (fa * fb > 0.0) {
    std::print("Error: Root must be bracketed for False Position\n");
    return;
  }

  std::print("\n k {:>10}  {:>10}  {:>10}\n", "xmin", "xmax", "dx");

  for (int k = 1; k <= maxIter; k++) {

    // Correct False position formula: c = a - fa * (b - a) / (fb - fa)
    double c = a - fa * (b - a) / (fb - fa);
    double fc = f(c);
    double dx = b - a;

    if (fa * fc < 0.0) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }

    std::print("{:2}  {:10.6f}  {:10.6f}  {:10.6f}\n", k, a, b, dx);

    if (std::abs(fc) <= tol || dx <= tol) {
      std::print("Result: {:.6f}\n", c);
      return;
    }
  }

  std::print("False Position method failed to converge after {} iterations\n",
             maxIter);
}

void numericalRecipesFalsePosition(double x1, double x2, double tol,
                                   equation &f) {
  try {
    double root = rtflsp(f, x1, x2, tol);
    std::print("False Position Root found: {:.8f}\n", root);
    std::print("False Position Verification: f({:.8f}) = {:.8f}\n", root,
               f(root));
  } catch (const char *error) {
    std::print("False Position Error: {}\n", error);
  }
}

void solveRidders(double a, double b, double tol, equation &f,
                  int maxIter = 100) {
  double fa = f(a);
  double fb = f(b);

  if (fa * fb > 0.0) {
    std::print("Error: Root must be bracketed for Ridders method\n");
    return;
  }

  std::print("\n{:>5} {:>15} {:>15}\n", "k", "x", "dx");

  double xl = a, xh = b;
  double fl = fa, fh = fb;
  double ans = -9.99e99;

  for (int k = 1; k <= maxIter; k++) {
    double xm = 0.5 * (xl + xh);
    double fm = f(xm);

    double s = std::sqrt(fm * fm - fl * fh);
    if (s == 0.0) {
      std::print("{:5} {:15.6f} {:15.6g}\n", k, ans, 0.0);
      std::print("Result: {:.6f}\n", ans);
      return;
    }

    double xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);

    double dx = (k == 1) ? std::abs(xnew - xm) : std::abs(xnew - ans);

    std::print("{:5} {:15.6f} {:15.6g}\n", k, xnew, dx);

    if (k > 1 && std::abs(xnew - ans) <= tol) {
      std::print("Result: {:.6f}\n", xnew);
      return;
    }

    ans = xnew;
    double fnew = f(ans);

    if (std::abs(fnew) <= tol) {
      std::print("Result: {:.6f}\n", ans);
      return;
    }

    if (fm * fnew < 0.0) {
      xl = xm;
      fl = fm;
      xh = ans;
      fh = fnew;
    } else if (fl * fnew < 0.0) {
      xh = ans;
      fh = fnew;
    } else if (fh * fnew < 0.0) {
      xl = ans;
      fl = fnew;
    }

    if (std::abs(xh - xl) <= tol) {
      std::print("Result: {:.6f}\n", ans);
      return;
    }
  }

  std::print("Ridders method failed to converge after {} iterations\n",
             maxIter);
}

void numericalRecipesRidders(double x1, double x2, double tol, equation &f) {
  try {
    double root = zriddr(f, x1, x2, tol);
    std::print("Ridders Root found: {:.8f}\n", root);
    std::print("Ridders Verification: f({:.8f}) = {:.8f}\n", root, f(root));
  } catch (const char *error) {
    std::print("Ridders Error: {}\n", error);
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

  std::print("\n=== Newton-Raphson Method (library) ===\n");
  numericalRecipesNewton(a, b, tol, fWithDerivative);

  std::print("\n---------------------------------------------------------------"
             "----\n");

  std::print("\n=== Secant Method (with table) ===\n");
  solveSecant(a, tol, f);

  std::print("\n=== Secant Method (library) ===\n");
  numericalRecipesSecant(a, b, tol, f);

  std::print("\n---------------------------------------------------------------"
             "----\n");

  std::print("\n=== False Position Method (with table) ===\n");
  solveFalsePosition(a, b, tol, f);

  std::print("\n=== False Position Method (library) ===\n");
  numericalRecipesFalsePosition(a, b, tol, f);

  std::print("\n---------------------------------------------------------------"
             "----\n");

  std::print("\n=== Ridders Method (with table) ===\n");
  solveRidders(a, b, tol, f);

  std::print("\n=== Ridders Method (library) ===\n");
  numericalRecipesRidders(a, b, tol, f);

  return 0;
}
