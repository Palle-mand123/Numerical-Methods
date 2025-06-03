#include "ludcmp.h"
#include "nr3.h"
#include "quadrature.h"
#include "utilities.h"
#include <cmath>
#include <iostream>
#include <print>
#include <vector>

struct Func1 {
  double operator()(double x) { return cos(pow(x, 2)) * exp(-x); }
};

void midPointRec(const double a, const double b, const int n) {
  Func1 f;
  Midpnt<Func1> integrator(f, a, b);
  std::vector<Doub> A_h;
  double order_esimate = 0.0;
    // Print header
  std::cout << std::left << std::setw(10) << "i" << std::setw(15) << "A(h_i)"
            << std::setw(25) << "A(h_(i-1)) - A(h_i)" << std::setw(20)
            << "alp^k" << std::setw(20) << "Rich-error" << std::setw(20)
            << "f-calculations" << std::setw(15) << "order-estimate"
            << std::endl;

  for (int i = 1; i <= n; i++) {
    Doub A_hi = integrator.next();
    A_h.push_back(A_hi);

    std::cout << std::left << std::setw(10) << i;
    std::cout << std::fixed << std::setprecision(6) << std::setw(15)
              << A_h[i - 1];

    if (i >= 2) {
      std::cout << std::setw(25) << (A_h[i - 2] - A_h[i - 1]);
    } else {
      std::cout << std::setw(25) << "-";
    }

    if (i >= 3) {
      Doub alpha = (A_h[i - 3] - A_h[i - 2]) / (A_h[i - 2] - A_h[i - 1]);
      std::cout << std::setw(20) << alpha;
      Doub error = A_h[i - 1] + ((A_h[i - 1] - A_h[i - 2]) / (pow(2,2)-1));
      std::cout << std::setw(20)
                << error;
      int f_calculations = 4 << (i - 3);
      std::cout << std::setw(20) << f_calculations;
      if(i>3)
      {
      order_esimate = log(fabs((A_h[i - 2]-A_h[i - 3])/(A_h[i - 1]-A_h[i-2]))/log(2.0));
      }
      std::cout << std::setw(20) << order_esimate;
    } else {
      std::cout << std::setw(20) << "-";
      std::cout << std::setw(20) << "-";
      std::cout << std::setw(20) << "-";
      std::cout << std::setw(20) << "-";

    }

    std::cout << std::endl;
  }
}

/*
void midPointRec(const double a, const double b, const int n) {
  Func1 f;
  std::vector<Doub> A_h;
  Midpnt<Func1> integrator(f,a,b);
  double result;

      // Print header
    std::cout << std::left
              << std::setw(10) << "i"
              << std::setw(15) << "A(h_i)"
              << std::setw(30) << "A(h_(i-1)) - A(h_i)"
              << std::setw(15) << "alp^k"
              << std::setw(15) << "Rich-error"
              << std::setw(20) << "f-calculations"
              << std::setw(15) << "order-estimate"
              << std::endl;

  for (int i = 1; i <= n; i++) {
    result = integrator.next();
    A_h.push_back(result);
  }

  for (auto i : A_h) {
    std::cout << i << std::endl;
  }
}
*/

void trapezoidal(const double a, const double b, const int n) {
  Func1 f;
  Trapzd<Func1> integrator(f, a, b);
  double result;

  for (int i = 1; i <= n; i++) {
    result = integrator.next();
    std::cout << "step " << i << ": " << result << std::endl;
  }
}

void simpson(const double a, const double b, const int n) {
  Func1 f;
  double result = qsimp(f, a, b);

  std::cout << "result:" << result << std::endl;
}

int main() {

  std::cout << "\n"
            << "Midpoint" << std::endl;
  const double a = 0.0;
  const double b = 1.0;
  const int n = 16;

  midPointRec(a, b, n);
  std::cout << "\n"
            << "Trapezoidal" << std::endl;
  trapezoidal(a, b, n);

  std::cout << "\n"
            << "Simpson" << std::endl;
  simpson(a, b, n);

  return 0;
}
