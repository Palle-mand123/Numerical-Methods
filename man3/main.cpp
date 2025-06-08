#include "ludcmp.h"
#include "nr3.h"
#include "quadrature.h"
#include "utilities.h"
#include <cmath>
#include <derule.h>
#include <functional>
#include <iostream>
#include <vector>

struct equation {
  double operator()(double x) { return (cos(pow(x, 3)) * exp(-x)) / sqrt(x); }
};

// Midpoint integration method
double midPoint(double a, double b, int n,
                std::function<double(double)> integral) {
  double h = (b - a) / n;
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    double x_mid = a + (i + 0.5) * h;
    sum += integral(x_mid);
  }
  return h * sum;
}

void problemII(double a, double b, std::function<double(double)> integral,
               int iterations, int order, double tolerance = pow(10, -3)) {

  std::vector<std::vector<double>> A(iterations + 1,
                                     std::vector<double>(iterations + 1, 0.0));
  vector<double> differences(iterations + 1, 0.0);

  for (int i = 1; i <= iterations; i++) {
    int n = pow(2, i) + 1;
    A[i][1] = midPoint(a, b, n, integral);
  }

  cout << "\n" << endl;
  cout << setw(5) << "i" << setw(15) << "A(h_i)" << setw(25)
       << "A(h_(i-1)) - A(h_i)" << setw(10) << "alp^k" << setw(17)
       << "Rich-error" << setw(20) << "f-calculations" << setw(23)
       << "order-estimate \n"
       << endl;

  // First row
  cout << setw(5) << 1 << setw(16) << fixed << setprecision(6) << A[1][1]
       << endl;

  // Second row
  differences[2] = A[1][1] - A[2][1];
  cout << setw(5) << 2 << setw(16) << fixed << setprecision(6) << A[2][1]
       << setw(19) << fixed << setprecision(6) << differences[2] << endl;

  // Remaining rows
  for (int i = 3; i <= iterations; i++) {
    // Differences for current row
    differences[i] = A[i - 1][1] - A[i][1];

    // (alp^k)
    double alpha = 0.0;
    if (fabs(differences[i - 1]) > 1e-20) {
      alpha = differences[i - 1] / differences[i];
    }

    // Richardson error estimate
    double richError = 0.0;
    if (i > 2) {
      richError = (A[i][1] - A[i - 1][1]) / (pow(2, order) - 1);
    }

    // Order estimate
    double order_estimate = 0.0;
    if (i > 3 && fabs(differences[i - 1]) > 1e-20) {
      order_estimate =
          log(fabs(differences[i - 1] / differences[i])) / log(2.0);
    }

    // Function calculations
    int f_cal = pow(2, i - 1);

    cout << setw(5) << i << setw(16) << fixed << setprecision(6) << A[i][1]
         << setw(19) << fixed << setprecision(6) << differences[i] << setw(16)
         << fixed << setprecision(6) << alpha << setw(15) << fixed
         << setprecision(6) << richError << setw(15) << f_cal << setw(24)
         << fixed << setprecision(6) << order_estimate << endl;

    if (fabs(richError) < tolerance && i > 3) {
      cout << "\nConverged to desired accuracy." << endl;
      break;
    }
  }
}

struct DErule_func {
  int callCount = 0;

  double operator()(double x, double dx) {
    callCount++;
    return ((cos(pow(x, 3)) * exp(-x)) / sqrt(x));
  }

  void resetCount() { callCount = 0; }
  int getCount() const { return callCount; }
};

void problemIII(int a, int b, int iterations, double tolerance = pow(10, -3)) {
  DErule_func func;
  DErule<DErule_func> de(func, a, b);
  double prev = 0.0;
  cout << setw(5) << "\nApproximations" << setw(20) << "f-calculations"
       << "\n"
       << endl;
  for (int i = 0; i < iterations; i++) {
    func.resetCount();
    double current = de.next();
    cout << setw(11) << fixed << setprecision(6) << de.next() << setw(17)
         << func.getCount() << endl;
    if (i > 0) {
      if (fabs(current - prev) < tolerance) {
        cout << "\nconvereged towards final result"
             << "\n"
             << endl;
        break;
      }
    }
    prev = current;
  }
}

int main() {
  equation integral;
  double a = 0.0;
  double b = 4.0;
  int order = 2;
  int iterations = 16;

  problemII(a, b, integral, iterations, order);

  problemIII(a, b, iterations);

  return 0;
}