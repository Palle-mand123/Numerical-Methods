#include "ludcmp.h"
#include "nr3.h"
#include "quadrature.h"
#include "utilities.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <print>
#include <vector>

double midPointRec(double (*f)(double), double a_low, double b_high,
                   int iterations = 100) {

  double n = (b_high - a_low) / iterations;
  double results = 0.0;

  for (int i = 0; i < iterations; i++) {
    double x1 = a_low + i * n;
    double x2 = x1 + n;
    results += f((x1 + x2) / 2) * n;
  }
  return results;
}

double trapezoidal(double (*f)(double), const double a_low, const double b_high,
                   int iterations = 100) {
  double n = (b_high - a_low) / iterations;
  double results = 0.5 * f(a_low);

  for (int i = 1; i < iterations; i++) {
    double x = a_low + i * n;
    results += f(x);
  }

  results += 0.5 * f(b_high);
  results *= n;
  return results;
}

double simpson(double (*f)(double), const double a_low, const double b_high,
               int iterations = 100) {

  if (iterations % 2 != 0) {
    std::cerr << "Error: subIntervals must be even for Simpson's rule."
              << std::endl;
    ++iterations;
  }

  double n = (b_high - a_low) / iterations;
  double results = f(a_low) + f(b_high);

  for (int i = 1; i < iterations; i++) {
    double x = a_low + i * n;
    results += f(x) * (i % 2 == 0 ? 2 : 4);
  }

  results *= n / 3;
  return results;
}

std::vector<double> A_h_diff(std::vector<double> &A_h) {
  std::vector<double> diff = {NAN};
  for (size_t i = 1; i < A_h.size(); i++) {
    diff.push_back(A_h[i - 1] - A_h[i]);
  }
  return diff;
}

std::vector<double> alphaK(std::vector<double> &A_h) {
  std::vector<double> alphaK = {NAN, NAN};
  for (size_t i = 2; i < A_h.size(); ++i) {
    double A1 = A_h[i - 2];
    double A2 = A_h[i - 1];
    double A3 = A_h[i];
    double A_i = (A1 - A2) / (A2 - A3);
    alphaK.push_back(A_i);
  }
  return alphaK;
}

std::vector<double> richError(std::vector<double> &A_h,
                              const double alphaK_exp) {
  std::vector<double> richError = {NAN};
  for (size_t i = 1; i < A_h.size(); ++i) {
    double A1 = A_h[i - 1];
    double A2 = A_h[i];
    double R_i = (A2 - A1) / (alphaK_exp - 1.0);
    richError.push_back(R_i);
  }
  return richError;
}

double richError_current(const std::vector<double> &A_h,
                         const double alphaK_exp) {
  double error = std::numeric_limits<double>::max();
  if (A_h.size() < 2) {
    return error;
  } else {
    const double A1 = A_h[A_h.size() - 2];
    const double A2 = A_h[A_h.size() - 1];

    error = (A2 - A1) / (alphaK_exp - 1);
    return error;
  }
}

std::vector<double> order_estimates(const std::vector<double> &alphaK) {
  std::vector<double> order_est = {NAN, NAN};
  for (size_t i = 2; i < alphaK.size(); ++i) {
    double order = log2(abs(alphaK[i]));
    order_est.push_back(order);
  }
  return order_est;
}

std::pair<std::vector<double>, std::vector<int>> computeMethod(
    std::function<double(double (*f)(double), double, double, int)> method,
    double (*f)(double), double a_low, double b_high, double accuracy = 1e-10,
    int max_iterations = 10000) {
  std::vector<double> A_h;
  std::vector<int> f_computations;
  double exp_order = 0.0;
  double error = std::numeric_limits<double>::max();
  int iterations = 2;
  while (abs(error) > accuracy && iterations < max_iterations) {
    A_h.push_back(method(f, a_low, b_high, iterations));
    f_computations.push_back(iterations);
    exp_order = 2.0;
    error = richError_current(A_h, pow(2, exp_order));
    iterations *= 2;
  }
  return {A_h, f_computations};
}

double func_1(double x) { return cos(x * x) * exp(-x); }

int main() {
  double a_low = 0.0;
  double b_high = 1.0;

  auto result = computeMethod(midPointRec, func_1, a_low, b_high);
  std::vector<double> A_h = result.first;
  std::vector<double> diff = A_h_diff(A_h);
  std::vector<double> alpha = alphaK(A_h);
  std::vector<double> rich_error = richError(A_h, pow(2, 2));
  std::vector<double> order_est = order_estimates(alpha);
  std::vector<int> f_computations = result.second;

  std::print("\n{:>3} {:>12} {:>18} {:>12} {:>12} {:>12} {:>10}\n", "i",
             "A(h_i)", "A(h_{i-1}) - A(h_i)", "alpha^k", "Rich error",
             "Order est.", "f comps");

  for (size_t i = 0; i < A_h.size(); ++i) {
    std::print("{:3} {:12.6f} {:18.6f} {:12.6f} {:12.6f} {:12.6f} {:10d}\n",
               i + 1, A_h[i], diff[i], alpha[i], rich_error[i], order_est[i],
               f_computations[i]);
  }

  return 0;
}
