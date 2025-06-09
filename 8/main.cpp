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
    double (*f)(double), double a_low, double b_high, double exp_order,
    double accuracy = 1e-10, int max_iterations = 10000) {
  std::vector<double> A_h;
  std::vector<int> f_computations;
  double error = std::numeric_limits<double>::max();
  int iterations = 2;
  while (abs(error) > accuracy && iterations < max_iterations) {
    A_h.push_back(method(f, a_low, b_high, iterations));
    f_computations.push_back(iterations);
    error = richError_current(A_h, pow(2, exp_order));
    iterations *= 2;
  }
  return {A_h, f_computations};
}

double func_1(double x) { return cos(x * x) * exp(-x); }

int main() {
  double a_low = 0.0;
  double b_high = 1.0;

  std::print("\n--------------------------------------------\n");
  std::print("Using Midpoint rule\n");

  auto mid_result = computeMethod(midPointRec, func_1, a_low, b_high, 2);
  std::vector<double> mid_A_h = mid_result.first;
  std::vector<double> mid_diff = A_h_diff(mid_A_h);
  std::vector<double> mid_alpha = alphaK(mid_A_h);
  std::vector<double> mid_rich_error = richError(mid_A_h, pow(2, 2));
  std::vector<double> mid_order_est = order_estimates(mid_alpha);
  std::vector<int> mid_f_computations = mid_result.second;

  std::print("\n{:>3} {:>12} {:>18} {:>12} {:>12} {:>12} {:>10}\n", "i",
             "A(h_i)", "A(h_{i-1}) - A(h_i)", "alpha^k", "Rich error",
             "Order est.", "f comps");

  for (size_t i = 0; i < mid_A_h.size(); ++i) {
    std::print("{:3} {:12.6f} {:18.6f} {:12.6f} {:12.6f} {:12.6f} {:10d}\n",
               i + 1, mid_A_h[i], mid_diff[i], mid_alpha[i], mid_rich_error[i],
               mid_order_est[i], mid_f_computations[i]);
  }

  std::print("\n--------------------------------------------\n");
  std::print("Using Simpson's rule\n");

  auto sim_result = computeMethod(simpson, func_1, a_low, b_high, 4);
  std::vector<double> sim_A_h = sim_result.first;
  std::vector<double> sim_diff = A_h_diff(sim_A_h);
  std::vector<double> sim_alpha = alphaK(sim_A_h);
  std::vector<double> sim_rich_error = richError(sim_A_h, pow(2, 2));
  std::vector<double> sim_order_est = order_estimates(sim_alpha);
  std::vector<int> sim_f_computations = sim_result.second;

  std::print("\n{:>3} {:>12} {:>18} {:>12} {:>12} {:>12} {:>10}\n", "i",
             "A(h_i)", "A(h_{i-1}) - A(h_i)", "alpha^k", "Rich error",
             "Order est.", "f comps");

  for (size_t i = 0; i < sim_A_h.size(); ++i) {
    std::print("{:3} {:12.6f} {:18.6f} {:12.6f} {:12.6f} {:12.6f} {:10d}\n",
               i + 1, sim_A_h[i], sim_diff[i], sim_alpha[i], sim_rich_error[i],
               sim_order_est[i], sim_f_computations[i]);
  }

  std::print("\n--------------------------------------------\n");
  std::print("Using Trapezoidal rule\n");

  auto trap_result = computeMethod(trapezoidal, func_1, a_low, b_high, 2);
  std::vector<double> trap_A_h = trap_result.first;
  std::vector<double> trap_diff = A_h_diff(trap_A_h);
  std::vector<double> trap_alpha = alphaK(trap_A_h);
  std::vector<double> trap_rich_error = richError(trap_A_h, pow(2, 2));
  std::vector<double> trap_order_est = order_estimates(trap_alpha);
  std::vector<int> trap_f_computations = trap_result.second;

  std::print("\n{:>3} {:>12} {:>18} {:>12} {:>12} {:>12} {:>10}\n", "i",
             "A(h_i)", "A(h_{i-1}) - A(h_i)", "alpha^k", "Rich error",
             "Order est.", "f comps");

  for (size_t i = 0; i < trap_A_h.size(); ++i) {
    std::print("{:3} {:12.6f} {:18.6f} {:12.6f} {:12.6f} {:12.6f} {:10d}\n",
               i + 1, trap_A_h[i], trap_diff[i], trap_alpha[i],
               trap_rich_error[i], trap_order_est[i], trap_f_computations[i]);
  }

  return 0;
}
