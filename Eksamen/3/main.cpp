#include "nr3.h"
#include "utilities.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <print>
#include <vector>

double amax = 4.0;
double vdes = 25.0;
double D0 = 50.0;
double T = 1.5;
double acom = 2.0;
double t0 = -10.0;
static size_t f_comps_current = 0;

double XF(double t) { return 250 + 15 * t - 5 * sqrt(1 + pow(t, 2)); }
double XF_prime(double t) { return 15 - (5 * t) / sqrt(1 + pow(t, 2)); }

VecDoub derivs(const Doub t, VecDoub_I &x) {
  VecDoub_O dxdt(2);
  dxdt[0] = x[1];
  dxdt[1] = amax * (1 - pow((x[1] / vdes), 4) -
                    pow((D0 + MAX(0, x[1] * T + (x[1] * (x[1] - XF_prime(t))) /
                                                    (2 * acom))) /
                            (XF(t) - x[0]),
                        2));

  return dxdt;
}

VecDoub second_order_runge_kuttea_method(double low, double high, int steps,
                                         VecDoub_I &x,
                                         VecDoub derivs(const Doub x,
                                                        VecDoub_I &y)) {
  const double h = (high - low) / (double)steps;
  VecDoub x_n = x;
  for (double t_n = low; t_n < high; t_n += h) {
    f_comps_current += 2;
    auto k1 = h * derivs(t_n, x_n);
    auto k2 = h * derivs(t_n + 0.5 * h, x_n + 0.5 * k1);
    auto y_n_next = x_n + k2;
    x_n = y_n_next;
  }
  return x_n;
}

int main() {
  std::print("\n--------------------Problem II------------------------\n");
  double vals[2] = {0.0, 15.0};
  VecDoub x1(2, vals);
  VecDoub result = derivs(t0, x1);

  util::print(result);

  std::print("\n--------------------Problem III------------------------\n");

  VecDoub mid_result =
      second_order_runge_kuttea_method(-10, 10.0, 80, x1, derivs);

  util::print(mid_result);

  return 0;
}