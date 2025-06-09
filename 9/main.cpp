
#include "nr3.h"
#include "quadrature.h"
#include <iostream>
#include <vector>

template <class T> struct DErule : Quadrature {
  Doub a, b, hmax, s;
  T &func;

  DErule(T &funcc, const Doub aa, const Doub bb, const Doub hmaxx = 3.7)
      : func(funcc), a(aa), b(bb), hmax(hmaxx) {
    n = 0;
  }

  Doub next() {
    Doub del, fact, q, sum, t, twoh;
    Int it, j;
    n++;

    if (n == 1) {
      fact = 0.25;
      s = hmax * 2.0 * (b - a) * fact * func(0.5 * (b + a), 0.5 * (b - a));
      return s;

    } else {
      for (it = 1, j = 1; j < n - 1; j++)
        it <<= 1;
      twoh = hmax / it;
      t = 0.5 * twoh;
      for (sum = 0.0, j = 0; j < it; j++) {
        q = exp(-2.0 * sinh(t));
        del = (b - a) * q / (1.0 + q);
        fact = q / SQR(1.0 + q) * cosh(t);
        sum += fact * (func(a + del, del) + func(b - del, del));
        t += twoh;
      }
      s = 0.5 * s + (b - a) * twoh * sum;
      return s;
    }
  }
};

double func_1(double x, double delta) {
  (void)delta;

  return cos(x * x) * exp(-x);
}

std::vector<double> A_h_diff(std::vector<double> &s_values) {
  std::vector<double> diff = {NAN};
  for (size_t i = 1; i < s_values.size(); i++) {
    diff.push_back(s_values[i - 1] - s_values[i]);
  }
  return diff;
}

int main() {
  auto DE_rule = DErule(func_1, 0.0, 1.0);
  double diff = std::numeric_limits<double>::max();
  std::vector<double> s_values;
  double prev = 0;
  while (abs(diff) > 1e-10) {
    double cur = DE_rule.next();
    s_values.push_back(cur);
    diff = cur - prev;
    prev = cur;
  }
  std::vector<double> diff_values = A_h_diff(s_values);

  std::print("\n i {:>10}  {:>10}\n", "A(h_i)", "A(h(i-1)) - A(h_i)");

  for (size_t i = 0; i < s_values.size(); ++i) {
    std::print("{:2}  {:10.6f}  {:10.6f}\n", i + 1, s_values[i],
               diff_values[i]);
  }
}
