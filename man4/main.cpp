#include "nr3.h"
#include "utilities.h"
#include <cmath>
#include <iomanip>
#include <iostream>

#include "ludcmp.h"
#include "qrdcmp.h"

template <class T>
void lnsrch(VecDoub_I &xold, const Doub fold, VecDoub_I &g, VecDoub_IO &p,
            VecDoub_O &x, Doub &f, const Doub stpmax, Bool &check, T &func) {
  const Doub ALF = 1.0e-4, TOLX = numeric_limits<Doub>::epsilon();
  Doub a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
  Doub rhs1, rhs2, slope = 0.0, sum = 0.0, temp, test, tmplam;
  Int i, n = xold.size();
  check = false;
  for (i = 0; i < n; i++)
    sum += p[i] * p[i];
  sum = sqrt(sum);
  if (sum > stpmax)
    for (i = 0; i < n; i++)
      p[i] *= stpmax / sum;
  for (i = 0; i < n; i++)
    slope += g[i] * p[i];
  if (slope >= 0.0)
    throw("Roundoff problem in lnsrch.");
  test = 0.0;
  for (i = 0; i < n; i++) {
    temp = abs(p[i]) / MAX(abs(xold[i]), 1.0);
    if (temp > test)
      test = temp;
  }
  alamin = TOLX / test;
  alam = 1.0;
  for (;;) {
    for (i = 0; i < n; i++)
      x[i] = xold[i] + alam * p[i];
    f = func(x);
    if (alam < alamin) {
      for (i = 0; i < n; i++)
        x[i] = xold[i];
      check = true;
      return;
    } else if (f <= fold + ALF * alam * slope)
      return;
    else {
      if (alam == 1.0)
        tmplam = -slope / (2.0 * (f - fold - slope));
      else {
        rhs1 = f - fold - alam * slope;
        rhs2 = f2 - fold - alam2 * slope;
        a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
        b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) /
            (alam - alam2);
        if (a == 0.0)
          tmplam = -slope / (2.0 * b);
        else {
          disc = b * b - 3.0 * a * slope;
          if (disc < 0.0)
            tmplam = 0.5 * alam;
          else if (b <= 0.0)
            tmplam = (-b + sqrt(disc)) / (3.0 * a);
          else
            tmplam = -slope / (b + sqrt(disc));
        }
        if (tmplam > 0.5 * alam)
          tmplam = 0.5 * alam;
      }
    }
    alam2 = alam;
    f2 = f;
    alam = MAX(tmplam, 0.1 * alam);
  }
}
template <class T> struct NRfdjac {
  const Doub EPS;
  T &func;
  NRfdjac(T &funcc) : EPS(1.0e-8), func(funcc) {}
  MatDoub operator()(VecDoub_I &x, VecDoub_I &fvec) {
    Int n = x.size();
    MatDoub df(n, n);
    VecDoub xh = x;
    for (Int j = 0; j < n; j++) {
      Doub temp = xh[j];
      Doub h = EPS * abs(temp);
      if (h == 0.0)
        h = EPS;
      xh[j] = temp + h;
      h = xh[j] - temp;
      VecDoub f = func(xh);
      xh[j] = temp;
      for (Int i = 0; i < n; i++)
        df[i][j] = (f[i] - fvec[i]) / h;
    }
    return df;
  }
};
template <class T> struct NRfmin {
  VecDoub fvec;
  T &func;
  Int n;
  NRfmin(T &funcc) : func(funcc) {}
  Doub operator()(VecDoub_I &x) {
    n = x.size();
    Doub sum = 0;
    fvec = func(x);
    for (Int i = 0; i < n; i++)
      sum += SQR(fvec[i]);
    return 0.5 * sum;
  }
};
template <class T> void newt(VecDoub_IO &x, Bool &check, T &vecfunc) {
  const Int MAXITS = 200;
  const Doub TOLF = 1.0e-8, TOLMIN = 1.0e-12, STPMX = 100.0;
  const Doub TOLX = numeric_limits<Doub>::epsilon();
  Int i, j, its, n = x.size();
  Doub den, f, fold, stpmax, sum, temp, test;
  VecDoub g(n), p(n), xold(n);
  MatDoub fjac(n, n);
  NRfmin<T> fmin(vecfunc);
  NRfdjac<T> fdjac(vecfunc);
  VecDoub &fvec = fmin.fvec;
  f = fmin(x);
  test = 0.0;
  for (i = 0; i < n; i++)
    if (abs(fvec[i]) > test)
      test = abs(fvec[i]);
  if (test < 0.01 * TOLF) {
    check = false;
    return;
  }
  sum = 0.0;
  for (i = 0; i < n; i++)
    sum += SQR(x[i]);
  stpmax = STPMX * MAX(sqrt(sum), Doub(n));
  for (its = 0; its < MAXITS; its++) {
    fjac = fdjac(x, fvec);
    for (i = 0; i < n; i++) {
      sum = 0.0;
      for (j = 0; j < n; j++)
        sum += fjac[j][i] * fvec[j];
      g[i] = sum;
    }
    for (i = 0; i < n; i++)
      xold[i] = x[i];
    fold = f;
    for (i = 0; i < n; i++)
      p[i] = -fvec[i];
    LUdcmp alu(fjac);
    alu.solve(p, p);
    lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);
    test = 0.0;
    for (i = 0; i < n; i++)
      if (abs(fvec[i]) > test)
        test = abs(fvec[i]);
    if (test < TOLF) {
      check = false;
      return;
    }
    if (check) {
      test = 0.0;
      den = MAX(f, 0.5 * n);
      for (i = 0; i < n; i++) {
        temp = abs(g[i]) * MAX(abs(x[i]), 1.0) / den;
        if (temp > test)
          test = temp;
      }
      check = (test < TOLMIN);
      return;
    }
    test = 0.0;
    for (i = 0; i < n; i++) {
      temp = (abs(x[i] - xold[i])) / MAX(abs(x[i]), 1.0);
      if (temp > test)
        test = temp;
    }
    if (test < TOLX)
      return;
  }
  throw("MAXITS exceeded in newt");
}
template <class T> void broydn(VecDoub_IO &x, Bool &check, T &vecfunc) {
  const Int MAXITS = 200;
  const Doub EPS = numeric_limits<Doub>::epsilon();
  const Doub TOLF = 1.0e-8, TOLX = EPS, STPMX = 100.0, TOLMIN = 1.0e-12;
  Bool restrt, skip;
  Int i, its, j, n = x.size();
  Doub den, f, fold, stpmax, sum, temp, test;
  VecDoub fvcold(n), g(n), p(n), s(n), t(n), w(n), xold(n);
  QRdcmp *qr;
  NRfmin<T> fmin(vecfunc);
  NRfdjac<T> fdjac(vecfunc);
  VecDoub &fvec = fmin.fvec;
  f = fmin(x);
  test = 0.0;
  for (i = 0; i < n; i++)
    if (abs(fvec[i]) > test)
      test = abs(fvec[i]);
  if (test < 0.01 * TOLF) {
    check = false;
    return;
  }
  for (sum = 0.0, i = 0; i < n; i++)
    sum += SQR(x[i]);
  stpmax = STPMX * MAX(sqrt(sum), Doub(n));
  restrt = true;
  for (its = 1; its <= MAXITS; its++) {
    if (restrt) {
      qr = new QRdcmp(fdjac(x, fvec));
      if (qr->sing) {
        MatDoub one(n, n, 0.0);
        for (i = 0; i < n; i++)
          one[i][i] = 1.0;
        delete qr;
        qr = new QRdcmp(one);
      }
    } else {
      for (i = 0; i < n; i++)
        s[i] = x[i] - xold[i];
      for (i = 0; i < n; i++) {
        for (sum = 0.0, j = i; j < n; j++)
          sum += qr->r[i][j] * s[j];
        t[i] = sum;
      }
      skip = true;
      for (i = 0; i < n; i++) {
        for (sum = 0.0, j = 0; j < n; j++)
          sum += qr->qt[j][i] * t[j];
        w[i] = fvec[i] - fvcold[i] - sum;
        if (abs(w[i]) >= EPS * (abs(fvec[i]) + abs(fvcold[i])))
          skip = false;
        else
          w[i] = 0.0;
      }
      if (!skip) {
        qr->qtmult(w, t);
        for (den = 0.0, i = 0; i < n; i++)
          den += SQR(s[i]);
        for (i = 0; i < n; i++)
          s[i] /= den;
        qr->update(t, s);
        if (qr->sing)
          throw("singular update in broydn");
      }
    }
    qr->qtmult(fvec, p);
    for (i = 0; i < n; i++)
      p[i] = -p[i];
    for (i = n - 1; i >= 0; i--) {
      for (sum = 0.0, j = 0; j <= i; j++)
        sum -= qr->r[j][i] * p[j];
      g[i] = sum;
    }
    for (i = 0; i < n; i++) {
      xold[i] = x[i];
      fvcold[i] = fvec[i];
    }
    fold = f;
    qr->rsolve(p, p);
    Doub slope = 0.0;
    for (i = 0; i < n; i++)
      slope += g[i] * p[i];
    if (slope >= 0.0) {
      restrt = true;
      continue;
    }
    lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);
    test = 0.0;
    for (i = 0; i < n; i++)
      if (abs(fvec[i]) > test)
        test = abs(fvec[i]);
    if (test < TOLF) {
      check = false;
      delete qr;
      return;
    }
    if (check) {
      if (restrt) {
        delete qr;
        return;
      } else {
        test = 0.0;
        den = MAX(f, 0.5 * n);
        for (i = 0; i < n; i++) {
          temp = abs(g[i]) * MAX(abs(x[i]), 1.0) / den;
          if (temp > test)
            test = temp;
        }
        if (test < TOLMIN) {
          delete qr;
          return;
        } else
          restrt = true;
      }
    } else {
      restrt = false;
      test = 0.0;
      for (i = 0; i < n; i++) {
        temp = (abs(x[i] - xold[i])) / MAX(abs(x[i]), 1.0);
        if (temp > test)
          test = temp;
      }
      if (test < TOLX) {
        delete qr;
        return;
      }
    }
  }
  throw("MAXITS exceeded in broydn");
}

void problemI(double t, double v1, double v2, double v3) {
  double v1_prime = std::exp(-t) * std::cos(v2) + v3 * v3 - v1;
  double v2_prime = std::cos(v3 * v3) - v2;
  double v3_prime = std::cos(t) * std::exp(-v1 * v1) - v3;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "v1'(0) = " << v1_prime << std::endl;
  std::cout << "v2'(0) = " << v2_prime << std::endl;
  std::cout << "v3'(0) = " << v3_prime << std::endl;
}

VecDoub Trapezoidal_method(double low, double high, int steps, VecDoub_I y,
                           VecDoub derivs(const Doub x, VecDoub_I &y)) {

  const double h = (high - low) / (double)steps;

  VecDoub y_n = y;

  for (double x_n = low; x_n < high; x_n += h) {
    auto y_n_next = y_n + h * derivs(x_n, y_n); // Euler step
    auto func = [&](VecDoub y_guess) {
      return y_guess - y_n -
             h / 2 * (derivs(x_n, y_n) + derivs(x_n + h, y_guess));
    };
    bool check;
    newt(y_n_next, check, func);
    y_n = y_n_next;
  }
}

VecDoub derivs(const Doub t, VecDoub_I &v) {
  VecDoub dvdt(3);
  dvdt[0] = std::exp(-t) * std::cos(v[1]) + v[2] * v[2] - v[0];
  dvdt[1] = std::cos(v[2] * v[2]) - v[1];
  dvdt[2] = std::cos(t) * std::exp(-v[0] * v[0]) - v[2];
  return dvdt;
}

void problemII() {
  VecDoub y0(3);
  y0[0] = 1.0;
  y0[1] = 2.0;
  y0[2] = 3.0;

  double low = 0.0;
  double high = 5.0;
  int N_values[] = {50, 100, 200, 400, 800};

  std::cout << std::fixed << std::setprecision(6);
  for (int i = 0; i < 5; i++) {
    int N = N_values[i];
    VecDoub result = Trapezoidal_method(low, high, N, y0, derivs);
    std::cout << "N = " << N << ": ";
    std::cout << "v1(5) = " << result[0] << ", ";
    std::cout << "v2(5) = " << result[1] << ", ";
    std::cout << "v3(5) = " << result[2] << std::endl;
  }
}

int main() {
  problemI(0.0, 1.0, 2.0, 3.0);
  VecDoub y0(3);
  y0[0] = 5;
  y0[1] = 5;
  y0[2] = 5;
  VecDoub result = Trapezoidal_method(0, 5, 50, y0, derivs);
  std::cout << "v1(5) = " << result[0] << ", ";
  std::cout << "v2(5) = " << result[1] << ", ";
  std::cout << "v3(5) = " << result[2] << std::endl;
  return 0;
}
