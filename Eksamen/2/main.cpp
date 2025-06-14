#include "ludcmp.h"
#include "nr3.h"
#include "qrdcmp.h"
#include "utilities.h"
#include <cassert>
#include <cmath>
#include <print>
#include <vector>
template <class T> void newt(VecDoub_IO &x, Bool &check, T &vecfunc);

void printRoots(const std::vector<VecDoub> &roots_x,
                const std::vector<double> &lambda_vals) {
  std::vector<VecDoub> dx_k;
  dx_k.push_back(VecDoub(4));
  for (size_t i = 1; i < roots_x.size(); i++) {
    dx_k.push_back(roots_x[i] - roots_x[i - 1]);
  }

  std::vector<double> convergence;
  convergence.push_back(NAN);
  convergence.push_back(NAN);
  for (size_t i = 2; i < dx_k.size(); i++) {
    convergence.push_back(util::norm(dx_k[i]) /
                          (pow(util::norm(dx_k[i - 1]), 2)));
  }

  std::vector<double> errorX_k;
  for (size_t i = 0; i < dx_k.size(); i++) {
    errorX_k.push_back(convergence.at(i) * pow(util::norm(dx_k[i]), 2));
  }

  std::vector<double> backtracking;
  for (size_t i = 0; i < lambda_vals.size(); i++) {
    if (lambda_vals[i] < 1) {
      backtracking.push_back(lambda_vals[i]);
    } else {
      backtracking.push_back(NAN);
    }
  }

  std::print(
      "\n k  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}  {:>10}\n",
      "x0", "x1", "x2", "x3", "dx_k", "lambda", "convergence", " error");

  for (size_t i = 0; i < 7; ++i) {
    std::print("{:2}  {:10.6f}  {:10.6f}  {:10.6f}  {:10.6f}  {:10.6f}  "
               "{:10.6f}  {:10.6f}  {:10.6f}\n",
               i + 1, roots_x[i][0], roots_x[i][1], roots_x[i][2],
               roots_x[i][3], dx_k[i][0], backtracking[i], convergence[i],
               errorX_k[i]);
  }
}

VecDoub vecfunc(VecDoub_I x) {
  assert(x.size() == 4);

  VecDoub f_cal(4);

  f_cal[0] = 3 * x[0] + x[1] * sin(x[2]) - cos(x[0]) + cos(x[1] * x[1]) + 4.2;
  f_cal[1] = 3 * x[1] + x[0] * x[2] * x[3] + sin(x[1]) - 5.1;
  f_cal[2] = -pow(x[1], 2) + x[2] * x[3] * x[3] + 3 * x[2] + 5.2;
  f_cal[3] =
      x[0] + 3 * x[3] + sin(pow(x[2], 2) * pow(x[3], 2)) + cos(x[1]) - 2.3;

  return f_cal;
}

int main() {

  std::print("\n-------------------- Problem I ------------------------\n");
  Doub x_vals[] = {-0.7, 1.2, 2.3, -4.1};
  auto solutions = vecfunc(VecDoub(4, x_vals));
  std::cout << "\n" << std::endl;
  util::print(solutions);

  std::print(
      "\n--------------------Problem II & III------------------------\n");
  Doub x_guess[] = {0, 0, 0, 0};
  VecDoub_IO x(4, x_guess);
  auto functions = vecfunc;
  bool check;
  newt(x, check, functions);
  std::cout << "\n" << std::endl;
}

//----------------------- roots_multidim.h -----------------------

template <class T>
void lnsrch(VecDoub_I &xold, const Doub fold, VecDoub_I &g, VecDoub_IO &p,
            VecDoub_O &x, Doub &f, const Doub stpmax, Bool &check, T &func,
            std::vector<double> &lambda_vals) {
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
      lambda_vals.push_back(alam);
      return;
    } else if (f <= fold + ALF * alam * slope) {
      lambda_vals.push_back(alam);
      return;
    } else {
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
  std::vector<VecDoub> roots_x;
  std::vector<double> lambda_vals;
  Int i, j, its, n = x.size();
  Doub den, f, fold, stpmax, sum, temp, test;
  VecDoub g(n), p(n), xold(n);
  MatDoub fjac(n, n);
  NRfmin<T> fmin(vecfunc);
  NRfdjac<T> fdjac(vecfunc);
  VecDoub &fvec = fmin.fvec;
  f = fmin(x);
  // roots_x.push_back(x);
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
    lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin, lambda_vals);
    test = 0.0;
    for (i = 0; i < n; i++)
      if (abs(fvec[i]) > test)
        test = abs(fvec[i]);
    if (test < TOLF) {
      check = false;
      roots_x.push_back(x);
      printRoots(roots_x, lambda_vals);
      ;
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
      roots_x.push_back(x);
      printRoots(roots_x, lambda_vals);

      return;
    }
    test = 0.0;
    for (i = 0; i < n; i++) {
      temp = (abs(x[i] - xold[i])) / MAX(abs(x[i]), 1.0);
      if (temp > test)
        test = temp;
    }
    if (test < TOLX) {
      roots_x.push_back(x);
      printRoots(roots_x, lambda_vals);

      return;
    }
    roots_x.push_back(x);
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