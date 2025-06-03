#ifndef _RANPT_H_
#define _RANPT_H_
#include "nr3.h"
#include "ran.h"
void ranpt(VecDoub_O &pt, VecDoub_I &regn) {
  static const int RANSEED = 5331;
  static Ran ran(RANSEED);
  Int j, n = pt.size();
  for (j = 0; j < n; j++)
    pt[j] = regn[j] + (regn[n + j] - regn[j]) * ran.doub();
}
#endif //
