#ifndef _INTERP_LINEAR_H_
#define _INTERP_LINEAR_H_
#include "interp_1d.h"
#include "nr3.h"
struct Linear_interp : Base_interp {
  Linear_interp(VecDoub_I &xv, VecDoub_I &yv) : Base_interp(xv, &yv[0], 2) {}
  Doub rawinterp(Int j, Doub x) {
    if (xx[j] == xx[j + 1])
      return yy[j];
    else
      return yy[j] + ((x - xx[j]) / (xx[j + 1] - xx[j])) * (yy[j + 1] - yy[j]);
  }
};
#endif //
