#ifndef SIMON_UTIL_DR_H
#define SIMON_UTIL_DR_H

#include <math.h>

inline double dR2(double eta1, double phi1, double eta2, double phi2) {
  /*
    * Compute delta R^2 between (eta1, phi1) and (eta2, phi2)
    */
  double deta = eta2 - eta1;
  double dphi = phi2 - phi1;
  if (dphi > M_PI) {
    dphi = 2 * M_PI - dphi;
  } else if (dphi < -M_PI) {
    dphi = 2 * M_PI + dphi;  //correct up to sign
  }
  return deta * deta + dphi * dphi;
}

#endif
