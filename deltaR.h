#ifndef SIMON_UTIL_DR_H
#define SIMON_UTIL_DR_H

#include <math.h>
#include "SRothman/armadillo-12.2.0/include/armadillo"
#include "jets.h"

inline double deltaphi(double phi1, double phi2){
    double dphi = phi2 - phi1;
    if(dphi>M_PI){
        dphi -= 2*M_PI;
    } else if(dphi<-M_PI){
        dphi += 2*M_PI;
    }
    return dphi;
}

inline arma::vec deltaphi(const arma::vec& phi1, const arma::vec& phi2){
    arma::vec dphi(phi1.n_elem, arma::fill::none);

    for(unsigned i=0; i<dphi.n_elem; i++){
        dphi(i) = deltaphi(phi1(i), phi2(i));
    }

    return dphi;
}

inline double dR2(double eta1, double phi1, double eta2, double phi2) {
  /*
    * Compute delta R^2 between (eta1, phi1) and (eta2, phi2)
    */
  double deta = eta2 - eta1;
  double dphi = deltaphi(phi1, phi2);

  return deta * deta + dphi * dphi;
}

inline double dR(double eta1, double phi1, double eta2, double phi2) {
  /*
    * Compute delta R between (eta1, phi1) and (eta2, phi2)
    */
  return sqrt(dR2(eta1, phi1, eta2, phi2));
}

inline double dR2(const particle& p1, const particle& p2) {
  /*
    * Compute delta R^2 between p1 and p2
    */
  return dR2(p1.eta, p1.phi, p2.eta, p2.phi);
}

inline double dR(const particle& p1, const particle& p2) {
  /*
    * Compute delta R between p1 and p2
    */
  return dR(p1.eta, p1.phi, p2.eta, p2.phi);
}

#endif
