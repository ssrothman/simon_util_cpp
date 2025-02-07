#ifndef SIMON_UTIL_DR_H
#define SIMON_UTIL_DR_H

#include <math.h>
#include "simon_jet.h"

inline double deltaEta(double eta1, double eta2){
    return eta2 - eta1;
}

inline double deltaPhi(double phi1, double phi2){
    double dphi = phi2 - phi1;
    if(dphi>M_PI){
        dphi -= 2*M_PI;
    } else if(dphi<-M_PI){
        dphi += 2*M_PI;
    }
    return dphi;
}

inline Eigen::VectorXd deltaPhi(const Eigen::VectorXd& phi1, const Eigen::VectorXd& phi2){
    Eigen::VectorXd dphi(phi1.rows());

    for(unsigned i=0; i<dphi.rows(); i++){
        dphi(i) = deltaPhi(phi1(i), phi2(i));
    }

    return dphi;
}

inline double deltaR2(double eta1, double phi1, double eta2, double phi2) {
  /*
    * Compute delta R^2 between (eta1, phi1) and (eta2, phi2)
    */
  double deta = eta2 - eta1;
  double dphi = deltaPhi(phi1, phi2);

  return deta * deta + dphi * dphi;
}

inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
  /*
    * Compute delta R between (eta1, phi1) and (eta2, phi2)
    */
  return sqrt(deltaR2(eta1, phi1, eta2, phi2));
}

inline double deltaR2(const simon_particle& p1, const simon_particle& p2) {
  /*
    * Compute delta R^2 between p1 and p2
    */
  return deltaR2(p1.eta, p1.phi, p2.eta, p2.phi);
}

inline double deltaR(const simon_particle& p1, const simon_particle& p2) {
  /*
    * Compute delta R between p1 and p2
    */
  return deltaR(p1.eta, p1.phi, p2.eta, p2.phi);
}

#endif
