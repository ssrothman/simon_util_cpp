#ifndef SIMONTOOLS_COMMON_H
#define SIMONTOOLS_COMMON_H

#include <vector>
#include <math.h>
#include "SRothman/armadillo-12.2.0/include/armadillo"
//#include <boost/math/ccmath/ccmath.hpp>

struct particle{
    double pt, eta, phi;
    double dpt, deta, dphi;
    unsigned pdgid; //absolute value
    int charge;

    particle(double pt, double eta, double phi,
             double dpt, double deta, double dphi,
             unsigned pdgid, int charge):
        pt(pt), eta(eta), phi(phi),
        dpt(dpt), deta(deta), dphi(dphi),
        pdgid(pdgid), charge(charge) {}

    particle() :
        pt(0), eta(0), phi(0),
        dpt(0), deta(0), dphi(0), 
        pdgid(0), charge(0) {}
};


struct jet{
    double pt, eta, phi;
    unsigned nPart;
    std::vector<particle> particles;
    double sumpt;
    unsigned iJet;

    inline arma::vec ptvec() const{
        arma::vec ans(nPart, arma::fill::none);
        for(unsigned i=0; i<nPart; ++i){
            ans(i) = particles[i].pt;
        }
        return ans;
    }

    inline arma::vec etavec() const{
        arma::vec ans(nPart, arma::fill::none);
        for(unsigned i=0; i<nPart; ++i){
            ans(i) = particles[i].eta;
        }
        return ans;
    }

    inline arma::vec phivec() const{
        arma::vec ans(nPart, arma::fill::none);
        for(unsigned i=0; i<nPart; ++i){
            ans(i) = particles[i].phi;
        }
        return ans;
    }

    inline arma::vec dptvec() const{
        arma::vec ans(nPart, arma::fill::none);
        for(unsigned i=0; i<nPart; ++i){
            ans(i) = particles[i].dpt;
        }
        return ans;
    }

    inline arma::vec detavec() const{
        arma::vec ans(nPart, arma::fill::none);
        for(unsigned i=0; i<nPart; ++i){
            ans(i) = particles[i].deta;
        }
        return ans;
    }

    inline arma::vec dphivec() const{
        arma::vec ans(nPart, arma::fill::none);
        for(unsigned i=0; i<nPart; ++i){
            ans(i) = particles[i].dphi;
        }
        return ans;
    }
};

//constexpr double inv_sqrt_2pi = 1/boost::math::ccmath::sqrt(M_PI);
constexpr double inv_sqrt_2pi = 0.3989422804; 
inline double normal_pdf(double x, double m, double s){
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

#endif
