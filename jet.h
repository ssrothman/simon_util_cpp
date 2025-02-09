#ifndef SIMONTOOLS_COMMON_H
#define SIMONTOOLS_COMMON_H

#include <vector>
#include <math.h>
#include <cmath>
#include <Eigen/Dense>

namespace simon{
    struct particle{
        double pt, eta, phi;
        double dpt, deta, dphi;
        unsigned pdgid; //absolute value
        int charge;

        //extras
        double vtx_x, vtx_y, vtx_z;
        double dxy, dz;
        int fromPV;
        double puppiweight;

        particle(double pt, double eta, double phi):
            pt(pt), eta(eta), phi(phi),
            dpt(-1), deta(-1), dphi(-1),
            pdgid(0), charge(0),
            vtx_x(9999), vtx_y(9999), vtx_z(9999),
            dxy(9999), dz(9999),
            fromPV(9999), puppiweight(9999) {}

        particle(double pt, double eta, double phi,
                 double dpt, double deta, double dphi,
                 unsigned pdgid, int charge):
            pt(pt), eta(eta), phi(phi),
            dpt(dpt), deta(deta), dphi(dphi),
            pdgid(pdgid), charge(charge),
            vtx_x(9999), vtx_y(9999), vtx_z(9999),
            dxy(9999), dz(9999),
            fromPV(9999), puppiweight(9999) {}

        particle(double pt, double eta, double phi,
                 unsigned pdgid, int charge,
                 double vtx_x, double vtx_y, double vtx_z,
                 double dxy, double dz,
                 int fromPV, double puppiweight):
            pt(pt), eta(eta), phi(phi),
            dpt(-1), deta(-1), dphi(-1),
            pdgid(pdgid), charge(charge),
            vtx_x(vtx_x), vtx_y(vtx_y), vtx_z(vtx_z),
            dxy(dxy), dz(dz),
            fromPV(fromPV), puppiweight(puppiweight) {}

        particle() :
            pt(0), eta(0), phi(0),
            dpt(0), deta(0), dphi(0), 
            pdgid(0), charge(0),
            vtx_x(9999), vtx_y(9999), vtx_z(9999),
            dxy(9999), dz(9999),
            fromPV(9999), puppiweight(9999) {}
    };


    struct jet{
        double pt, sumpt, rawpt;

        double eta, phi;

        double mass;

        unsigned nPart;
        std::vector<particle> particles;
        unsigned iJet;

        //extras
        unsigned nEM0, nHAD0, nHADCH, nELE, nMU;
        double jecfactor;
        std::vector<unsigned> iCHS;
        unsigned iSV;

        jet() :
            pt(0), sumpt(0), rawpt(0),
            eta(0), phi(0), 
            mass(0),
            nPart(0), particles(),
            iJet(9999),
            nEM0(0), nHAD0(0), nHADCH(0), 
            nELE(0), nMU(0), jecfactor(9999), 
            iCHS(), iSV(9999) {}


        inline Eigen::VectorXd ptvec() const{
            Eigen::VectorXd ans(nPart);
            for(unsigned i=0; i<nPart; ++i){
                ans[i] = particles[i].pt;
            }
            return ans;
        }

        inline Eigen::VectorXd etavec() const{
            Eigen::VectorXd ans(nPart);
            for(unsigned i=0; i<nPart; ++i){
                ans[i] = particles[i].eta;
            }
            return ans;
        }

        inline Eigen::VectorXd phivec() const{
            Eigen::VectorXd ans(nPart);
            for(unsigned i=0; i<nPart; ++i){
                ans[i] = particles[i].phi;
            }
            return ans;
        }

        inline Eigen::VectorXd dptvec() const{
            Eigen::VectorXd ans(nPart);
            for(unsigned i=0; i<nPart; ++i){
                ans[i] = particles[i].dpt;
            }
            return ans;
        }

        inline Eigen::VectorXd detavec() const{
            Eigen::VectorXd ans(nPart);
            for(unsigned i=0; i<nPart; ++i){
                ans[i] = particles[i].deta;
            }
            return ans;
        }

        inline Eigen::VectorXd dphivec() const{
            Eigen::VectorXd ans(nPart);
            for(unsigned i=0; i<nPart; ++i){
                ans[i] = particles[i].dphi;
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
};

#endif
