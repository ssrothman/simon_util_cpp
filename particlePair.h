#ifndef SROTHMAN_SIMONTOOLS_PARTICLEPAIR_H
#define SROTHMAN_SIMONTOOLS_PARTICLEPAIR_H

#include "histutil.h"

namespace simon{
    template <bool distances_squared>
    struct pairentry_resolved {
        double deta;
        double dphi;
        double floatDR;

        pairentry_resolved() noexcept :
            deta(0), dphi(0), floatDR(0) {}

        pairentry_resolved(double deta, double dphi) noexcept :
                deta(deta), dphi(dphi) {
            if constexpr(distances_squared){
                floatDR = deta*deta + dphi*dphi;
            } else {
                floatDR = std::sqrt(deta*deta + dphi*dphi);
            }
        }

        pairentry_resolved(double deta, double dphi, double floatDR) noexcept :
                deta(deta), dphi(dphi), floatDR(floatDR) {}

        template <typename T>
        pairentry_resolved(const T& p1, const T& p2) noexcept {
            deta = simon::deltaEta(p1.eta, p2.eta); 
            dphi = simon::deltaPhi(p1.phi, p2.phi);
            if constexpr(distances_squared){
                floatDR = deta*deta + dphi*dphi;
            } else {
                floatDR = std::sqrt(deta*deta + dphi*dphi);
            }
        }
    };

    struct pairentry_projected {
        unsigned dR_index;

        pairentry_projected() noexcept : dR_index(0) {}

        pairentry_projected(unsigned dR_index) noexcept : dR_index(dR_index) {}

        template <typename AXIS, bool distances_squared>
        pairentry_projected(const pairentry_resolved<distances_squared>& resolved,
                           const AXIS& axis) noexcept {
            if constexpr(distances_squared){
                dR_index = simon::getIndex(resolved.floatDR, axis);
            } else {
                dR_index = simon::getIndex(std::sqrt(resolved.floatDR), axis);
            }
        }

        template <typename AXIS>
        pairentry_projected(double floatDR, const AXIS& axis) noexcept:
            dR_index(simon::getIndex(floatDR, axis)) {}

        template <typename AXIS, typename T>
        pairentry_projected(const T& p1, const T& p2,
                            const AXIS& axis) noexcept:
            dR_index(simon::getIndex(simon::deltaR(
                            p1.eta, p1.phi,
                            p2.eta, p2.phi), 
                        axis)) {}
    };

    template <bool distances_squared>
    inline double angle_between(
            const pairentry_resolved<distances_squared>& pair1,
            const pairentry_resolved<distances_squared>& pair2) noexcept {

        if(pair1.floatDR == 0 || pair2.floatDR == 0){
            return 0;
        }
        double dot = pair1.deta*pair2.deta + pair1.dphi*pair2.dphi;
        
        double cosphi;
        if constexpr(distances_squared){
            cosphi = dot / std::sqrt(pair1.floatDR * pair2.floatDR);
        } else {
            cosphi = dot / (pair1.floatDR * pair2.floatDR);
        }

        if(cosphi > 1){
            if (cosphi > 1.001){
                printf("uh oh.... cosphi-1 = %0.2g\n", cosphi-1);
                assert(false);
            }
            cosphi = 1;
        } else if (cosphi < -1){
            if (cosphi < -1.001){
                printf("uh oh.... cosphi+1 = %0.2g\n", cosphi+1);
                assert(false);
            }
            cosphi = -1;
        }

        double angle = std::acos(cosphi);
        if (angle > M_PI/2){
            angle = M_PI - angle;
        }
        return angle;
    }
};

#endif
