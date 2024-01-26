#ifndef SIMONTOOLS_VTXCUTS_H
#define SIMONTOOLS_VTXCUTS_H

#include "SRothman/SimonTools/src/jets.h"

struct vtxCuts {
    int fromPVcut;
    double puppiCut;
    double maxDZ;
    double maxDXY;

    vtxCuts() = default;

    vtxCuts(int fromPVcut_, double puppiCut_, 
            double maxDZ_, double maxDXY_)
        : fromPVcut(fromPVcut_), 
        puppiCut(puppiCut_), 
        maxDZ(maxDZ_), 
        maxDXY(maxDXY_) {}

    vtxCuts(const edm::ParameterSet& conf)
        : fromPVcut(conf.getParameter<int>("fromPVcut")),
        puppiCut(conf.getParameter<double>("puppiCut")),
        maxDZ(conf.getParameter<double>("maxDZ")),
        maxDXY(conf.getParameter<double>("maxDXY")) {}

    bool pass(int fromPV, double puppiWt,
              double dz, double dxy, int charge) const{
        if(charge != 0){
            return fromPV >= fromPVcut
                && puppiWt >= puppiCut
                && std::abs(dz) <= maxDZ
                && std::abs(dxy) <= maxDXY;
        } else {
            return fromPV >= fromPVcut
                && puppiWt >= puppiCut;
        }
    }

    template <typename P>
    bool pass(const P* const partptr) const {
        int fromPV;
        double puppiWt;
        double dz;
        double dxy;
        int charge = partptr -> charge();
        if constexpr (std::is_same_v<P, pat::PackedCandidate>){
            fromPV = partptr->fromPV();
            puppiWt = partptr->puppiWeight();
            dz = partptr->dz();
            dxy = partptr->dxy();
        } else {
            fromPV = true;
            puppiWt = 1.0;
            dz = 0.0;
            dxy = 0.0;
        }
        return pass(fromPV, puppiWt, dz, dxy, charge);
    }

    template <typename P>
    bool pass(const P& part) const {
        return pass(&part);
    }

    bool pass(const particle& part) const {
        return pass(part.fromPV, part.puppiweight, 
                part.dz, part.dxy, part.charge);
    }

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc){
        iDesc.add<int>("fromPVcut", 0);
        iDesc.add<double>("puppiCut", 0.0);
        iDesc.add<double>("maxDZ", 0.0);
        iDesc.add<double>("maxDXY", 0.0);
    }
};

#endif
