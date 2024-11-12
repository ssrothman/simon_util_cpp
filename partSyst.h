#ifndef SIMONTOOLS_PARTSYST_H
#define SIMONTOOLS_PARTSYST_H

#include <vector>
#include <limits.h>
#include <random>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "jets.h"
#include "isID.h"
#include "etaRegion.h"
#include "deltaR.h"

class partSyst{
public:
    std::default_random_engine rng;
    std::uniform_real_distribution<double> rand;
    std::normal_distribution<double> normal;

    enum SYSTEMATIC{
        NOM,
        EM0_UP,
        EM0_DN,
        HAD0_UP,
        HAD0_DN,
        CH_UP,
        CH_DN,
        TRK_EFF,
        TRK_EFF_DR,
        RECO_EFF,
    };

    static SYSTEMATIC getSystEnum(std::string name){
        if(name == "EM0_UP") return EM0_UP;
        if(name == "EM0_DN") return EM0_DN;
        if(name == "HAD0_UP") return HAD0_UP;
        if(name == "HAD0_DN") return HAD0_DN;
        if(name == "CH_UP") return CH_UP;
        if(name == "CH_DN") return CH_DN;
        if(name == "TRK_EFF") return TRK_EFF;
        if(name == "TRK_EFF_DR") return TRK_EFF_DR;
        if(name == "RECO_EFF") return RECO_EFF;
        return NOM;
    }
    /*
     * Detector energy scales:
     *    EM0scale: amount by which to vary photon energy scale
     *    HAD0scale: amount by which to vary hadron energy scale
     *    CHscale: amount by which to vary charged hadron energy scale
     */
    double EM0scale;
    double HAD0scale;
    double CHscale;

    /*
     * Tracking efficiency systematic
     *   trkDropProb: probability of dropping a track
     *   trkDropSmear: amount by which to smear the energy of a dropped track
     */
    double trkDropProb; 
    double trkDropSmear;

    /*
     * DR-dependent tracking efficiency systematic
     *  DRtrkDropProbs: vector of probabilities of dropping a track as a function of DR
     *  DRtrkDropEdges: vector of edges for DR bins
     */
    std::vector<double> DRtrkDropProbs;
    std::vector<double> DRtrkDropEdges;

    /*
     * Reconstruction efficiency systematic
     *  pDropEM0: probability of dropping a photon as a function of pT
     *  pDropHAD0: probability of dropping a hadron as a function of pT
     *  pDropELE: probability of dropping an electron as a function of pT
     *  pDropMU: probability of dropping a muon as a function of pT
     *  pDropHADCH: probability of dropping a charged hadron as a function of pT
     *  pDropEdges: vector of edges for pT bins
     */
    std::vector<double> pDropEM0;
    std::vector<double> pDropHAD0;
    std::vector<double> pDropELE;
    std::vector<double> pDropMU;
    std::vector<double> pDropHADCH;
    std::vector<double> pDropEdges;

    partSyst():
        rng(0),
        rand(0.0, 1.0),
        normal(0.0, 1.0),

        EM0scale(0.03),
        HAD0scale(0.05),
        CHscale(0.01),

        trkDropProb(0.03),
        trkDropSmear(0.10),

        DRtrkDropProbs({0.0}),
        DRtrkDropEdges({0.0, std::numeric_limits<double>::max()}),

        pDropEM0({0.0}),
        pDropHAD0({0.0}),
        pDropELE({0.0}),
        pDropMU({0.0}),
        pDropHADCH({0.0}),
        pDropEdges({0.0, std::numeric_limits<double>::max()})
    {};

    partSyst(const double EM0scale,
                 const double HAD0scale,
                 const double CHscale,

                 const double trkDropProb,
                 const double trkDropSmear,

                 const std::vector<double> DRtrkDropProbs,
                 const std::vector<double> DRtrkDropEdges,

                 const std::vector<double> pDropEM0,
                 const std::vector<double> pDropHAD0,
                 const std::vector<double> pDropELE,
                 const std::vector<double> pDropMU,
                 const std::vector<double> pDropHADCH,
                 const std::vector<double> pDropEdges):
        rng(0),
        rand(0.0, 1.0),
        normal(0.0, 1.0),

        EM0scale(EM0scale),
        HAD0scale(HAD0scale),
        CHscale(CHscale),

        trkDropProb(trkDropProb),
        trkDropSmear(trkDropSmear),

        DRtrkDropProbs(DRtrkDropProbs),
        DRtrkDropEdges(DRtrkDropEdges),

        pDropEM0(pDropEM0),
        pDropHAD0(pDropHAD0),
        pDropELE(pDropELE),
        pDropMU(pDropMU),
        pDropHADCH(pDropHADCH),
        pDropEdges(pDropEdges){};

    partSyst(const edm::ParameterSet& conf):
        rng(0),
        rand(0.0, 1.0),
        normal(0.0, 1.0),

        EM0scale(conf.getParameter<double>("EM0scale")),
        HAD0scale(conf.getParameter<double>("HAD0scale")),
        CHscale(conf.getParameter<double>("CHscale")),

        trkDropProb(conf.getParameter<double>("trkDropProb")),
        trkDropSmear(conf.getParameter<double>("trkDropSmear")),

        DRtrkDropProbs(conf.getParameter<std::vector<double>>("DRtrkDropProbs")),
        DRtrkDropEdges(conf.getParameter<std::vector<double>>("DRtrkDropEdges")),

        pDropEM0(conf.getParameter<std::vector<double>>("pDropEM0")),
        pDropHAD0(conf.getParameter<std::vector<double>>("pDropHAD0")),
        pDropELE(conf.getParameter<std::vector<double>>("pDropELE")),
        pDropMU(conf.getParameter<std::vector<double>>("pDropMU")),
        pDropHADCH(conf.getParameter<std::vector<double>>("pDropHADCH")),
        pDropEdges(conf.getParameter<std::vector<double>>("pDropEdges")){};

    static void fillPSetDescription(edm::ParameterSetDescription& desc){
        desc.add<double>("EM0scale", 0.03);
        desc.add<double>("HAD0scale", 0.05);
        desc.add<double>("CHscale", 0.01);

        desc.add<double>("trkDropProb", 0.03);
        desc.add<double>("trkDropSmear", 0.10);

        desc.add<std::vector<double>>("DRtrkDropProbs", {0.0});
        desc.add<std::vector<double>>("DRtrkDropEdges", {0.0, std::numeric_limits<double>::max()});

        desc.add<std::vector<double>>("pDropEM0", {0.0});
        desc.add<std::vector<double>>("pDropHAD0", {0.0});
        desc.add<std::vector<double>>("pDropELE", {0.0});
        desc.add<std::vector<double>>("pDropMU", {0.0});
        desc.add<std::vector<double>>("pDropHADCH", {0.0});
        desc.add<std::vector<double>>("pDropEdges", {0.0, std::numeric_limits<double>::max()});
    }

    bool applySystematic(const SYSTEMATIC sys, particle& part, double jetEta, double jetPhi) {
        switch(sys){
            case EM0_UP:
                if(isEM0(part)){
                    part.pt *= (1.0 + EM0scale);
                }
                return true;
            case EM0_DN:
                if(isEM0(part)){
                    part.pt *= (1.0 - EM0scale);
                }
                return true;
            case HAD0_UP:
                if(isHAD0(part)){
                    part.pt *= (1.0 + HAD0scale);
                }
                return true;
            case HAD0_DN:
                if(isHAD0(part)){
                    part.pt *= (1.0 - HAD0scale);
                }
                return true;
            case CH_UP:
                if(part.charge !=0){
                    part.pt *= (1.0 + CHscale);
                }
                return true;
            case CH_DN:
                if(part.charge !=0){
                    part.pt *= (1.0 - CHscale);
                }
                return true;
            case TRK_EFF:
               {if(part.charge !=0){
                    if(rand(rng) < trkDropProb){
                        if(isMU(part)){
                            return false;
                        } else if(isELE(part)){
                            part.pdgid = 22;
                        } else {
                            part.pdgid = 130;
                        }
                        part.charge = 0;
                        double smear = normal(rng) * trkDropSmear + 1.0;
                        part.pt *= smear;
                    }
                }
                return true;}
            case TRK_EFF_DR:
               {if(part.charge !=0){
                    double deltaR = dR(jetEta, jetPhi, part.eta, part.phi);
                    int bin = getEtaRegion(deltaR, DRtrkDropEdges);
                    double prob = DRtrkDropProbs.at(bin);
                    if(rand(rng) < prob){
                        if(isMU(part)){
                            return false;
                        } else if(isELE(part)){
                            part.pdgid = 22;
                        } else {
                            part.pdgid = 130;
                        }
                        part.charge = 0;
                        double smear = normal(rng) * trkDropSmear + 1.0;
                        part.pt *= smear;
                    }
                }
                return true;}
            case RECO_EFF:
               {int bin = getEtaRegion(part.pt, pDropEdges);
                double prob=0;
                if(isEM0(part)){
                    prob = pDropEM0.at(bin);
                } else if(isHAD0(part)){
                    prob = pDropHAD0.at(bin);
                } else if(isELE(part)){
                    prob = pDropELE.at(bin);
                } else if(isMU(part)){
                    prob = pDropMU.at(bin);
                } else if(isHADCH(part)){
                    prob = pDropHADCH.at(bin);
                }
                if(rand(rng) < prob){
                    return false;
                }
                return true;}
            default:
                return true;
        }
    }
};

#endif
