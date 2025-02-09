#ifndef SIMONTOOLS_PARTICLETHRESHOLDS_H
#define SIMONTOOLS_PARTICLETHRESHOLDS_H

#include "jet.h"
#include "isID.h"
#include "etaRegion.h"

namespace simon{
    struct particleThresholds{
        std::vector<double> EM0thresholds, HAD0thresholds,
            HADCHthresholds, ELEthresholds, MUthresholds;

        std::vector<double> etaRegions;

        particleThresholds() = default;
        particleThresholds(const std::vector<double>& EM0thresh,
                           const std::vector<double>& HAD0thresh,
                           const std::vector<double>& HADCHthresh,
                           const std::vector<double>& ELEthresh,
                           const std::vector<double>& MUthresh,
                           const std::vector<double>& etaRegions):
            EM0thresholds(EM0thresh),
            HAD0thresholds(HAD0thresh),
            HADCHthresholds(HADCHthresh),
            ELEthresholds(ELEthresh),
            MUthresholds(MUthresh),
            etaRegions(etaRegions){}

        particleThresholds(const edm::ParameterSet& conf):
            EM0thresholds(conf.getParameter<std::vector<double>>("EM0thresholds")),
            HAD0thresholds(conf.getParameter<std::vector<double>>("HAD0thresholds")),
            HADCHthresholds(conf.getParameter<std::vector<double>>("HADCHthresholds")),
            ELEthresholds(conf.getParameter<std::vector<double>>("ELEthresholds")),
            MUthresholds(conf.getParameter<std::vector<double>>("MUthresholds")),
            etaRegions(conf.getParameter<std::vector<double>>("etaRegions")){}

        template <typename P>
        double getThreshold(const P& part, int region) const {
            if(region<0 || region>=(int)etaRegions.size()-1){
                throw std::runtime_error(
                        "particleThresholds::getThreshold: "
                        "particle eta out of range");
            } else if(isEM0(part)){
                return EM0thresholds[region];
            } else if(isHAD0(part)){
                return HAD0thresholds[region];
            } else if(isHADCH(part)){
                return HADCHthresholds[region];
            } else if(isELE(part)){
                return ELEthresholds[region];
            } else if(isMU(part)){
                return MUthresholds[region];
            } else {
                /*throw std::runtime_error(
                        "particleThresholds::getThreshold: "
                        "particle flavor not recognized");*/
                printf("WARNING: particleThresholds::getThreshold: "
                        "particle flavor not recognized\n");
                return 0;
                
            }
        }

        template <typename P>
        double getThreshold(const P& part) const{
            int region = getEtaRegion(part.eta(), etaRegions);
            if(region > (int)etaRegions.size()-1){
                throw std::runtime_error(
                        "particleThresholds::getThreshold: "
                        "particle eta out of range " + std::to_string(part.eta()));
            }
            return getThreshold(part, region);
        }

        double getThreshold(const particle& part) const{
            int region = getEtaRegion(part.eta, etaRegions);
            if(region > (int)etaRegions.size()-1){
                throw std::runtime_error(
                        "particleThresholds::getThreshold: "
                        "particle eta out of range " + std::to_string(part.eta));
            }
            return getThreshold(part, region);
        }

        template <typename P>
        double getThreshold(const P* const part) const{
            return getThreshold(*part);
        }

        static void fillPSetDescription(edm::ParameterSetDescription& desc){
            desc.add<std::vector<double>>("EM0thresholds", {0.0, 0.0, 0.0});
            desc.add<std::vector<double>>("HAD0thresholds", {0.0, 0.0, 0.0});
            desc.add<std::vector<double>>("HADCHthresholds", {0.0, 0.0, 0.0});
            desc.add<std::vector<double>>("ELEthresholds", {0.0, 0.0, 0.0});
            desc.add<std::vector<double>>("MUthresholds", {0.0, 0.0, 0.0});
            desc.add<std::vector<double>>("etaRegions", {0.0, 1.0, 2.0, 3.0});
        }
    };
};

#endif

