#ifndef SIMONTOOLS_PART_SELECTOR_H
#define SIMONTOOLS_PART_SELECTOR_H

#include <vector>
#include <limits.h>
#include <random>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "jet.h"
#include "isID.h"
#include "etaRegion.h"
#include "deltaR.h"

namespace simon{
    class syst_parameters{
    public:
        double EM0scale, HAD0scale, CHscale;
        double trkDropProb, trkDropSmear;
        std::vector<double> EM0thresholds, HAD0thresholds, ELEthresholds, MUthresholds, HADCHthresholds;
        int minFromPV;
        double minPuppiWt, maxDZ, maxDXY;

        syst_parameters(const edm::ParameterSet& conf):
            EM0scale(conf.getParameter<double>("EM0scale")),
            HAD0scale(conf.getParameter<double>("HAD0scale")),
            CHscale(conf.getParameter<double>("CHscale")),

            trkDropProb(conf.getParameter<double>("trkDropProb")),
            trkDropSmear(conf.getParameter<double>("trkDropSmear")),

            EM0thresholds(conf.getParameter<std::vector<double>>("EM0thresholds")),
            HAD0thresholds(conf.getParameter<std::vector<double>>("HAD0thresholds")),
            ELEthresholds(conf.getParameter<std::vector<double>>("ELEthresholds")),
            MUthresholds(conf.getParameter<std::vector<double>>("MUthresholds")),
            HADCHthresholds(conf.getParameter<std::vector<double>>("HADCHthresholds")),

            minFromPV(conf.getParameter<int>("minFromPV")),
            minPuppiWt(conf.getParameter<double>("minPuppiWt")),
            maxDZ(conf.getParameter<double>("maxDZ")),
            maxDXY(conf.getParameter<double>("maxDXY")) {}

        static void fillPSetDescription(edm::ParameterSetDescription& desc){
            desc.add<double>("EM0scale");
            desc.add<double>("HAD0scale");
            desc.add<double>("CHscale");

            desc.add<double>("trkDropProb");
            desc.add<double>("trkDropSmear");

            desc.add<std::vector<double>>("EM0thresholds");
            desc.add<std::vector<double>>("HAD0thresholds");
            desc.add<std::vector<double>>("ELEthresholds");
            desc.add<std::vector<double>>("MUthresholds");
            desc.add<std::vector<double>>("HADCHthresholds");

            desc.add<int>("minFromPV");
            desc.add<double>("minPuppiWt");
            desc.add<double>("maxDZ");
            desc.add<double>("maxDXY");
        }
    };

    class syst_settings{
    public:
        std::string EM0scale, HAD0scale, CHscale;
        std::string trkDrop;
        std::string EM0threshold, HAD0threshold;
        std::string ELEthreshold, MUthreshold, HADCHthreshold;
        std::string requireVertex;

        bool applyPuppi, onlyCharged;

        syst_settings(const edm::ParameterSet& conf):
            EM0scale(conf.getParameter<std::string>("EM0scale")),
            HAD0scale(conf.getParameter<std::string>("HAD0scale")),
            CHscale(conf.getParameter<std::string>("CHscale")),
            trkDrop(conf.getParameter<std::string>("trkDrop")),
            EM0threshold(conf.getParameter<std::string>("EM0threshold")),
            HAD0threshold(conf.getParameter<std::string>("HAD0threshold")),
            ELEthreshold(conf.getParameter<std::string>("ELEthreshold")),
            MUthreshold(conf.getParameter<std::string>("MUthreshold")),
            HADCHthreshold(conf.getParameter<std::string>("HADCHthreshold")),
            requireVertex(conf.getParameter<std::string>("requireVertex")),
            applyPuppi(conf.getParameter<bool>("applyPuppi")),
            onlyCharged(conf.getParameter<bool>("onlyCharged")) {}

        static void fillPSetDescription(edm::ParameterSetDescription& desc){
            desc.add<std::string>("EM0scale");
            desc.add<std::string>("HAD0scale");
            desc.add<std::string>("CHscale");
            desc.add<std::string>("trkDrop");
            desc.add<std::string>("EM0threshold");
            desc.add<std::string>("HAD0threshold");
            desc.add<std::string>("ELEthreshold");
            desc.add<std::string>("MUthreshold");
            desc.add<std::string>("HADCHthreshold");
            desc.add<std::string>("requireVertex");

            desc.add<bool>("applyPuppi");
            desc.add<bool>("onlyCharged");
        }
    };

    class particleSelector{
    public:
        std::default_random_engine rng;
        std::uniform_real_distribution<double> rand;
        std::normal_distribution<double> normal;

        /*
         */
        double EM0factor;
        double HAD0factor;
        double CHfactor;

        /*
         * Tracking efficiency systematic
         *   trkDropProb: probability of dropping a track
         *   trkDropSmear: amount by which to smear the 
         *                 energy of a dropped track
         */
        double trkDropProb; 
        double trkDropSmear;

        /*
         * Reconstruction thresholds:
         */
        double EM0threshold;
        double HAD0threshold;
        double ELEthreshold;
        double MUthreshold;
        double HADCHthreshold;

        /*
         * Vertexing selections
         */
        int minFromPV;
        double minPuppiWt;
        double maxDZ;
        double maxDXY;

        /*
         * Config booleans
         */
        bool applyPuppi;
        bool onlyCharged;

        double DN_NOM_UP(const std::string& str){
            if(str == "NOM"){
                return 0.0;
            } else if(str == "UP"){
                return 1.0;
            } else if(str == "DN"){
                return -1.0;
            } else {
                throw std::runtime_error("Invalid systematic variation: " + str);
            }
        }

        bool ON_OFF(const std::string& str){
            if(str == "ON"){
                return true;
            } else if(str == "OFF"){
                return false;
            } else {
                throw std::runtime_error("Invalid boolean: " + str);
            }
        }

        particleSelector(const edm::ParameterSet& conf):
                rng(0),
                rand(0.0, 1.0),
                normal(0.0, 1.0) {
        
            syst_parameters params(conf.getParameter<edm::ParameterSet>("parameters"));
            syst_settings settings(conf.getParameter<edm::ParameterSet>("settings"));

            EM0factor = 1 + DN_NOM_UP(settings.EM0scale) * params.EM0scale;
            HAD0factor = 1 + DN_NOM_UP(settings.HAD0scale) * params.HAD0scale;
            CHfactor = 1 + DN_NOM_UP(settings.CHscale) * params.CHscale;

            trkDropProb = ON_OFF(settings.trkDrop) ? params.trkDropProb : 0.0;
            trkDropSmear = ON_OFF(settings.trkDrop) ? params.trkDropSmear : 0.0;

            EM0threshold = params.EM0thresholds[1+DN_NOM_UP(settings.EM0threshold)];
            HAD0threshold = params.HAD0thresholds[1+DN_NOM_UP(settings.HAD0threshold)];
            ELEthreshold = params.ELEthresholds[1+DN_NOM_UP(settings.ELEthreshold)];
            MUthreshold = params.MUthresholds[1+DN_NOM_UP(settings.MUthreshold)];
            HADCHthreshold = params.HADCHthresholds[1+DN_NOM_UP(settings.HADCHthreshold)];

            minFromPV = ON_OFF(settings.requireVertex) ? params.minFromPV : 0;
            minPuppiWt = ON_OFF(settings.requireVertex) ? params.minPuppiWt : 0.0;
            maxDZ = ON_OFF(settings.requireVertex) ? params.maxDZ : 999.0;
            maxDXY = ON_OFF(settings.requireVertex) ? params.maxDXY : 999.0;

            applyPuppi = settings.applyPuppi;
            onlyCharged = settings.onlyCharged;
        }

        static void fillPSetDescription(edm::ParameterSetDescription& desc){
            edm::ParameterSetDescription params;
            syst_parameters::fillPSetDescription(params);
            desc.add<edm::ParameterSetDescription>("parameters", params);

            edm::ParameterSetDescription settings;
            syst_settings::fillPSetDescription(settings);
            desc.add<edm::ParameterSetDescription>("settings", settings);
        }

        template <typename T>
        bool makeParticle(const T* const partptr, particle& nextpart){
            nextpart.pt = 0;

            //lookup particle properties
            int fromPV;
            double puppiWeight, dxy, dz;
            if constexpr (std::is_same_v<T, pat::PackedCandidate>){
                fromPV = partptr->fromPV();
                puppiWeight = partptr->puppiWeight();
                dxy = partptr->dxy();
                dz = partptr->dz();
            } else {
                fromPV = true;
                puppiWeight = 1.0;
                dxy = 0;
                dz = 0;
            }

            double nextpt = partptr->pt();
            if(applyPuppi){
                nextpt *= puppiWeight;
            }

            unsigned pdgid = std::abs(partptr->pdgId());

            //reject neutrinos early
            if(pdgid == 12 || pdgid == 14 || pdgid == 16){
                return false;
            }

            //build the particle
            nextpart.pt = nextpt;
            nextpart.eta = partptr->eta();
            nextpart.phi = partptr->phi();
            nextpart.pdgid = std::abs(partptr->pdgId());
            nextpart.charge = partptr->charge();
            nextpart.dxy = dxy;
            nextpart.dz = dz;
            nextpart.fromPV = fromPV;
            nextpart.vtx_x = partptr->vertex().x();
            nextpart.vtx_y = partptr->vertex().y();
            nextpart.vtx_z = partptr->vertex().z();
            nextpart.puppiweight = puppiWeight;

            //do track dropping
            if(nextpart.charge !=0 && rand(rng) < trkDropProb){
                nextpart.pt *= trkDropSmear*normal(rng) + 1;
                nextpart.charge = 0;
                if(nextpart.pdgid == 11){
                    nextpart.pdgid = 22;
                } else {
                    nextpart.pdgid = 130;
                }
            }

            //apply energy scale systematics
            if(nextpart.pdgid == 22){
                nextpart.pt *= EM0factor;
            } else if(nextpart.charge == 0){
                nextpart.pt *= HAD0factor;
            } else {
                nextpart.pt *= CHfactor;
            }

            if (onlyCharged && partptr->charge() == 0){
                return false;
            }

            //apply thresholds
            if(pdgid==11 && nextpt < ELEthreshold){
                return false;
            }
            if(pdgid==13 && nextpt < MUthreshold){
                return false;
            }
            if(pdgid==22 && nextpt < EM0threshold){
                return false;
            }
            if(pdgid >=100 && partptr->charge() ==0 && nextpt < HAD0threshold){
                return false;
            }
            if(pdgid >=100 && partptr->charge() !=0 && nextpt < HADCHthreshold){
                return false;
            }

            //apply vertexing cuts
            if(std::abs(dz) < maxDZ || std::abs(dxy) < maxDXY 
                    || puppiWeight < minPuppiWt || fromPV < minFromPV){
                return false;
            }

            return true;
        }

        template <typename COLLECTION, typename F=std::function<bool(const edm::Ptr<reco::Candidate>&)>>
        void buildJet(const COLLECTION& parts, jet& result,
                      F* filter = nullptr){
            result.rawpt = 0;
            result.sumpt = 0;

            for (const auto& part : parts){
                if (filter && !((*filter)(part))){
                    continue;
                }
                
                const auto* recoptr = dynamic_cast<const pat::PackedCandidate*>(part.get());
                const auto* genptr = dynamic_cast<const pat::PackedGenParticle*>(part.get());
                const auto* genptr2 = dynamic_cast<const reco::GenParticle*>(part.get());
        
                particle nextpart;
                bool passed;
                if(recoptr){
                    passed = makeParticle(recoptr, nextpart);
                } else if(genptr){
                    passed = makeParticle(genptr, nextpart);
                } else if(genptr2){
                    passed = makeParticle(genptr2, nextpart);
                } else {
                    throw std::runtime_error("Unknown particle type");
                }

                result.rawpt += nextpart.pt;
                if (passed){
                    result.sumpt += nextpart.pt;
                    result.particles.push_back(nextpart);
                    ++result.nPart;

                    if(simon::isELE(nextpart)){
                        ++result.nELE;
                    } else if(simon::isMU(nextpart)){
                        ++result.nMU;
                    } else if(simon::isEM0(nextpart)){
                        ++result.nEM0;
                    } else if(simon::isHADCH(nextpart)){
                        ++result.nHADCH;
                    } else if(simon::isHAD0(nextpart)){
                        ++result.nHAD0;
                    } else {
                        printf("WARNING: Unknown pdgid\n");
                        printf("pdgid: %u charge: %d\n", std::abs(part->pdgId()), part->charge());
                        printf("pt: %f eta: %f phi: %f\n", part->pt(), part->eta(), part->phi());
                        printf("\n");
                    }
                }
            }

            if (result.nPart == 0){
                result.particles.emplace_back(
                    -1, 0, 0, 
                    -1, -1, -1,
                    13, 1
                );
                ++result.nPart;
            }
        }
    };
};

#endif
