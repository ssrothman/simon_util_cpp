#ifndef SIMONTOOLS_LUND_H
#define SIMONTOOLS_LUND_H

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "Math/Boost.h"
#include "Math/DisplacementVector3D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzVector.h"
#include "jet.h"

namespace simon {

/*
    Struct for 1->23 splitting, where 1 is the mother and 2,3 are the daughters.
    Ordered such that pt2 > pt3
*/
struct SingleSplittingInfo {
    math::PtEtaPhiMLorentzVector p1, p2, p3;
    double deltaR, z, kt;
    int pdgId1, pdgId2, pdgId3;

    // angle p2 and p3 make in the (y, phi) plane
    double psi_type1() const {
        double deltaPhi = reco::deltaPhi(p2.phi(), p3.phi());
        double deltaY = p2.Rapidity() - p3.Rapidity();
        double psi = atan2(deltaY, deltaPhi);
        if (psi < -M_PI) {
            psi += 2*M_PI;
        }
        if (psi > M_PI) {
            psi -= 2*M_PI;
        }
        return psi;
    }

    //angle p2 and p3 make in the (y, phi) plane
    //after boosting into the rest frame of the reference vector
    double psi_type2(math::PtEtaPhiMLorentzVector ref) const {
        ROOT::Math::Boost boostToRest(ref.BoostToCM());
        math::PtEtaPhiMLorentzVector p2_boosted = boostToRest(p2);
        math::PtEtaPhiMLorentzVector p3_boosted = boostToRest(p3);
        
        double deltaPhi = reco::deltaPhi(p2_boosted.phi(), p3_boosted.phi());
        double deltaY = p2_boosted.Rapidity() - p3_boosted.Rapidity();
        double psi = atan2(deltaY, deltaPhi);
        if (psi < -M_PI) {
            psi += 2*M_PI;
        }
        if (psi > M_PI) {
            psi -= 2*M_PI;
        }
        return psi;
    }

    // normalized cross product of the spatial coordinates of p2 and p3
    ROOT::Math::XYZVector psi_type3() const {
        ROOT::Math::XYZVector p2_vec(p2.Px(), p2.Py(), p2.Pz());
        ROOT::Math::XYZVector p3_vec(p3.Px(), p3.Py(), p3.Pz());
        ROOT::Math::XYZVector cross = p2_vec.Cross(p3_vec);
        
        // if the cross product is very small, return a default vector to avoid numerical issues
        if (cross.Mag2() < 1e-4) {
            return ROOT::Math::XYZVector(0,0,0);
        } else {
            return cross.Unit();
        }
    }

    // normalized cross product of the spatial coordinates of p2 and p3
    // after boosting into the rest frame of the reference vector
    ROOT::Math::XYZVector psi_type4(math::PtEtaPhiMLorentzVector ref) const {
        ROOT::Math::Boost boostToRest(ref.BoostToCM());
        math::PtEtaPhiMLorentzVector p2_boosted = boostToRest(p2);
        math::PtEtaPhiMLorentzVector p3_boosted = boostToRest(p3);
        ROOT::Math::XYZVector p2_vec(p2_boosted.Px(), p2_boosted.Py(), p2_boosted.Pz());
        ROOT::Math::XYZVector p3_vec(p3_boosted.Px(), p3_boosted.Py(), p3_boosted.Pz());
        ROOT::Math::XYZVector cross = p2_vec.Cross(p3_vec);

        // if the cross product is very small, return a default vector to avoid numerical issues
        if (cross.Mag2() < 1e-4) {
            return ROOT::Math::XYZVector(0,0,0);
        } else {
            return cross.Unit();
        }
    }

    // constructor from three genParticle objects
    SingleSplittingInfo(const reco::GenParticle& genPart1, const reco::GenParticle& genPart2, const reco::GenParticle& genPart3) :
        p1(genPart1.pt(), genPart1.eta(), genPart1.phi(), genPart1.mass()),
        p2(genPart2.pt(), genPart2.eta(), genPart2.phi(), genPart2.mass()),
        p3(genPart3.pt(), genPart3.eta(), genPart3.phi(), genPart3.mass()),
        pdgId1(genPart1.pdgId()),
        pdgId2(genPart2.pdgId()),
        pdgId3(genPart3.pdgId()) {
            deltaR = reco::deltaR(p2.eta(), p2.phi(), p3.eta(), p3.phi());
            z = p3.pt() / (p2.pt() + p3.pt());
            kt = p3.pt() * deltaR;
        }

    // constructor from three fastjet::PseudoJet objects
    SingleSplittingInfo(const fastjet::PseudoJet& pseudoJet1, const fastjet::PseudoJet& pseudoJet2, const fastjet::PseudoJet& pseudoJet3) :
        p1(pseudoJet1.pt(), pseudoJet1.eta(), pseudoJet1.phi(), pseudoJet1.m()),
        p2(pseudoJet2.pt(), pseudoJet2.eta(), pseudoJet2.phi(), pseudoJet2.m()),
        p3(pseudoJet3.pt(), pseudoJet3.eta(), pseudoJet3.phi(), pseudoJet3.m()),
        pdgId1(0), pdgId2(0), pdgId3(0) {
            deltaR = reco::deltaR(p2.eta(), p2.phi(), p3.eta(), p3.phi());
            z = p3.pt() / (p2.pt() + p3.pt());
            kt = p3.pt() * deltaR;
        }

    //empty constructor
    SingleSplittingInfo() :
        p1(0,0,0,0),
        p2(0,0,0,0),
        p3(0,0,0,0),
        deltaR(0),
        z(0),
        kt(0),
        pdgId1(0),
        pdgId2(0),
        pdgId3(0) {}
};

/*
Represents 1->23 and subsequent 4->56
sorted such that 
pt2 > pt3 and pt5 > pt6
*/
struct DoubleSplittingInfo{
    SingleSplittingInfo split123;
    SingleSplittingInfo split456;

    double deltaPsi_type1() const {
        return reco::deltaPhi(split123.psi_type1(), split456.psi_type1());
    }

    double deltaPsi_type2() const {
        return deltaPsi_type2(split456.p1);
    }

    double deltaPsi_type2(math::PtEtaPhiMLorentzVector ref) const {
        return reco::deltaPhi(split123.psi_type2(ref), split456.psi_type2(ref));
    }

    double deltaPsi_type3() const {
        ROOT::Math::XYZVector cross1 = split123.psi_type3();
        ROOT::Math::XYZVector cross2 = split456.psi_type3();
        double dot = cross1.Dot(cross2);
        return acos(dot);
    }

    double deltaPsi_type4(math::PtEtaPhiMLorentzVector ref) const {
        ROOT::Math::XYZVector cross1 = split123.psi_type4(ref);
        ROOT::Math::XYZVector cross2 = split456.psi_type4(ref);
        double dot = cross1.Dot(cross2);
        return acos(dot);
    }

    double deltaPsi_type4() const {
        return deltaPsi_type4(split456.p1);
    }

    DoubleSplittingInfo(const reco::GenParticle& genPart1, const reco::GenParticle& genPart2, const reco::GenParticle& genPart3,
                        const reco::GenParticle& genPart4, const reco::GenParticle& genPart5, const reco::GenParticle& genPart6) :
        split123(genPart1, genPart2, genPart3),
        split456(genPart4, genPart5, genPart6) {}

    DoubleSplittingInfo(const fastjet::PseudoJet& pseudoJet1, const fastjet::PseudoJet& pseudoJet2, const fastjet::PseudoJet& pseudoJet3,
                        const fastjet::PseudoJet& pseudoJet4, const fastjet::PseudoJet& pseudoJet5, const fastjet::PseudoJet& pseudoJet6) :
        split123(pseudoJet1, pseudoJet2, pseudoJet3),
        split456(pseudoJet4, pseudoJet5, pseudoJet6) {}

    DoubleSplittingInfo() :
        split123(),
        split456() {}
};

inline void LundDeclustered(const simon::jet & j, 
                            const bool hardSide_,
                            double zcut1,
                            double zcut2,
                            std::vector<simon::DoubleSplittingInfo>& result){
    // Create constituents from particles
    std::vector<fastjet::PseudoJet> constituents;
    constituents.reserve(j.nPart);
    
    for(const auto& particle : j.particles) {
        double px = particle.pt * cos(particle.phi);
        double py = particle.pt * sin(particle.phi);
        double pz = particle.pt * sinh(particle.eta);
        double E = sqrt(px*px + py*py + pz*pz);
        
        constituents.push_back(fastjet::PseudoJet(px, py, pz, E));
    }

    // Build jet from constituents
    // Radius arbitrarily large to ensure all constituents are clustered into one jet
    fastjet::JetDefinition jd(fastjet::antikt_algorithm, fastjet::JetDefinition::max_allowable_R);
    fastjet::ClusterSequence cs(constituents, jd);
    
    // Double-check that we get exactly one jet from the clustering
    // Otherwise something has gone wrong
    if (cs.inclusive_jets().size() != 1){
        printf("There were %lu constituents\n", constituents.size());
        printf("Resulting in %lu jets\n", cs.inclusive_jets().size());
        printf("The 0th constituent has pt %f\n", constituents[0].pt());
        throw std::runtime_error("Expected exactly one jet from clustering, got " + std::to_string(cs.inclusive_jets().size()));
    }

    // Get the jet
    fastjet::PseudoJet jet = cs.inclusive_jets()[0];

    //recluster with C/A
    fastjet::Recluster rc(fastjet::cambridge_algorithm);
    fastjet::PseudoJet jCA = rc.result(jet);

    // walk the tree to find the highest-kt splitting
    fastjet::PseudoJet jj, j1, j2;

    fastjet::PseudoJet j_first; 
    double kt_best = -1.0;

    jj = jCA;;
    while (jj.has_parents(j1, j2)) {
        SingleSplittingInfo current_split(jj, j1, j2);
        if (current_split.z > zcut1 && current_split.kt > kt_best) {
            kt_best = current_split.kt;
            j_first = jj;
        }
        // follow harder branch
        jj = j1;
    }

    if (kt_best < 0.0) {
        // no splittings found
        return;
    }

    // now find the next splitting down the hard/soft branch
    j_first.has_parents(j1, j2);
    if (hardSide_){
        jj = j1;
    } else {
        jj = j2;
    }

    fastjet::PseudoJet j_second;
    double kt_second = -1.0;

    while (jj.has_parents(j1, j2)) {
        SingleSplittingInfo current_split(jj, j1, j2);
        if (current_split.z > zcut2 && current_split.kt > kt_second) {
            kt_second = current_split.kt;
            j_second = jj;
        }

        jj = j1;
    }
    if (kt_second < 0.0) {
        // no second splitting found
        return;
    }

    fastjet::PseudoJet j_first1, j_first2;
    fastjet::PseudoJet j_second1, j_second2;
    j_first.has_parents(j_first1, j_first2);
    j_second.has_parents(j_second1, j_second2);

    result.emplace_back(
        j_first, j_first1, j_first2,
        j_second, j_second1, j_second2
    );
}

} // namespace simon

#endif