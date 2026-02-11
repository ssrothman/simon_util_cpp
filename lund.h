#ifndef SIMONTOOLS_LUND_H
#define SIMONTOOLS_LUND_H

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"

#include "jet.h"

namespace simon {

// Code below from F. Dreyer
// https://github.com/rappoccio/fastjet-tutorial/blob/master/lund-jet-example/lund.hh
// by way of Jennifer Roloff 
// https://gitlab.com/jroloff/leshouches2023/-/blob/main/MC_LESHOUCHES23.cc?ref_type=heads#L701
struct Declustering {
    // the (sub)jet, and its two declustering parts, for subsequent use
    fastjet::PseudoJet jj, j1, j2;
    // variables of the (sub)jet about to be declustered
    double pt, m;
    // properties of the declustering; NB kt is the _relative_ kt
    // = pt2 * delta_R (ordering is pt2 < pt1)
    double pt1, pt2, delta_R, z, kt, varphi;
};

/*
 This struct holds information about the hardest splitting 
 in the jet C/A declustering sequence, 
 as defined for the primary Lund plane.

 The diagram looks like
 
         4
        /
       /
      2  
     / \ 
    /   \
   1     5
    \ 
     \
      3
*/
struct SplittingInfo { 
    double pt1, pt2, pt3, pt4, pt5;
    double delta_R23, z23, kt23, phi23;
    double delta_R45, z45, kt45, phi45;
    double deltaPsi;
};

inline std::vector<Declustering>  jet_declusterings(const fastjet::PseudoJet & jet_in){
    fastjet::Recluster rc(fastjet::cambridge_algorithm);
    fastjet::PseudoJet j = rc.result(jet_in);

    std::vector<Declustering > result;
    fastjet::PseudoJet jj, j1, j2;
    jj = j;

    while (jj.has_parents(j1,j2)) {
        Declustering declust;

        // make sure j1 is always harder branch 
        if (j1.pt2() < j2.pt2()) {
            fastjet::PseudoJet jTemp;
            jTemp = j1;
            j1 = j2;
            j2 = jTemp;
        }

        // store the subjets themselves
        declust.jj   = jj;
        declust.j1   = j1;
        declust.j2   = j2;

        // get info about the jet 
        declust.pt   = jj.pt();
        declust.m    = jj.m();

        // collect info about the declustering
        declust.pt1     = j1.pt();
        declust.pt2     = j2.pt();
        declust.delta_R = j1.delta_R(j2);
        declust.z       = declust.pt2 / (declust.pt1 + declust.pt2);
        declust.kt      = j2.pt() * declust.delta_R;

        // this is now phi along the jet axis, defined in a
        // long. boost. inv. way
        declust.varphi = atan2(j1.rap()-j2.rap(), j1.delta_phi_to(j2));

        // add it to our result
        result.push_back(declust);

        // follow harder branch
        jj = j1;
    }

    return result;
}

inline void hardest_splitting_info(const simon::jet & j,
                                   struct SplittingInfo & info) {

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
    
    // Perform C/A reclustering
    auto declust = jet_declusterings(jet);

    if (declust.size() < 2){
        info.pt1 = -10;
        info.pt2 = -10;
        info.pt3 = -10;
        info.pt4 = -10;
        info.pt5 = -10;
        info.delta_R23 = -10;
        info.z23 = -10;
        info.kt23 = -10;
        info.phi23 = -10;
        info.delta_R45 = -10;
        info.z45 = -10;
        info.kt45 = -10;
        info.phi45 = -10;
        info.deltaPsi = -10;
        return;
    } else {
        const auto& declus123 = declust[0];
        const auto& declus245 = declust[1];

        info.pt1 = declus123.pt;
        info.pt2 = declus123.pt1;
        info.pt3 = declus123.pt2;

        info.pt4 = declus245.pt1;
        info.pt5 = declus245.pt2;

        info.delta_R23 = declus123.delta_R;
        info.z23 = declus123.z;
        info.kt23 = declus123.kt;
        info.phi23 = declus123.varphi;

        info.delta_R45 = declus245.delta_R;
        info.z45 = declus245.z;
        info.kt45 = declus245.kt;
        info.phi45 = declus245.varphi;

        info.deltaPsi = info.phi23 - info.phi45;
        if(info.deltaPsi < 0) {
            info.deltaPsi += 2*M_PI;
        }
        if(info.deltaPsi > 2*M_PI){
            info.deltaPsi -= 2*M_PI;
        }
    }    
}

} // namespace simon

#endif