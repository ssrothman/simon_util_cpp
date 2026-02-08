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

inline double deltaPsi(std::vector< Declustering> clust){
    if(clust.size() < 2){
        return -10;
    }

    double varPhi1 = clust[0].varphi;
    double varPhi2 = clust[1].varphi;
    double deltaVarPhi = varPhi1 - varPhi2;
    
    if(deltaVarPhi < 0) {
        deltaVarPhi = 2*M_PI - deltaVarPhi;
    }
    if(deltaVarPhi > 2*M_PI){
        deltaVarPhi = deltaVarPhi - 2*M_PI;
    } 

    return deltaVarPhi;
}

inline double get_dpsi(const simon::jet & j){
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
    
    // Cluster constituents - keep ClusterSequence alive during declusterings
    fastjet::JetDefinition jd(fastjet::cambridge_algorithm, 1000.0);
    fastjet::ClusterSequence cs(constituents, jd);
    if (cs.inclusive_jets().size() != 1){
        printf("There were %lu constituents\n", constituents.size());
        printf("Resulting in %lu jets\n", cs.inclusive_jets().size());
        printf("The 0th constituent has pt %f\n", constituents[0].pt());
        throw std::runtime_error("Expected exactly one jet from clustering, got " + std::to_string(cs.inclusive_jets().size()));
    }
    fastjet::PseudoJet jet = cs.inclusive_jets()[0];
    
    // Compute declusterings while ClusterSequence is alive
    auto declust = jet_declusterings(jet);
    return deltaPsi(declust);
}

} // namespace simon

#endif