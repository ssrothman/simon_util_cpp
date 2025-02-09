#include "computeJetMass.h"
#include <Math/Vector4D.h>

void simon::computeJetMass(jet& thejet){
    ROOT::Math::PtEtaPhiMVector jet_p4;
    for (const auto& part: thejet.particles){
        jet_p4 += ROOT::Math::PtEtaPhiMVector(part.pt, part.eta, part.phi, 0);
    }
    thejet.mass = jet_p4.M();
}
