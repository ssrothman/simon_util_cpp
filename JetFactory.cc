#include "JetFactory.h"

JetFactory::JetFactory(){
    initialize();
}

void JetFactory::initialize(){
    rng.seed(0);
    norm = std::normal_distribution<double>(0, 0.4);
    gamma = std::gamma_distribution<double>(2, 2);
    unif = std::uniform_real_distribution<double>(0, 1);
}

void JetFactory::makeJet(simon::jet& J, const int nPart){
    J.nPart = nPart;

    J.sumpt = 0;
    J.pt = 0;
    J.rawpt = 0;

    J.eta = 0;
    J.phi = 0;

    J.particles.clear();
    for(int i=0; i<nPart; ++i){
        double pt = gamma(rng);
        double phi = norm(rng);
        double eta = norm(rng);
        J.particles.emplace_back(pt, phi, eta);
        J.sumpt += pt;
        J.pt += pt;
        J.rawpt += pt;
    }
}

void JetFactory::makeTransferJet(const simon::jet& J1, simon::jet& J2, Eigen::MatrixXd& tmat){
    J2.nPart = 0;

    J2.sumpt = 0;
    J2.pt = 0;
    J2.rawpt = 0;

    J2.eta = J1.eta;
    J2.phi = J1.phi;

    tmat.resize(J1.nPart, J1.nPart);
    for(unsigned i=0; i<J1.nPart; ++i){
        for(unsigned j=0; j<J1.nPart; ++j){
            tmat(i, j) = 0;
        }
    }

    for (unsigned i=0; i<J1.nPart; ++i){
        double pt = J1.particles[i].pt;
        double phi = J1.particles[i].phi;
        double eta = J1.particles[i].eta;

        ++J2.nPart;
        J2.particles.emplace_back(pt, eta, phi);
        J2.sumpt += pt;
        J2.pt += pt;
        J2.rawpt += pt;

        tmat(i, i) = 1;
    }
}

void JetFactory::makeMatchedVec(const simon::jet& J, std::vector<bool>& matched, double prob){
    matched.clear();
    for(unsigned i=0; i<J.nPart; ++i){
        matched.push_back(unif(rng) < prob);
    }
}
