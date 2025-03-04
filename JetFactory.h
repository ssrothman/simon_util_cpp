#ifndef SROTHMAN_SIMONTOOLS_JETFACTORY_H
#define SROTHMAN_SIMONTOOLS_JETFACTORY_H

#include "jet.h"
#include <Eigen/Dense>
#include <random>

class JetFactory {
public:
    JetFactory();

    void initialize();

    void makeJet(simon::jet& J, const int nPart);

    void makeTransferJet(const simon::jet& J1, simon::jet& J2, Eigen::MatrixXd& tmat);

    void makeMatchedVec(const simon::jet& J, std::vector<bool>& matched, double prob);

private:
    std::default_random_engine rng;
    std::normal_distribution<double> norm;
    std::gamma_distribution<double> gamma;
    std::uniform_real_distribution<double> unif;
};

#endif
