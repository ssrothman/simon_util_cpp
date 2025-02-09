#ifndef SIMONTOOLS_PRINTPART_H
#define SIMONTOOLS_PRINTPART_H

#include <stdio.h>
#include "jet.h"

namespace simon{
    inline void printPart(const particle& part){
        printf("particle: (%0.2f, %0.2f, %0.2f)\n", part.pt, part.eta, part.phi);
        printf("\tcharge: %d\n", part.charge);
        printf("\tpdgid: %d\n", part.pdgid);
        printf("\tvtx: (%0.2f, %0.2f, %0.2f)\n", part.vtx_x, part.vtx_y, part.vtx_z);
        printf("\tdxy: %0.2f\n", part.dxy);
        printf("\tdz: %0.2f\n", part.dz);
        printf("\tfromPV: %d\n", part.fromPV);
        printf("\tpuppiweight: %0.2f\n", part.puppiweight);
    }
};

#endif
