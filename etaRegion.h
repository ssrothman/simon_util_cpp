#ifndef SIMON_TOOLS_ETA_REGION_H
#define SIMON_TOOLS_ETA_REGION_H

namespace simon{
    inline int getEtaRegion(const double& eta, const std::vector<double>& boundaries){
        for(unsigned i=0; i<boundaries.size(); ++i){
            if(std::abs(eta) < boundaries[i]){
                return i-1;
            }
        }
        return boundaries.size();
    }
};

#endif
