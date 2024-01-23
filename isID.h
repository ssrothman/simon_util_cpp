#ifndef SIMONTOOLS_ISID_H
#define SIMONTOOLS_ISID_H

template <typename P>
inline bool isEM0(const P& partptr){
    unsigned pdgid = std::abs(partptr->pdgId());
    return pdgid==22;
}

template <typename P>
inline bool isELE(const P& partptr){
    unsigned pdgid = std::abs(partptr->pdgId());
    return pdgid==11;
}

template <typename P>
inline bool isMU(const P& partptr){
    unsigned pdgid = std::abs(partptr->pdgId());
    return pdgid==13;
}

template <typename P>
inline bool isHADCH(const P& partptr){
    unsigned pdgid = std::abs(partptr->pdgId());
    return pdgid>=100 && partptr->charge()!=0;
}

template <typename P>
inline bool isHAD0(const P& partptr){
    unsigned pdgid = std::abs(partptr->pdgId());
    return pdgid>=100 && partptr->charge()==0;
}

#include "jets.h"
inline bool isEM0(const particle& part){
    return part.pdgid==22;
}

inline bool isELE(const particle& part){
    return part.pdgid==11;
}

inline bool isMU(const particle& part){
    return part.pdgid==13;
}

inline bool isHADCH(const particle& part){
    return part.pdgid>=100 && part.charge!=0;
}

inline bool isHAD0(const particle& part){
    return part.pdgid>=100 && part.charge==0;
}

#endif
