#ifndef SIMON_UTIL_UTIL_H
#define SIMON_UTIL_UTIL_H

#include <vector>
#include <string>
#include <cstdarg>

template <typename T>
void printVec(std::vector<T> v){
    for(unsigned i=0; i<v.size(); ++i){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

template <typename T>
bool uniform(std::vector<T> v){
    if (v.size() < 2) 
        return true;

    for(unsigned i=1; i<v.size(); ++i){
        if (v[i] != v[0]){
            return false;
        }
    }
    return true;
}

template <typename T>
inline T square(const T& x){
    return x*x;
}

// requires at least C++11
//from https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf/49812018#49812018
const std::string vformat(const char * const zcFormat, ...);

#endif
