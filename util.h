#ifndef SIMON_UTIL_UTIL_H
#define SIMON_UTIL_UTIL_H

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

#endif
