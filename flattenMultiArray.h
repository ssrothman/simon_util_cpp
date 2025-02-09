#ifndef SIMONTOOLS_FLATTENMULTIARRAY_H
#define SIMONTOOLS_FLATTENMULTIARRAY_H

#include <boost/multi_array.hpp>
#include <vector>

namespace simon{
    template <typename T1, unsigned long Ndim, typename T2>
    size_t flattenMultiArray(const boost::multi_array<T1, Ndim>& arr, std::vector<T2>& flat){
        size_t N = arr.num_elements();
        flat.reserve(flat.size() + N);
        for(auto it=arr.data(); it!=arr.data()+N; ++it){
            flat.emplace_back(*it);
        }
        return N;
    }
};

#endif
