#ifndef SIMONTOOLS_COPYMULTIARRAY_H
#define SIMONTOOLS_COPYMULTIARRAY_H

#include <boost/multi_array.hpp>

namespace simon{
    template <typename T1, unsigned long Ndim, typename T2>
    void copyMultiArray(const boost::multi_array<T1, Ndim>& source,
                        boost::multi_array<T2, Ndim>& dest){
        std::vector<size_t> shape(Ndim);
        const long unsigned int* shapearr = source.shape();
        for(unsigned i=0; i<=Ndim; ++i){
            shape[i] = shapearr[i];
        }
        dest.resize(shape);
        std::copy(source.data(), source.data() + source.num_elements(), dest.data());
    }

    template <typename T1, typename T2>
    void copyMultiArray(const std::vector<T1>& source, 
                        std::vector<T2>& dest){
        dest.resize(source.size());
        std::copy(source.begin(), source.end(), dest.begin());
    }
};

#endif
