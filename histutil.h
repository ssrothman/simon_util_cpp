#ifndef SROTHMAN_SIMONTOOLS_HISTUTIL_H
#define SROTHMAN_SIMONTOOLS_HISTUTIL_H

#include <boost/histogram.hpp>

template <typename T>
unsigned getIndex(const T& val, const boost::histogram::axis::variable<T>& ax){
    return static_cast<unsigned>(ax.index(val) + 1);
}

#endif
