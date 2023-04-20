#ifndef SIMON_UTIL_ITERATING_H
#define SIMON_UTIL_ITERATING_H

#include <vector>
#include <iostream>

template <typename T>
void printOrd(const std::vector<T> ord) {
  std::cout << "(";
  unsigned int i = 0;
  for (i = 0; i < ord.size() - 1; ++i) {
    std::cout << ord[i] << ", ";
  }
  std::cout << ord[ord.size() - 1] << ")";
}

inline std::vector<unsigned> ord0_full(unsigned dim){
    return std::vector<unsigned>(dim, 0u);
}

inline std::vector<unsigned> ord0_nodiag(unsigned dim){
    std::vector<unsigned> ans(dim);
    for(unsigned i=0; i<dim; ++i){
        ans[i] = i;
    }
    return ans;
}

template <typename T, typename R>
bool iterate_awkward(const std::vector<T>& dims, std::vector<R>& ord){
    for(int dim=dims.size()-1; dim>=0; --dim){
        if(ord[dim] < dims[dim]-1){
            ++ord[dim];
            break;
        } else {
            ord[dim] = 0;
            if(dim==0){
                std::fill(ord.begin(), ord.end(), 0);
                return false;
            }
        }
    }
    return true;
}

template <typename T, typename R>
bool iterate_nodiag(const T dims, std::vector<R>& ordinates, const T nParts) {
  // iterate over dimensions in reverse...
  for (int dim = dims - 1; dim >= 0; --dim) {
    if (ordinates[dim] < nParts - dims + dim) {
      ++ordinates[dim];
      for (int d2 = dim + 1; d2 < dims; ++d2) {
        ordinates[d2] = ordinates[d2 - 1] + 1;
      }
      return true;
    }
    if (dim > 0) {
      ordinates[dim] = ordinates[dim - 1] + 1;
    } else {
      ordinates[dim] = 0;
      if(dim==0){
          std::fill(ordinates.begin(), ordinates.end(), 0);
          return false;
      }
    }
  }
  return true;
}

template <typename T, typename R>
bool iterate_wdiag(const T dims, std::vector<R>& ordinates, const T nParts) {
  // iterate over dimensions in reverse...
  for (int dim = dims - 1; dim >= 0; --dim) {
    if (ordinates[dim] < nParts-1) {
      ++ordinates[dim];
      for (int d2 = dim + 1; d2 < dims; ++d2) {
        ordinates[d2] = ordinates[d2 - 1];
      }
      return true;
    }
    if (dim > 0) {
      ordinates[dim] = ordinates[dim - 1];
    } else {
      ordinates[dim] = 0;
      if(dim==0){
          std::fill(ordinates.begin(), ordinates.end(), 0);
          return false;
      }
    }
  }
  return true;
}



size_t getIndex(const std::vector<int>& ordinates, const int nPart);

template <typename T>
bool iterate_full(const T dims, std::vector<T>& ordinates, const T nParts) {
  // iterate over dimensions in reverse...
  for (int dim = dims - 1; dim >= 0; --dim) {
    if (ordinates[dim] < nParts-1) {
      ++ordinates[dim];
      break;
    } else {
      ordinates[dim] = 0;
      if(dim==0){
          std::fill(ordinates.begin(), ordinates.end(), 0);
          return false;
      }
    }
  }
  return true;
}

#endif
