#ifndef VECND_H
#define VECND_H

#include <vector>
#include "combinatorics.h"
#include "iterating.h"
#include <algorithm>
#include <cstddef>
#include <stdio.h>

namespace vecND{

enum type{
    FULL = 0,
    SYM = 1,
    NODIAG = 2
};

template <class valuetype, enum type vectype>
class vecND{
  public:

    vecND(unsigned N, unsigned dim) : dim_(dim), N_(N),
                                      size_(getsize_()),
                                      vec_(size_){}
    
    vecND(unsigned N, unsigned dim, valuetype fillval) : dim_(dim), N_(N),
                                                 size_(getsize_()),
                                                 vec_(size_, fillval){}


    template <typename I>
    valuetype& at(const std::vector<I>& ord){
      if(ord.size() != dim_){
        throw std::invalid_argument("trying to index vecND with the wrong dimensionality");
      }
      size_t idx = idx_(ord);
      printf("\t%lu\n", idx);
      return vec_.at(idx);
    }

    //only should be used if you know what you're doing
    valuetype& at(const size_t i){
        return vec_.at(i);
    }

    std::vector<unsigned> ord0(){
        std::vector<unsigned> ans(dim_, 0u);
        if constexpr(vectype == type::NODIAG){
            for(unsigned i=0; i<dim_; ++i){
                ans[i] = i;
            }
        }
        return ans;
    }

  unsigned nPart() const {
    return N_;
  }

  unsigned dim() const {
      return dim_;
  }

  size_t size() const{
      return size_;
  }

  bool iterate(std::vector<unsigned>& ord){
        if constexpr(vectype==type::FULL){
            return iterate_full(dim_, ord, N_);
        } else if constexpr(vectype == type::SYM){
            return iterate_wdiag(dim_, ord, N_);
        } else if constexpr(vectype == type::NODIAG){
            return iterate_nodiag(dim_, ord, N_);
        }
  }

  private:
    const unsigned dim_, N_;
    const size_t size_;
    std::vector<valuetype> vec_;

    size_t getsize_() const{
        if constexpr(vectype==type::FULL){
            return intPow(N_, dim_);
        } else if constexpr(vectype==type::SYM){
            return choose(N_ + dim_ -1, dim_);
        } else if constexpr(vectype==type::NODIAG){
            return choose(N_, dim_);
        }
    }

    //NB assumes ord is the right size
    //sorts a copy of ord
    template <typename I>
    size_t idx_(const std::vector<I>& ord){
        if constexpr(vectype==type::FULL){
            size_t result=0;
            for(size_t i=0; i<dim_; ++i){
                result += ord[i] * intPow<size_t>(N_, dim_-i-1);
            }
            return result;
        } else {
            throw std::logic_error("Not Implemented!");
            std::vector<I> sorted(ord);
            std::sort(sorted.begin(), sorted.end());

            if constexpr(vectype==type::SYM){
                return 0;
            } else if constexpr(vectype==type::NODIAG){
                return 0;
            }
        }
    }
};

};
#endif

