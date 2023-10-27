#ifndef SIMON_UTIL_VECND_H
#define SIMON_UTIL_VECND_H

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

    vecND() : dim_(0), N_(0), size_(0), vec_() {}

    vecND(unsigned N, unsigned dim) : dim_(dim), N_(N),
                                      size_(getsize_()),
                                      vec_(size_){}
    
    vecND(unsigned N, unsigned dim, valuetype fillval) : dim_(dim), N_(N),
                                                 size_(getsize_()),
                                                 vec_(size_, fillval){}


    template <typename I, typename R=unsigned>
    valuetype& at(const std::vector<I>& ord, R* idxout = nullptr){
      if(ord.size() != dim_){
        throw std::invalid_argument("trying to index vecND with the wrong dimensionality");
      }
      size_t idx = idx_(ord);
      if(idxout){
          *idxout = idx;
      }
      return vec_.at(idx);
    }

    template <typename I, typename R=unsigned>
    const valuetype& at(const std::vector<I>& ord, R* idxout = nullptr) const{
      if(ord.size() != dim_){
        throw std::invalid_argument("trying to index vecND with the wrong dimensionality");
      }
      size_t idx = idx_(ord);
      if(idxout){
          *idxout = idx;
      }
      return vec_.at(idx);
    }

    //only should be used if you know what you're doing
    template <typename I>
    valuetype& at(const I i){
        return vec_.at(i);
    }

    template <typename I>
    const valuetype& at(const I i) const{
        return vec_.at(i);
    }

    std::vector<unsigned> ord0() const{
        if constexpr(vectype == type::NODIAG){
            return ord0_nodiag(dim_);
        } else {
            return ord0_full(dim_);
        }
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

  bool iterate(std::vector<unsigned>& ord) const{
        if constexpr(vectype==type::FULL){
            return iterate_full(dim_, ord, N_);
        } else if constexpr(vectype == type::SYM){
            return iterate_wdiag(dim_, ord, N_);
        } else if constexpr(vectype == type::NODIAG){
            return iterate_nodiag(dim_, ord, N_);
        }
  }

  const std::vector<valuetype>& data() const{
    return vec_;
  }

  private:
    unsigned dim_, N_;
    size_t size_;
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
    size_t idx_(const std::vector<I>& ord) const{
        if constexpr(vectype==type::FULL){
            size_t result=0;
            for(size_t i=0; i<dim_; ++i){
                result += ord[i] * intPow<size_t>(N_, dim_-i-1);
            }
            return result;
        } else {
            std::vector<I> sorted(ord);
            std::sort(sorted.begin(), sorted.end());

            if constexpr(vectype==type::NODIAG){
                if(dim_==2){
                    return size_ - tri(N_-sorted[0]-1) + (sorted[1] - sorted[0] - 1);
                }
            }
            throw std::logic_error("Not Implemented!");
        }
    }
};


typedef vecND<double, type::FULL> fullvec;;
typedef vecND<double, type::SYM> symvec;
typedef vecND<double, type::NODIAG> nodiagvec;
typedef vecND<unsigned, type::FULL> ufullvec;
typedef vecND<unsigned, type::SYM> usymvec;
typedef vecND<unsigned, type::NODIAG> unodiagvec;
typedef vecND<int, type::FULL>iufullvec;
typedef vecND<int, type::SYM> isymvec;
typedef vecND<int, type::NODIAG> inodiagvec;

};

template <typename Value>
class ArbitraryMatrix {
    public:
        ArbitraryMatrix() : dims_(), values_() {setoffsets();}

        ArbitraryMatrix(const std::vector<unsigned>& dims, const Value& fillval){
            dims_ = dims;
            size_ = 1;
            for(auto& d : dims_){
                size_ *= d;
            }
            values_ = std::vector<Value>(size_, fillval);

            setoffsets();
        }

        ArbitraryMatrix(const std::vector<unsigned>& dims) : 
            ArbitraryMatrix(dims, 0) {}

        ArbitraryMatrix(unsigned dimSize, unsigned nDim, const Value& fillval){
            dims_ = std::vector<unsigned>(nDim, dimSize);
            size_ = intPow(dimSize, nDim);
            values_ = std::vector<Value>(size_, fillval);
            setoffsets();
        }

        ArbitraryMatrix(unsigned dimSize, unsigned nDim) : 
            ArbitraryMatrix(dimSize, nDim, 0) {}

        const Value& at(const std::vector<unsigned>& ord) const {
            size_t idx = get_idx(ord);
            return values_.at(idx);
        }

        Value& at(const std::vector<unsigned>& ord){
            size_t idx = get_idx(ord);
            return values_.at(idx);
        }

        std::vector<unsigned> ord0() const{
            return std::vector<unsigned>(dims_.size(), 0);
        }

        bool iterate(std::vector<unsigned>& ord) const{
            for(int i=ord.size()-1; i>=0; --i){
                if(ord[i] < dims_[i]-1){
                    ++ord[i];
                    return true;
                } else {
                    ord[i] = 0;
                }
            }
            return false;
        }

        size_t get_idx(const std::vector<unsigned>& ord) const {
            size_t ans=0;
            for(unsigned i=0; i<ord.size(); ++i){
                ans += ord[i] * offsets_[i];
            }
            return ans;
        }

        std::vector<Value>& data(){
            return values_;
        }

        const std::vector<Value>& data() const{
            return values_;
        }
    private:
        std::vector<unsigned> dims_;
        std::vector<size_t> offsets_;
        std::vector<Value> values_;
        size_t size_;

        void setoffsets(){
            offsets_.resize(dims_.size(), 1);
            size_t offset = 1;
            for(int i=dims_.size()-1; i>=0; --i){
                offsets_[i] = offset;
                offset *= dims_[i];
            }
        }
};

#endif

