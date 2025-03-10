#ifndef SIMON_UTIL_COMBINATORICS_H
#define SIMON_UTIL_COMBINATORICS_H

#include <vector>
#include <cstddef>

namespace simon {
    struct comp{
        std::vector<unsigned> composition;
        unsigned factor;

        comp(){}

        comp(std::vector<unsigned> composition, unsigned factor):
            composition(composition), factor(factor) {}
    };

    //NB assumes ord is sorted
    size_t getNodiagIdx(const std::vector<unsigned>& ord,
                        const unsigned N, const unsigned dim,
                        const unsigned istart=0, const unsigned off=0);

    typedef std::vector<std::vector<comp>> comp_t;

    int fact(int n);

    void fillCompositions(const int n, comp_t& out);

    template <typename T, typename R>
    T intPow(const T a, const R b) {
      //mine, very stupid
      //should be upgraded to to square multiply
      //also, assumes b > 0
      T result = 1;
      for (unsigned i = 0; i < (unsigned)b; ++i) {
        result *= a;
      }
      return result;
    }

    size_t choose(int n, int k);

    size_t simp(const size_t nPart, const size_t order);

    template <typename T>
    size_t tri(const T N){
        return N*(N+1)/2;
    }
}

#endif
