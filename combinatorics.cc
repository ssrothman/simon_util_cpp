#include "combinatorics.h"

#include <stdio.h>
#include <math.h>
#include <algorithm>

#include <boost/container_hash/hash.hpp>

int simon::fact(int n) {
  int result = 1;
  for (int i = 2; i < n + 1; ++i) {
    result *= i;
  }
  return result;
}

static unsigned getFactor(const unsigned N, const std::vector<unsigned>& composition){
    unsigned result = simon::fact(N);
    for(const unsigned& i : composition){
        result /= simon::fact(i);
    }
    return result;
}

void simon::fillCompositions(const int n, comp_t& out) {
  out.clear();
  for (int i = 0; i < n; ++i) {
    std::vector<comp> next;
    out.push_back(next);
  }

  std::vector<int> q(n + 1);
  int* a = q.data();

  int k, y, x, l;
  for (int q = 0; q < n + 1; ++q) {
    a[q] = 0;
  }
  k = 1;
  y = n - 1;
  while (k != 0) {
    x = a[k - 1] + 1;
    k -= 1;
    while (2 * x <= y) {
      a[k] = x;
      y -= x;
      k += 1;
    }
    l = k + 1;
    while (x <= y) {
      a[k] = x;
      a[l] = y;
      //yield a[:k + 2];
      do {
        comp next;
        next.composition.reserve(k + 1);
        for (int q = 0; q < k + 2; ++q) {
          next.composition.push_back(a[q]);
        }
        next.factor = getFactor(n, next.composition);
        out[k + 1].push_back(next);
      } while (std::next_permutation(a, a + k + 2));
      x += 1;
      y -= 1;
    }
    a[k] = x + y;
    y = x + y - 1;
    //yield a[:k + 1];
    do {
      comp next;
      next.composition.reserve(k);
      for (int q = 0; q < k + 1; ++q) {
        next.composition.push_back(a[q]);
      }
      next.factor = getFactor(n, next.composition);
      out[k].push_back(next);
    } while (std::next_permutation(a, a + k + 1));
  }
}

size_t simon::choose(int n, int k) {
    /*static std::unordered_map<std::pair<int, int>, size_t, boost::hash<std::pair<int,int>>> cache;*/
  if (k > n)
    return 0;
  if (k * 2 > n)
    k = n - k;
  if (k == 0)
    return 1;

  /*std::pair<int, int> key(n, k);
  auto it = cache.find(key);
  if(it != cache.end()){
    return it->second;
  }*/

  size_t result = n;
  for (int i = 2; i <= k; ++i) {
    result *= (n - i + 1);
    result /= i;
  }

  /*cache[key] = result;*/
  return result;
}

size_t simon::simp(const size_t nPart, const size_t order){
  size_t result = 1;
  for(unsigned i=0; i<(unsigned)order; ++i){
    result *= (nPart - i);
  }
  result/=fact(order);
  printf("simp(%lu, %lu) = %lu\n", nPart, order, result);
  return result;
}

size_t simon::getNodiagIdx(const std::vector<unsigned>& ord,
                    const unsigned N, const unsigned dim,
                    const unsigned istart, const unsigned off){
    if(dim == 0){
        return 0;
    } else if(dim == 1){
        return ord[istart]-off;
    } else {//if(dim == 2){
        size_t s1 = choose(N-off, dim) - choose(N-ord[istart], dim);
        size_t s2 = getNodiagIdx(ord, N, dim-1, istart+1, ord[istart]+1);
        return s1+s2;
    }
}
