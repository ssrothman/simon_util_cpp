#include <stdio.h>
#include <iostream>
#include <vector>
#include "vecND.h"
#include "iterating.h"

int main(){
    unsigned nPart=10;
    unsigned nDim=5;

    vecND::vecND<float, vecND::type::NODIAG> vec(nPart, nDim, 0.f);
    std::vector<unsigned> ord = vec.ord0();

    printf("vec parameters: (%u, %u, %lu)\n", vec.nPart(), vec.dim(), vec.size());

    for(size_t i=0, loop=true; loop; loop=vec.iterate(ord), ++i){
        printf("%lu: ", i);
        printOrd(ord);
        printf("\n");
        printf("getIdx() returns %lu\n", getNodiagIdx(ord, nPart, nDim));
    }
    return 0;
}
