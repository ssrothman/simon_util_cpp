#include <stdio.h>
#include <iostream>
#include <vector>
#include "vecND.h"
#include "iterating.h"

int main(){
    unsigned nPart=5;
    unsigned nDim=2;

    std::vector<unsigned> nadj{5, 3, 4};

    std::vector<unsigned> ord = ord0_full(3);

    do{
        printOrd(ord);
        printf("\n");
    } while(iterate_awkward(nadj, ord));
    return 0;
}
