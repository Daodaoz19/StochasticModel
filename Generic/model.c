/*
model.c file provied by user
*/
#include "model.h" 
#include <stdlib.h>


void initial(long ***species_population) {
    for (int j = 0; j < W; j++) {
        (*species_population)[0][j] = 2; 
        (*species_population)[1][j] = 1; 
    }
}

double prop3(int numbin, long **species_population) {

    double prop = 0;
    for (int j = 0; j < numbin; j++) {
        prop +=  species_population[0][j] * ( species_population[0][j] - 1) * species_population[1][j];
    }
    double rateConstant = 1 / (W * W);
    prop *= rateConstant * numbin * numbin;
    return prop;
}

//PropensityFunc propensityFuncs[M] = {NULL, NULL, NULL, prop3};


