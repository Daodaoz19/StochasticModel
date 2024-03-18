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

double prop3(double *a, long **species_population) { 
    *a = 0; 
    for (int i = 0; i < W; i++) {
        long x = species_population[0][W - 1]; 
        long y = species_population[1][W - 1];
        long temp = x * (x - 1) * y;
        *a += temp;
    }
    double rateConstant = 1 / (W * W);
    *a *= rateConstant * W * W;
    return *a;
}

void stateChangeMinus(int i, int j, long **species_population) {
    species_population[i][j] --;
}

void stateChangeAdd(int i, int j, long **species_population) {
    species_population[i][j] ++;
}

//PropensityFunc propensityFuncs[M] = {NULL, NULL, NULL, prop3};


