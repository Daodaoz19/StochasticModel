/*
Create file contain initial population for each species base on initial()
*/
#include <stdio.h>
#include <stdlib.h> 
#include "model.h"


int main() {

    long **species_population;
    species_population = (long **)malloc(N * sizeof(long *));
    for (int n = 0; n < N; n++) {
        species_population[n] = (long *)malloc(W * sizeof(long));
    }

    initial(&species_population);

    FILE *file = fopen("init_population", "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < W; j++) {
            fprintf(file, "%ld\t", species_population[i][j]); 
        }
        fprintf(file, "\n"); 
    }

    fclose(file); 

    for (int n = 0; n < N; n++) {
        free(species_population[n]);
    }
    free(species_population);

    return 0;
}


