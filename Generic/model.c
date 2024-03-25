/*
model.c file provied by user
*/
#include "model.h" 
#include <stdlib.h>
#include <stdio.h>


void initial(long ***species_population) {
    for (int j = 0; j < W; j++) {
        (*species_population)[0][j] = 2; 
        (*species_population)[1][j] = 1; 
    }
}

void prop3(double **a,long **species_population,int rulenum) {  
    for (int i = 0; i < W; i++) {
        long x = species_population[0][i]; 
        long y = species_population[1][i];
        a[rulenum][i] = x * (x - 1) * y; 
       } 
}

//reaction i will hapen in j bin
void stateChange(int i,int j, long **species_population,long *tot) {
   switch(i) {
    case 1:
        species_population[0][j]++;
        tot[0]++; 
        break;
    case 2:
        species_population[0][j]--;
        tot[0]--; 
        break;
    case 3:
        species_population[1][j]++;
        tot[1]++; 
        break;
    case 4:
        species_population[0][j]++;
        species_population[1][j]--;
        tot[0]++; 
        tot[1]--; 
        break;
    default:
        exit(1);
        break;
    }
}


//PropensityFunc propensityFuncs[M] = {NULL, NULL, NULL, ,



















