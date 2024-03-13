/*
model.h file provied by user
*/
#ifndef MODEL_H
#define MODEL_H

#define N 2 
#define M 4 
#define W 100 
#define L 10.0
#define S 100 

extern double diffusionRates[N];
extern int reactionTypes[M];
extern double rateConstants[M];

#define DIFFUSION_RATE {5e-4,(5e-4 * 10000)}
#define REACTION_TYPE {1, 2, 1, 3}
#define RATE_CONSTANTS {1, 2, 3, 1.0/W*W}

void initial(long ***species_population);
double prop3(int numBins, long  **species_population);
//extern PropensityFunc propensityFuncs[M];

#endif