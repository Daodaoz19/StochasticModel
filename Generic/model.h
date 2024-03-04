/*
model.h
*/
#ifndef MODEL_H
#define MODEL_H

#define N 2 
#define M 4 
#define W 100 
#define L 10

extern double diffusionRates[N];
extern int reactionTypes[M];
extern double rateConstants[M];
extern int sepciesRatio[N];

#define DIFFUSION_RATE {5e-4,(5e-4 * 10000)}
#define REACTION_TYPE {1, 2, 1, 4}
#define RATE_CONSTANTS {1, 2, 3, 1.0/W*W}
#define SPECIES_RATIO {2,1}

#endif