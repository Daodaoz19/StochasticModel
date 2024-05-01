#ifndef MODEL_PopZ_H
#define MODEL_PopZ_H
#include <stdlib.h> 

typedef void (*PropensityFunc)(double **, long **, int);
#define n_species 12
#define m_reaction 7
#define w_bin 100 
#define L 10.0

struct Species{
    int index;
    char name[100];
    int initialValue;
    double diffusionRate;
};

struct Reaction{
    int type; 
    double rateConstant;
    int* reactant; 
    PropensityFunc calculatePropensity;
};


struct Species* createSpeciesArray();
struct Reaction* createReactionArray();
void prop3(double **a,long **species_population,int rulenum);
void stateChange(int i,int j, long **species_population,long *tot);
#endif