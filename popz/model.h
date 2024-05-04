#ifndef MODEL_PopZ_H
#define MODEL_PopZ_H
#include <stdlib.h> 

typedef void (*PropensityFunc)(double **, long **, int);
#define n_species 2
#define m_reaction 7
#define w_bin 100 
#define L 10.0
#define NUMofRUNS 1
#define SimulationTime 1
#define Mode Auto//Choose from Auto or Manual 

struct Species{
    int index;
    char name[100];
    int initialValue;
    double diffusionRate;
};

struct Reaction{
    int type; 
    double rateConstant;
    int reactant; 
    PropensityFunc calculatePropensity;
};


struct Species* createSpeciesArray();
struct Reaction* createReactionArray();
void prop3(double **a,long **species_population,int rulenum);
void prop5(double **a,long **species_population,int rulenum);
void reactionChange(int i_reaction,int j_bin, long **species_population,long *total_population);
#endif