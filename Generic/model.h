#ifndef MODEL_H
#define MODEL_H
#include <stdlib.h> 

typedef void (*PropensityFunc)(double **, long **, int);
#define n_species 2 
#define m_reaction 4 
#define w_bin 100 
#define L 10.0
#define NUMofRUNS 1
#define SimulationTime 10
#define Mode Auto //Choose from Auto or Manual 

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
void reactionChange(int i_reaction,int j_bin, long **species_population,long *total_population, double **a, double *total_a);
void diffusionChange(int ruleIndex, int target_bin, long **population, double **a, double *total_a);

#endif