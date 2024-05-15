#ifndef MODEL_H
#define MODEL_H
#include <stdlib.h> 

typedef void (*PropensityFunc)(double **, long **, int);
#define n_species 2 //Number of Species
#define m_reaction 4  //Number of Reactions
#define w_bin 100 //Total Number of Bins
#define L 10.0 //Interval Length
#define NUMofRUNS 1  //Number of Simulation
#define SimulationTime 10 //Total Simulation time
#define OutputPoints 100 //Output length for populationfile
#define Mode Auto//Choose from Auto or Manual 

/*
Struct definition for Species
*/
struct Species{
    int index; // Index of the species
    char name[100]; // Name of the species
    int initialValue;// Initial population value of the species
    double diffusionRate;// Diffusion rate of the species
};

/*
Struct definition for Reaction
*/
struct Reaction{
    int type; // Reaction type of reaction
    double rateConstant;// Reaction Rate of reaction
    int reactant;  // Index of reactant of reaction (only for reaction type 1)
    PropensityFunc calculatePropensity;// Function name to calculate propensity value (only for reaction type 2)
};

// Function prototypes
struct Species* createSpeciesArray();
struct Reaction* createReactionArray();
void prop3(double **a,long **species_population,int rulenum);
void reactionChange(int i_reaction,int j_bin, long **species_population,long *total_population);
void reaction_propensityChange(int i_reaction,int j_bin,long **species_population, long total_population[], double **a, double total_a[], double jumpRate[]);
void diffusionChange(int speciesIndex, int bin, long **population, double **a, double *total_a,long *total_population,double jumpRate[]);
#endif