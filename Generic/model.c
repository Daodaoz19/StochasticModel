#include "model.h"
#include <string.h>
#include <stdio.h>

struct Species* createSpeciesArray() {
    struct Species* species  = (struct Species*)malloc(n_species * sizeof(struct Species));
    for (size_t i = 0; i < n_species; ++i) {
        species[i].index = i; //index
        strcpy(species[i].name, "NAME"); // Default type
        species[i].initialValue = 0; // Default type
        species[i].diffusionRate = 0; // Default type
    }
    return species;
}

struct Reaction* createReactionArray() {
    struct Reaction* reaction= (struct Reaction*)malloc(m_reaction * sizeof(struct Reaction));
    for (size_t i = 0; i < m_reaction; ++i) {
        reaction[i].type = 1; // Default type
        reaction[i].rateConstant = 0; // Default type
        reaction[i].reactant = 0;// Default type
        reaction[i].calculatePropensity =NULL;// Default type
    }
    return reaction;
}

void initializeSimulation(struct Species** species, struct Reaction** reactions) {
    *species = createSpeciesArray(n_species);
    *reactions = createReactionArray(m_reaction);
    
    strcpy((*species)[0].name, "X");
    (*species)[0].initialValue = 2;
    (*species)[0].diffusionRate = 5e-4;

    strcpy((*species)[1].name, "Y");
    (*species)[1].initialValue = 1;
    (*species)[1].diffusionRate = (5e-4 * 10000);

    // A -> X
    (*reactions)[0].type = 0;
    (*reactions)[0].rateConstant = 1.0;

    
    // X -> C
    (*reactions)[1].type = 1;
    (*reactions)[1].rateConstant = 2.0;
    (*reactions)[1].reactant= 0;

    //B -> Y
    (*reactions)[2].type = 0;
    (*reactions)[2].rateConstant = 3.0;


    //2X + Y -> 3X
    (*reactions)[3].type = 2;
    (*reactions)[3].rateConstant = 1.0/w_bin*w_bin;
    (*reactions)[3].calculatePropensity = prop3;
   
}

void prop3(double **a,long **species_population,int rulenum) {  
    for (int i = 0; i < w_bin; i++) {
        long x = species_population[0][i]; 
        long y = species_population[1][i];
        a[rulenum][i] = x * (x - 1) * y * 1.0/(w_bin*w_bin)* w_bin*w_bin; 
       } 

}

void reactionChange(int i_reaction,int j_bin, long **species_population,long *total_population, double **a, double *total_a) {
   switch(i_reaction) {
    case 0:// A -> X
        species_population[0][j_bin]++;
        total_population[0]++; 

        total_a[i_reaction] =1.0*10.0;//RateConstant * L
        a[i_reaction][j_bin] =1.0*10.0;

        //A -> X fire, we need to update the propensity for reaction that involve X as a reactant:
        //X -> C 
        total_a[i_reaction+1] =2.0*total_population[0];//RateConstant * reactant.population
        a[i_reaction+1][j_bin] =2.0*total_population[0];
        //2X + Y -> 3X
        total_a[i_reaction+3] =200;//hardcoded dummy value
        a[i_reaction+3][j_bin] =200;
        break;
    case 1:// X -> C
        species_population[0][j_bin]--;
        total_population[0]--; 
        
        total_a[i_reaction] =2.0*total_population[0];//RateConstant * reactant.population
        a[i_reaction][j_bin] =2.0*total_population[0];
        
        //X -> C fire, we need to update the propensity for reaction that involve X as a reactant:
        //2X + Y -> 3X
        total_a[i_reaction+2] =200;//hardcoded dummy value
        a[i_reaction+2][j_bin] =200;
        break;
    case 2://B -> Y
        species_population[1][j_bin]++;
        total_population[1]++; 
 
        total_a[i_reaction] =3.0*10.0;//RateConstant * L
        a[i_reaction][j_bin] =3.0*10.0;

        //B -> Y fire, we need to update the propensity for reaction that involve Y as a reactant:
        //2X + Y -> 3X
        total_a[i_reaction+1] =200;//hardcoded dummy value
        a[i_reaction+1][j_bin] =200;
        break;
    case 3://2X + Y -> 3X
        species_population[0][j_bin]++;
        species_population[1][j_bin]--;
        total_population[0]++; 
        total_population[1]--; 

        total_a[i_reaction] =200;//hardcoded dummy value
        a[i_reaction][j_bin] =200;

        //2X + Y -> 3X fire, we need to update the propensity for reaction that involve X or Y as a reactant:
        // X -> C
        total_a[i_reaction-2] =2.0*total_population[0];//RateConstant * reactant.population
        a[i_reaction-2][j_bin] =2.0*total_population[0];


        break;
    default:
        break;
    }
}

void diffusionChange(int speciesIndex, int bin, long **population, double **a, double *total_a) {
    if (speciesIndex==0){//update propensity value for x diffuse
    total_a[speciesIndex ] =5e-4 * population[0][bin];
    a[speciesIndex][bin] = 5e-4 * population[0][bin];
    }else{//update propensity value for y diffuse
    total_a[speciesIndex ] =5 * population[1][bin];
    a[speciesIndex][bin] = 5 * population[1][bin];
    }   
}
