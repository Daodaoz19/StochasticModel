#include "model.h"
#include <string.h>


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
        reaction[i].reactant = NULL;// Default type
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
    (*reactions)[1].reactant = (int*)malloc((*reactions)[1].type * sizeof(int));
    (*reactions)[1].reactant[0]= 0;

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
        a[rulenum][i] = x * (x - 1) * y; 
       } 
}

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
