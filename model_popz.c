#include "model_PopZ.h"
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
    
    strcpy((*species)[0].name, "PopZm");
    (*species)[0].diffusionRate = 835;

    strcpy((*species)[1].name, "PopZp");
    (*species)[1].diffusionRate = 0.0005;


    // Reactions

    // Null -> PopZm
    (*reactions)[0].type = 0;
    (*reactions)[0].rateConstant = 150; 

    // PopZm -> Null
    (*reactions)[1].type = 1;
    (*reactions)[1].rateConstant = 0.02;
    (*reactions)[1].reactant[0]= 0;

    //PopZp -> PopZm
    (*reactions)[2].type = 1;
    (*reactions)[2].rateConstant = 0.15; 
    (*reactions)[2].reactant[1]= 0;

    //PopZm + 2PopZp -> 3PopZp (with enzyme of PodJLt)
    (*reactions)[3].type = 2;
    (*reactions)[3].calculatePropensity = prop3;

    //PopZm -> PopZp
    (*reactions)[4].type = 1;
    (*reactions)[4].rateConstant = 40
    (*reactions)[4].reactant[0]= 0;

    //PopZm + 2PopZp -> 3PopZp
    (*reactions)[5].type = 2;
    (*reactions)[5].calculatePropensity = prop5;

    // PopZp -> Null
    (*reactions)[6].type = 1;
    (*reactions)[6].rateConstant = 0.02;
    (*reactions)[6].reactant[1]= 0;


}

void prop3(double **a,long **species_population,int rulenum) {  
    //PopZm + 2PopZp -> 3PopZp (with enzyme of PodJLt)
    //k * a_popzpodj * [PodJ] * [PopZm] * [PopZp]**2
    for (int i = 0; i < w_bin; i++) {
        long x = species_population[0][i]; 
        long y = species_population[1][i];
        a[rulenum][i] = x * (x - 1) * y * 0.0012 // Assume the population of [PodJ] is 1 (0.001uM)
       } 
}
void prop5(double **a,long **species_population,int rulenum) {  
    //PopZm + 2PopZp -> 3PopZp, k*[PopZm]*[PopZp]**2
    for (int i = 0; i < w_bin; i++) {
        long x = species_population[1][i]; 
        long y = species_population[0][i];
        a[rulenum][i] = 0.00001 * x * (x - 1) * y; 
       } 
       
}

void stateChange(int i,int j, long **species_population,long *tot) {
   switch(i) {
    case 0: // Null -> PopZm
        species_population[0][j]++;
        tot[0]++; 
        break;
    case 1: //PopZm -> Null
        species_population[0][j]--;
        tot[0]--; 
        break;
    case 2: //PopZp -> PopZm
        species_population[0][j]++;
        species_population[1][j]--;
        tot[0]++;
        tot[1]--; 
        break;
    case 3: //PopZm + 2PopZp -> 3PopZp (with enzyme of PodJLt)
        species_population[1][j]++;
        species_population[0][j]--;
        tot[1]++;
        tot[0]--
        break;
    case 4: //PopZm -> PopZp
        species_population[0][j]--;
        species_population[1][j]++;
        tot[0]--; 
        tot[1]++; 
        break;
    case 5: //PopZm + 2PopZp -> 3PopZp
        species_population[1][j]++;
        species_population[0][j]--;
        tot[1]++;
        tot[0]-- 
        break;
    case 6: // PopZp -> Null
        species_population[1][j]--;
        tot[1]--;  
        break;
    
    default:
        exit(1);
        break;
    }
}