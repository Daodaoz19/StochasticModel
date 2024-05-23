#include "model.h"
#include <string.h>
#include <stdio.h>

/*
Dynamically allocates memory for an array of Species structs 
Initializes each species parameter with default values.
*/
struct Species *createSpeciesArray()
{
    struct Species *species = (struct Species *)malloc(n_species * sizeof(struct Species));
    for (size_t i = 0; i < n_species; ++i)
    {
        species[i].index = i;            // Index of the species
        strcpy(species[i].name, "NAME"); // Name of the species
        species[i].initialValue = 0;     // Initial population value of the species
        species[i].diffusionRate = 0;    // Diffusion rate of the species
    }
    return species;
}

/*
Dynamically allocates memory for an array of Reaction structs 
Initializes each species parameter with default values.

Reaction type 0: 0 molecule reaction
Reaction type 1: 1 molecule reaction
Reaction type 2: 2 or more molecules reaction

*/
struct Reaction *createReactionArray()
{
    struct Reaction *reaction = (struct Reaction *)malloc(m_reaction * sizeof(struct Reaction));
    for (size_t i = 0; i < m_reaction; ++i)
    {
        reaction[i].type = 0;                   // Reaction type of reaction
        reaction[i].rateConstant = 0;           // Reaction Rate of reaction
        reaction[i].reactant = 0;               // Index of reactant of reaction (only for reaction type 1)
        reaction[i].calculatePropensity = NULL; // Function name to calculate propensity value (only for reaction type 2)
    }
    return reaction;
}

/*
Initializes the simulation by creating arrays for species and reactions
with specific values for each parameter in species and reactions base on specific model.
Here, we use Shnakenberg model as example:
A -> X
X -> C
B -> Y
2X + Y -> 3X
*/
void initializeSimulation(struct Species **species, struct Reaction **reactions)
{
    *species = createSpeciesArray(n_species);
    *reactions = createReactionArray(m_reaction);
     
    // parameter of species X
    strcpy((*species)[0].name, "X");
    (*species)[0].initialValue = 1;
    (*species)[0].diffusionRate = 5e-4;
    // parameter of species Y
    strcpy((*species)[1].name, "Y");
    (*species)[1].initialValue = 2;
    (*species)[1].diffusionRate = (5e-4 * 10000);

    // parameter of reaction: A -> X
    (*reactions)[0].type = 0;
    (*reactions)[0].rateConstant = 10.0;

    // parameter of reaction: X -> C 
    (*reactions)[1].type = 1;
    (*reactions)[1].rateConstant = 2.0;
    (*reactions)[1].reactant = 0;

    // parameter of reaction: B -> Y
    (*reactions)[2].type = 0;
    (*reactions)[2].rateConstant = 30.0;

    // parameter of reaction: 2X + Y -> 3X
    (*reactions)[3].type = 2;
    (*reactions)[3].calculatePropensity = prop3;
}

/*
Function to calculate propneisty value for  2X + Y -> 3X
*/
void prop3(double **a, long **species_population, int rulenum)
{
    for (int i = 0; i < w_bin; i++)
    {
        long x = species_population[0][i];// Population of species X
        long y = species_population[1][i];// Population of species Y
        a[rulenum][i] = x * (x - 1) * y;// Calculate propensity value
    }
}

/*
Update the population when each reaction occurs
*/
void reactionChange(int i_reaction, int j_bin, long **species_population, long *total_population)
{
    switch (i_reaction)
    {
    case 0: // A -> X
        species_population[0][j_bin]++; //Population of X in j_bin increase 
        total_population[0]++; //Total population of X increase
        break;
    case 1: // X -> C
        species_population[0][j_bin]--;//Population of X in j_bin decrease 
        total_population[0]--;//Total population of X decrease 
        break;
    case 2: // B -> Y
        species_population[1][j_bin]++;//Population of Y in j_bin increase 
        total_population[1]++;//Total population of Y increase
        break;
    case 3: // 2X + Y -> 3X
        species_population[0][j_bin]++;//Population of X in j_bin increase 
        species_population[1][j_bin]--;//Population of Y in j_bin decrease 
        total_population[0]++;//Total population of X increase
        total_population[1]--;//Total population of Y decrease
        break;
    default:
        break;
    }
}

/*

Update the propensity value when each reaction fires in Mannual Mode

total[0]: X jump
total[1]: y jump
total[2]: A -> X
total[3]: X -> C
total[4]: B -> Y
total[5]: 2X + Y -> 3X
*/
void reaction_propensityChange(int i_reaction, int j_bin, long **species_population, long total_population[], double **a, double total_a[], double jumpRate[])
{
    long x = species_population[0][j_bin];
    long y = species_population[1][j_bin];
    switch (i_reaction)
    {
    case 0: // A -> X
        total_a[0] = jumpRate[0] * total_population[0];
        total_a[3] = 2 * total_population[0];
        total_a[5]-=a[3][j_bin];
        a[3][j_bin]= x * (x - 1) * y;
         total_a[5] += a[3][j_bin];
        break;
    case 1: // X -> C
        total_a[0] = jumpRate[0] * total_population[0];
        total_a[3] = 2 * total_population[0];
        total_a[5]-=a[3][j_bin];
        a[3][j_bin]= x * (x - 1) * y;
         total_a[5] += a[3][j_bin];
        break;
    case 2: // B -> Y
        total_a[1] = jumpRate[1] * total_population[1];
        total_a[3] = 2 * total_population[0];
        total_a[5]-=a[3][j_bin];
        a[3][j_bin]= x * (x - 1) * y;
        total_a[5] += a[3][j_bin];
        break;
    case 3: // 2X + Y -> 3X
        total_a[0] = jumpRate[0] * total_population[0];
        total_a[1] = jumpRate[1] * total_population[1];
        total_a[3] = 2 * total_population[0];
        total_a[5]-=a[3][j_bin];
        a[3][j_bin]= x * (x - 1) * y;
        total_a[5] += a[3][j_bin];
        break;
    }
}

/*

Update the propensity value when diffusion fires in Mannual Mode

total[0]: X jump
total[1]: y jump
total[2]: A -> X
total[3]: X -> C
total[4]: B -> Y
total[5]: 2X + Y -> 3X
*/
void diffusionChange(int speciesIndex, int j_bin, long **species_population, double **a, double *total_a, long *total_population, double jumpRate[])
{
    long x = species_population[0][j_bin];
    long y = species_population[1][j_bin];
    switch (speciesIndex)
    {
    case 0://X jump

        total_a[0] = jumpRate[0] * total_population[0];

        total_a[3] = 2 * total_population[0];

        total_a[5]-=a[3][j_bin];
        a[3][j_bin]= x * (x - 1) * y;
        total_a[5] += a[3][j_bin];
        break;
    case 1://Y jump

        total_a[1] = jumpRate[1] * total_population[1];

        total_a[5]-=a[3][j_bin];
        a[3][j_bin]= x * (x - 1) * y;
        total_a[5] += a[3][j_bin];
        break;
    }
}