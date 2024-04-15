#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "model.h"

extern void initializeSimulation(struct Species **species, struct Reaction **reactions);
long **initializePopulation(struct Species *species);
int main(int argc, char *argv[])
{
    struct Species *species = createSpeciesArray(n_species);
    struct Reaction *reactions = createReactionArray(m_reaction);
    initializeSimulation(&species, &reactions);
    //Variable for model parameter
    double reaction_rate;
    int reaction_type,reactant_index;
    int RULENUM = n_species + m_reaction;
    
    int i, j;

    // Calculate interval length
    double h = L / w_bin;

    // Calculate the bidirectional jump rate for each species
    double jumpRate[n_species];
    for (i = 0; i < n_species; i++)
    {
        jumpRate[i] = 2 * species[i].diffusionRate / (h * h);
    }

    // Variable for SSA process
    double a0, total_a[RULENUM];
    double r1, r2, tau;
    int ruleIndex, target_bin, reaction_index;
    double r2a0, sum_a, residue, propensity_sum;
    long species_sum;
    
    long total_population[n_species];
    for (i = 0; i < n_species; i++)
    {
        total_population[i] = 0;
    }

    // Set up files for logging final populations and reactions
    long firings[RULENUM], location[w_bin];
    FILE *final_population = fopen("population", "w");
    FILE *total = fopen("total", "w");
    FILE *firefile = fopen("firings", "w");
    FILE *locationfile = fopen("location", "w");

    printf("Begin the model simulation ...\n");
    printf("BINNUM = %d\n", w_bin);

    srand(time(NULL));
    for (int real = 0; real < NUMofRUNS; real++)
    {
        //Set up initial population
        long **population = initializePopulation(species);
        //Calculate total population for each species
        for (i = 0; i < n_species; i++)
        {
            for (j = 0; j < w_bin; j++)
            {
                total_population[i] += population[i][j];
            }
        }
        //Allocate memory for the propensity array a[RULENUM][W_BIN]
        double **a;
        a = (double **)malloc(RULENUM * sizeof(double *));
        for (i = 0; i < RULENUM; i++)
        {
            a[i] = (double *)malloc(w_bin * sizeof(double));
        }
        //Initialize log file
        for (i = 0; i < RULENUM; i++)
        {
            firings[i] = 0;
        }
        for (i = 0; i < w_bin; i++)
        {
            location[i] = 0;
        }

        double timeTracker = 0.0;
        while (timeTracker < SimulationTime)
        {
            //Calculate propensity for diffusions
            for (i = 0; i < n_species; i++) 
            {
                total_a[i] = jumpRate[i] * (double)(total_population[i]);
            }
            //Calculate propensity for reactions base on reaction_type
            for (int run = 0; run < m_reaction; run++)
            {
                reaction_rate = reactions[run].rateConstant;
                reaction_type = reactions[run].type;
                if (reaction_type == 0)// For reaction type 0, compute propensity values by reaction rate * length
                {
                    total_a[n_species + run] = reaction_rate * L;
                }
                else if (reaction_type == 1)// For reaction type 1, compute propensity values by reaction rate * total population
                {
                    reaction_index = reactions[run].reactant;
                    total_a[n_species + run] = reaction_rate * total_population[reaction_index];
                }
                else// For reaction type 2, compute propensity values using a user-defined function
                {
                    total_a[n_species + run] = 0;
                    reactions[run].calculatePropensity(a, population, run);
                    for (int i = 0; i < w_bin; i++)
                    {
                        total_a[n_species + run] += a[run][i];
                    }
                }
            }

            //calculate a0
            a0 = 0;
            for (i = 0; i < RULENUM; i++)
            {
                a0 += total_a[i];
            }
            //Calculate the time when the next reaction occurs
            r1 = (double)rand() / RAND_MAX;
            tau = -1.0 / a0 * log(r1);
            timeTracker += tau;

            //Find the index of the next reaction occurs
            r2 = (double)rand() / RAND_MAX;
            r2a0 = r2 * a0;
            sum_a = total_a[0];
    
            ruleIndex = 0;//Index of the firing reaction
            while (sum_a < r2a0)
            {
                ruleIndex++;
                sum_a += total_a[ruleIndex];
            }
            sum_a -= total_a[ruleIndex];
            residue = r2a0 - sum_a;
            firings[ruleIndex]++;

            //Find the bin index where the jump occurs
            if (ruleIndex < n_species)
            {
                target_bin = 0;
                species_sum = population[ruleIndex][0];
                residue = residue / jumpRate[ruleIndex];
                while (species_sum < residue)
                {
                    target_bin++; //find the bin that the reaction occur
                    species_sum += population[ruleIndex][target_bin];
                }
                residue = species_sum - residue;
                if (0.5 * population[ruleIndex][target_bin] > residue)// jump to right
                { 
                    if (target_bin < w_bin - 1)
                    {
                        population[ruleIndex][target_bin]--;
                        population[ruleIndex][target_bin + 1]++;

                    } 
                    //if it's the last bin on the right, do nothing
                }
                else// jump to left
                { 
                    if (target_bin > 0)
                    {
                        population[ruleIndex][target_bin]--;
                        population[ruleIndex][target_bin - 1]++;
                    } 
                    //if it's the first bin on the left, do nothing
                }
            }
            else //Find the bin index where the reaction occurs
            {
                reaction_index = ruleIndex - n_species;
                if (reactions[reaction_index].type == 0)
                {
                    target_bin = floor(((residue) / total_a[ruleIndex]) * w_bin); //find the target bin
                    location[target_bin]++;
                }
                else if (reactions[reaction_index].type == 1)
                {
                    target_bin = 0;
                    reactant_index = reactions[reaction_index].reactant;
                    species_sum = population[reactant_index][0];

                    residue = residue * total_population[reactant_index] / total_a[ruleIndex];
                    //partial sum of population exceed the population residue
                    while (species_sum < residue)
                    {
                        target_bin++;//find the target bin
                        species_sum += population[reactant_index][target_bin];
                    }
                    location[target_bin]++;
                }
                else
                {
                    target_bin = 0;
                    propensity_sum = a[reaction_index][0];
                    //partial sum of propensity exceeds the propensity residue
                    while (propensity_sum < residue)
                    {
                        target_bin++;//find the target bin
                        propensity_sum += a[reaction_index][target_bin];
                    }
                    location[target_bin]++;
                }
                reactionChange(reaction_index, target_bin, population, total_population);
            }
        }
        // end of simulation, collect the statistics
        for (i = 0; i < n_species; i++)
        {
            fprintf(total, "%ld\t", total_population[i]);
        }
          fprintf(total, "\n");

        for (i = 0; i < n_species; i++)
        {
            for (j = 0; j < w_bin; j++)
            {
                fprintf(final_population, "%ld\t", population[i][j]);
            }
            fprintf(final_population, "\n");
        }
        
        for (i = 0; i < RULENUM; i++)
        {
            fprintf(firefile, "%ld\t", firings[i]);
        }
        fprintf(firefile, "\n");

        for (i = 0; i < w_bin; i++)
        {
            fprintf(locationfile, "%ld\t", location[i]);
        }
        fprintf(locationfile, "\n");

    }

    free(species);
    free(reactions);
    fclose(final_population);
    fclose(total);
    fclose(firefile);
    fclose(locationfile); 
    printf("End of simulation ...\n");
    return 0;
}

long **initializePopulation(struct Species *species)
{   //Allocate memory for the population array population[n_species][W_BIN]
    long **population = (long **)malloc(n_species * sizeof(long *));
    for (int i = 0; i < n_species; i++)
    {
        population[i] = (long *)malloc(w_bin * sizeof(long));
        for (int j = 0; j < w_bin; j++)
        { // Set the initial population for each bin by the user-defined initial value for each species
            population[i][j] = species[i].initialValue;
        }
    }
    return population;
}