#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "model.h"

enum Mode{ Auto, Manual };
extern void initializeSimulation(struct Species **species, struct Reaction **reactions);
long **initializePopulation(struct Species *species);
void calculatePropensities(double total_a[], long total_population[], struct Reaction *reactions, struct Species *species, double **a, long **population, double jumpRate[]);
int main(int argc, char *argv[])
{
    struct Species *species = createSpeciesArray(n_species);
    struct Reaction *reactions = createReactionArray(m_reaction);
    initializeSimulation(&species, &reactions);
    
    int i, j;

    //Set up initial population
    long **population;
    char* inputFilename;
    if (argc ==1) {// Use default initial values from model.c
        population = initializePopulation(species);
        printf("starting simulation with initial population from initialValue\n");
    } else {//Use inital population from user provied file 
       FILE *inputfile = fopen(argv[1], "r");
        if (!inputfile) {
            fprintf(stderr, "Error opening file %s\n", argv[1]);
            return EXIT_FAILURE;
        }

        population = (long **)malloc(n_species * sizeof(long *));
        if (!population) {
            fprintf(stderr, "Memory allocation failed\n");
            fclose(inputfile);
            return EXIT_FAILURE;
        }

        char buffer[1000];  
        for (i = 0; i < n_species; i++) {
            population[i] = (long *)malloc(w_bin * sizeof(long));
            char *token = strtok(buffer, ",");  
            for (j = 0; token != NULL && j < w_bin; j++) {
                population[i][j] = strtol(token, NULL, 10);
                token = strtok(NULL, ",");
            }
        }
        fclose(inputfile);
        printf("Starting simulation with initial population from %s\n", argv[1]);
    }

    //Variable for model parameter
    double reaction_rate;
    int reaction_type,reactant_index;
    int RULENUM = n_species + m_reaction;
    
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
    FILE *logfile = fopen("log", "w");

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

    calculatePropensities(total_a, total_population, reactions, species, a, population, jumpRate); 
    
    srand(time(NULL));
    //srand(2001);
    for (int real = 0; real < NUMofRUNS; real++)
    {
       
        //Initialize output file
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
            //Calculate propensity value for each reaction 
            if(Mode == Auto) {
               calculatePropensities(total_a, total_population, reactions, species, a, population, jumpRate); 
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
 
            fprintf(logfile, "at time %f, reaction %d fire, a0 is %f\n", timeTracker, ruleIndex, a0);
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
                    target_bin++; // Increment to find the target bin
                    species_sum += population[ruleIndex][target_bin];
                }
                residue = species_sum - residue;
                if (0.5 * population[ruleIndex][target_bin] > residue)// jump to right
                { 
                    if (target_bin < w_bin - 1)
                    {
                        population[ruleIndex][target_bin]--;
                        population[ruleIndex][target_bin + 1]++;
                        // if(Mode == Manual){
                        // diffusionChange(ruleIndex,target_bin, population, a,total_a);
                        // diffusionChange(ruleIndex,target_bin+1, population, a,total_a);
                        // }
                       
                    } 
                    //if it's the last bin on the right, do nothing
                }
                else// jump to left
                { 
                    if (target_bin > 0)
                    {
            
                        population[ruleIndex][target_bin]--;
                        population[ruleIndex][target_bin - 1]++;
                        //  if(Mode == Manual){
                        // diffusionChange(ruleIndex,target_bin, population, a,total_a);
                        // diffusionChange(ruleIndex,target_bin-1, population, a,total_a);
                        //  }
                    } 
                    //if it's the first bin on the left, do nothing
                }
                
            }
            else //Find the bin index where the reaction occurs
            {
                reaction_index = ruleIndex - n_species;
                switch (reactions[reaction_index].type) {
                    case 0:
                        // For reaction type 0, find the target bin based on the residue location
                        target_bin = floor((residue / total_a[ruleIndex]) * w_bin);
                        location[target_bin]++;
                        break;
                    case 1:
                        // For reaction type 1, calculate target bin based on population distribution
                        target_bin = 0;
                        reactant_index = reactions[reaction_index].reactant;
                        species_sum = population[reactant_index][0];
                        residue = residue * total_population[reactant_index] / total_a[ruleIndex];
                        while (species_sum < residue) {
                            target_bin++; // Increment to find the target bin
                            species_sum += population[reactant_index][target_bin];
                        }
                        location[target_bin]++;
                        break;
                    case 2:
                        // For reaction type 2, calculate based on a propensity function
                        target_bin = 0;
                        propensity_sum = a[reaction_index][0];
                        while (propensity_sum < residue) {
                            target_bin++; // Increment to find the target bin
                            propensity_sum += a[reaction_index][target_bin];
                        }
                        location[target_bin]++;
                        break;
                }
                reactionChange(reaction_index, target_bin, population, total_population);
                // if(Mode == Manual){
                // reaction_propensityChange(reaction_index, target_bin,population,total_population, a,total_a);   
                // }
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
    fclose(logfile);
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

void calculatePropensities(double total_a[], long total_population[], struct Reaction *reactions, struct Species *species, double **a, long **population, double jumpRate[]) {
    int i, run;

   //Calculate propensity for diffusions based on jump rates calculated
    for (i = 0; i < n_species; i++) {
        total_a[i] = jumpRate[i] * (double)(total_population[i]);
    }

    // Calculate propensity for reactions based on reaction type
    for (run = 0; run < m_reaction; run++) {
        double reaction_rate = reactions[run].rateConstant;
        int reaction_type = reactions[run].type;
        switch (reaction_type) {
            case 0: // For reaction type 0, compute propensity values by reaction rate * length
                total_a[n_species + run] = reaction_rate * L;
                break;
            case 1: // For reaction type 1, compute propensity values by reaction rate * total population of the reactant
                int reactant_index = reactions[run].reactant;
                total_a[n_species + run] = reaction_rate * total_population[reactant_index];
                break;
            case 2: // For reaction type 2, compute propensity values using a user-defined function
                total_a[n_species + run] = 0;
                if (reactions[run].calculatePropensity != NULL) {
                    reactions[run].calculatePropensity(a, population, run);

                    for (i = 0; i < w_bin; i++) {
                        total_a[n_species + run] += a[run][i];
                    }
                }
                break;
        }
    }
}

