#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "model.h"

/*
Defult value 
*/
int NUMofRUNS = 1; 
int SimulationTime = 10; 

extern void initializeSimulation(struct Species** species, struct Reaction** reactions);
long** initializePopulation(struct Species* species);
int main(int argc, char *argv[])
{
    //  if (argc > 1) {
    //         NUMofRUNS = atoi(argv[0]);
    //  }
    //     if (argc > 2) {
    //         SimulationTime = atoi(argv[1]);
    //     }
    //  else {
    //     fprintf(stderr, "Usage: %s [NUMofRUNS] [SimulationTime]\n", argv[0]);
    //     return 1; 
    // }

    struct Species* species = createSpeciesArray(n_species);
    struct Reaction* reactions = createReactionArray(m_reaction);
    initializeSimulation(&species, &reactions);

    double h = L / w_bin;

    double jumpRate[n_species];
    for (int i = 0; i < n_species; i++)
    {
        jumpRate[i] = 2*species[i].diffusionRate / (h * h);
    }

    const int RULENUM = n_species + m_reaction;
    double a0,total_a[RULENUM];
    double r1, r2, tau;
    int ita;
    double r2a0, sum_a, sum_x;


    long firings[RULENUM], location[w_bin];
    FILE *total = fopen("total", "w");
    FILE *firefile = fopen("firings", "w");
    FILE *locationfile = fopen("location", "w"); 

    printf("Begin the model simulation ...\n");
    printf("BINNUM = %d\n", w_bin);

    srand(time(NULL)); 
   
    for (int real = 0; real < NUMofRUNS; real++)
    {
        /*
        Initialize total population for each species
        */
        long tot[n_species];
        for (int n = 0; n < n_species; n++)
        {
            tot[n] = 0;
        }

        long** population = initializePopulation(species);
        
        for (int n = 0; n < n_species; n++) {
        for (int w = 0; w < w_bin; w++) { 
            tot[n] += population[n][w]; 
        }}

        double **a; 
        a = (double **)malloc(RULENUM * sizeof(double *));
        for (int n = 0; n < RULENUM; n++)
        {    
            a[n] = (double *)malloc(w_bin * sizeof(double));
        }

        double timeTracker = 0.0; 
        while (timeTracker < SimulationTime)
        {
           
            for(int n =0; n<n_species;n++){
                total_a[n] = jumpRate[n] * (double)(tot[n]);  
                double propensityPerBin = total_a[n]/w_bin; 
                    for(int i = 0; i < w_bin; i++) {
                        a[n][i] = propensityPerBin;
                    }
            }
             for (int run = 0; run < m_reaction; run++)
            {
                double rate = reactions[run].rateConstant;
                int type =reactions[run].type;
                 
                if (type == 0 || type == 1) {
                  for (int i = 0; i < w_bin; i++) {
                    if (type == 0) {
                        total_a[n_species + run] = rate * L; 
                        double propensityPerBin = total_a[n_species + run]/w_bin; 
                        a[n_species+run][i] = propensityPerBin;
                    } else if (type == 1) {
                        int reactantIndex = reactions[run].reactant[0];
                        total_a[n_species + run] = rate*tot[reactantIndex]; 
                        double propensityPerBin = total_a[n_species + run]/w_bin; 
                        a[n_species+run][i] = propensityPerBin;
                    }
                   }
                } else if (type == 2) {
                     //total_a[n_species +run] = 0; 
                     reactions[run].calculatePropensity(a, population, run);
                     for(int i=0;i<w_bin;i++){
                       total_a[n_species +run]+=a[run][i];
                        double propensityPerBin = total_a[n_species + run]/w_bin; 
                        a[n_species+run][i] = propensityPerBin;
                     }
                     
                }     
                
            }
            

            // printf("total_a is %f\n",total_a[0]);
            // printf("total_a is %f\n",total_a[1]);
            // printf("total_a is %f\n",total_a[2]);
            // printf("total_a is %f\n",total_a[3]);
            // printf("total_a is %f\n",total_a[4]);
            // printf("total_a is %f\n",total_a[5]);
            // printf("\n");
            // printf("a[0][bin] is %f\n",a[0][10]);
            // printf("a[1][bin] is %f\n",a[1][10]);
            // printf("a[2][bin] is %f\n",a[2][10]);
            // printf("a[3][bin] is %f\n",a[3][10]);
            // printf("a[4][bin] is %f\n",a[4][10]);
            // printf("a[5][bin] is %f\n",a[5][10]);
      
            a0=0;
            for(int i =0;i<RULENUM;i++){
                a0 += total_a[i];  
            }
            r1 = (double)rand() / RAND_MAX;
            tau = -1.0/a0*log(r1);
            timeTracker += tau; 
             
            r2 = (double)rand() / RAND_MAX;
            r2a0 = r2 * a0; 
            sum_a = total_a[0];
            /*
            select which reaction will happen
            */
            ita = 0;//index of reaction
              while (sum_a < r2a0)
            {
                ita++;
                sum_a += total_a[ita];
            }
            sum_a -= total_a[ita]; 
            printf("At time %f, reaction %d fire\n",timeTracker,ita);
            double residue = r2a0 - sum_a;
            firings[ita]++;
          
            int target=0;//index of bin

            /*
            ----------------------------
                                        |
                                        |
                                       \ / r2a0 =90010, ita =1
                    
            reaction           |__X__|_____Y______|___1__|__2__|__3__|__4__|
            index              0     1      2      3     4     5     6
            propensity         0    200         90000   100   200   300   400
    
            */
            if (ita < n_species) 
            {   
                long species_sum = population[ita][0]; 
                residue = residue/jumpRate[ita];
                while(species_sum<residue){  
                    target++;
                    species_sum +=population[ita][target];  
                    residue = residue-species_sum;
                }   
                if(0.5*population[ita][target] > residue){
                    population[ita][target]--;
                    population[ita][target + 1]++;
                    //printf("jump to right\n");
                }
                else{
                    population[ita][target]--;
                    population[ita][target - 1]++;
                    //printf("jump to left\n");
                }
                //break;
                //propensityChange(ita,target,a[ita][target]);
            }
            else{ 
               int reactionType = ita -n_species;
               if (reactions[reactionType].type == 0) {
                    int target = floor(((residue)/total_a[ita])*w_bin);
                    location[target]++;
                }else if (reactions[reactionType].type == 1) {
                    /*
                    for reaction type 1, we want to find the target bin contain the maximum value of species
                    population.
                    */
                    int type =reactions[ita].type; //type of species that reacts
                    int boundary_bin = floor(((residue)/total_a[ita])*w_bin);
                    long species_sum = population[type][0]; 
                    //search for largest population bin within the boundary
                    for(int i=0;i<boundary_bin;i++){
                    if(species_sum < population[type][i])
                    {
                        species_sum = population[type][i]; 
                        target = i; }
                    }
                    location[target]++;
                }else if (reactions[reactionType].type == 2) {
                    /*
                    for reaction type 2, we want to find the target bin contain the maximum value of species
                    population. We propobally need modeler to input the index of species involve in this reaction
                    */
                    location[target]++;
                }
                //after finding the bin, apply statechange to the reaction base on modeler provied function
                stateChange(ita, target, population, tot);
            }





        }


     
        

  






     for (size_t i = 0; i < n_species; ++i) {
        free(population[i]);
         }
     free(population);
    
    }
    free(species);
    free(reactions);

    return 0;

}


long** initializePopulation(struct Species* species) {
    long** population = (long**)malloc(n_species * sizeof(long*));
    for (int i = 0; i < n_species; ++i) {
        population[i] = (long*)malloc(w_bin * sizeof(long));
        for (int j = 0; j < w_bin; ++j) {
            population[i][j] = species[i].initialValue;
        }
    }
    return population;
}