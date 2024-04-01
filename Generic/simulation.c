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
        jumpRate[i] = species[i].diffusionRate / (h * h);
    }

    const int RULENUM = n_species + m_reaction;
    double a0,total_a[RULENUM];
    double r1, r2, tau;
    int ita;
    double sum, sum_a, sum_x;


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
                total_a[n] = 2*jumpRate[n] * (double)(tot[n]);  
            }
             for (int run = 0; run < m_reaction; run++)
            {
                double rate = reactions[run].rateConstant;
                int type =reactions[run].type;
                if (type == 0)
                {
                   total_a[n_species + run] = rate * L; 

                } 
                else if (type == 1)
                {
                    int reactantIndex = reactions[run].reactant[0];
                    total_a[n_species + run] = rate*tot[reactantIndex]; 
                }
                else if (type == 2)
                {    
                     total_a[run] = 0; 
                     reactions[run].calculatePropensity(a, population, run);
                     for(int i=0;i<w_bin;i++){
                       total_a[n_species +run]+=a[run][i];
                     }
                       //printf("a[run][i] is %f\n",a[5][99]) ;
                }
                 
                
            }

            a0=0;
            for(int i =0;i<RULENUM;i++){
                a0 += total_a[i];  
            }
        
            do
            {
                r1 = (double)rand() / RAND_MAX;
            } while (r1 <= 0 || r1 >= 1);

            tau = -1.0/a0*log(r1);
           
            timeTracker += tau;  
           
            r2 = (double)rand() / RAND_MAX;
            sum = r2 * a0; 
            sum_a = total_a[0];
            /*
            select which reaction will happen
            */
            ita = 0;//index of reaction
              while (sum_a < sum)
            {
                ita++;
                sum_a += total_a[ita];
            }
            //sum_a -= total_a[ita];
            firings[ita]++;
          
            
            int target;//index of bin
            if (ita < n_species) 
            {     
                int totP = 0; 
                for (int i = 0; i < w_bin; i++) {
                    totP += population[ita][i];
                }
           
                target = -1; 
                double bin_ra0;
                int bin_sum = 0;
                if (r2 < 0.5) {// to left
                    bin_ra0 = r2 * (totP - population[ita][0]);
                    target = w_bin - 1;
                    bin_sum = population[ita][target];
                    while (bin_sum < bin_ra0 && target > 0) {
                        target--;
                        bin_sum += population[ita][target];
                    }
                     printf("target for left is %d\n",target);
                } else {
                    bin_ra0 = r2 * (totP - population[ita][w_bin - 1]);
                    target = 0;
                    bin_sum = population[ita][target];
                    while (bin_sum < bin_ra0 && target < w_bin - 2) {
                        target++;
                        bin_sum += population[ita][target];
                    }
                    printf("target for right is %d\n",target);
                }

                if (target != -1 && target > 0 && r2 < 0.5) { 
                    population[ita][target]--;
                    population[ita][target - 1]++;
                   
                } else if (target != -1 && target < w_bin - 1 && r2 >= 0.5) { 
                    population[ita][target]--;
                    population[ita][target + 1]++;
                  
                } else {
                    break;
                    //printf("Error\n");
                }

            }
            else{ 
               ita = ita -n_species;
               if (reactions[ita].type == 1) {
                    r2 = (sum - sum_a)/total_a[ita];
                    //i = floor(r2 * w_bin / a0); 
                    location[target]++;
                }else if (reactions[ita].type == 2) {
                    r2 = (sum - sum_a)/total_a[ita];
                    //i = 0;
                    sum_x = population[(ita - n_species) / 2][0];
                    sum = (r2 - (sum_a - total_a[ita])) / total_a[ita] * tot[(ita - n_species) / 2];
                    while (sum_x < sum && target < w_bin - 1) {
                        target++;
                        sum_x += population[(ita - n_species) / 2][target];
                    }
                    location[target]++;
                }else if (reactions[ita].type == 3) {
                    //i = 0;
                    sum_x = population[(ita - n_species) / 2 - 1][0] * (population[(ita - n_species) / 2 - 1][0] - 1) * population[(ita - n_species) / 2][0];
                    sum = (r2 - (sum_a - total_a[ita])) / total_a[ita] * reactions[ita].rateConstant * w_bin * w_bin;
                    while (sum_x < sum && target < w_bin - 1) {
                        target++;
                        sum_x += population[(ita - n_species) / 2 - 1][target] * (population[(ita - n_species) / 2 - 1][target] - 1) * population[(ita - n_species) / 2][target];
                    }
                    location[target]++;
                }
                //after finding the bin, apply statechange to the reaction base on user provied function
                stateChange(ita, ita, population, tot);
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