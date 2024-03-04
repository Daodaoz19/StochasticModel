#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "model.h"

/*
These value should be input in command prompy by user
*/
//#define BINNUM 100
#define NUMofRUNS 1 
#define SimulationTime 10

double diffusionRates[N] = DIFFUSION_RATE;
int reactionTypes[M] = REACTION_TYPE;
double rateConstants[M] = RATE_CONSTANTS;
int sepciesRatio[N] = SPECIES_RATIO;
int BINNUM =100; //this value should also input in command prompy by user
int main(int argc, char *argv[])
{
    int pm = BINNUM / BINNUM;
    
    double w = 100; //scale factor, it is a fixed number?
    /*
    calculate h
    */
    double h[N];
    for (int i = 0; i < N; i++)
    {
        h[i] = (double)L / BINNUM;
    }
 
    /*
    calculate jump rate
    */
    double jumpRate[N];
    for (int i = 0; i < N; i++)
    {
        jumpRate[i] = diffusionRates[i] / (h[i] * h[i]);
    }
   
    const int RULENUM = 2 * N + M;
    double a[RULENUM], a0;
    double r1, r2, tau;
    int ita;
    double sum, sum_a, sum_x;

    /*
    for output
    */
    long firings[RULENUM], location[BINNUM];
    FILE *population = fopen("population", "w");
    FILE *total = fopen("total", "w");
    FILE *firefile = fopen("firings", "w");
    FILE *locationfile = fopen("location", "w"); 

    printf("Begin the model simulation ...\n");
    printf("XBINNUM = %d, YBINNUM = %d\n", BINNUM, BINNUM);

    for (int real = 0; real < NUMofRUNS; real++)
    {
        
        srand(time(NULL));
        /*
        set total population for each species
        */
        long tot[N];
        for (int n = 0; n < N; n++)
        {
            tot[n] = sepciesRatio[n] * w;
        }

        /*
        set population value for each species in each bin
        */
        long **space;
        space = (long **)malloc(N * sizeof(long *));
        for (int n = 0; n < N; n++)
        {    
            space[n] = (long *)malloc(BINNUM * sizeof(long));
            for (int m = 0; m < BINNUM; m++)
            {
                space[n][m] = sepciesRatio[n] * w / BINNUM;
            }
        }

        //printf("space is %ld",space[0][0]);
        /*
        for writing output
        */
        for (int i = 0; i < RULENUM; i++)
            firings[i] = 0;
        for (int i = 0; i < BINNUM; i++)
            location[i] = 0;

        double timeTracker = 0.0; 
        //printf("a0 is %f and %f \n", jumpRate[0],(tot[0] - space[0][0]));
        while (timeTracker < SimulationTime)
        {
            /*
            calculate propensity for diffusion for each species
            */
            for (int n = 0; n < N; n++)
            {
                a[2 * n] = jumpRate[n] * (tot[n] - space[n][0]);
                a[2 * n + 1] = jumpRate[n] * (tot[n] - space[n][BINNUM - 1]);
            }
            /*
            calculate propneisty for each reaction channel
            */
            for (int run = 0; run < M; run++)
            {
                if (reactionTypes[run] == 1)
                {
                    a[2 * N + run] = rateConstants[run] * w;

                }
                else if (reactionTypes[run] == 2)
                {
                    a[2 * N + run] = rateConstants[run] * tot[0];
                   
                }
                else if (reactionTypes[run] == 3)
                {
                    ;
                }
                else if (reactionTypes[run] == 4)
                {
                    a[2 * N + run] = 0;
                    for (int i = 0; i < BINNUM; i++)
                    {          
                        long temp = space[N - 2][i] * (space[N - 2][i] - 1) * space[N-1][i / pm];
                        //printf("%ld and %ld and %ld\n",space[N - 2][i],(space[N - 2][i] - 1),space[N-1][i / pm]);
                        a[2 * N + run] += temp; 
                       
                    }
                   
                     a[2 * N + run] *= rateConstants[run] * BINNUM * BINNUM;
                }
    
            }
            // printf("a0 is %f\n",a[0]);
            // printf("a1 is %f\n",a[1]);
            // printf("a2 is %f\n",a[2]);
            // printf("a3 is %f\n",a[3]);
            // printf("a4 is %f\n",a[4]);
            // printf("a5 is %f\n",a[5]);
            // printf("a6 is %f\n",a[6]);
            // printf("a7 is %f\n",a[7]);

            /*
            sum total propensity
            */
            a0 = 0;
            for (int i = 0; i < RULENUM; i++){
                  a0 += a[i];
            } 
        
            do
            {
                r1 = (double)rand() / RAND_MAX;
            } while (r1 <= 0 || r1 >= 1);

            tau = -1.0/a0*log(r1);
           
            timeTracker += tau;
            
            r2 = (double)rand() / RAND_MAX;

            sum = r2 * a0;
            sum_a = a[0];
            ita = 0;
            while (sum_a < sum)
            {
                ita++;
                sum_a += a[ita];
            }
            sum_a -= a[ita];
            firings[ita]++;

            int i;
            if (ita < 2 * N)
            {    
                if (ita % 2 == 0) //even number is the species diffuse to left
                {
                    r2 = (sum - sum_a) / a[ita];
                    sum = (tot[ita / 2] - space[(ita / 2)][0]) * r2;
                    sum_x = space[(ita / 2)][1];
                    i = 1;
                    while (sum_x < sum)
                    {
                        i++;
                        sum_x += space[(ita / 2)][i];
                    }
                    space[(ita / 2)][i]--;
                    space[(ita / 2)][i - 1]++;
                    location[i]++;
                  
                }

                else// odd number is species diffuse to right
                { 
                    i = 0;
                    r2 = (sum - sum_a) / a[ita];
                    sum = (tot[ita / 2] - space[ita / 2][BINNUM - 1]) * r2;
                    sum_x = space[ita / 2][0];

                    while (sum_x < sum)
                    {
                        i++;
                        sum_x += space[ita / 2][i];
                    }
                    space[ita / 2][i]--;
                    space[ita / 2][i + 1]++;
                    location[i]++;
                
                }
            }
            else
            {
                if (reactionTypes[ita - 2 * N] == 1)
                { 
                    r2 = (sum - sum_a) / a[ita];
                    i = r2 * BINNUM; // Find the index of the bin for the new X to appear at
                    space[(ita - 2 * N) / 2][i]++;
                    tot[(ita - 2 * N) / 2]++;
                    location[i]++;
                  
                }
                if (reactionTypes[ita - 2 * N] == 2)
                { 
                    r2 = (sum - sum_a) / a[ita];
                    sum_x = space[(ita - 2 * N) / 2][0];
                    sum = r2 * tot[(ita - 2 * N) / 2];
                    i = 0;
                    while (sum_x < sum)
                    {
                        i++;
                        sum_x += space[(ita - 2 * N) / 2][i];
                    }

                    space[(ita - 2 * N) / 2][i]--;
                    tot[(ita - 2 * N) / 2]--;
                    location[i]++;
                   
                }//rule type 3 omit for now
                if (reactionTypes[ita - 2 * N] == 4)
                {
                    sum = (sum - sum_a) / (rateConstants[ita - 2 * N] * BINNUM * BINNUM);
                    sum_x = space[(ita - 2 * N) / 2 - 1][0] * (space[(ita - 2 * N) / 2 - 1][0] - 1) * space[(ita - 2 * N) / 2][0];
                    i = 0;
                    while (sum_x < sum)
                    {
                        i++;
                        sum_x += space[(ita - 2 * N) / 2 - 1][i] * (space[(ita - 2 * N) / 2 - 1][i] - 1) * space[(ita - 2 * N) / 2][i / pm];
                    }
                    space[(ita - 2 * N) / 2 - 1][i]++;
                    tot[(ita - 2 * N) / 2 - 1]++;
                    space[(ita - 2 * N) / 2][i / pm]--;
                    tot[(ita - 2 * N) / 2]--;
                    location[i]++;
                  
                }
            }
        }

        for (int i = 0; i < N; i++)
        {
            if (i = N - 1)
            { // the last species
                fprintf(total, "%ld\n", tot[i]);
            }
            else
            {
                fprintf(total, "%ld\t", tot[i]);
            }
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < BINNUM; j++)
            {
                fprintf(population, "%ld\t", space[i][j]);
            }
        }
        fprintf(population, "\n");

        for (int i = 0; i < RULENUM; i++)
            fprintf(firefile, "%ld\t", firings[i]);
        fprintf(firefile, "\n");
        for (int i = 0; i < BINNUM; i++)
            fprintf(locationfile, "%ld\t", location[i]);
        fprintf(locationfile, "\n");

         for (int n = 0; n < N; n++) {
        free(space[n]); 
    }
    free(space); 
    }

    printf("End of simulation ...\n");

    fclose(population);
    fclose(total);
    fclose(firefile);
    fclose(locationfile);
    
   
    return 0;
}