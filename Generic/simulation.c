/*
main simulation part handle everything
Usage: ./program <input_filename> [NUMofRUNS] [SimulationTime]
*/
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

double diffusionRates[N] = DIFFUSION_RATE;
int reactionTypes[M] = REACTION_TYPE;
double rateConstants[M] = RATE_CONSTANTS;


int main(int argc, char *argv[])
{
    char* inputFilename;
    if (argc > 1) {
        inputFilename = argv[1];
        if (argc > 2) {
            NUMofRUNS = atoi(argv[2]);
        }
        if (argc > 3) {
            SimulationTime = atoi(argv[3]);
        }
    } else {
        fprintf(stderr, "Usage: %s <input_filename> [NUMofRUNS] [SimulationTime]\n", argv[0]);
        return 1; 
    }
 
    /*
    calculate h
    */
    double h[N];
    for (int i = 0; i < N; i++)
    {
        h[i] = L / W;
    }
 
    /*
    calculate jump rate
    */
    double jumpRate[N];
    for (int i = 0; i < N; i++)
    {
        jumpRate[i] = diffusionRates[i] / (h[i] * h[i]);
    }
   
    
    const int RULENUM = N + M;
    double a[RULENUM], a0;
    double r1, r2, tau;
    int ita;
    double sum, sum_a, sum_x;


    long firings[RULENUM], location[W];
    FILE *population = fopen("population", "w");
    FILE *total = fopen("total", "w");
    FILE *firefile = fopen("firings", "w");
    FILE *locationfile = fopen("location", "w"); 

    printf("Begin the model simulation ...\n");
    printf("BINNUM = %d\n", W);

    srand(time(NULL)); 
    for (int real = 0; real < NUMofRUNS; real++)
    {
        /*
        set total population for each species
        */
        long tot[N];
        for (int n = 0; n < N; n++)
        {
            tot[n] =  S;
        }

        /*
        set population value for each species in each bin from the input file
        */
        FILE *inputfile;
        inputfile = fopen(inputFilename, "r");

        long **species_population;
        species_population = (long **)malloc(N * sizeof(long *));
        for (int n = 0; n < N; n++)
        {    
            species_population[n] = (long *)malloc(W * sizeof(long));
        }
        
        for (int n = 0; n < N; n++) {
        for (int m = 0; m < W; m++) {
            fscanf(inputfile, "%ld", &species_population[n][m]);
        }
    }
        /*
        for writing output
        */
        for (int i = 0; i < RULENUM; i++)
            firings[i] = 0;
        for (int i = 0; i < W; i++)
            location[i] = 0;

        double timeTracker = 0.0; 
        while (timeTracker < SimulationTime)
        {
            /*
            calculate propensity for diffusion for each species
            */
              for (int n = 0; n < N; n++)
            {
                a[n] = jumpRate[n] * (tot[n]);
            
            }
            /*
            calculate propneisty for each reaction channel
            */
            for (int run = 0; run < M; run++)
            {
                if (reactionTypes[run] == 1)
                {
                    a[N + run] = rateConstants[run] * S;

                }
                else if (reactionTypes[run] == 2)
                {
                    a[N + run] = rateConstants[run] * tot[0];
                   
                }
                else if (reactionTypes[run] == 3)
                {

                    a[N + run]=prop3(W, species_population);
                   
                }
    
            }
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

            int i;//index of bin
            if (ita < N)
            {    
                if (ita % 2 == 0) //even number is the species diffuse to left
                {
                    r2 = (sum - sum_a) / a[ita];
                    sum = (tot[ita / 2] - species_population[(ita / 2)][0]) * r2;
                    sum_x = species_population[(ita / 2)][1];
                    i = 1;
                    while (sum_x < sum)
                    {
                        i++;
                        sum_x += species_population[(ita / 2)][i];
                    }
                    species_population[(ita / 2)][i]--;
                    species_population[(ita / 2)][i - 1]++;
                    location[i]++;  
                  
                }

                else// odd number is species diffuse to right
                { 
                    r2 = (sum - sum_a) / a[ita];
                    sum = (tot[ita / 2+1] - species_population[(ita / 2+1)][0]) * r2;
                    sum_x = species_population[(ita / 2+1)][1];
                    i = 1;
                    while (sum_x < sum)
                    {
                        i++;
                        sum_x += species_population[(ita / 2+1)][i];
                    }
                    species_population[(ita / 2+1)][i]--;
                    species_population[(ita / 2+1)][i - 1]++;
                    location[i]++;

                    i = 0;
                    r2 = (sum - sum_a) / a[ita];
                    sum = (tot[ita / 2+1] - species_population[ita / 2+1][W - 1]) * r2;
                    sum_x = species_population[ita / 2+1][0];

                    while (sum_x < sum)
                    {
                        i++;
                        sum_x += species_population[ita / 2+1][i];
                    }
                    species_population[ita / 2+1][i]--;
                    species_population[ita / 2+1][i + 1]++;
                    location[i]++;
                
                }
            }
            else
            {
                if (reactionTypes[ita - N] == 1)
                { 
                    r2 = (sum - sum_a) / a[ita];
                    i = r2 * W; 
                    species_population[(ita - N) / 2][i]++;
                    tot[(ita - N) / 2]++;
                    location[i]++;
                  
                }
                if (reactionTypes[ita - N] == 2)
                { 
                    r2 = (sum - sum_a) / a[ita];
                    sum_x = species_population[(ita - N) / 2][0];
                    sum = r2 * tot[(ita - N) / 2];
                    i = 0;
                    while (sum_x < sum)
                    {
                        i++;
                        sum_x += species_population[(ita -  N) / 2][i];
                    }

                    species_population[(ita -  N) / 2][i]--;
                    tot[(ita -N) / 2]--;
                    location[i]++;
                   
                }
                if (reactionTypes[ita - N] == 3)
                {
                    sum = (sum - sum_a) / (rateConstants[ita -  N] * W * W);
                    sum_x = species_population[(ita - N) / 2 - 1][0] * (species_population[(ita - N) / 2 - 1][0] - 1) * species_population[(ita - N) / 2][0];
                    i = 0;
                    while (sum_x < sum)
                    {
                        i++;
                        sum_x += species_population[(ita -  N) / 2 - 1][i] * (species_population[(ita - N) / 2 - 1][i] - 1) * species_population[(ita - N) / 2][i];
                    }
                    species_population[(ita -  N) / 2 - 1][i]++;
                    tot[(ita - N) / 2 - 1]++;
                    species_population[(ita - N) / 2][i]--;
                    tot[(ita - N) / 2]--;
                    location[i]++;
                  
                }
            }
        }

        for (int i = 0; i < N; i++)
        {
            if (i == N - 1)
            { 
                fprintf(total, "%ld\n", tot[i]);
            }
            else
            {
                fprintf(total, "%ld\t", tot[i]);
            }
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < W; j++)
            {
                fprintf(population, "%ld\t", species_population[i][j]);
            }
            fprintf(population, "\n");
        }
      

        for (int i = 0; i < RULENUM; i++)
            fprintf(firefile, "%ld\t", firings[i]);
        fprintf(firefile, "\n");
        for (int i = 0; i < W; i++)
            fprintf(locationfile, "%ld\t", location[i]);
        fprintf(locationfile, "\n");

         for (int n = 0; n < N; n++) {
        free(species_population[n]); 
    }
    free(species_population); 
    fclose(inputfile);
    }

    printf("End of simulation ...\n");

    fclose(population);
    fclose(total);
    fclose(firefile);
    fclose(locationfile);
    
   
    return 0;
}

