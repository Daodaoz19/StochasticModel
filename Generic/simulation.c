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
    double a[W][RULENUM], a0[W];
    double r1, r2[W], tau;
    int ita[W];
    double sum[W], sum_a[W], sum_x[W];


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
        Initialize total population for each species
        */
        long tot[N];
        for (int n = 0; n < N; n++)
        {
            tot[n] = 0;
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
            tot[n] += species_population[n][m]; 
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
              for (int j = 0; j < W; j++)
            {
                for(int n =0; n<N;n++){
                   a[j][n] = jumpRate[n] * (tot[n]);  
                }
            }
            /*
            calculate propneisty for each reaction channel
            */
            int Type2_index = 0;
            for (int j = 0; j < W; j++)
            {
            for (int run = 0; run < M; run++)
            {
                if (reactionTypes[run] == 1)
                {
                    a[j][N+run] = rateConstants[run] * S;

                }
                else if (reactionTypes[run] == 2)
                {
                    a[j][N+run] = rateConstants[run] * tot[Type2_index];
                    Type2_index+=1;
                }
                else if (reactionTypes[run] == 3)
                {
                   a[j][N+run]  = prop3(&a[j][N+run], species_population);
                   
                }//printf("a[j][N+run] %d\nis ",a[j][N+run]);
            }
            }
            /*
            sum total propensity
            */
            for(int i =0;i<W;i++){
                a0[i] = 0;  
            }

            for(int i =0;i<W;i++){
            for (int j = 0; j < M-1; j++){ 
                    a0[i] += a[i][j]; 
            
            } 
           }
        
            

            for(int i =0; i<W;i++){
                do
            {
                r1 = (double)rand() / RAND_MAX;
            } while (r1 <= 0 || r1 >= 1);
                tau = -1.0/a0[i]*log(r1); 

            }
            timeTracker += tau;

            for(int i=0;i<W;i++){
                  r2[i] = (double)rand() / RAND_MAX;
                  sum[i] = r2[i]*a0[i];
                  sum_a[i] = r2[i] * a[i][0];
                  ita[i] = 0;
                        
                    while (sum_a[i]< sum[i])
                    {
                        ita[i]++;
                        sum_a[i] +=a[i][ita[i]];
                    }
                  sum_a[i] -= a[i][ita[i]];  
                  firings[ita[i]]++;

                    int q;//index of bin
            if (ita[i] < N)
            {    
                if (ita[i] % 2 == 0) //even number is the species diffuse to left
                {
                    r2[i] = (sum - sum_a) / a[i][ita[i]];
                    sum[i] = (tot[ita[i] / 2] - species_population[(ita[i] / 2)][0]) * r2[i];
                    sum_x[i] = species_population[(ita[i] / 2)][1];
                    q = 1;
                    while (sum_x < sum)
                    {
                        q++;
                        sum_x[i] += species_population[(ita[i] / 2)][q];
                    }
                    species_population[(ita[i] / 2)][q]--;
                    species_population[(ita[i] / 2)][q - 1]++;
                    location[q]++;  
                  
                }

                else// odd number is species diffuse to right
                { 
                    r2[i] = (sum[i] - sum_a[i]) / a[i][ita[i]];
                    sum[i] = (tot[ita[i] / 2+1] - species_population[(ita[i] / 2+1)][0]) * r2[i];
                    sum_x[i] = species_population[(ita[i] / 2+1)][1];
                    q = 1;
                    while (sum_x[i] < sum[i])
                    {
                        q++;
                        sum_x[i] += species_population[(ita[i] / 2+1)][q];
                    }
                    species_population[(ita[i] / 2+1)][q]--;
                    species_population[(ita[i] / 2+1)][q - 1]++;
                    location[q]++;

                    q = 0;
                    r2[i] = (sum[i] - sum_a[i]) / a[i][ita[i]];
                    sum[i] = (tot[ita[i] / 2+1] - species_population[ita[i] / 2+1][W - 1]) * r2[i];
                    sum_x[i] = species_population[ita[i] / 2+1][0];

                    while (sum_x[i] < sum[i])
                    {
                        q++;
                        sum_x[i] += species_population[ita[i] / 2+1][q];
                    }
                    species_population[ita[i] / 2+1][q]--;
                    species_population[ita[i] / 2+1][q + 1]++;
                    location[q]++;
                
                }
            }
            else
            {
                if (reactionTypes[ita[i] - N] == 1)
                { 
                    r2[i] = (sum[i] - sum_a[i]) / a[i][ita[i]];
                    q = r2[i] * W; 
                    stateChangeAdd((ita[i] - N) / 2,q,species_population);
                    tot[(ita[i] - N) / 2]++;
                    location[q]++;
                  
                }
                if (reactionTypes[ita[i] - N] == 2)
                { 
                    r2[i] = (sum[i] - sum_a[i]) / a[i][ita[i]];
                    sum_x[i] = species_population[(ita[i] - N) / 2][0];
                    sum[i] = r2[i] * tot[(ita[i] - N) / 2];
                    q = 0;
                    while (sum_x[i] < sum[i])
                    {
                        q++;
                        sum_x[i] += species_population[(ita[i] -  N) / 2][q];
                    }

                    stateChangeMinus((ita[i] - N) / 2,q,species_population);
                    tot[(ita[i] -N) / 2]--;
                    location[q]++;
                   
                }
                if (reactionTypes[ita[i] - N] == 3)
                {
                    sum[i] = (sum[i] - sum_a[i]) / (rateConstants[ita[i] -  N] * W * W);
                    sum_x[i] = species_population[(ita[i] - N) / 2 - 1][0] * (species_population[(ita[i] - N) / 2 - 1][0] - 1) * species_population[(ita[i] - N) / 2][0];
                    q = 0;
                    while (sum_x < sum)
                    {
                        q++;
                        sum_x[i] += species_population[(ita[i] -  N) / 2 - 1][q] * (species_population[(ita[i] - N) / 2 - 1][q] - 1) * species_population[(ita[i] - N) / 2][q];
                    }
                    stateChangeAdd((ita[i] - N) / 2 -1,q,species_population);
                    tot[(ita[i] - N) / 2 - 1]++;
                    stateChangeMinus((ita[i] - N) / 2,q,species_population);
                    tot[(ita[i] - N) / 2]--;
                    location[q]++;
                  
                }
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

