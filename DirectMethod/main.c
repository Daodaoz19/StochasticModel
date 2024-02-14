#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "propensity.h"
#include "ssa.h" 
#include "constant.h"

typedef void (*PropensityFunc)(const double*, double*);
void directMethod(int stoich_matrix[MAX_REACTIONS][MAX_SPECIES], PropensityFunc propfunc, double tspan[2], double* x0, double* t_out, double* x_out);

int main() {
    clock_t begin = clock();
    srand(time(NULL));

     FILE* file1 = fopen("schlogl.txt", "w");

    if (file1 == NULL) {
        perror("Error opening file");
        return 1;
    }

    double* t;
    double* x;
    
    t = (double*)malloc(MAX_TIMESTEPS * sizeof(double));
    x = (double*)malloc(MAX_TIMESTEPS * MAX_SPECIES * sizeof(double));
    
    for (int sim = 0; sim < 10000; sim++) { //sim: number of simulation
        
        directMethod(stoich_matrix, Schlogl_propensity, tspan, x0, t, x);

        for (int i = 0; i < MAX_TIMESTEPS; i++) {
            fprintf(file1, "%f %f\t", t[i], x[i * MAX_SPECIES]);
            fprintf(file1, "\n");
        }

    }

    fclose(file1);
    free(t);
    free(x);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;    
    printf("Time elapsed is %f seconds\n", time_spent);
    
    return 0;
}

void directMethod(int stoich_matrix[MAX_REACTIONS][MAX_SPECIES], PropensityFunc propfunc, double tspan[2], double* x0, double* t_out, double* x_out) {

    double T = tspan[0];
    double X[MAX_SPECIES];
    memcpy(X, x0, MAX_SPECIES * sizeof(double));

    int rxn_count = 0;

    /*
    (t_out)[rxn_count] = T;
    memcpy(x_out + rxn_count * MAX_SPECIES, X, MAX_SPECIES * sizeof(double));
    */
    
    double outputsize = (tspan[1] - tspan[0])/MAX_TIMESTEPS;
    double outputtime = outputsize + tspan[0];

    while (T < tspan[1]) {
        double a[MAX_REACTIONS];
        propfunc(X, a);

        double a0 = a[0] + a[1] + a[2] + a[3];
        if (a0 <= 0) break;

        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;
        double tau = -log(r1) / a0;

        double cumulativeSum = 0.0;
        int mu = -1;
        for (int i = 0; i < MAX_REACTIONS; i++) {
            cumulativeSum += a[i];
            if (cumulativeSum >= r2 * a0) {
                mu = i;
                break;
            }
        }


        if (T + tau > outputtime) {
//          printf("output: %f %f\n", *t_out, *x_out);
            (t_out)[rxn_count] = outputtime;
            memcpy(x_out + rxn_count * MAX_SPECIES, X, MAX_SPECIES * sizeof(double));
            rxn_count++;
            outputtime += outputsize;
        }


        T += tau;
        for (int i = 0; i < MAX_SPECIES; i++) {
            X[i] += stoich_matrix[mu][i];
        }
    }
    
}