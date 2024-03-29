#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define MAX_SPECIES 10
#define MAX_REACTIONS 10
#define MAX_TIMESTEPS 10000

typedef void (*PropensityFunc)(const double*, double*);

int directMethod(int stoich_matrix[MAX_REACTIONS][MAX_SPECIES], PropensityFunc propfunc, double tspan[2], double* x0, double** t_out, double** x_out);
void Lotka_propensity(const double* y, double* a);
void Schlogl_propensity(const double* x, double* a);

int main() {
    clock_t begin = clock();
    srand(time(NULL));
    /**
     * lotka ssa 
    */
    // double tspan[2] = {0,50};
    // double x0[MAX_SPECIES] = {200, 200};
    // int stoich_matrix[MAX_REACTIONS][MAX_SPECIES] = {
    //     {1, 0},
    //     {-1, 1},
    //     {0, -1}
    // };  
    /**
     * schlogl ssa
    */
    double tspan[2] = {0, 10};
    double x0[MAX_SPECIES] = {250};
    int stoich_matrix[MAX_REACTIONS][MAX_SPECIES] = {
    {1},
    {-1},
    {1},
    {-1}
    };

 
    FILE* file1 = fopen("X.txt", "w");
    FILE* file2 = fopen("Y.txt", "w");
    FILE* file3 = fopen("phasePlot.txt", "w");

    if (file1 == NULL || file2 == NULL || file3 == NULL) {
        perror("Error opening file");
        return 1;
    }

    for (int sim = 0; sim < 1; sim++) { //sim: number of simulation
        double* t;
        double* x;
        int rxn_count = directMethod(stoich_matrix, Schlogl_propensity, tspan, x0, &t, &x);

    // for (int i = 0; i < MAX_TIMESTEPS && t[i] <= tspan[1]; i++) {
        
    //     printf("Time: %f, Species 1: %f, Species 2: %f\n", t[i], x[i * MAX_SPECIES], x[i * MAX_SPECIES + 1]);
    // }
        for (int i = 0;i <= rxn_count; i++) {
            
            fprintf(file1, "%f\t%f\n", t[i], x[i * MAX_SPECIES]);           
            fprintf(file2, "%f\t%f\n", t[i], x[i * MAX_SPECIES+1]);
            fprintf(file3, "%f\t%f\n", x[i * MAX_SPECIES], x[i * MAX_SPECIES + 1]);
        }

        free(t);
        free(x);
    }

    fclose(file1);
    fclose(file2);
    fclose(file3);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;    
    printf("Time elapsed is %f seconds\n", time_spent);
    
    return 0;
}

void Lotka_propensity(const double* y, double* a) {
    double k1 = 1.0;
    double k2 = 0.003;
    double k3 = 1.0;
    
    a[0] = k1 * y[0];
    a[1] = k2 * y[0] * y[1];
    a[2] = k3 * y[1];
}

void Schlogl_propensity(const double* x, double* a) {
    double C1 = 3e-7, C2 = 1e-4, C3 = 1e-3, C4 = 3.5;
    double N1 = 1e5, N2 = 2e5;

    a[0] = C1 / 2 * N1 * x[0] * (x[0] - 1);
    a[1] = C2 / 6 * x[0] * (x[0] - 1) * (x[0] - 2);
    a[2] = C3 * N2;
    a[3] = C4 * x[0];
}



int directMethod(int stoich_matrix[MAX_REACTIONS][MAX_SPECIES], PropensityFunc propfunc, double tspan[2], double* x0, double** t_out, double** x_out) {

    double T = tspan[0];
    double X[MAX_SPECIES];
    memcpy(X, x0, MAX_SPECIES * sizeof(double));

    int rxn_count = 0;
    *t_out = (double*)malloc(MAX_TIMESTEPS * sizeof(double));
    *x_out = (double*)malloc(MAX_TIMESTEPS * MAX_SPECIES * sizeof(double));

    (*t_out)[rxn_count] = T;
    memcpy(*x_out + rxn_count * MAX_SPECIES, X, MAX_SPECIES * sizeof(double));

    double outputsize = (tspan[1] - tspan[0]) / 1000.0;
    double outputtime = outputsize + tspan[0];

    while (T < tspan[1]) {
        double a[MAX_REACTIONS];
        propfunc(X, a);

        double a0 = a[0] + a[1] + a[2];
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

        if (mu == -1) continue;

        if (T + tau > outputtime) {
            rxn_count++;
            (*t_out)[rxn_count] = T;
            memcpy(*x_out + rxn_count * MAX_SPECIES, X, MAX_SPECIES * sizeof(double));
            outputtime += outputsize;
        }

        T += tau;
        for (int i = 0; i < MAX_SPECIES; i++) {
            X[i] += stoich_matrix[mu][i];
        }  
      //printf("outputtime %f\n",outputtime);
    }
   

    if (T >= tspan[1]) {
        rxn_count++;
        (*t_out)[rxn_count] = T;
        memcpy(*x_out + rxn_count * MAX_SPECIES, X, MAX_SPECIES * sizeof(double));
    }
    return  rxn_count;
    //printf("%d",rxn_count);
}