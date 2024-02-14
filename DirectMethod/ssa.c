#include "ssa.h" 
#include "constant.h"

double tspan[2] = {0, 20};
double x0[MAX_SPECIES] = {250}; 
int stoich_matrix[MAX_REACTIONS][MAX_SPECIES] = {
    {1},
    {-1},
    {1},
    {-1}
}; 
