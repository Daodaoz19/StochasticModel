#include "ssa.h" // Include the header file for declarations
#include "constant.h"
// Definitions of global variables
double tspan[2] = {0, 20};
double x0[MAX_SPECIES] = {250}; // Ensure MAX_SPECIES is defined or included from a common header
int stoich_matrix[MAX_REACTIONS][MAX_SPECIES] = {
    {1},
    {-1},
    {1},
    {-1}
}; // Ensure MAX_REACTIONS and MAX_SPECIES are defined or included
