#include "propensity.h"

void Schlogl_propensity(const double* x, double* a) {
    double C1 = 3e-7, C2 = 1e-4, C3 = 1e-3, C4 = 3.5;
    double N1 = 1e5, N2 = 2e5;

    a[0] = C1 / 2 * N1 * x[0] * (x[0] - 1);
    a[1] = C2 / 6 * x[0] * (x[0] - 1) * (x[0] - 2);
    a[2] = C3 * N2;
    a[3] = C4 * x[0];
}

void Lotka_propensity(const double* y, double* a) {
    double k1 = 1.0;
    double k2 = 0.003;
    double k3 = 1.0;
    
    a[0] = k1 * y[0];
    a[1] = k2 * y[0] * y[1];
    a[2] = k3 * y[1];
}
