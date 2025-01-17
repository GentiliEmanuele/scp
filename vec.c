#include "vec.h"
#include <stdlib.h>
#include <stdio.h>

inline double *d_const(int n, double c) {
    double *v = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        v[i] = c;
    }
    return v;
    
}

inline double *d_random(int n) {
    double *v = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        v[i] = (double)rand() / RAND_MAX;
    }
    return v;
    
}

inline int d_veceq(double *u, double *v, int n, double eps) {
    for (int i = 0; i < n; i++) {
        if (fabs(u[i] - v[i]) > eps) {
            return 0;
        }
    }
    return 1;
}

