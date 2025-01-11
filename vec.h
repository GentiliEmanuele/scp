#ifndef __VEC_H__
#define __VEC_H__

#include <math.h>
#include <stdlib.h>

#define d_ones(n) d_const(n, 1.0)
#define i_ones(n) d_const(n, 1)
#define d_zeros(n) d_const(n, 0.0)
#define i_zeros(n) i_const(n, 0)

inline double *d_const(int n, double c) {
    double *v = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        v[i] = c;
    }
    return v;
    
}

inline int *i_const(int n, int c) {
    int *v = malloc(n * sizeof(int));
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

inline int *i_random(int n) {
    int *v = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        v[i] = rand();
    }
    return v;
}

inline int d_veceq(double *u, double *v, int n, double eps) {
    for (int i = 0; i < n; i++) {
        if (fabs(u[i] - v[i]) > eps)
            return 0;
    }
    return 1;
}

inline int i_veceq(int *u, int *v, int n) {
    for (int i = 0; i < n; i++) {
        if (u[i] != v[i])
            return 0;
    }
    return 1;
}

#endif