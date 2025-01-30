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
        v[i] = (double)rand() / RAND_MAX + 1;
    }
    return v;
    
}

inline int d_veceq(double *u, double *v, int n, double eps) {
    for (int i = 0; i < n; i++) {
        if (fabs(u[i] - v[i]) > eps) {
            printf("%d: %f != %f\n", u[i], v[i]);
            return 0;
        }
    }
    return 1;
}

int read_vector(double *vector, int n, const char *path) {
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        printf("cannot open file: %s\n", path);
        return 1;
    }
    for (int i = 0; i < n; i++) {
        fscanf(f, "%lg", &vector[i]);
    }
    return 0;
}

void print_vec(double *v, int n) {
    for (int i = 0; i < n; i++)
    {
        printf("%d %f\n", i, v[i]);
    }
}
