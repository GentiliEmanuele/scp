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
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        double d = fabs(u[i] - v[i]);
        if (d > eps) {
            s += d;
        }
    }
    s /= n;
    if (s > 0.0)
        printf("error=%f\n", s);
    return s > 0.0 ? 0 : 1;
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

float std_dev(float *samples, float avg, int n) {
    float dev = 0.0;
    for (int i = 0; i < n; i++) {
        dev += (samples[i] - avg) * (samples[i] - avg);
    }
    return dev / n;
}