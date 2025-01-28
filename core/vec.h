#ifndef __VEC_H__
#define __VEC_H__

#ifdef __cplusplus
extern "C"
{
#endif

#include <math.h>
#include <stdlib.h>

#define d_ones(n) d_const(n, 1.0)
#define d_zeros(n) d_const(n, 0.0)

double *d_const(int n, double c);

double *d_random(int n);

int d_veceq(double *u, double *v, int n, double eps);

int read_vector(double *vector, int n, const char *path);

void print_vec(double *v, int n);

#ifdef __cplusplus
} // extern "C"
#endif
#endif
