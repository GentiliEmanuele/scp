#ifndef __VEC_H__
#define __VEC_H__

#include <math.h>
#include <stdlib.h>

#define d_ones(n) d_const(n, 1.0)
#define d_zeros(n) d_const(n, 0.0)

double *d_const(int n, double c);

double *d_random(int n);

int d_veceq(double *u, double *v, int n, double eps);

#endif