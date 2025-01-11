#ifndef __VEC_H__
#define __VEC_H__

#include <math.h>
#include <stdlib.h>

#define d_ones(n) d_const(n, 1.0)
#define i_ones(n) d_const(n, 1)
#define d_zeros(n) d_const(n, 0.0)
#define i_zeros(n) i_const(n, 0)

double *d_const(int n, double c);

int *i_const(int n, int c);

double *d_random(int n);

int *i_random(int n);

int d_veceq(double *u, double *v, int n, double eps);

int i_veceq(int *u, int *v, int n);


#endif