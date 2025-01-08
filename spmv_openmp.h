#ifndef __SPMV_OPENMP_H__
#define __SPMV_OPENMP_H__
#include "csr.h"

int d_spmv_csr_par(double *res, struct csr *csr, double *v, int n);
int d_spmv_csr_par_slow(double *res, struct csr *csr, double *v, int n);
int i_spmv_csr_par(int *res, struct csr *csr, int *v, int n);
#endif
