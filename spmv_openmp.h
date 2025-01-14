#ifndef __SPMV_OPENMP_H__
#define __SPMV_OPENMP_H__
#include "csr.h"
#include "hll.h"

int d_spmv_csr_par(double *res, struct csr *csr, double *v, int n);
int i_spmv_csr_par(int *res, struct csr *csr, int *v, int n);
int d_spmv_hll_par(double *res, struct hll *hll, double *v, int n);
int i_spmv_hll_par(int *res, struct hll *hll, int *v, int n);
int d_spmv_hll_par2(double *res, struct hll *hll, double *v, int n);

#endif
