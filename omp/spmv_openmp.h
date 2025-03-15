#ifndef __SPMV_OPENMP_H__
#define __SPMV_OPENMP_H__

#ifdef __cplusplus
extern "C"
{
#endif

#include "csr.h"
#include "hll.h"

int spmv_csr_par(double *res, struct csr *csr, double *v, int n);
int spmv_csr_par_v2(double *res, struct csr *csr, double *v, int n);
int spmv_hll_par(double *res, struct hll *hll, double *v, int n);
int spmv_hll_par_v2(double *res, struct hll *hll, double *v, int n);

#ifdef __cplusplus
} // extern "C"
#endif
#endif
