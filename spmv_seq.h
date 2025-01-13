#ifndef SPMV_SEQ_H
#define SPMV_SEQ_H
#include "csr.h"
#include "hll.h"

int d_spmv_csr_seq(double *res, struct csr *csr, double *v, int n);
int i_spmv_csr_seq(int *res, struct csr *csr, int *v, int n);
#endif