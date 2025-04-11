#ifndef SPMV_SEQ_H
#define SPMV_SEQ_H
#include "csr.h"
#include "hll.h"

int spmv_csr_seq(double *res, struct csr *csr, double *v, int n);
int spmv_hll_seq(double *res, struct hll *hll, double *v, int n);
#endif