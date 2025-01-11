#include "csr.h"
#include "hll.h"
#include "omp_test.h"
#include "spmv_openmp.h"
#include "spmv_seq.h"
#include "utils.h"
#include "vec.h"

#define EPS (1e-6)

int test_csr(const char *file) {
    struct MatrixMarket mm;
    if (read_mtx(file, &mm)) {
        return 1;
    }

    struct csr sm;
    if (csr_init(&sm, &mm)) {
        return 1;
    }

    int n = sm.num_cols;
    double *v = d_random(n);
    double *p = d_zeros(n);
    double *s = d_zeros(n);
    if (d_spmv_csr_par(p, &sm, v, n)) {
        return 1;
    }
    if (d_spmv_csr_seq(p, &sm, v, n)) {
        return 1;
    }
    if (!d_veceq(p, s, n, EPS)) {
        printf("(csr) test failed for matrix: %s\n", file);
        return 1;
    }
    return 0;
}

int test_hll(const char *file, int hack_size) {
    struct MatrixMarket mm;
    if (read_mtx(file, &mm)) {
        return 1;
    }

    struct csr csr;
    if (csr_init(&csr, &mm)) {
        return 1;
    }

    struct hll sm;
    if (hll_init(&sm, hack_size, &mm)) {
        return 1;
    }

    int n = sm.num_cols;
    double *v = d_random(n);
    double *p = d_zeros(n);
    double *s = d_zeros(n);
    if (d_spmv_hll_par(p, &sm, v, n)) {
        return 1;
    }
    if (d_spmv_csr_seq(p, &csr, v, n)) {
        return 1;
    }
    if (!d_veceq(p, s, n, EPS)) {
        printf("(csr) test failed for matrix: %s\n", file);
        return 1;
    }
    return 0;
}