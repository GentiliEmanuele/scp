#include "csr.h"
#include "fmt.h"
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
        mtx_cleanup(&mm);
        return 1;
    }

    // write_csr_mtx(&sm, &mm);
    int m = sm.num_rows;
    int n = sm.num_cols;
    double *v = d_random(n);
    double *p = d_zeros(m);
    double *s = d_zeros(m);
    if (d_spmv_csr_par(p, &sm, v, n)) {
        mtx_cleanup(&mm);
        csr_cleanup(&sm);
        return 1;
    }
    if (d_spmv_csr_seq(s, &sm, v, n)) {
        mtx_cleanup(&mm);
        csr_cleanup(&sm);
        return 1;
    }
    if (!d_veceq(p, s, n, EPS)) {
        mtx_cleanup(&mm);
        csr_cleanup(&sm);
        printf("(csr) test failed for matrix: %s\n", file);
        return 1;
    }
    mtx_cleanup(&mm);
    csr_cleanup(&sm);
    return 0;
}

int test_hll(const char *file, int hack_size) {
    struct MatrixMarket mm;
    if (read_mtx(file, &mm)) {
        return 1;
    }

    struct csr csr;
    if (csr_init(&csr, &mm)) {
        mtx_cleanup(&mm);
        return 1;
    }

    struct hll sm;
    if (hll_init(&sm, hack_size, &mm)) {
        csr_cleanup(&csr);
        mtx_cleanup(&mm);
        return 1;
    }

    int m = sm.num_rows;
    int n = sm.num_cols;
    double *v = d_random(n);
    double *p = d_zeros(m);
    double *s = d_zeros(m);
    if (d_spmv_hll_par(p, &sm, v, n)) {
        csr_cleanup(&csr);
        mtx_cleanup(&mm);
        hll_cleanup(&sm);
        return 1;
    }
    if (d_spmv_csr_seq(s, &csr, v, n)) {
        csr_cleanup(&csr);
        mtx_cleanup(&mm);
        hll_cleanup(&sm);
        return 1;
    }
    if (!d_veceq(p, s, n, EPS)) {
        printf("(hll) test failed for matrix: %s\n", file);
        csr_cleanup(&csr);
        mtx_cleanup(&mm);
        hll_cleanup(&sm);
        return 1;
    }
    write_hll(&sm, &mm);
    csr_cleanup(&csr);
    mtx_cleanup(&mm);
    hll_cleanup(&sm);
    return 0;
}