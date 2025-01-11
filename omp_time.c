#include "omp_time.h"
#include "csr.h"
#include "spmv_openmp.h"
#include "utils.h"
#include "vec.h"
#include <omp.h>
#include <stdlib.h>
#include <time.h>

void omp_time_csr(const char *file, int num_runs, int num_threads) {
    struct MatrixMarket mm;
    if (read_mtx(file, &mm)) {
        return;
    }

    struct csr sm;
    if (csr_init(&sm, &mm)) {
        return;
    }
    
    int n = sm.num_cols;
    double *r = d_zeros(n);
    srand(42);
    double *v = d_random(n);
    
    omp_set_num_threads(num_threads);
    double sum = 0;
    for (int i = 0; i < num_runs; i++) {
        clock_t start = clock();
        if (d_spmv_csr_par(r, &sm, v, n)) {
            printf("warning: couldn't complete sparse matrix-vector product of run %d\n", i);
        }
        clock_t end = clock();
        sum += (end - start) / (double)CLOCKS_PER_SEC;  
    }

    double mean_time = sum / num_runs;
    double flops = (2 * mm.nz) / (double)mean_time;
    printf("mean time: %f s\n", mean_time);
    printf("MEGA FLOPS=%f\n", flops * 1e-6);
}

void omp_time_hll(const char *file, int hack_size, int num_runs, int num_threads) {
    struct MatrixMarket mm;
    if (read_mtx(file, &mm)) {
        return;
    }

    struct hll sm;
    if (hll_init(&sm, hack_size, &mm)) {
        return;
    }
    
    int n = sm.num_cols;
    double *r = d_zeros(n);
    srand(42);
    double *v = d_random(n);
    
    omp_set_num_threads(num_threads);
    double sum = 0;
    for (int i = 0; i < num_runs; i++) {
        clock_t start = clock();
        if (d_spmv_hll_par(r, &sm, v, n)) {
            printf("warning: couldn't complete sparse matrix-vector product of run %d\n", i);
        }
        clock_t end = clock();
        sum += (end - start) / (double)CLOCKS_PER_SEC;  
    }

    double mean_time = sum / num_runs;
    double flops = (2 * mm.nz) / (double)mean_time;
    printf("mean time: %f s\n", mean_time);
    printf("MEGA FLOPS=%f\n", flops * 1e-6);
}