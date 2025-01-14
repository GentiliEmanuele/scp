#include "omp_time.h"
#include "csr.h"
#include "spmv_openmp.h"
#include "utils.h"
#include "vec.h"
#include <omp.h>
#include <stdlib.h>
#include <time.h>

/**
 * @param file (in)                 path to the matrix
 * @param num_runs (in)             number of times the measurement must be performed
 * @param num_thread (in)           number of thread used in the executions
 * @param time_measurement (out)    struct that contain the measurement
 * */  
int omp_time_csr(const char *file, int num_runs, int num_threads, time_measurement_t *time_measurement) {
    struct MatrixMarket mm;
    if (read_mtx(file, &mm)) {
        return 1;
    }

    struct csr sm;
    if (csr_init(&sm, &mm)) {
        mtx_cleanup(&mm);
        return 1;
    }
    
    int m = sm.num_rows;
    int n = sm.num_cols;
    double *r = d_zeros(m);
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
    time_measurement -> mean_time = sum / num_runs;
    time_measurement -> flops = (2 * mm.nz) / (double)time_measurement->mean_time;
    time_measurement->num_threads = num_threads;
    time_measurement->num_runs = num_runs;
    csr_cleanup(&sm);
    mtx_cleanup(&mm);
}

/**
 * @param file (in)                 path to the matrix
 * @param hack_size                 number of rows of each hack
 * @param num_runs (in)             number of times the measurement must be performed
 * @param num_thread (in)           number of thread used in the executions
 * @param time_measurement (out)    struct that contain the measurement
 * */  
int omp_time_hll(const char *file, int hack_size, int num_runs, int num_threads, time_measurement_t *time_measurement) {
    struct MatrixMarket mm;
    if (read_mtx(file, &mm)) {
        return 1;
    }

    struct hll sm;
    if (hll_init(&sm, hack_size, &mm)) {
        mtx_cleanup(&mm);
        return 1;
    }

    int m = sm.num_rows;
    int n = sm.num_cols;
    double *r = d_zeros(m);
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
    time_measurement -> mean_time = sum / num_runs;
    time_measurement -> flops = (2 * mm.nz) / (double)time_measurement->mean_time;
    time_measurement->num_threads = num_threads;
    time_measurement->num_runs = num_runs;
    mtx_cleanup(&mm);
    hll_cleanup(&sm);
}