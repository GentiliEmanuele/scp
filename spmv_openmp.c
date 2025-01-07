#include "spmv_openmp.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

static struct component *tmp_res;
#pragma omp threadprivate(tmp_res)

int i_spmv_csr_par(int *res, struct csr *csr, int *v, int n) {
    int *data = (int*)csr->data;
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }

    int **results = malloc(omp_get_num_threads() * sizeof(int*));

    int tid;
#pragma omp parallel private(tid)
{
    tid = omp_get_thread_num();
    printf("Thread id %d\n", tid);
    if (tid == 0) {
        printf("number of threads: %d\n", omp_get_num_threads());
    }
}
    int j;
    #pragma omp for private(j) schedule(static)
    // int *tmp = (int *)tmp_res;
    for (int i = 0; i < csr->num_rows; ++i) {
        res[i] = 0;
        for (j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            //tmp[i] += data[j] * v[csr->col_index[j]];
        }
    }
    return 0;
}

struct component {
    int pos;
    union {
        int i_val;
        double d_val;
    };
};

int __init_component_arrays(struct component ***out, int num_threads, int chunk_size) {
    struct component **arrays = malloc(num_threads * sizeof(struct component*));
    if (arrays == NULL) {
        return 1;
    }
    for (int i = 0; i < num_threads; ++i) {
        arrays[i] = malloc((chunk_size + 1) * sizeof(struct component));
        if (arrays[i] == NULL) {
            printf("malloc failed while allocating space for temporary results\n");
            for (int j = 0; j < i; ++j) {
                free(arrays[j]);
            }
            return 1;
        }
    }
    *out = arrays;
    return 0;
}

int d_spmv_csr_par(double *res, struct csr *csr, double *v, int n) {
    double *data = (double*)csr->data;
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns");
        return 1;
    }
    int num_threads;

    struct component **results;
#pragma omp parallel
{
    if (omp_get_thread_num() == 0) {
        num_threads = omp_get_num_threads();
        printf("number of threads: %d\n", omp_get_num_threads());
    }
}

if (__init_component_arrays(&results, num_threads,  csr->num_cols / num_threads)) {
    printf("malloc failed while allocating space for temporary results\n");
    return 1;
}
    
#pragma omp parallel
{
    tmp_res = results[omp_get_thread_num()];
    int k = 0;
    struct component *p_res = tmp_res;
    int j;
    #pragma omp for private(j) schedule(static)
    for (int i = 0; i < csr->num_rows; ++i) {
        p_res[k].pos = i;
        p_res[k].d_val = 0.0;
        for (j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            p_res[k].d_val += data[j] * v[csr->col_index[j]];
        }
        ++k;
    }
    p_res[k].pos = -1;
}

    for (int i = 0; i < num_threads; ++i) {
        struct component *p = results[i];
        for (int j = 0; p[j].pos > 0; ++j) {
            res[p[j].pos] = p[j].d_val;
        }
    }

    return 0;
}