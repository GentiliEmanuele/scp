#include "spmv_openmp.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

struct component {
    int pos;
    union {
        int i_val;
        double d_val;
    };
};

static int init_component_arrays(struct component ***out, int num_threads, int chunk_size) {
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
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }

    int num_threads;
    struct component **results;
#pragma omp parallel
{
    if (omp_get_thread_num() == 0) {
        num_threads = omp_get_num_threads();
        // printf("executing sparse matrix vector product with %d threads\n", omp_get_num_threads());
    }
}

    int chunk_size = ceil(((double)csr->num_rows) / num_threads);
    if (init_component_arrays(&results, num_threads, chunk_size)) {
        printf("malloc failed while allocating space for temporary results\n");
        return 1;
    }
    
#pragma omp parallel
{
    struct component *tmp_res = results[omp_get_thread_num()];
    int k = 0;
    #pragma omp for schedule(static)
    for (int i = 0; i < csr->num_rows; ++i) {
        tmp_res[k].pos = i;
        tmp_res[k].d_val = 0.0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            tmp_res[k].d_val += data[j] * v[csr->col_index[j]];
        }
        ++k;
    }
    tmp_res[k].pos = -1;
}

    for (int i = 0; i < num_threads; ++i) {
        struct component *p = results[i];
        for (int j = 0; j < chunk_size && p[j].pos >= 0; ++j) {
            res[p[j].pos] = p[j].d_val;
        }
    }
    return 0;
}

int i_spmv_csr_par(int *res, struct csr *csr, int *v, int n) {
    int *data = (int*)csr->data;
    if (n != csr->num_cols) {
        printf("matrix and vector must have the same number of columns\n");
        return 1;
    }

    int num_threads;
    struct component **results;
#pragma omp parallel
{
    if (omp_get_thread_num() == 0) {
        num_threads = omp_get_num_threads();
        // printf("executing sparse matrix vector product with %d threads\n", omp_get_num_threads());
    }
}

    int chunk_size = ceil(((double)csr->num_rows) / num_threads);
    if (init_component_arrays(&results, num_threads, chunk_size)) {
        printf("malloc failed while allocating space for temporary results\n");
        return 1;
    }

    #pragma omp parallel
{
    struct component *tmp_res = results[omp_get_thread_num()];
    int k = 0;
    #pragma omp for schedule(static)
    for (int i = 0; i < csr->num_rows; ++i) {
        tmp_res[k].pos = i;
        tmp_res[k].i_val = 0;
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            tmp_res[k].i_val += data[j] * v[csr->col_index[j]];
        }
        ++k;
    }
    tmp_res[k].pos = -1;
}

    for (int i = 0; i < num_threads; ++i) {
        struct component *p = results[i];
        for (int j = 0; j < chunk_size && p[j].pos >= 0; ++j) {
            res[p[j].pos] = p[j].i_val;
        }
    }
    for (int i = 0; i < num_threads; i++) {
        free(results[i]);
    }
    
    return 0;
}