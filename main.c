#include "mmio.h"
#include "csr.h"
#include "hll.h"
#include "utils.h"
#include "spmv_seq.h"
#include "spmv_openmp.h"
#include <errno.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void d_write_mrk_mtx(struct MatrixMarket *mm) {
    for (int i = 0; i < mm->nz; ++i) {
        printf("%d %d %f\n", mm->rows[i], mm->cols[i], ((double*)mm->data)[i]);
    }
}

void i_write_mrk_mtx(struct MatrixMarket *mm) {
    for (int i = 0; i < mm->nz; ++i) {
        printf("%d %d %d\n", mm->rows[i], mm->cols[i], ((int*)mm->data)[i]);
    }
}

void write_csr_mtx(struct csr *csr, struct MatrixMarket *mm) {
    for (int i = 0; i < csr->num_rows; ++i) {
        int num_cols = csr->row_pointer[i+1] - csr->row_pointer[i];
        printf("row no. = %d has %d non-zero elements\n", i+1, num_cols);
    }
    for (int i = 0; i < csr->num_rows; ++i) {
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            int col = csr->col_index[j];
            if (mm_is_integer(mm->typecode)) {
                int v = ((int*)csr->data)[j];
                printf("%d %d %d\n", i, col, v);
            } else {
                double v = ((double*)csr->data)[j];
                printf("%d %d %f\n", i, col, v);
            }
        }
    }
}

void write_hll(struct hll *hll, struct MatrixMarket *mm) {
    printf("offsets:\n");
    for (int i = 0; i < hll->offsets_num; i++) {
        printf( i < hll->offsets_num - 1 ? "%d, " : "%d", hll->offsets[i]);
    }
    printf("\n");
    printf("indices:\n");
    int offp = 0;
    for (int i = 0; i < hll->data_num; ++i) {
        if (offp < hll->offsets_num && i == hll->offsets[offp]) {
            offp++;
            printf("|");
        }
        printf("%d ", hll->col_index[i]);
    }
    printf("|\n");
    printf("values:\n");
    offp = 0;
    for (int i = 0; i < hll->data_num; ++i) {
        if (offp < hll->offsets_num && i == hll->offsets[offp]) {
            offp++;
            printf("|");
        }
        if (mm_is_integer(mm->typecode)) {
            printf("%d ", ((int*)hll->data)[i]);
        } else {
            printf("%f ", ((double*)hll->data)[i]);
        }
    }
    printf("|\n");
}

int d_csr_test(struct csr sm, int num_iterations, struct MatrixMarket mm, int verbose) {
    double *r = malloc(sm.num_cols * sizeof(double));
    double *s = malloc(sm.num_cols * sizeof(double));
    double *v = d_random(sm.num_cols);
    for (int i = 0; i < sm.num_cols; i++) {
        v[i] = 1.0;
    }
    double sum_of_times = 0;

    for (int i = 0; i < num_iterations; i++) {
        clock_t start = clock();
        if (d_spmv_csr_seq(r, &sm, v, sm.num_cols)) {
            return 1;
        }
        clock_t end = clock();
        sum_of_times += (end - start) / (double)CLOCKS_PER_SEC;  
    }

    double mean_time = sum_of_times / num_iterations;
    double flops = (2 * mm.nz) / (double)mean_time;
    printf("mean time: %f s\n", mean_time);
    printf("MEGA FLOPS=%f\n", flops * 1e-6);

    if (d_spmv_csr_seq(s, &sm, v, sm.num_cols)) {
        return 1;
    }
    if (d_veceq(r, s, sm.num_cols, 1e-6)) {
        printf("test failed\n");
    }
    if (verbose) {
        for (int i = 0; i < sm.num_cols; ++i) {
            printf("%f ", r[i]);
        }
        printf("\n");
    }
    return 0;
}

int i_csr_test(struct csr sm, int num_iterations, struct MatrixMarket mm, int verbose) {
    int *r = malloc(sm.num_cols * sizeof(int));
    int *s = malloc(sm.num_cols * sizeof(int));
    int *v = i_random(sm.num_cols);
    for (int i = 0; i < sm.num_cols; i++) {
        v[i] = 1;
    }

    double sum_of_times = 0;

    for (int i = 0; i < num_iterations; i++) {
        clock_t start = clock();
        if (i_spmv_csr_seq(r, &sm, v, sm.num_cols)) {
            return 1;
        }
        clock_t end = clock();
        sum_of_times += (end - start) / (double)CLOCKS_PER_SEC;  
    }

    double mean_time = sum_of_times / num_iterations;
    double flops = (2 * mm.nz) / (double)mean_time;
    printf("mean time: %f s\n", mean_time);
    printf("MEGA FLOPS=%f\n", flops * 1e-6);

    if (i_spmv_csr_seq(s, &sm, v, sm.num_cols)) {
        return 1;
    }
    if (i_veceq(r, s, sm.num_cols)) {
        printf("test failed\n");
    }
    if(verbose) {
        for (int i = 0; i < sm.num_cols; ++i) {
            printf("%d ", r[i]);
        }
        printf("\n");   
    }
    return 0;
}

int d_hll_test(struct hll sm, int num_iterations, struct MatrixMarket mm, int verbose) {
    double *r = malloc(sm.num_cols * sizeof(double));
    double *v = d_random(sm.num_cols);
    for (int i = 0; i < sm.num_cols; i++) {
        v[i] = 1.0;
    }
    double sum_of_times = 0;

    for (int i = 0; i < num_iterations; i++) {
        clock_t start = clock();
        if (d_spmv_hll_seq(r, &sm, v, sm.num_cols)) {
            return 1;
        }
        clock_t end = clock();
        sum_of_times += (end - start) / (double)CLOCKS_PER_SEC;  
    }

    double mean_time = sum_of_times / num_iterations;
    double flops = (2 * mm.nz) / (double)mean_time;
    printf("mean time: %f s\n", mean_time);
    printf("MEGA FLOPS=%f\n", flops * 1e-6);

    if (verbose) {
        for (int i = 0; i < sm.num_cols; ++i) {
            printf("%f ", r[i]);
        }
        printf("\n");
    }
    return 0;
}

void print_hll(struct hll hll) {
    double *data = hll.data;
    for (int i = 0; i < hll.data_num; i++) {
        printf("% d : %f \t", i, data[i]);
    }
    printf("\n");
}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        printf("see usage: program num_threads num_runs matrix\n");
        return 1;
    }
    long num_threads = strtol(argv[--argc], NULL, 10);
    if (num_threads == 0) {
        printf("an error occurred while parsing string %s to long\n", argv[argc]);
        return 1;
    }
    long num_iterations = strtol(argv[--argc], NULL, 10);
    if (num_iterations == 0) {
        printf("an error occured while parsing string %s to long", argv[argc]);
        return 1;
    }
    struct MatrixMarket mm;
    if (read_mtx(argv[--argc], &mm)) {
        printf("cannot read matrix: %s\n", argv[argc]);
        return 1;
    }
    printf("matrix has %d rows and %d cols and number of non-zeros %d\n", mm.num_rows, mm.num_cols, mm.nz);
  
    d_write_mrk_mtx(&mm);
    struct hll sm;
    if (hll_init(&sm, 2, &mm)) {
        printf("cannot read matrix into HLL format\n");
        return 1;
    }

    if (d_hll_test(sm, num_iterations, mm, 1)) {
        printf("An error occured while test the funtion for double matrix\n");
        return 1;
    }


/*
    struct csr sm;
    if (csr_init(&sm, &mm)) { 
        printf("cannot read matrix into CSR format\n");
        return 1;
    }
    
    write_csr_mtx(&sm, &mm);
    omp_set_num_threads(num_threads);
    srand(42);
    if (mm_is_real(mm.typecode)) {
        if (d_csr_test(sm, num_iterations, mm, 0)) {
            printf("An error occured while test the funtion for double matrix");
            return 1;
        }
    } else {
        if (i_csr_test(sm, num_iterations, mm, 0)) {
            printf("An error occured while test the funtion for int matrix");
            return 1;
        }
    }
    */
    return 0;
}
