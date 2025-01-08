#include "mmio.h"
#include "csr.h"
#include "hll.h"
#include "utils.h"
#include "spmv_seq.h"
#include "spmv_openmp.h"
#include <errno.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void write_csr_mtx(struct csr *csr, struct MatrixMarket *mm) {
    for (int i = 0; i < csr->num_rows; ++i) {
        int num_cols = csr->row_pointer[i+1] - csr->row_pointer[i];
        printf("row no. = %d has %d non-zero elements\n", i+1, num_cols);
    }
    for (int i = 0; i < csr->num_rows; ++i) {
        for (int j = csr->row_pointer[i]; j < csr->row_pointer[i+1]; ++j) {
            int col = csr->col_index[j]+1;
            if (mm_is_integer(mm->typecode)) {
                int v = ((int*)csr->data)[j];
                printf("%d %d %d\n", i+1, col, v);
            } else {
                double v = ((double*)csr->data)[j];
                printf("%d %d %f\n", i+1, col, v);
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

int main(int argc, char *argv[])
{
    if (--argc != 3) {
        printf("see usage: program num_threads num_runs matrix\n");
        return -1;
    }
    long num_threads = strtol(argv[argc--], NULL, 10);
    if (num_threads == 0) {
        printf("an error (%s) occurred while parsing num_threads: %s\n", argv[argc+1]);
        return 1;
    }
    long num_iterations = strtol(argv[argc--], NULL, 10);
    if (num_iterations == 0) {
        printf("The error (%s) occured while convert string to long", strerror);
        return -1;
    }
    struct MatrixMarket mm;
    if (read_mtx(argv[argc--], &mm)) {
        printf("Read error!");
        return 1;
    }
    printf("matrix has %d rows and %d cols\n", mm.num_rows, mm.num_cols);
    
    struct csr sm;
    if (csr_init(&sm, &mm)) { 
        printf("cannot read matrix into CSR format\n");
        return 1;
    }

    omp_set_num_threads(num_threads);
    srand(42);
    double *r = malloc(sm.num_cols * sizeof(double));
    double *v = d_random(sm.num_cols);
    double sum_of_times = 0;
    for (int i = 0; i < num_iterations; i++) {
        clock_t start = clock();
        if (d_spmv_csr_par(r, &sm, v, sm.num_cols)) {
            return 1;
        }
        clock_t end = clock();
        sum_of_times += (end - start) / (double)CLOCKS_PER_SEC;  
    }
    double mean_time = sum_of_times / num_iterations;
    double flops = (2 * mm.nz) / (double)mean_time;
    printf("mean time: %f s\n", mean_time);
    printf("MEGA FLOPS=%f\n", flops * 1e-6);
    return 0;
}
