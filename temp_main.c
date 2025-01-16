#include "csr.h"
#include "hll.h"
#include "utils.h"
#include "vec.h"
#include "spmv_openmp.h"
#include <stdio.h>
#include <string.h>

int read_vector(double *vector, int n, const char *path) {
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        printf("cannot open file: %s\n", path);
        return 1;
    }
    for (int i = 0; i < n; i++) {
        fscanf(f, "%f", &vector[i]);
    }
    return 0;
}

int test_hll(const char *file, const char *vfile, const char *rfile, int hack_size) {
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
    double *v = d_zeros(m);
    if (read_vector(v, m, vfile)) {
        mtx_cleanup(&mm);
        return 0;
    }
    double *p = d_zeros(m);
    if (d_spmv_hll_par(p, &sm, v, n)) {
        mtx_cleanup(&mm);
        hll_cleanup(&sm);
        return 1;
    }
    double *r = d_zeros(m);
    read_vector(r, m, rfile);
    int err = d_veceq(r, p, m, 1e-6) ? 0 : 1;
    mtx_cleanup(&mm);
    hll_cleanup(&sm);
    return 0;
}

int test_csr(const char *file, const char *vfile, const char *rfile) {
    struct MatrixMarket mm;
    if (read_mtx(file, &mm)) {
        return 0;
    }

    struct csr sm;
    if (csr_init(&sm, &mm)) {
        mtx_cleanup(&mm);
        return 1;
    }

    int m = sm.num_rows;
    int n = sm.num_cols;
    double *v = d_zeros(m);
    if (read_vector(v, m, vfile)) {
        mtx_cleanup(&mm);
        return 0;
    }
    double *p = d_zeros(m);
    if (d_spmv_csr_par(p, &sm, v, n)) {
        mtx_cleanup(&mm);
        csr_cleanup(&sm);
        return 1;
    }
    double *r = d_zeros(m);
    read_vector(r, m, rfile);
    int err = d_veceq(r, p, m, 1e-6) ? 0 : 1;
    mtx_cleanup(&mm);
    csr_cleanup(&sm);
    return 0;
}

int main(int argc, char **argv) {
    if (--argc != 1) {
        printf("see usage: program matrices_file\n");
        return 1;
    }
    FILE* fp = fopen(argv[1], "r");
    if (!fp) {
        printf("cannot open matrices file: %s\n", argv[1]);
        return 1;
    }
    char line[1024];
    char rfile[1024];
    char mfile[1024];
    char vfile[1024];
    while (fgets(line, 1024, fp) != NULL) {
        int n = strlen(line);
        line[--n] = 0;
        sprintf(mfile, "data/%s", line);
        sprintf(rfile, "test/output/%s.result", line);
        sprintf(vfile, "test/output/%s.vector", line);
        if (test_csr(mfile, vfile, rfile)) {
            printf("test failed for %s\n", line);
            return 1;
        }
    }
    return 0;
}