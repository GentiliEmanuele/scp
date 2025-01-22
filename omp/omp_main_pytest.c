#include "csr.h"
#include "hll.h"
#include "utils.h"
#include "vec.h"
#include "spmv_openmp.h"
#include "spmv_seq.h"
#include <stdio.h>
#include <string.h>

enum {
    SERIAL = 0,
    CSR,
    HLL
};

static int test_hll(const char *file, const char *vfile, const char *rfile, int hack_size) {
    struct MatrixMarket mm;
    if (read_mtx(file, &mm)) {
        return 0;
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
    if (spmv_hll_par(p, &sm, v, n)) {
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

static int test_csr(const char *file, const char *vfile, const char *rfile) {
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
    if (spmv_csr_par(p, &sm, v, n)) {
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

static int test_csr_seq(const char *file, const char *vfile, const char *rfile) {
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
    if (spmv_csr_seq(p, &sm, v, n)) {
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
    --argc;
    if (argc != 3 && argc != 4) {
        printf("see usage: program matrix_path matrices_file (csr|hll|serial) {hack_size}\n");
        return 1;
    }
    int ttype;
    if (!strcmp(argv[3], "csr")) {
        ttype = CSR;
    } else if (!strcmp(argv[3], "hll")) {
        ttype = HLL;
    } else if (!strcmp(argv[3], "serial")) {
        ttype = SERIAL;
    } else {
        printf("expected one of csr, hll or serial but got %s\n", argv[3]);
        return 1;
    }
    int hack_size = 32;
    if (ttype == HLL) {
        if (argc == 3) {
            printf("warning: hack_size not provided, defaulting to 32\n");
        } else {
            hack_size = strtol(argv[4], NULL, 10);
            if (hack_size <= 0) {
                printf("invalid hack_size provided %s\n", argv[4]);
                return 1;
            }
        }
    }
    FILE* fp = fopen(argv[2], "r");
    if (!fp) {
        printf("cannot open matrices file: %s\n", argv[2]);
        return 1;
    }
    char line[1024];
    char rfile[1024];
    char mfile[1024];
    char vfile[1024];
    int total = 0;
    int failed = 0;
    while (fgets(line, 1024, fp) != NULL) {
        int n = strlen(line);
        line[--n] = 0;
        sprintf(mfile, "%s/data/%s", argv[1], line);
        sprintf(rfile, "%s/test/output/%s.result", argv[1], line);
        sprintf(vfile, "%s/test/output/%s.vector", argv[1], line);
        ++total;
        switch (ttype)
        {
        case SERIAL:
            if (test_csr_seq(mfile, vfile, rfile)) {
                printf("test failed for %s\n", line);
                ++failed;
            }
            break;
        case CSR:
            if (test_csr(mfile, vfile, rfile)) {
                printf("test failed for %s\n", line);
                ++failed;
            }
            break;
        case HLL:
            if (test_hll(mfile, vfile, rfile, 32)) {
                printf("test failed for %s\n", line);
                ++failed;
            }
            break;
        default:
            printf("cannot determine which kind of test it's supposed to run\n");
            return 1;
        }
        
    }
    printf("test passed: %d/%d\n", total-failed, total);
    return 0;
}