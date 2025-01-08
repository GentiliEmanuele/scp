#include "utils.h"
#include "mmio.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>

static inline char *get_format_string(struct MatrixMarket *m) {
    if (mm_is_pattern(m->typecode)) {
        return "%d %d\n";
    } else {
        if (mm_is_integer(m->typecode)) {
            return "%d %d %ld\n";
        } else {
            return "%d %d %lg\n";
        }
    }
}

static inline void iswap(int *v1, int i, int j) {
    int t = v1[i];
    v1[i] = v1[j];
    v1[j] = t;
}

static inline void dswap(double *v1, int i, int j) {
    double t = v1[i];
    v1[i] = v1[j];
    v1[j] = t;
}

static inline void vswap(void *v, int i, int j, struct MatrixMarket *mm) {
    if (mm_is_pattern(mm->typecode) || mm_is_real(mm->typecode)) {
        dswap((double*)v, i, j);
    } else {
        iswap((int*)v, i, j);
    }
}

static inline int __lt(int *rows, int *cols, int i, int j) {
    return rows[i] < rows[j] || rows[i] == rows[j] && cols[i] < cols[j];
}

static int partition(int *rows, int *cols, void *values, int lo, int hi, struct MatrixMarket *mm) {
    int pivot = hi;
    int i = lo;
    for (int j = lo; j < hi; ++j) {
        if (__lt(rows, cols, j, pivot)) {
            iswap(rows, i, j);
            iswap(cols, i, j);
            vswap(values, i, j, mm);
            i += 1;
        }
    }
    iswap(rows, i, hi);
    iswap(cols, i, hi);
    vswap(values, i, hi, mm);
    return i;
}

static void aux_sort(int *rows, int *cols, void *values, int lo, int hi, struct MatrixMarket *mm) {
    if (lo >= hi || lo < 0) {
        return;
    }
    int p = partition(rows, cols, values, lo, hi, mm);
    aux_sort(rows, cols, values, lo, p - 1, mm);
    aux_sort(rows, cols, values, p + 1, hi, mm);
}

void sort(int *rows, int *cols, void *values, int n, struct MatrixMarket *mm) {
    aux_sort(rows, cols, values, 0, n - 1, mm);
}

static void unpack(struct MatrixMarket *mm, int real_nz) {
    void *values = malloc(real_nz * get_element_size(mm));
    int *rows = malloc(real_nz * sizeof(int));
    int *cols = malloc(real_nz * sizeof(int));
    int p = 0;
    int *packed_rows = mm->rows;
    int *packed_cols = mm->cols;
    void *packed_values = mm->data;
    for (int k = 0; k < mm->nz; ++k) {
        int i = packed_rows[k];
        int j = packed_cols[k];
        rows[p] = i;
        cols[p] = j;
        if (mm_is_integer(mm->typecode)) {
            ((int*)values)[p] = ((int*)packed_values)[k];
        } else if (mm_is_real(mm->typecode)) {
            ((double*)values)[p] = ((double*)packed_values)[k];
        }
        ++p;
        if (i != j) {
            rows[p] = j;
            cols[p] = i;
            if (mm_is_integer(mm->typecode)) {
                ((int*)values)[p] = ((int*)packed_values)[k];
            } else if (mm_is_real(mm->typecode)) {
                ((double*)values)[p] = ((double*)packed_values)[k];
            }
            ++p;
        }
    }
    sort(rows, cols, values, real_nz, mm);
    free(mm->rows);
    free(mm->cols);
    free(mm->data);
    mm->rows = rows;
    mm->cols = cols;
    mm->data = values;
    mm->nz = real_nz;
}

static int parse_rows(FILE *f, struct MatrixMarket *mm) {
    int *rows = malloc(mm->nz * sizeof(int));
    int *cols = malloc(mm->nz * sizeof(int));
    void *values = malloc(mm->nz * get_element_size(mm));
    char *fmt = get_format_string(mm);
    int real_nz = 0;
    for (int i = 0; i < mm->nz; i++) {
        if (!mm_is_pattern(mm->typecode)) {
            if (mm_is_real(mm->typecode)) {
                fscanf(f, fmt, &cols[i], &rows[i], &((double*)values)[i]);
            } else {
                fscanf(f, fmt, &cols[i], &rows[i], &((int*)values)[i]);
            }
        } else {
            fscanf(f, fmt, &cols[i], &rows[i]);
            ((double*)values)[i] = 1.0;
        }
        rows[i]--;  /* adjust from 1-based to 0-based */
        cols[i]--;
        real_nz++;
        if (mm_is_symmetric(mm->typecode) && rows[i] != cols[i]) {
            real_nz++;
        }
    }
    mm->rows = rows;
    mm->cols = cols;
    mm->data = values;
    if (mm->nz != real_nz && mm_is_symmetric(mm->typecode)) {
        printf("parsing a symmetric matrix, real number of non-zero values goes from %d to %d\n", mm->nz, real_nz);
        unpack(mm, real_nz);
    } else if (!mm_is_symmetric(mm->typecode) && mm->nz != real_nz) {
        printf("expected %d non-zeros values but got %d\n", mm->nz, real_nz);
        return 1;
    }
    return 0;
}

/**
 * read_mtx -- read a file in MatrixMarket format
 * @param path -- path of the file to read
 * @param mm   -- output struct
 * @return 1 if the file cannot be opened or the matrix is in a format not supported,
 * 0 otherwise.
 */
int read_mtx(const char *path, struct MatrixMarket *mm) {
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        printf("cannot open file %s\n", path);
        return 1;
    }
    
    if (mm_read_banner(f, &(mm->typecode)) != 0) {
        printf("could not process Matrix Market banner\n");
        return 1;
    }

    if (mm_is_complex(mm->typecode) && mm_is_matrix(mm->typecode) && mm_is_sparse(mm->typecode)) {
        printf("sorry, this application does not support ");
        printf("matrix market type: [%s]\n", mm_typecode_to_str(mm->typecode));
        return 1;
    }

    int M, N, nz;
    int ret_code =  mm_read_mtx_crd_size(f, &M, &N, &nz);
    if (ret_code) {
        return 1;
    }

    mm->num_rows = M;
    mm->num_cols = N;
    mm->nz = nz;
    int ir = parse_rows(f, mm);
    fclose(f);
    return ir;
}

double *d_random(int n) {
    double *v = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        v[i] = (double)rand() / RAND_MAX;
    }
    return v;
    
}

int d_veceq(double *u, double *v, int n, double eps) {
    for (int i = 0; i < n; i++) {
        if (fabs(u[i] - v[i]) > eps)
            return 1;
    }
    return 0;
}