#include "utils.h"
#include "mmio.h"
#include "omp_time.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#define IS_ZERO(x) (fabs(x) == 0.0)

struct vec3d {
    int row;
    int col;
    int pos;
};

static inline int lt(struct vec3d *a, struct vec3d *b) {
    if ((a->row < b->row) || (a->row == b-> row) && (a->col < b->col)){
        return -1;
    } else if ((a->row == b->row) && (a -> col == b -> col)) {
        return 0;
    } else {
        return 1;
    }
}

static inline int qcmp(const void *a, const void *b) {
    return lt((struct vec3d*)a, (struct vec3d*)b);
}

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

static int readline(FILE *f, const char *fmt, int *row, int *col, void *val, int vp, struct MatrixMarket *mm) {
    int not_zero = 1;
    if (!mm_is_pattern(mm->typecode)) {
        if (mm_is_real(mm->typecode)) {
            fscanf(f, fmt, col, row, &((double*)val)[vp]);
            not_zero = !IS_ZERO(((double*)val)[vp]);
        } else {
            fscanf(f, fmt, col, row, &((int*)val)[vp]);
            not_zero = ((int*)val)[vp] != 0;
        }
    } else {
        fscanf(f, fmt, col, row);
        ((double*)val)[vp] = 1.0;
    }
    // Adjust from one-based to zero-based
    *col -= 1;
    *row -= 1;
    return not_zero;
}

static int parse_rows_sy(FILE *f, struct MatrixMarket *mm) {
    struct vec3d *items = malloc(2 * mm->nz * sizeof(struct vec3d));
    void *data = malloc(mm->nz * get_element_size(mm));
    char *fmt = get_format_string(mm);
    int real_nz = 0;
    int count_zero = 0;
    for (int k = 0; k < mm->nz; ++k) {
        if (readline(f, fmt, &items[real_nz].row, &items[real_nz].col, data, k, mm)) {
            items[real_nz].pos = k;
            int i = items[real_nz].row;
            int j = items[real_nz].col;
            ++real_nz;
            if (i != j) {
                items[real_nz].row = j;
                items[real_nz].col = i;
                items[real_nz].pos = k;
                ++real_nz;
            }
        } else count_zero ++;
    }
    qsort(items, real_nz, sizeof(struct vec3d), qcmp);
    mm->rows = malloc(real_nz * sizeof(int));
    mm->cols = malloc(real_nz * sizeof(int));
    int element_size = get_element_size(mm);
    mm->data = malloc(real_nz * element_size);
    for (int i = 0; i < real_nz; ++i) {
        mm->rows[i] = items[i].row;
        mm->cols[i] = items[i].col;
        if (mm_is_real(mm->typecode) || mm_is_pattern(mm->typecode)) {
            memcpy(&((double*)mm->data)[i], &((double*)data)[items[i].pos], element_size);
        } else {
            memcpy(&((int*)mm->data)[i], &((int*)data)[items[i].pos], element_size);            
        }
    }
    free(data);
    printf("symmetric matrix: number of non-zeros goes from %d to %d (explicit zeros = %d)\n", mm->nz, real_nz, count_zero);
    mm->nz = real_nz;
    return 0;
}


static int parse_rows_ns(FILE *f, struct MatrixMarket *mm) {
    int *rows = malloc(mm->nz * sizeof(int));
    int *cols = malloc(mm->nz * sizeof(int));
    void *values = malloc(mm->nz * get_element_size(mm));
    char *fmt = get_format_string(mm);
    int real_nz = 0;
    int explicit_zeros = 0;
    for (int i = 0; i < mm->nz; i++) {
        if (readline(f, fmt, &rows[i], &cols[i], values, i, mm)) {
            ++real_nz;
        } else {
            ++explicit_zeros;
        }
    }
    mm->rows = rows;
    mm->cols = cols;
    mm->data = values;
    mm->nz = real_nz;
    printf("non symmetric matrix: number of non-zeros %d (explicit zeros = %d)\n", mm->nz, explicit_zeros);
    return 0;
}

static int parse_rows(FILE *f, struct MatrixMarket *mm) { 
    if (mm_is_symmetric(mm->typecode)) {
        return parse_rows_sy(f, mm);
    } else {
        return parse_rows_ns(f, mm);
    }
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
    printf("matrix %s:\n", path);
    int ir = parse_rows(f, mm);
    fclose(f);
    printf("matrix has %d rows and %d cols and number of non-zeros %d\n", mm->num_rows, mm->num_cols, mm->nz);
    return ir;
}

/**
 * read_and_test for csr format-- read a file in txt with the name of all matrices to use
 * @param path (in)                 path of the file to read
 * @param num_runs (in)             number of times the measurement must be performed
 * @param num_thread (in)           number of thread used in the executions
 * @param out_path (in)             path to write results
 */
void read_and_measure_csr(char *path, int num_runs, int num_thread, char *out_path) {
    if (path == NULL) {
        printf("Please pass the path to the file with matrices name\n");
        return;
    }
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        printf("cannot open file %s\n", path);
        return;
    }
    size_t len = 0;
    size_t read;
    char *line = NULL;
    time_measurement_t time_measurement;
    FILE *results = fopen(out_path, "w");
    if (results == NULL) {
        printf("cannot open file %s\n", path);
        return;   
    }
    while ((read = getline(&line, &len, f)) != -1) {
        printf("Get time for the matrix %s", line);
        size_t last_idx = strlen(line) - 1;
        if( line[last_idx] == '\n' ) {
            line[last_idx] = '\0';
        }
        omp_time_csr(line, num_runs, num_thread, &time_measurement);
        fprintf(results, "Matrix name: %s \t MFLOPS %f \t mean_time %f\n", line, time_measurement.flops, time_measurement.mean_time);
    }
    printf("\n");
    fclose(f);
    if (line)
        free(line);
}

/**
 * read_and_test for hll format-- read a file in txt with the name of all matrices to use
 * @param path (in)                 path of the file to read
 * @param hack_size (in)            number of rows of each hack
 * @param num_runs (in)             number of times the measurement must be performed
 * @param num_thread (in)           number of thread used in the executions
 * @param out_path (in)             path to write results
 */
void read_and_measure_hll(char *path, int hack_size, int num_runs, int num_thread, char *out_path) {
    if (path == NULL) {
        printf("Please pass the path to the file with matrices name\n");
        return;
    }
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        printf("cannot open file %s\n", path);
        return;
    }
    size_t len = 0;
    size_t read;
    char *line = NULL;
    time_measurement_t time_measurement;
    FILE *results = fopen(out_path, "w");
    if (results == NULL) {
        printf("cannot open file %s\n", path);
        return;   
    }
    while ((read = getline(&line, &len, f)) != -1) {
        printf("Get time for the matrix %s", line);
        size_t last_idx = strlen(line) - 1;
        if( line[last_idx] == '\n' ) {
            line[last_idx] = '\0';
        }
        omp_time_hll(line, hack_size, num_runs, num_thread, &time_measurement);
        fprintf(results, "Matrix name: %s \t MFLOPS %f \t mean_time %f\n", line, time_measurement.flops, time_measurement.mean_time);
    }
    printf("\n");
    fclose(f);
    if (line)
        free(line);
}
