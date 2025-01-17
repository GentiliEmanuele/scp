#ifndef UTILS_H
#define UTILS_H
#include "mmio.h"

#ifdef SCP_VERBOSE
#define PATH_BUFFER_SZ 48
#endif
struct MatrixMarket {
    MM_typecode typecode;
    double *data;
    int    *rows; 
    int    *cols;
    int    num_rows;
    int    num_cols;
    int    nz;
#ifdef SCP_VERBOSE
    char  *path;
    char  __path_buffer[PATH_BUFFER_SZ];
#endif
};

int read_mtx(const char *path, struct MatrixMarket *mm);
void read_and_measure_csr(char *path, int num_runs, int num_thread, char *out_path);
void read_and_measure_hll(char *path, int hack_size, int num_runs, int num_thread, char *out_path);
void mtx_cleanup(struct MatrixMarket *mm);
#endif