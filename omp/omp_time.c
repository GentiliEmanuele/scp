#include "omp_time.h"
#include "csr.h"
#include "spmv_openmp.h"
#include "utils.h"
#include "vec.h"
#include <omp.h>
#include <stdlib.h>
#include <string.h>
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
    double *samples = malloc(num_runs * sizeof(double));
    if (!samples) {
        printf("out of memory\n");
        csr_cleanup(&sm);
        mtx_cleanup(&mm);
        return 1;
    }
    double sum = 0;
    for (int i = 0; i < num_runs; i++) {
        clock_t start = clock();
        if (spmv_csr_par(r, &sm, v, n)) {
            printf("warning: couldn't complete sparse matrix-vector product of run %d\n", i);
        }
        clock_t end = clock();
        sum += (end - start) / (double)CLOCKS_PER_SEC;  
    }
    time_measurement->mean_time = sum / num_runs;
    time_measurement->flops = (2 * mm.nz) / (double)time_measurement->mean_time;
    time_measurement->std_dev = std_devl(samples, time_measurement->mean_time, num_runs);
    time_measurement->num_threads = num_threads;
    time_measurement->num_runs = num_runs;
    csr_cleanup(&sm);
    mtx_cleanup(&mm);
    free(samples);
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
    double *samples = malloc(num_runs * sizeof(double));
    if (!samples) {
        printf("out of memory\n");
        mtx_cleanup(&mm);
        hll_cleanup(&sm);
        return 1;
    }
    double sum = 0;
    for (int i = 0; i < num_runs; i++) {
        clock_t start = clock();
        if (spmv_hll_par(r, &sm, v, n)) {
            printf("warning: couldn't complete sparse matrix-vector product of run %d\n", i);
        }
        clock_t end = clock();
        double t = (end - start) / (double)CLOCKS_PER_SEC;
        samples[i] = t;
        sum += t;
    }
    time_measurement->mean_time = sum / num_runs;
    time_measurement->flops = (2 * mm.nz) / (double)time_measurement->mean_time;
    time_measurement->std_dev = std_devl(samples, time_measurement->mean_time, num_runs);
    time_measurement->num_threads = num_threads;
    time_measurement->num_runs = num_runs;
    mtx_cleanup(&mm);
    hll_cleanup(&sm);
    free(samples);
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
    fprintf(results, "File,FLOPS,mean_time [s],std_dev,num_threads,num_runs\n");
    while ((read = getline(&line, &len, f)) != -1) {
        printf("Get time for the matrix %s", line);
        size_t last_idx = strlen(line) - 1;
        if( line[last_idx] == '\n' ) {
            line[last_idx] = '\0';
        }
        omp_time_csr(line, num_runs, num_thread, &time_measurement);
        fprintf(results, "%s,%f,%f,%f,%d,%d\n",
            line,
            time_measurement.flops,
            time_measurement.mean_time,
            time_measurement.std_dev,
            time_measurement.num_threads,
            time_measurement.num_runs);
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
    fprintf(results, "File,FLOPS,mean_time [s],std_dev,num_threads,num_runs\n");
    while ((read = getline(&line, &len, f)) != -1) {
        printf("Get time for the matrix %s", line);
        size_t last_idx = strlen(line) - 1;
        if( line[last_idx] == '\n' ) {
            line[last_idx] = '\0';
        }
        omp_time_hll(line, hack_size, num_runs, num_thread, &time_measurement);
        fprintf(results, "%s,%f,%f,%f,%d,%d\n",
            line,
            time_measurement.flops,
            time_measurement.mean_time,
            time_measurement.std_dev,
            time_measurement.num_threads,
            time_measurement.num_runs);
    }
    printf("\n");
    fclose(f);
    if (line)
        free(line);
}