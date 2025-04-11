#include "seq_time.h"
#include "utils.h"
#include "csr.h"
#include "vec.h"
#include "spmv_seq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int seq_time_csr(const char *file, int num_runs, time_measurement_seq_t *time_measurement) {
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
    double *samples = malloc(num_runs * sizeof(double));
    if (!samples) {
        printf("out of memory\n");
        csr_cleanup(&sm);
        mtx_cleanup(&mm);
        return 1;
    }
    double max = -1;
    double min = 1e10;
    double sum = 0;
    double start, end, execution_time;
    for (int i = 0; i < num_runs; i++) {
        start = omp_get_wtime();
        if (spmv_csr_seq(r, &sm, v, n)) {
            printf("warning: couldn't complete sparse matrix-vector product of run %d\n", i);
        }
        end = omp_get_wtime();
        execution_time = end - start;
        if (execution_time < min) {
            min = execution_time;
        }
        if (execution_time > max) {
            max = execution_time;
        }
        samples[i] = execution_time;
        sum += execution_time;  
    }
    time_measurement->mean_time = sum / num_runs;
    time_measurement->flops = (2 * mm.nz) / time_measurement->mean_time;
    time_measurement->std_dev = std_devl(samples, time_measurement->mean_time, num_runs);
    time_measurement->num_runs = num_runs;
    time_measurement->min = min;
    time_measurement->max = max;
    csr_cleanup(&sm);
    mtx_cleanup(&mm);
    free(samples);
}

void read_and_measure_seq_csr(char *path, int num_runs, char *out_path) {
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
    time_measurement_seq_t time_measurement;
    FILE *results = fopen(out_path, "w");
    if (results == NULL) {
        printf("cannot open file %s\n", path);
        return;   
    }
    fprintf(results, "File,FLOPS,mean_time [s],std_dev,min,max,num_threads,num_runs\n");
    while ((read = getline(&line, &len, f)) != -1) {
        printf("Get time for the matrix %s", line);
        size_t last_idx = strlen(line) - 1;
        if( line[last_idx] == '\n' ) {
            line[last_idx] = '\0';
        }
        seq_time_csr(line, num_runs, &time_measurement);
        fprintf(results, "%s,%f,%f,%f,%f,%f,%d\n",
            line,
            time_measurement.flops,
            time_measurement.mean_time,
            time_measurement.std_dev,
            time_measurement.min,
            time_measurement.max,
            time_measurement.num_runs);
    }
    printf("\n");
    fclose(f);
    if (line)
        free(line);
}