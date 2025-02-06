#ifndef SCPA_OMP_TEST_H
#define SCPA_OMP_TEST_H

typedef struct time_measurement {
    double mean_time;
    double flops;
    double std_dev;
    int num_threads;
    int num_runs;
} time_measurement_t;

int omp_time_csr(const char *file, int num_runs, int num_threads, time_measurement_t *time_measurement);
int omp_time_hll(const char *file, int hack_size, int num_runs, int num_threads,  time_measurement_t *time_measurement);
void read_and_measure_csr(char *path, int num_runs, int num_thread, char *out_path);
void read_and_measure_hll(char *path, int hack_size, int num_runs, int num_thread, char *out_path);
#endif