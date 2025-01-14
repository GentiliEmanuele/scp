#ifndef SCPA_OMP_TEST_H
#define SCPA_OMP_TEST_H

typedef struct time_measurement {
    double mean_time;
    double flops;
} time_measurement_t;

int omp_time_csr(const char *file, int num_runs, int num_threads, time_measurement_t *time_measurement);
int omp_time_hll(const char *file, int hack_size, int num_runs, int num_threads,  time_measurement_t *time_measurement);

#endif