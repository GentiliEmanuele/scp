#ifndef SCPA_OMP_TEST_H
#define SCPA_OMP_TEST_H

void omp_time_csr(const char *file, int num_runs, int num_threads);
void omp_time_hll(const char *file, int hack_size, int num_runs, int num_threads);

#endif