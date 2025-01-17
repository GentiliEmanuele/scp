#include <stdio.h>
#include <stdlib.h>
#include "omp_time.h"
#include "utils.h"

int main(int argc, char **argv) {
    if (argc != 5) {
        printf("usage:\nprogram matrices_list num_threads num_runs result_file_path\n");
        return 1;
    }
    char *out_path = argv[--argc];
    int num_iterations = atoi(argv[--argc]);
    int num_threads = atoi(argv[--argc]);
    char *path = argv[--argc];
    time_measurement_t time_measurement;
    read_and_measure_csr(path, num_iterations, num_threads, out_path);
    return 0;
}