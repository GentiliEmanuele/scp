#include <stdio.h>
#include <stdlib.h>
#include "omp_time.h"

// this main take as parameters (in order) path, num_thread and num_iterations
int main(int argc, char **argv) {
    if (argc != 4) {
        printf("Please pass as arguments: path to the matrix, num_thread and num_iterations\n");
        return 1;
    }
    int num_iterations = atoi(argv[--argc]);
    int num_threads = atoi(argv[--argc]);
    char *path = argv[--argc];
    time_measurement_t time_measurement;
    omp_time_csr(path, num_iterations, num_threads, &time_measurement);
    return 0;
}