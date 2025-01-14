#include <stdio.h>
#include <stdlib.h>
#include "omp_time.h"
#include "utils.h"

// this main take as parameters (in order) path to the file with the matrices path, num_thread and num_iterations
int main(int argc, char **argv) {
    if (argc != 4) {
        printf("Please pass as arguments: path to the matrix, num_thread and num_iterations\n");
        return 1;
    }
    int num_iterations = atoi(argv[--argc]);
    int num_threads = atoi(argv[--argc]);
    char *path = argv[--argc];
    time_measurement_t time_measurement;
    read_and_measure_csr(path, num_iterations, num_threads);
    return 0;
}