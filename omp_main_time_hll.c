#include <stdio.h>
#include <stdlib.h>
#include "omp_time.h"
#include "utils.h"

int main(int argc, char **argv) {
    if (argc != 6) {
        printf("usage:\nprogram matrices_list num_runs num_threads hack_size result_file_path\n");
        return 1;
    }
    char *out_path = argv[--argc];
    int num_iterations = atoi(argv[--argc]);
    int num_thread = atoi(argv[--argc]);
    int hack_size = atoi(argv[--argc]);
    char *path = argv[--argc];
    read_and_measure_hll(path, hack_size, num_iterations, num_thread, out_path);
    return 0;
}