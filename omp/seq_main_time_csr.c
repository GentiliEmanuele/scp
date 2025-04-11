#include <stdio.h>
#include <stdlib.h>
#include "seq_time.h"

int main(int argc, char **argv) {
    if (argc != 4) {
        printf("usage:\nprogram matrices_list num_runs result_file_path\n");
        return 1;
    }
    char *out_path = argv[--argc];
    int num_iterations = atoi(argv[--argc]);
    char *path = argv[--argc];
    read_and_measure_seq_csr(path, num_iterations, out_path);
    return 0;
}