#include <stdio.h>
#include "seq_time.h"
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 5) {
        printf("usage:\nprogram matrices_list num_runs hack_size result_file_path\n");
        return 1;
    }
    char *out_path = argv[4];
    int num_iterations = atoi(argv[2]);
    int hack_size = atoi(argv[3]);
    char *path = argv[1];
    read_and_measure_seq_hll(path, hack_size, num_iterations, out_path);
    return 0;
}

