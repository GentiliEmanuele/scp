#include <stdio.h>
#include <stdlib.h>
#include "omp_time.h"
#include "utils.h"

int main(int argc, char **argv) {
    if (argc != 5) {
        printf("Please pass the path to the txt with the matrices path, the hack_size, the number of the thread and the number of the iterations as arguments \n");
        return 1;
    }
    int num_iterations = atoi(argv[--argc]);
    int num_thread = atoi(argv[--argc]);
    int hack_size = atoi(argv[--argc]);
    char *path = argv[--argc];
    read_and_measure_hll(path, hack_size, num_iterations, num_thread);
    return 0;
}