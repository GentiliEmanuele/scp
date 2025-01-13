#include <stdio.h>
#include <stdlib.h>
#include "omp_time.h"

int main(int argc, char **argv) {
    if (argc != 5) {
        printf("Please pass the path to the matrix, the hack_size, the number of the thread and the number of the iterations as arguments \n");
        return 1;
    }
    int num_iterations = atoi(argv[--argc]);
    int num_thread = atoi(argv[--argc]);
    int hack_Size = atoi(argv[--argc]);
    char *path = argv[--argc];
    time_measurement_t time_measurement;
    omp_time_hll(path, hack_Size, num_iterations, num_thread, &time_measurement);
    return 0;
}