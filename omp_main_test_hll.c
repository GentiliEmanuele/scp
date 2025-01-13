#include <stdio.h>
#include <stdlib.h>
#include "omp_test.h"

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Please pass the path to the matrix and hack_size as arguments \n");
        return 1;
    }
    int hack_size = atoi(argv[--argc]);
    char *path = argv[--argc];
    test_hll(path, hack_size);
    return 0;
}