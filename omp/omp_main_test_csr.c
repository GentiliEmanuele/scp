#include <stdio.h>
#include "omp_test.h"

// this main function requires only path to the matrix
int main(int argc, char **argv) {
    if (argc != 2) {
        printf("Please pass the path to the matrix as argument \n");
        return 1;
    }
    test_csr(argv[--argc]);
    return 0;
}