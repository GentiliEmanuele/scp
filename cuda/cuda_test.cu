#include "cuda_test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NOTEST 0
#define CSR    1
#define HLL    2

int parse_test_type(char *s) {
    if (!strcmp(s, "csr")) {
        return CSR;
    } else if (!strcmp(s, "hll")) {
        return HLL;
    } else {
        return NOTEST;
    }
}

int main(int argc, char *argv[]) {
    --argc;
    if (argc != 2 && argc != 3) {
        printf("see usage: program matrix [hll|csr] {hack_size}\n");
        return -1;
    }

    int test_type = parse_test_type(argv[2]);
    if (test_type == NOTEST) {
        printf("expected one of csr or hll but got %s\n", argv[2]);
        return -1;
    }

    int hack_size = 0;
    if (test_type == HLL) {
        if (argc == 3) {
            hack_size = atoi(argv[3]);
            if (hack_size == 0) {
                printf("An error occurred while converting hack_size\n");
                return -1;
            }
        } else {
            hack_size = 32;
            printf("no hack_size specified, default value (32) set\n");
        }
    }
    
    printf("matrix %s\n", argv[1]);
    if (test_type == CSR) {
        csr_test(argv[1]);
    } else if (test_type == HLL) {
        hll_test(argv[1], hack_size);
    }
}
