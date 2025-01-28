#include "cuda_time.h"
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
    if (argc != 4 && argc != 5) {
        printf("see usage: program output matrices_list runs_num [hll|csr] {hack_size}\n");
        return -1;
    }

    int runs_num = atoi(argv[3]);
    if (runs_num == 0) {
        printf("An error occurred while converting hack_size\n");
        return -1;
    }
    int test_type = parse_test_type(argv[4]);
    if (test_type == NOTEST) {
        printf("expected one of csr or hll but got %s\n", argv[4]);
        return -1;
    }

    int hack_size = 0;
    if (test_type == HLL) {
        if (argc == 5) {
            hack_size = atoi(argv[5]);
            if (hack_size == 0) {
                printf("An error occurred while converting hack_size\n");
                return -1;
            }
        } else {
            hack_size = 32;
            printf("no hack_size specified, default value (32) set\n");
        }
    }
    
    char filename[1024];
    if (hack_size) {
        sprintf(filename, "%s_%d_%d.csv", argv[1], runs_num, hack_size);
    } else {
        sprintf(filename, "%s_%d.csv", argv[1], runs_num);
    }
    FILE* off = fopen(filename, "w");
    if (off == NULL) {
        printf("cannot open file %s\n", filename);
        return -1;
    }
    FILE* iff = fopen(argv[2], "r");
    if (iff == NULL) {
        printf("cannot open file %s\n", argv[2]);
        fclose(off);
        return -1;
    }
    fprintf(off, "matrix,time,flops\n");
    INIT_TIME_INFO(ti);
    char line[1024];
    while (fgets(line, 1024, iff) != NULL) {
        int n = strlen(line);
        if (line[--n] == '\n') {
            line[n] = 0;
        }
        printf("matrix %s\n", line);
        ti.flops = 0.0;
        ti.millis = 0.0;
        if (test_type == CSR) {
            csr_time(line, runs_num, &ti);
        } else if (test_type == HLL) {
            hll_time(line, runs_num, hack_size, &ti);
        }
        fprintf(off, "\"%s\",%f,%f\n", line, ti.millis, ti.flops);
    }
    fclose(iff);
    fclose(off);
}
