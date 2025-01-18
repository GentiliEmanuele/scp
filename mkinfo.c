#include "csr.h"
#include "hll.h"
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_LINE 255

static char* mm_name(char *path) {
    int n = strlen(path);
    char *pp = &path[n];
    while (*pp != '.') {
        --pp;
    }
    *pp = 0;
    while (pp != path && *pp != '/') {
        --pp;
    }
    return *pp == '/' ? pp + 1 : pp;
}

int main(int argc, char *argv[])
{
    if (--argc != 2) {
        printf("see usage: program matrices_file output_file\n");
        return -1;
    }
    struct MatrixMarket mm;
    char line[MAX_LINE + 1];
    FILE* off = fopen(argv[2], "w");
    FILE* iff = fopen(argv[1], "r");
    if (!iff) {
        printf("cannot open file %s\n", argv[1]);
    }
    if (!off) {
        printf("cannot create file %s\n", argv[2]);
        return -1;
    }
    fprintf(off, "Matrix name,M,N,NZ,csr_size,hll_size(32),hll_size(64),hll_size(128),hll_size(160)\n");
    while (fgets(line, MAX_LINE, iff)) {
        int n = strlen(line);
        if (line[--n] == '\n')
            line[n] = 0;
        if (read_mtx(line, &mm)) {
            continue;
        }
        struct csr csr;
        if (csr_init(&csr, &mm)) {
            mtx_cleanup(&mm);
            continue;
        }
        size_t csr_size = csr_get_size(&csr, mm.nz);
        csr_cleanup(&csr);
        size_t hll_size[4] = {0};
        for (int i = 0; i < 4; i++) {
            struct hll hll;
            if (hll_init(&hll, 32 * (i+1), &mm)) {
                continue;
            }
            hll_size[i] = hll_get_size(&hll);
            hll_cleanup(&hll);
        }
        fprintf(off, "%s,%d,%d,%d,%lu,%lu,%lu,%lu,%lu\n",
            mm_name(line), mm.num_rows, mm.num_cols, mm.nz,
            csr_size,
            hll_size[0], hll_size[1], hll_size[2], hll_size[3]);
        mtx_cleanup(&mm);
    }
    fclose(iff);
    fclose(off);
    return 0;
}
