/* 
*   Matrix Market I/O example program
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   Usage:  a.out [filename] > output
*
*       
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing 
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "csr.h"
#include "utils.h"

int main(int argc, char *argv[])
{
    double *val;
    int *I, *J, M, N, nz;
    if (read_mtx(argv[--argc], &val, &I, &J, &M, &N, &nz)) {
        printf("Error");
        return 1;
    }
    
    double data[256];
    int row_pointers[256];
    int col_index[256];
    struct csr *csr= malloc(sizeof(struct csr *));
    csr -> col_index = col_index;
    csr -> row_pointer = row_pointers;
    csr -> data = data;
    csr -> num_rows = 9;
    csr -> num_cols = 9;
    /*
    struct csr csr = {
        .col_index = col_index,
        .row_pointer = row_pointers,
        .data = data,
        .num_rows = 9,
        .num_cols = 9
    };
    */
    if (csr_dinit(csr, val, J, I, nz)) {
        printf("Error!\n");
    }
    printf("==================================\n");
    for (int i = 0; i < nz; ++i) {
        printf("%d %20.19g\n", csr->col_index[i], data[i]);
    }
    for (int i = 0; i < 10; i++) {
        printf("i=%d rp=%d \n", i, csr->row_pointer[i]);
    }  
    printf("ok\n");
	return 0;
}

