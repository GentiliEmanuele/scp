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
    struct MatrixMarket mm;
    if (read_mtx(argv[--argc], &mm)) {
        printf("Error");
        return 1;
    }
    
    double data[256];
    int row_pointers[256];
    int col_index[256];
    struct csr csr = {
        .col_index = col_index,
        .row_pointer = row_pointers,
        .data = data,
        .num_rows = 9,
        .num_cols = 9
    };
    if (csr_dinit(&csr, &mm)) { 
        printf("Error!\n");
    }
    mm_write_banner(stdout, mm.typecode);
    mm_write_mtx_crd_size(stdout, mm.num_rows, mm.num_cols, mm.nz);
	return 0;
}

