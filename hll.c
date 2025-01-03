#include "hll.h"
#include <stdlib.h>
#include <string.h>

static inline int __max(int *v, int lo, int hi_exclusive) {
    int max = -1;
    for (int i = lo; i < hi_exclusive; ++i) {
        if (max < v[i]) {
            max = v[i];
        }
    }
    return max;
}

/**
  * __get_nonzeros calculates the number of nonzero items for each row in the
  * matrix mm and saves them in nzr.
  * @param mm  (in)  matrix
  * @param nzr (out) vector thats holds the number of nonzeros for each row of mm -i.e.
  * nzr[i] holds the number of nonzeros for the i-th row in mm
  */
void __get_nonzeros(struct MatrixMarket *mm, int *nzr) {
    int curr_nz = 0;
    int prev = mm->rows[0];
    int k = 0;
    for (int i = 0; i < mm->nz; ++i) {
        int curr = mm->rows[i];
        if (prev != curr) {
            prev = curr;
            nzr[k++] = curr_nz;
            curr_nz = 0;
        }
        curr_nz++;
    }
    // the last row is skipped in the loop
    nzr[k++] = curr_nz;
}

int __get_data_size(int *nzr, int hack_size, int num_rows) {
    int data_size = 0;
    int start = 0;
    while (start + hack_size < num_rows) {
        data_size += hack_size * __max(nzr, start, start + hack_size);
        start += hack_size;
    }
    return data_size + (num_rows - start) * __max(nzr, start, num_rows);
}

static void inline ifill(int *v, int start, int n, int val) {
    for (int i = 0; i < n; ++i) {
        v[start + i] = val;
    }
}

static void inline dfill(double *v, int start, int n, double val) {
    for (int i = 0; i < n; ++i) {
        v[start + i] = val;
    }
}

static inline void __vcopy(void *dest, int dstart, void *src, int sstart, int n, int itmsz) {
    memcpy(dest + (dstart * itmsz), src + (sstart * itmsz), n * itmsz);
}

/**
 * __fill_padding add n values of padding to data and col_index if the row has less non-zero
 * values of the maximum number of non-zeros in the hack:
 * - 0 to data
 * - the last index of column in the row to col_index
 * @param mm (in)   it is needed to determine if the matrix contains integers or doubles
 * @param hll (out) wip hll representation of hll 
 * @param hllp (in) starting index of hll->data and hll->col_index where to start writing padding
 * @param n (in)    number of values of padding to write
 */
static inline void __fill_padding(struct MatrixMarket *mm, struct hll *hll, int hllp, int n) {
    if (n > 0) {
        ifill(hll->col_index, hllp, n, hll->col_index[hllp-1]);
        if (mm_is_integer(mm->typecode)) {
            ifill((int*)hll->data, hllp, n, 0);
        } else {
            dfill((double*)hll->data, hllp, n, 0.0);
        }
    }
}

/**
 * __read_row reads a row of a sparse matrix in MatrixMarket format and produces
 * the equivalent HLL representation of that row.
 * @param mm (in)       input matrix in MatrixMarket format
 * @param mmp (in/out)  pointer to the next matrix component to read
 * @param hll (in/out)  hll representation of mm 
 * @param hllp (hllp)   pointer to the free area of data and col_index where to save
 *                      matrix components and column indices respectively
 * @param maxnzr (in)   max number of non-zeros in the hack
 */
void __read_row(struct MatrixMarket *mm, int *mmp, struct hll *hll, int *hllp, int maxnzr) {
    // Legge tutta la riga della matrice:
    // n è il numero di elementi non zero della riga
    int row = mm->rows[*mmp];
    int n = 0;
    for (int i = *mmp; i < mm->nz && row == mm->rows[i]; ++i) {
        ++n;
    }
    
    __vcopy(hll->data, *hllp, mm->data, *mmp, n, get_element_size(mm));
    __vcopy(hll->col_index, *hllp, mm->cols, *mmp, n, sizeof(int));
    *hllp += n;
    *mmp += n;

    __fill_padding(mm, hll, *hllp, maxnzr - n);
    *hllp += maxnzr - n;
}

/**
 * @param hll (out)      HLL representation of matrix mm
 * @param hack_size (in) number of rows of each hack
 * @param mm (in)        input matrix in MatrixMarket format
 * */  
int hll_init(struct hll *hll, int hack_size, struct MatrixMarket *mm) {
    // nzr is the array that keeps track of how many non-zeros each
    // row has -i.e. nzr[i] is the number of non-zeros of the i-th row of mm
    int *nzr = malloc(mm->num_rows * sizeof(int));
    if (nzr == NULL) {
        return 1;
    }
    __get_nonzeros(mm, nzr);

    int hacks_num = mm->num_rows / hack_size;
    if (mm->num_rows % hack_size) {
        hacks_num++;
    }
    hll->hacks_num = hacks_num;

    hll->offsets = malloc((hacks_num + 1) * sizeof(int));
    if (hll->offsets == NULL) {
        free(nzr);
        return 1;
    }
    int size = __get_data_size(nzr, hack_size, mm->num_rows);
    hll->data = malloc(size * get_element_size(mm));
    if (hll->data == NULL) {
        free(nzr);
        free(hll->offsets);
        hll->offsets = NULL;
        return 1;
    }
    hll->col_index = malloc(size * sizeof(int));
    if (hll->col_index == NULL) {
        free(nzr);
        free(hll->offsets);
        hll->offsets = NULL;
        free(hll->data);
        hll->data = NULL;
        return 1;
    }
    
    int hllp = 0;
    int offp = 0;
    int mmp = 0;
    int lo = 0; // current number of rows read
    hll->offsets[offp++] = 0;
    while (hacks_num-- > 0) {
        // the last hack may have an hack_size smaller than the other
        // because the number of rows of a matrix may not be a multiple o hack_size
        // if it is not the case mm->num_rows - lo is equal to hack_size, otherwise
        // it is equal to mm->num_rows % hack_size
        if (hacks_num == 0) {
            hack_size = mm->num_rows - lo;
        }
        // maxznr is the maximum number of non-zeros in the current hack -i.e.
        // maxznr = max(nzr[lo], nzr[lo+1], ..., nzr[lo+hack_size-1])
        int maxnzr = __max(nzr, lo, lo + hack_size);
        for (int i = 0; i < hack_size; ++i) {
            __read_row(mm, &mmp, hll, &hllp, maxnzr);
        }
        hll->offsets[offp++] = hllp;
        lo += hack_size;
    }
    hll->offsets_num = offp;
    if (hllp != size) {
        printf("elements read mismatch with elements expected\n");
        free(nzr);
        free(hll->offsets);
        hll->offsets = NULL;
        free(hll->data);
        hll->data = NULL;
        free(hll->col_index);
        hll->col_index = NULL;
        return 1;
    } else {
        hll->data_num = hllp;
        free(nzr);
        return 0;
    }
}