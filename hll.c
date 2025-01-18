#include "hll.h"
#include <stdlib.h>
#include <string.h>

static inline int max(int *v, int lo, int hi_exclusive) {
    int max = -1;
    for (int i = lo; i < hi_exclusive; ++i) {
        if (max < v[i]) {
            max = v[i];
        }
    }
    return max;
}

/**
  * get_nonzeros calculates the number of nonzero items for each row in the
  * matrix mm and saves them in nzr.
  * @param mm  (in)  matrix
  * @param nzr (out) vector thats holds the number of nonzeros for each row of mm -i.e.
  * nzr[i] holds the number of nonzeros for the i-th row in mm
  */
static void get_nonzeros(struct MatrixMarket *mm, int *nzr) {
    memset(nzr, 0, sizeof(int) * mm->num_rows);
    for (int i = 0; i < mm->nz; ++i) {
        nzr[mm->rows[i]]++; 
    }
}

static int get_data_size(int *nzr, int hack_size, int num_rows) {
    int data_size = 0;
    int start = 0;
    while (start + hack_size < num_rows) {
        data_size += hack_size * max(nzr, start, start + hack_size);
        start += hack_size;
    }
    return data_size + (num_rows - start) * max(nzr, start, num_rows);
}

static inline void ifill(int *v, int start, int n, int val) {
    for (int i = 0; i < n; ++i) {
        v[start + i] = val;
    }
}

static inline void dfill(double *v, int start, int n, double val) {
    for (int i = 0; i < n; ++i) {
        v[start + i] = val;
    }
}

/**
 * read_row reads a row of a sparse matrix in MatrixMarket format which has
 * at least a non-zero component and produces
 * the equivalent HLL representation of that row.
 * @param mm (in)       input matrix in MatrixMarket format
 * @param mmp (in/out)  pointer to the next matrix component to read
 * @param hll (in/out)  hll representation of mm 
 * @param hllp (hllp)   pointer to the free area of data and col_index where to save
 *                      matrix components and column indices respectively
 * @param maxnzr (in)   max number of non-zeros in the hack
 */
static void read_row(struct MatrixMarket *mm, int *mmp, struct hll *hll, int *hllp, int maxnzr) {
    // Legge tutta la riga della matrice:
    // n è il numero di elementi non zero della riga
    int row = mm->rows[*mmp];
    int n = 0;
    for (int i = *mmp; i < mm->nz && row == mm->rows[i]; ++i) {
        ++n;
    }
    
    memcpy(&hll->data[*hllp], &mm->data[*mmp], n * sizeof(double));
    memcpy(&hll->col_index[*hllp], &mm->cols[*mmp], n * sizeof(int));
    *hllp += n;
    *mmp += n;

    ifill(hll->col_index, *hllp, maxnzr - n, hll->col_index[*hllp - 1]);
    dfill(hll->data, *hllp, maxnzr - n, 0.0);
    *hllp += maxnzr - n;
}

static inline int get_hacks_num(int num_rows, int hack_size) {
    int hacks_num = num_rows / hack_size;
    if (num_rows % hack_size) {
        hacks_num++;
    }
    return hacks_num;
}

void add_null_row(struct MatrixMarket *mm, struct hll *hll, int *hllp, int maxnzr) {
    ifill(hll->col_index, *hllp, maxnzr, 0);
    dfill(hll->data, *hllp, maxnzr, 0.0);
    *hllp += maxnzr;
}

/**
 * @param hll (out)      HLL representation of matrix mm
 * @param hack_size (in) number of rows of each hack
 * @param mm (in)        input matrix in MatrixMarket format
 * */  
int hll_init(struct hll *hll, int hack_size, struct MatrixMarket *mm) {
    hll->num_rows = mm->num_rows;
    hll->num_cols = mm->num_cols;
    hll->hack_size = hack_size;
    // nzr is the array that keeps track of how many non-zeros each
    // row has -i.e. nzr[i] is the number of non-zeros of the i-th row of mm
    int *nzr = malloc(mm->num_rows * sizeof(int));
    if (nzr == NULL) {
        return 1;
    }
    get_nonzeros(mm, nzr);

    if (hack_size > hll->num_rows) {
        hack_size = hll->num_rows;
        hll->hack_size = hack_size;
#ifdef SCP_VERBOSE
        printf("warning (%s): hack size > number of rows, hack size set to number of rows\n", mm->path);
#endif
    }

    hll->hacks_num = get_hacks_num(hll->num_rows, hack_size);
    hll->offsets = malloc((hll->hacks_num + 1) * sizeof(int));
    if (hll->offsets == NULL) {
        free(nzr);
        return 1;
    }

    int size = get_data_size(nzr, hack_size, mm->num_rows);
    hll->data = malloc(size * sizeof(double));
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
        free(hll->data);
        hll->offsets = NULL;
        hll->data = NULL;
        return 1;
    }
    hll->max_nzr = malloc(sizeof(int) * hll->hacks_num);
    if (hll->max_nzr == NULL) {
        free(nzr);
        free(hll->offsets);
        free(hll->data);
        free(hll->col_index);
        hll->offsets = NULL;
        hll->data = NULL;
        hll->col_index = NULL;
        return 1;
    }
    int maxnzrp = 0;
    int hllp = 0;
    int offp = 0;
    int mmp = 0;
    int lo = 0; // current number of rows read
    hll->offsets[offp++] = 0;
    int hacks_num = hll->hacks_num;
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
        int maxnzr = max(nzr, lo, lo + hack_size);
        if (maxnzr > 0) {
            for (int i = 0; i < hack_size; ++i) {
                if (nzr[lo + i]) {
                    read_row(mm, &mmp, hll, &hllp, maxnzr);
                } else {
                    add_null_row(mm, hll, &hllp, maxnzr);
                }
            }
        }
        hll->max_nzr[maxnzrp++] = maxnzr;
        hll->offsets[offp++] = hllp;
        lo += hack_size;
    }
    hll->offsets_num = offp;
    if (hllp != size) {
        printf("elements read mismatch with elements expected\n");
        free(nzr);
        free(hll->offsets);
        free(hll->data);
        free(hll->col_index);
        free(hll->max_nzr);
        nzr = NULL;
        hll->offsets = NULL;
        hll->data = NULL;
        hll->col_index = NULL;
        hll->max_nzr = NULL;
        return 1;
    } else {
        hll->data_num = hllp;
        return 0;
    }
}

size_t hll_get_footprint(struct hll *hll) {
    return hll->data_num * sizeof(double) + hll->data_num * sizeof(int)
    + hll->hacks_num * sizeof(int) + hll->offsets_num * sizeof(int) + sizeof(struct hll);
}

void hll_cleanup(struct hll *hll) {
    free(hll->offsets);
    free(hll->col_index);
    free(hll->data);
    free(hll->max_nzr);
}