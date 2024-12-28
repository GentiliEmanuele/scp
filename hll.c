#include "hll.h"
#include <stdlib.h>

int __max(int *v, int lo, int hi_exclusive)
{
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
  * @param mm  (in) matrix
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

int __set_offsets(int *offsets, int hack_size, struct MatrixMarket *mm) {
    // TODO
    return 0;
}

int __get_data_size(int *nzr, int hack_size, int num_rows) {
    int hacks_num = num_rows / hack_size;
    int data_size = 0;
    int start = 0;
    for (int i = 0; i < hacks_num; i++) {
        data_size += hack_size * __max(nzr, start, i + hack_size);
        start += hack_size;
    }
    return data_size + (num_rows - start) * __max(nzr, start, num_rows);
}

void ifill(int *v, int start, int n, int val) {
    for (int i = start; i < n; ++i) {
        v[i] = val;
    }
}

void dfill(double *v, int start, int n, double val) {
    for (int i = start; i < n; ++i) {
        v[i] = val;
    }
}

/**
  * @param mm (in)        matrice
  * @param nzr (in)       numero di elementi non zero per ogni riga della matrice mm
  * @param data (out)        
  * @param indices (out)
  * @param offsets (out)
  * @param hack_size (in)
  **/
void read_ellpack(struct MatrixMarket *mm, struct hll *hll, int hack_size)
{
    // nzr[i] contains the number of non-zeros for the i-th row
    int *nzr = malloc(mm->num_rows * sizeof(int));
    __get_nonzeros(mm, nzr);
    int hacks_num = mm->num_rows / hack_size;
    int *offsets = malloc((hacks_num + 1) * sizeof(int));
    int data_size = __get_data_size(nzr, hack_size, hacks_num);
    void *data = malloc(data_size * get_element_size(mm));
    int *indices = malloc(data_size * sizeof(int));
    // indice della prima riga iniziale dell'hack
    // serve a __max per calcolare il massimo numero
    // di non zeri di un hack
    int lo = 0;
    // è il massimo numero di non zeri di un hack
    int maxnzr = __max(nzr, lo, lo + hack_size);
    // prev, curr demarcano rispettivamente l'indice di riga
    // dell'elemento precedente e quello successivo, servono per capire
    // se si sta leggendo righe dello stesso hack oppure hack differenti
    int prev = mm->rows[0];
    offsets[0] = 0;
    int offsets_top = 0;
    int data_top = 0;
    int mm_top = 0;
    int i = 0;
    while (lo + hack_size < mm->num_rows) {
    
    }
    hack_size = mm->num_rows - lo;
    while (i < mm->nz) {
        if (mm_is_integer(mm->typecode)) {
            ((int*)data)[data_top] = ((int*)mm->data)[mm_top];
        } else {
            ((double*)data)[data_top] = ((double*)mm->data)[mm_top];
        }
        indices[data_top] = mm->cols[mm_top];
        data_top++;
        mm_top++;
        i++;
    }
    maxnzr = __max(nzr, lo, lo + hack_size);
    int padding = maxnzr - nzr[prev];
    if (mm_is_integer(mm->typecode)) {
        ifill(data, data_top, padding, 0);
    } else {
        dfill(data, data_top, padding, 0);
    }
    ifill(indices, data_top, padding, indices[data_top-1]);
    data_top += padding; 
    offsets[offsets_top++] = data_top;
    hll->data = data;
    hll->col_index = indices;
    hll->offsets = offsets;   
}

/**
  * @param hack_size number of rows of each hack */
int __hll_init(struct hll *hll, struct MatrixMarket *mm, int hack_size)
{    
    read_ellpack(mm, hll, hack_size);
}


int hll_iinit(struct hll *hll, struct MatrixMarket *mm, int hack_size)
{
    return __hll_init(hll, mm, hack_size);

}

int hll_dinit(struct hll *hll, struct MatrixMarket *mm, int hack_size)
{
    return __hll_init(hll, mm, hack_size);
}