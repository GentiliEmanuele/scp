#include "vec.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
    double a[] = {1, 2, 3, 4};
    double b[] = {1, 2, 3, 4};

    if (!d_veceq(a, b, 4, 1e-6)) {
        printf("a != b\n");
    } else {
        printf("a == b\n");
    }
}