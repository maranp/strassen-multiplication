#ifndef UTIL_H
#define UTIL_H
#include "strass.h"
#include <time.h>

#define MEM(a, n, i, j)  a[(i) * (n) + (j)]

/* helper function to create a unit matrix for testing */
void unit_matrix(ele_type *a, size_t size);
/* helper function to create a random matrix for testing */
void random_matrix(ele_type *a, size_t size);

#endif // UTIL_H
