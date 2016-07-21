#include "util.h"
#define MEM(a, n, i, j)  a[(i) * (n) + (j)]

void unit_matrix(ele_type *a, size_t size) {
  size_t i;
  memset(a, 0, size * size * sizeof(ele_type));
  for (i = 0; i < size; i++) {
	MEM(a, size, i, i) = 1;
    }
}

void random_matrix(ele_type *a, size_t size) {
  size_t i, j;
  time_t t;
  /* set the seed for random numbers */
  srand((unsigned) time(&t));

  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      MEM(a, size, i, j) = rand();
    }
  }
}
