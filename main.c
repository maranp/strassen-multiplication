#include "strass.h"
#include "util.h"
#define N 8

int test(size_t);

#ifndef NDEBUG
void sample_matrix_test();
#endif

int main() {

#ifndef NDEBUG
  sample_matrix_test();
#endif

  if(!test(1)) {
    printf("test(1) failed\n");
  } else {
    printf("test(1) passed\n");
  }

  if(!test(4)) {
    printf("test(4) failed\n");
  } else {
    printf("test(4) passed\n");
  }

  if(!test(16)) {
    printf("test(16) failed\n");
  } else {
    printf("test(16) passed\n");
  }

  if(!test(32)) {
    printf("test(32) failed\n");
  } else {
    printf("test(32) passed\n");
  }

  if(!test(64)) {
    printf("test(64) failed\n");
  } else {
    printf("test(64) passed\n");
  }

  if(!test(512)) {
    printf("test(512) failed\n");
  } else {
    printf("test(512) passed\n");
  }

  if(!test(1024)) {
    printf("test(1024) failed\n");
  } else {
    printf("test(1024) passed\n");
  }

  if(!test(2048)) {
    printf("test(2048) failed\n");
  } else {
    printf("test(2048) passed\n");
  }

  return 0;
}

int test(size_t n) {
  ele_type *unit, *rand, *result;
  //size_t i, j;
  int rc;

  unit = malloc(n * n * sizeof(ele_type));
  rand = malloc(n * n * sizeof(ele_type));
  result = malloc(n * n * sizeof(ele_type));


  unit_matrix(unit, n);
  random_matrix(rand, n);

  strass(rand, unit, result, n);

  printf("resultant matrix:\n");
  rc = memcmp(rand, result, n * n * sizeof(ele_type));
  /* for (i = 0; i < n; i++) { */
  /*   for (j = 0; j < n; j++) { */
  /*     //printf("%d ", MEM(result, n, i, j)); */
  /*     if (MEM(rand, n, i, j) != MEM(result, n, i, j)) { */
  /* 	printf("error\n"); */
  /* 	return 0; */
  /*     } */
  /*   } */
  /*   //printf("\n"); */
  /* } */

  return rc == 0;
}

#ifndef NDEBUG
void sample_matrix_test() {
  /* ele_type a[] = {1, 1, 0, 1, */
  /* 	     0, 1, 0, 0, */
  /* 	     -9, 7, 1, 0, */
  /* 	     1, 1, 0, 1,}; */

  /* ele_type b[] = {1, 0, 0, 0, */
  /* 	     0, 1, 0, 0, */
  /* 	     0, 0, 1, 0, */
  /* 	     0, 0, 0, 1, }; */

  ele_type a[] = {1, 1, 0, 1,  1, 1, 0, 1,
  	     0, 1, 0, 0,  0, 1, 0, 0,
  	     -9, 7, 1, 0, -9, 7, 1, 0,
  	     1, 1, 0, 1,  1, 1, 0, 1,

             1, 1, 0, 1,  1, 1, 0, 1,
  	     0, 1, 0, 0,  0, 1, 0, 0,
  	     -9, 7, 1, 0, -9, 7, 1, 0,
  	     1, 1, 0, 1,  1, 1, 0, 1, };


  ele_type b[] = {1, 0, 0, 0, 0, 0, 0, 0,
  	     0, 1, 0, 0, 0, 0, 0, 0,
  	     0, 0, 1, 0, 0, 0, 0, 0,
  	     0, 0, 0, 1, 0, 0, 0, 0,
  	     0, 0, 0, 0, 1, 0, 0, 0,
  	     0, 0, 0, 0, 0, 1, 0, 0,
  	     0, 0, 0, 0, 0, 0, 1, 0,
  	     0, 0, 0, 0, 0, 0, 0, 1, };

  /* ele_type a[] = {1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1, */
  /* 	     0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0, */
  /* 	     -9, 7, 1, 0, -9, 7, 1, 0, -9, 7, 1, 0, -9, 7, 1, 0, */
  /* 	     1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1, */
	     
  /*            1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1, */
  /* 	     0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0, */
  /* 	     -9, 7, 1, 0, -9, 7, 1, 0, -9, 7, 1, 0, -9, 7, 1, 0, */
  /* 	     1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1, */
             
  /*            1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1, */
  /* 	     0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0, */
  /* 	     -9, 7, 1, 0, -9, 7, 1, 0, -9, 7, 1, 0, -9, 7, 1, 0, */
  /* 	     1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1, */
	     
  /* 	     1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1, */
  /* 	     0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0, */
  /* 	     -9, 7, 1, 0, -9, 7, 1, 0, -9, 7, 1, 0, -9, 7, 1, 0, */
  /* 	     1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1,  1, 1, 0, 1, }; */

  /* ele_type b[] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, */
  /* 	     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, }; */

  ele_type c[N * N] = {0};
  size_t i, j;

  //mult(a, b, c, 4);
  //memset(c, 0, N * N * sizeof(ele_type));
  strass(a, b, c, N);

  printf("resultant matrix:\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      printf("%ld ", MEM(c, N, i, j));
    }
    printf("\n");
  }
}
#endif
