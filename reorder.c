#include <stdio.h>
#include <assert.h>

#define X(i, j, n)  (((i / (n / 2)) * (n / 2)) +  ((i % (n / 4)) * 2) \
		     + (j / (n / 2)))

#define Y(i, j, n) (((j % (n / 2))  + ((n / 2) * ((i / (n / 4)) % 2))))

#define MEM(a, n, i, j)  a[(i) * (n) + (j)]

void reorder(const int * restrict a, int * restrict b, int size) {
  int i, j;
  int n;

  printf("reorder before\n");
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      printf("%2d ", MEM(a, size, i, j));
    }
    printf("\n");
  }

  n = size;
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      int I = X(i, j, n);
      int J = Y(i, j, n);
      //printf("I: %d, J: %d\n", I, J);
      assert(I >= 0 && I < size);
      assert(J >= 0 && J < size);
      MEM(b, size, i, j) = MEM(a, size, I, J);
    }
  }

  printf("\n");
  printf("reorder after\n");
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      printf("%2d ", MEM(b, size, i, j));
    }
    printf("\n");
  }

}

void print_2d(const char *s, const int *a, int size) {
  int i, j;
  printf("print_2d: %s\n", s);
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      printf("%2d ", MEM(a, size, i, j));
    }
    printf("\n");
  }
}


/* #define N 8 */
/* int main() { */
/*   int a[N * N]; */
/*   int b[N * N]; */
/*   int i, j; */
/*   int size = N; */

/*   for (i = 0; i < size; i++) { */
/*     for (j = 0; j < size; j++) { */
/*       MEM(a, size, i, j) = i * size + j; */
/*     } */
/*   } */

/*   reorder(a, b, N); */
/* } */
