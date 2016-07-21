#include "strass.h"

/*
   (a b) * (e f) = (r s)
   (c d)   (g h)   (t u)
   at each level of recursion, before subdividing the input matrices,
   reorder them such each submatrix are placed contiguously in reordered
   matrix.

   For example, reorder function will rearrange the matrix elements such that
   their corresponding position in new array are as shown below

   original              reordered
   0  1  2  3      --->  0  1  4  5
   4  5  6  7            2  3  6  7
   8  9  10 11           8  9  12 13
   12 13 14 15           10 11 14 15

   if a(the top left quadrant of original 4 X 4 matrix) occupied indices
       a = (0 1)
           (4 5),
  in the reordered matrix, a occupies 0, 1, 2, 3 position, hence contiguous.
  Similar explanation applies for other 3 quadrants.

  Now, the split submatrices are already contiguous, which makes the further
  operations cache oblivious.

  The macros X and Y does the mapping of original (i, j) to reordered (i, j)
  where the reordering is as explained above.
*/
#define X(i, j, n)  (((i / (n / 2)) * (n / 2)) +  ((i % (n / 4)) * 2) \
		     + (j / (n / 2)))

#define Y(i, j, n) (((j % (n / 2))  + ((n / 2) * ((i / (n / 4)) % 2))))

static void reorder(const ele_type * a, ele_type * b, size_t size) {
  size_t i, j;
  size_t n;

  n = size;
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      size_t I = X(i, j, n);
      size_t J = Y(i, j, n);
      //printf("I: %d, J: %d\n", I, J);
      assert(I >= 0 && I < size);
      assert(J >= 0 && J < size);
      MEM(b, size, i, j) = MEM(a, size, I, J);
    }
  }

}

#ifndef NDEBUG
static void print_2d(const char *s, const ele_type *a, size_t size) {
  size_t i, j;
  printf("print_2d: %s\n", s);
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      printf("%2d ", MEM(a, size, i, j));
    }
    printf("\n");
  }
}

static void print_arr(const ele_type *a, size_t n) {
  size_t i;
  for(i = 0; i < n; i++) {
    printf("%d ", a[i]);
  }
  printf("\n");
}

#endif

/* it's not worth the overhead of blocked and strassen multiplication
   when the order of matrix is small. Do the traditional ijk multiplication */
static ele_type * mult_ijk(const ele_type * a,
			   const ele_type * b,
			   ele_type * c, size_t size) {
  size_t i, j, k;
  ele_type sum = 0;
  memset(c, 0, sizeof(ele_type) * size * size);
  for(i = 0; i < size; i++) {
    for(j = 0; j < size; j++) {
      sum = 0;
      for(k = 0; k < size; k++) {
	sum += MEM(a, size, i, k) * MEM(b, size, k, j);
      }
      MEM(c, size, i, j) = sum;
    }
  }
  return c;
}


/* blocked multiplication.
   basically, a strip (of block size) of matrix A is multiplied
   with a block (of block size) of matrix B and resultant matrix's
   strip is updated at each iteration. This method helps to minimise cache
   misses due colum wise access of matrix B in traditional ijk algorithm
   the block size should be chosen such that the block of B can fit in the
   cache (depends on cache configuration)
 */
static ele_type * mult_block(const ele_type * a,
	  const ele_type * b,
	  ele_type * c, size_t size) {
  size_t i, j, k;
  size_t j0, k0;
  ele_type sum = 0;
  size_t en = BLOCK_SIZE * (size / BLOCK_SIZE);
  memset(c, 0, sizeof(ele_type) * size * size);
  for(k0 = 0; k0 < en; k0 += BLOCK_SIZE) {
    for(j0 = 0; j0 < en; j0 += BLOCK_SIZE) {
      for(i = 0; i < size; i++) {
	for(j = j0; j < j0 + BLOCK_SIZE; j++) {
	  sum = MEM(c, size, i, j);
	  for(k = k0; k < k0 + BLOCK_SIZE; k++) {
	    sum += MEM(a, size, i, k) * MEM(b, size, k, j);
	  }
	  MEM(c, size, i, j) = sum;
	}
      }
    }
  }
  return c;
}

static void set_quad(ele_type *M, const ele_type *m, size_t size,
	      size_t rl1, size_t rl2, size_t cl1, size_t cl2) {
  size_t i, j, cnt;
  assert(M && m);

  cnt = 0;
  for(i = rl1; i < rl2; i++) {
    for(j = cl1; j < cl2; j++) {
      MEM(M, size, i, j) = m[cnt++];
    }
  }

  return;
}

static ele_type * add(const char *s, const ele_type *a, const ele_type *b, ele_type *r, size_t size) {
  size_t i;
  assert(r);

  for (i = 0; i < size * size; i++) {
    r[i] = a[i] + b[i];
  }

  return r;
}

static ele_type * minus(const char *s, const ele_type *a, const ele_type *b, ele_type *r, size_t size) {
  size_t i;
  assert(r);

  for (i = 0; i < size * size; i++) {
    r[i] = a[i] - b[i];
  }

  return r;
}

ele_type * strass(const ele_type * A,
	    const ele_type * B,
	    ele_type * C, size_t size) {
  size_t i;
  ele_type *a = 0;
  ele_type *b = 0;
  ele_type *c = 0;
  ele_type *d = 0;
  ele_type *e = 0;
  ele_type *f = 0;
  ele_type *g = 0;
  ele_type *h = 0;
  ele_type *r = 0;
  ele_type *s = 0;
  ele_type *t = 0;
  ele_type *u = 0;
  ele_type *Aa[8] = {0};
  ele_type *Bb[8] = {0};
  ele_type *P[8] = {0};;

  ele_type *t1 = 0;
  ele_type *t2 = 0;
  ele_type *Atemp = 0;
  ele_type *Btemp = 0;

  if (size <= BLOCK_CUTOFF) {
    return mult_ijk(A, B, C, size);;
  }

  if (size <= RECURSION_CUTOFF) {
    return mult_block(A, B, C, size);;
  }

  /* all allocations */
  /* 0th element of A, B, P unused */
  for (i = 1; i < 8; i++) {
    posix_memalign((void **)&Aa[i], ALIGN, sizeof(ele_type) * size * size / 4);
    assert(Aa[i] != NULL);
    posix_memalign((void **)&Bb[i], ALIGN, sizeof(ele_type) * size * size / 4);
    assert(Bb[i] != NULL);
    posix_memalign((void **)&P[i], ALIGN, sizeof(ele_type) * size * size / 4);
    assert(P[i] != NULL);
  }
  posix_memalign((void **)&t1, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(t1 != NULL);
  posix_memalign((void **)&t2, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(t2 != NULL);

  posix_memalign((void **)&r, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(r != NULL);
  posix_memalign((void **)&s, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(s != NULL);
  posix_memalign((void **)&t, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(t != NULL);
  posix_memalign((void **)&u, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(u != NULL);

  /* reordered matrices to apply strassen operators */
  posix_memalign((void **)&Atemp, ALIGN, sizeof(ele_type) * size * size);
  assert(Atemp != NULL);
  posix_memalign((void **)&Btemp, ALIGN, sizeof(ele_type) * size * size);
  assert(Btemp != NULL);

  reorder(A, Atemp, size);
  reorder(B, Btemp, size);

  a = Atemp + size * size * 0 / 4;

  b = Atemp + size * size * 1 / 4;

  c = Atemp + size * size * 2 / 4;

  d = Atemp + size * size * 3 / 4;

  // Aa[0] unused
  /* A1 = a */
  /* A2 = a + b */
  /* A3 = c + d */
  /* A4 = d */
  /* A5 = a + d */
  /* A6 = b - d */
  /* A7 = a - c */

  Aa[1] = a;
  add("a + b", a, b, Aa[2], size / 2);
  add("c + d", c, d, Aa[3], size / 2);

  Aa[4] = d;
  add("a + d", a, d, Aa[5], size / 2);

  minus("b - d", b, d, Aa[6], size / 2);

  minus("a - c", a, c, Aa[7], size / 2);

  e = Btemp + size * size * 0 / 4;
  f = Btemp + size * size * 1 / 4;
  g = Btemp + size * size * 2 / 4;
  h = Btemp + size * size * 3 / 4;

  // Bb[0] unused
  /* B1 = f - h */
  /* B2 = h */
  /* B3 = e */
  /* B4 = g - e */
  /* B5 = e + h */
  /* B6 = g + h */
  /* B7 = e + f */
  minus("f - h", f, h, Bb[1], size / 2);
  Bb[2] = h;
  Bb[3] = e;

  minus("g - e", g, e, Bb[4], size / 2);

  add("e + h", e, h, Bb[5], size / 2);

  add("g + h", g, h, Bb[6], size / 2);

  add("e + f", e, f, Bb[7], size / 2);

  // P[0] unused
  for (i = 1; i <= 7; i++)  {
    strass(Aa[i], Bb[i], P[i], size / 2);
#ifndef NDEBUG
    printf("P[%d]: ", i); 
    print_arr(P[i], size * size / 4);
#endif
  }

  /* r = P5 + P4 - P2 + P6; */
  /* s = P1 + P2; */
  /* t = P3 + P4; */
  /* u = P5 + P1 - P3 - P7; */

  add("p5 + p4", P[5], P[4], t1, size / 2);

  minus("t1 - p[2]", t1, P[2], t2, size / 2);

  add("t2 + p6", t2, P[6], r, size / 2);
  
  add("p1 + p2", P[1], P[2], s, size / 2);

  add("p3 + p4", P[3], P[4], t, size / 2);

  add("p5 + p1", P[5], P[1], t1, size / 2);

  minus("t1 - p3", t1, P[3], t2, size / 2);

  minus("t2 - p7", t2, P[7], u, size / 2);

#ifndef NDEBUG
  printf("result:\n");
  printf("r:\n");
  for(i = 0; i < size * size / 4; i++) {
    printf("%d ", r[i]);
  }
  printf("\n");

  printf("s:\n");
  for(i = 0; i < size * size / 4; i++) {
    printf("%d ", s[i]);
  }
  printf("\n");

  printf("t:\n");
  for(i = 0; i < size * size / 4; i++) {
    printf("%d ", t[i]);
  }
  printf("\n");

  printf("u:\n");
  for(i = 0; i < size * size / 4; i++) {
    printf("%d ", u[i]);
  }
  printf("\n");
#endif

  /* merge r s t u */
  set_quad(C, r, size, 0, size / 2, 0, size / 2);
  set_quad(C, s, size, 0, size / 2, size / 2, size);
  set_quad(C, t, size, size / 2, size, 0, size / 2);
  set_quad(C, u, size, size / 2, size, size / 2, size);

  /* all frees */
  for (i = 0; i < 8; i++) {
    if (i != 1 && i != 4)
      free(Aa[i]);
    if (i != 2 && i != 3)
      free(Bb[i]);
    free(P[i]);
  }
  free(t1);
  free(t2);
  free(Atemp);
  free(Btemp);
  free(r);
  free(s);
  free(t);
  free(u);


  return C;
}
