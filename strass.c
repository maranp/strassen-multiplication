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
      assert(I < size);
      assert(J < size);
      MEM(b, size, i, j) = MEM(a, size, I, J);
    }
  }

}

/* it's not worth the overhead of blocked and strassen multiplication
   when the order of matrix is small. Do the traditional ijk multiplication */
static void ijk_multiply(const ele_type * a,
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
}


/* blocked multiplication.
   basically, a strip (of block size) of matrix A is multiplied
   with a block (of block size) of matrix B and resultant matrix's
   strip is updated at each iteration. This method helps to minimise cache
   misses that could otherwise occur due to columnwise access of matrix B in
   traditional ijk algorithm. The block size should be chosen such that the
   block of B can fit in the cache (depends on cache configuration)
 */
static void block_multiply(const ele_type * a,
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

/* r = a + b
   where r, a, b are 1d vectors of length 'size' */
static ele_type * add(const ele_type *a, const ele_type *b, ele_type *r, size_t size) {
  size_t i;

  assert(a && b && r);

  /* vectorization friendly. SIMD instructions are used at -O3 (or -ftree-vectorize) */
  for (i = 0; i < size * size; i++) {
    r[i] = a[i] + b[i];
  }

  return r;
}

/* r = a - b
   where r, a, b are 1D vectors of length 'size' */
static ele_type * minus(const ele_type *a, const ele_type *b, ele_type *r, size_t size) {
  size_t i;

  assert(a && b && r);

  /* vectorization friendly. SIMD instructions are used at -O3 (or -ftree-vectorize) */
  for (i = 0; i < size * size; i++) {
    r[i] = a[i] - b[i];
  }

  return r;
}

static void * strassen_multiply_wrapper(void *arg) {
  strassen_params_t *params = arg;
  strassen_multiply(params->aterm, params->bterm, params->pterm, params->size);

  return NULL;
}


/* strassen multiplication:
   mat_c = mat_a * mat_b
   where mat_a, mat_b and mat_c are square matrix of length size * size.
   Quite a long function but covers the complete strassen method
   in single function body */

void strassen_multiply(const ele_type * mat_a,
		       const ele_type * mat_b,
		       ele_type * mat_c, size_t size) {
  /* index variable */
  size_t i;
  /* a, b, c, d,
     e, f, g, h
     r, s, t, u
     are intermediate submatrices of length
     (size / 2) * (size / 2). */
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
 
  /* Aterm saves the intermediate results of operations on mat_a's 4 quadrants */
  ele_type *Aterm[8] = {0};

  /* Aterm saves the intermediate results of operations on mat_a's 4 quadrants */
  ele_type *Bterm[8] = {0};

  /* Pterm saves the intermediate results of products Aterm and Bterm. */
  /* The purpose of Aterm, Bterm and P will be more evident if one refers to 
     strassen algorithm description.
     (http://www.cs.mcgill.ca/~pnguyen/251F09/matrix-mult.pdf) */
  ele_type *P[8] = {0};;

  /* temporary vectors */
  ele_type *t1 = 0;
  ele_type *t2 = 0;

  /* temporaries to save reordered A and B matrices */
  ele_type *A_reordered = 0;
  ele_type *B_reordered = 0;

  int rc; /* track the return value of library functions */

  if (size <= BLOCK_CUTOFF) {
    return ijk_multiply(mat_a, mat_b, mat_c, size);;
  }

  if (size <= RECURSION_CUTOFF) {
    return block_multiply(mat_a, mat_b, mat_c, size);;
  }

  /* all allocations */
  /* 0th element of Aterm, Bterm, P unused */
  for (i = 1; i < 8; i++) {
    rc = posix_memalign((void **)&Aterm[i], ALIGN, sizeof(ele_type) * size * size / 4);
    assert(rc == 0 && Aterm[i] != NULL);
    rc = posix_memalign((void **)&Bterm[i], ALIGN, sizeof(ele_type) * size * size / 4);
    assert(rc == 0 && Bterm[i] != NULL);
    rc = posix_memalign((void **)&P[i], ALIGN, sizeof(ele_type) * size * size / 4);
    assert(rc == 0 && P[i] != NULL);
  }
  rc = posix_memalign((void **)&t1, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(rc == 0 && t1 != NULL);
  rc = posix_memalign((void **)&t2, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(rc == 0 && t2 != NULL);

  rc = posix_memalign((void **)&r, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(rc == 0 && r != NULL);
  rc = posix_memalign((void **)&s, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(rc == 0 && s != NULL);
  rc = posix_memalign((void **)&t, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(rc == 0 && t != NULL);
  rc = posix_memalign((void **)&u, ALIGN, sizeof(ele_type) * size * size / 4);
  assert(rc == 0 && u != NULL);

  /* reordered matrices to apply strassen operators */
  rc = posix_memalign((void **)&A_reordered, ALIGN, sizeof(ele_type) * size * size);
  assert(rc == 0 && A_reordered != NULL);
  rc = posix_memalign((void **)&B_reordered, ALIGN, sizeof(ele_type) * size * size);
  assert(rc == 0 && B_reordered != NULL);

  reorder(mat_a, A_reordered, size);
  reorder(mat_b, B_reordered, size);

  /* since mat_a is reordered. The submatrices (4 quadrants)
     are now contiguous in memory. So, designating the submatrices is just pointing to
     the coorect memory address where each quadrant starts.
     mat_a = (a b)
             (c d)
  */
  a = A_reordered + size * size * 0 / 4;
  b = A_reordered + size * size * 1 / 4;
  c = A_reordered + size * size * 2 / 4;
  d = A_reordered + size * size * 3 / 4;

  /* Aterm[0] unused */
  /* Aterm1 = a */
  Aterm[1] = a;

  /* Aterm2 = a + b */
  add(a, b, Aterm[2], size / 2);

  /* Aterm3 = c + d */
  add(c, d, Aterm[3], size / 2);

  /* Aterm4 = d */
  Aterm[4] = d;

  /* Aterm5 = a + d */
  add(a, d, Aterm[5], size / 2);

  /* Aterm6 = b - d */
  minus(b, d, Aterm[6], size / 2);

  /* Aterm7 = a - c */
  minus(a, c, Aterm[7], size / 2);

  /* since mat_b is reordered. The submatrices (4 quadrants)
     are now contiguous in memory. So, designating the submatrices is just pointing to
     the coorect memory address where each quadrant starts.
     mat_b = (e f)
             (g h)
  */
  e = B_reordered + size * size * 0 / 4;
  f = B_reordered + size * size * 1 / 4;
  g = B_reordered + size * size * 2 / 4;
  h = B_reordered + size * size * 3 / 4;

  // Bterm[0] unused
  /* Bterm1 = f - h */
  minus(f, h, Bterm[1], size / 2);

  /* Bterm2 = h */
  /* Bterm3 = e */
  Bterm[2] = h;
  Bterm[3] = e;

  /* Bterm4 = g - e */
  minus(g, e, Bterm[4], size / 2);

  /* Bterm5 = e + h */
  add(e, h, Bterm[5], size / 2);

  /* Bterm6 = g + h */
  add(g, h, Bterm[6], size / 2);

  /* Bterm7 = e + f */
  add(e, f, Bterm[7], size / 2);

  // P[0] unused
  pthread_t thread[8];
  strassen_params_t *strassen_params[8];
  for (i = 1; i <= 7; i++)  {
    /* P[i] = Aterm[i] * Bterm[i] */
    // strassen_multiply(Aterm[i], Bterm[i], P[i], size / 2);
    strassen_params[i] = malloc(sizeof(strassen_params_t));
    strassen_params[i]->aterm = Aterm[i];
    strassen_params[i]->bterm = Bterm[i];
    strassen_params[i]->pterm = P[i];
    strassen_params[i]->size = size / 2;

    pthread_create(&thread[i], NULL, strassen_multiply_wrapper, strassen_params[i]);
  }

  void *res;
  for (i = 1; i <= 7; i++)  {
    pthread_join(thread[i], &res);
    free(strassen_params[i]);
  }

  /* The resultant matrix's quadrants are
     (r s)
     (t u) which are evaluated as follows. */
  
  /* r = P5 + P4 - P2 + P6
     which is evaluated with temporaries as
     t1 = P5 + P4
     t2 = t1 - P2
     r = t2 + P6
  */
  add(P[5], P[4], t1, size / 2);
  minus(t1, P[2], t2, size / 2);
  add(t2, P[6], r, size / 2);
  
  /* s = P1 + P2; */
  add(P[1], P[2], s, size / 2);

  /* t = P3 + P4; */
  add(P[3], P[4], t, size / 2);

  /* u = P5 + P1 - P3 - P7
     which is evaluated with temporaries as
     t1 = P5 + P1
     t2 = t1 - P3
     u = t2 - P7
  */
  add(P[5], P[1], t1, size / 2);
  minus(t1, P[3], t2, size / 2);
  minus(t2, P[7], u, size / 2);

  /* merge r s t u into the resultant matrix mat_c.
     That is, mat_c = (r s)
                      (t u) */
  set_quad(mat_c, r, size, 0, size / 2, 0, size / 2);
  set_quad(mat_c, s, size, 0, size / 2, size / 2, size);
  set_quad(mat_c, t, size, size / 2, size, 0, size / 2);
  set_quad(mat_c, u, size, size / 2, size, size / 2, size);

  /* all frees */
  for (i = 0; i < 8; i++) {
    if (i != 1 && i != 4)
      free(Aterm[i]);
    if (i != 2 && i != 3)
      free(Bterm[i]);
    free(P[i]);
  }
  free(t1);
  free(t2);
  free(A_reordered);
  free(B_reordered);
  free(r);
  free(s);
  free(t);
  free(u);
}
