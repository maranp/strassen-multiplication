#ifndef STRASS_H
#define STRASS_H
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <assert.h>
#include "config.h"

/* temporary buffer for matrices are allocated with posix_memalign.
   Keeping the ALIGN parameter to cache line size of the processor
   could help to mimise cache misses. 64 is cache line size on the 
   processor tested.
*/
#define ALIGN 64

/* convenience macro to access 1d buffer as 2d matrix elements */
#define MEM(a, n, i, j)  a[(i) * (n) + (j)]


/* THE function that does the core work - strassen multiplication */
void strassen_multiply(const ele_type * A,
		       const ele_type * B,
		       ele_type * C, size_t size);

#endif // STRASS_H
