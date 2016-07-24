#ifndef CONFIG_H
#define CONFIG_H

/* RECURSION_CUTOFF: strassen multiplication is worse when N (the order of
   square matrix) is small. So, stop recursing (dividing the submatrix) at
   a point and switch to O(N^3) algorithm which give better performance.
   The cutoff level depends on the hardware configuration and the algorithm
   needs to profiled on the hardware to find the optimal level.
   Note: Could find in one literature that 128 is a suitable cutoff level
*/
#define RECURSION_CUTOFF 128

/* BLOCK_CUTOFF: The cutoff above which blocked matrix multiplication
   is worth.
   basically, a strip (of block size) of matrix A is multiplied
   with a block (of block size) of matrix B and resultant matrix's
   strip is updated at each iteration. This method helps to minimise cache
   misses due column wise access of matrix B in traditional ijk algorithm
   the block size should be chosen such that the block of B can fit in the
   cache (depends on cache configuration). BLOCK_CUTOFF as 100 and
   BLOCK_SIZE of 25 is recommended.
   Need profiling on particular processor to choose the best numbers
*/

#define BLOCK_CUTOFF 64
#define BLOCK_SIZE 32

/* the type of the matrix elements. change the type here to allow the
   algorithm to work on matrices of different types */
typedef long int ele_type;

#endif // CONFIG_H
