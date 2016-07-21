To build:
make

To run:
./strass

Testing:
running ./strass tests multiplication of matrices are various sizes
and compares the output with the expected output to declare PASS or FAIL.

Features:
1. At each recursion level of strassen multiplication, the matrix is reordered
   before splitting into submatrix. This helps to keep the submatrices contiguous
   and thus have better spatial locality

2. After a certain recusion cutoff level is reached, blocked multiplication is
   used to bypass strassen multiplication overhead (for lower order matrices)
   and at the same time take advantage of temporal locality of the blocks to
   reduce cache misses.

3. For lower sizes of matrix, multiplication is switched to traditional
   ijk multiplication as its not worth the overhead of blocked multiplication.

Limitations:
1. Only square matrices of size in the order of power of 2 are supported.
   Other sizes can be supported by the method of dynamic peeling, but a TODO.

2. Strassen multiplication is implemented using recursion (upto recusion cutoff
    level) with a number of temporary buffers created at each level. More careful
    examination could help to reduce the amout of temporary used.
