**Copyright (c) 2024, Cismaru-Diana Iuliana (331CA / 2023 - 2024)**

# Optimization of Matrix Multiplication
### ASC - Assignment #3

## Description
> The assignment consists of implementing three solutions for the matrix expression
> C = (A<sup>t</sup> x B + B x A) x B<sup>t</sup>, where all matrices are square matrices
> of size N x N, and matrix A is **upper triangular**.

In all 3 variants, the matrices are allocated linearly in memory. The execution was
done on the `haswell` partition, using `sbatch`, and the compilation was done using
`gcc-8.5.0`.

## 1. BLAS Implementation
  In this part, I used functions from the BLAS (Basic Linear Algebra Subprograms)
library to perform various operations with matrices, but also with vectors.

The functions used are:
  - `dgemm` - for multiplying two general matrices
  - `dtrmm` - for multiplying a matrix with a triangular matrix
  - `dcopy` - for copying the contents of one vector to another
  - `daxpy` - for adding two vectors

Firstly, I performed the multiplication of matrices **A<sup>t</sup>** si **B**, storing
the result in the variable *A<sup>t</sup>B*. I used the `dtrmm` function because matrix
**A** is triangular. Given that in the desired calculation matrix **A** appears in its
transposed form, I used the **CblasTrans** parameter. Additionally, I used the 
**CblasUpper** parameter to specify that matrix **A** is upper triangular and the
**CblasLeft** parameter to specify that matrix **A** represents the left side of the
multiplication. `dtrmm` causes the multiplication result to be stored in the right-hand
matrix, which is why I allocated a copy of matrix **B** into **A<sup>t</sup>B**, using
the `dcopy` function.

Similarly, to perform the multiplication **B x A**, I used the `dtrmm` function, but
without the need to transpose any matrix, and the **CblasRight** parameter indicates
that matrix **A** represents the right side of the multiplication. The result is stored
in **BA**.

Given that the matrices resulting from **A<sup>t</sup>B** and **BA** are stored linearly
in memory, I was able to use the `daxpy` function to add them. The result is stored in
the second term of the addition, i.e. **BA**.

Finally, I multiplied the previous result with the transposed matrix **B**, using the
`dgemm` function. The final result is stored in matrix **C**.

## 2. Unoptimized Implementation
   The unoptimized version involves a brute calculation of each operation, without using
advanced tricks to optimize memory access.

   To begin with, the transposition of the upper triangular matrix **A** is done in the
`transpose_upper_triangular` function, avoiding unnecessary operations by allocating
the matrix **A<sup>t</sup>** through `calloc` and by copying only the elements above
the diagonal from matrix **A**.

   The addition of two matrices `add_matrices` involves the sum of the corresponding
elements from the two matrices, storing the result in a new matrix.

   The multiplication of matrices is done in the classic way, using 3 nested loops,
in which the sum of the products of the elements on each line and column is calculated,
according to the formula `c[i][j] = a[i][k] * b[k][j]`, with the mention that the matrices
are allocated linearly in memory, the formula becoming `c[i * N + j] = a[i * N + k] * b[k * N + j]`. However, the 3 necessary multiplications have obtained distinct
implementations, being differentiated by the type of matrices involved. Thus, we have the following cases:
   - **A<sup>t</sup> x B** - the function `multiply_lower_triangular` is used, because,
by transposing, matrix **A<sup>t</sup>** becomes a lower triangular matrix. Thus, the
multiplication with the elements below the main diagonal of matrix **A** is avoided, by
modifying the limits of the inner loop, in which **k** will go only up to **i** inclusive
- **B x A** - the function `multiply_upper_triangular` is used, analogous to the previous
case, but **k** will go only up to **j** inclusive
- **(A<sup>t</sup> x B + B x A) x B<sup>t</sup>** - the function `multiply_with_transpose`
is used, which uses the classic version, but in which matrix **B** is transposed on the spot: thus, instead of accessing the elements of matrix **B** as `B[k][j]`, they are accessed as `B[j][k]`.

## 3. Optimized Implementation
   The idea behind all optimizations starts from the code sequence presented in laboratory
9, through which memory access is optimized through pointers. Thus, vectorial access is
avoided by dereferencing:

```
for (i = 0; i < N; i++) {
	double *orig_pa = &a[i][0];
	for (j = 0; j < N; j++) {
		double *pa = orig_pa;
		double *pb = &b[0][j];
		register double sum = 0;
		for (k = 0; k < N; k++) {
			sum += *pa * *pb;
			pa++;
			pb += N;
		}
		c[i][j] = sum;
	}
}
```
My implementation adapts the code above to matrices that are allocated linearly in memory
and maintains a pointer even for the result matrix. Since pointers are used repeatedly
in calculations from for loops, they are declared as **register**, to send a signal to
the compiler to keep them in the processor registers for much faster access. However, to
avoid exhausting the CPU registers, the declarations were made at the beginning of the
functions, not inside the loops. This change brought significant improvements in the
execution time of the program.

In addition, in most of the implemented functions, to avoid repeating the calculation
`i * N`, an auxiliary register `index` was used, to which `N` was added at each iteration
of the outer for loop.

## Memory Access
In all 3 variants, I allocated memory with `calloc`, to ensure that all elements of the
matrices are initialized to 0. I also freed the memory used for matrices with `free`, to
avoid memory leaks.

To check that there are no memory access problems, the 3 executables were run with the
following command:
```
valgrind --tool=memcheck --leak-check=full ./tema3_varianta ../input/input_valgrind > ../memory/version.memory 2>&1
```

The results of the command are available in the `memory/` directory, highlighting the fact
that there are no memory access problems or memory leaks.

## Cache Access
To analyze information related to cache access, I used the following command on all 3
executables:
```
valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes ./tema3_varianta ../input/input_valgrind > ../cache/version.cache 2>&1
```

### Results obtained with `cachegrind`
- `I refs` - Intruction references
	* tema3_blas: 248,609,027
	* tema3_neopt: 3,626,867,963
	* tema3_opt_m: 1,765,290,607
- `I1 misses` - L1 instruction cache misses
	* tema3_blas: 23,240
	* tema3_neopt: 1,672
	* tema3_opt_m: 1,643
- `LLi misses` - L2 instruction cache misses
	* tema3_blas: 3,526
	* tema3_neopt: 1,592
	* tema3_opt_m: 1,568
- `D refs` - Data references
	* tema3_blas: 94,660,230
	* tema3_neopt: 1,944,113,346
	* tema3_opt_m: 722,999,257
- `D1 misses` - L1 data cache misses
	* tema3_blas: 1,754,629
	* tema3_neopt: 52,454,225
	* tema3_opt_m: 52,454,401
- `LLd misses` - L2 data cache misses
	* tema3_blas: 140,181
	* tema3_neopt: 133,278
	* tema3_opt_m: 133,273
- `LL refs` - Last level cache references
	* tema3_blas: 1,777,869
	* tema3_neopt: 52,455,897
	* tema3_opt_m: 52,456,044
- `LL misses` - Last level cache misses
	* tema3_blas: 143,707
	* tema3_neopt: 134,870
	* tema3_opt_m: 134,841
- `Branches` - Branch instructions
	* tema3_blas: 5,876,858
	* tema3_neopt: 133,573,337
	* tema3_opt_m: 133,573,272
- `Mispredicts` - Branch mispredictions
	* tema3_blas: 67,881
	* tema3_neopt: 503,137
	* tema3_opt_m: 503,108

#### Values Analysis
> The **BLAS** variant is the most efficient in terms of instruction references and
> data accesses, obtaining the lowest values for most metrics, except for cache misses.

> In **I refs** and **D refs** it can be observed that the unoptimized version executes
> by far the most instructions and accesses the most data. On the other hand, the 
> optimized version brings significant improvements in this regard, having a number of
> instruction references and data accesses almost 3 times lower, due to the fact that
> I used pointers and registers, no longer needing unnecessary memory accesses.
> While these variants bring numbers in the billions, the BLAS version has values in
> the hundreds of millions, which demonstrates the efficiency of using BLAS functions.

> For all other metrics, the optimized version has almost the same values as the 
> unoptimized version, the winner being, again, the BLAS version, which also uses
> algorithms with few and predictable branches, thus reducing the number of branches
> and mispredictions.

## Comparative Performance Analysis
For performance comparison, I ran the 3 variants on the `haswell` partition, using
5 tests with different sizes for N: 400, 800, 1000, 1200, 1400. The results are as
follows:
![grafic](https://i.ibb.co/m4sLT55/Screenshot-from-2024-05-18-22-44-32.png)

The results obtained are the expected ones: the BLAS version is the fastest, followed
by the optimized version, and the slowest is the unoptimized version.

The execution time increases almost linearly on the BLAS version, managing to remain
below 1 second even for N = 1400. On the other hand, the increases in the execution
times of the other two variants are more abrupt, having a similar slope, reaching
almost 11 seconds on the optimized implementation and 16 seconds on the unoptimized
one. We can say that the manual optimization has significantly improved performance,
but not at the level of the BLAS implementation.

## Resources:
* [BLAS Quick Reference Guide](https://www.netlib.org/blas/blasqr.pdf)
* [Laborator 9 - Tehnici de Optimizare de Cod - Inmultirea Matricelor](https://ocw.cs.pub.ro/courses/asc/laboratoare/09)
* [Tema 3 - ASC](https://ocw.cs.pub.ro/courses/asc/teme/tema3)
* [Picture URL](https://i.ibb.co/)
* [Graphic Maker](https://my.visme.co/)