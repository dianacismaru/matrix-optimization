/*
 * Tema 3 ASC
 * 2024 Spring
 */
#include "utils.h"
#include "helper.h"
#include <string.h>
#include <cblas.h>

double* my_solver(int N, double *A, double *B) {
	// Create a copy of B in AtB, to store the result of At * B
	double *AtB = allocate_matrix(N);
	cblas_dcopy(N * N, B, 1, AtB, 1);
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
				N, N, 1, A, N, AtB, N);

	// Create a copy of B in BA, to store the result of B * A
	double *BA = allocate_matrix(N);
	cblas_dcopy(N * N, B, 1, BA, 1);
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
				N, N, 1, A, N, BA, N);

	// Add the vectors AtB and BA and store the result in BA
	cblas_daxpy(N * N, 1, AtB, 1, BA, 1);

	// Multiply BA (currently containing At x B + B x A) with Bt
	double *C = allocate_matrix(N);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1, BA, 
				N, B, N, 0, C, N);

	// Free the allocated memory
	free(AtB);
	free(BA);

	return C;
}
