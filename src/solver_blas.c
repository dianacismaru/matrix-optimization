/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <string.h>
#include <cblas.h>

double* allocate_matrix(int N) {
	double *A = (double*) calloc(N * N, sizeof(double));
	if (!A) {
		printf("Memory allocation failed\n");
		exit(1);
	}
	return A;
}

double* add_matrices(int N, double *A, double *B, double *res) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			res[i * N + j] = A[i * N + j] + B[i * N + j];
		}
	}

	return res;
}

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
	// DGEMM - C = alpha * A * B + beta * C, matrice normale
	// DTRMM - B = alpha * A * B, matrice triunghiulare
	// At * B, A superior triangular

	// Copiez in AtB matricea B, pentru a putea stoca si rezultatul tot aici
	double *AtB = allocate_matrix(N);
	memcpy(AtB, B, N * N * sizeof(double));
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
				N, N, 1, A, N, AtB, N);

	double *BA = allocate_matrix(N);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, B, N, A,
				N, 0, BA, N);

	double *paranthesis = add_matrices(N, AtB, BA, allocate_matrix(N));

	double *result = allocate_matrix(N);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1, paranthesis, 
				N, B, N, 0, result, N);

	free(AtB);
	free(BA);
	free(paranthesis);

	return result;
}
