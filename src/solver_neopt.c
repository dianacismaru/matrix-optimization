/*
 * Tema 3 ASC
 * 2024 Spring
 */
#include "utils.h"
#include "helper.h"

// Multiply a normal matrix with a transposed one
double* multiply_with_transpose(int N, double *A, double *B, double *res) {	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0;
			for (int k = 0; k < N; k++) {
				// Access B's elements in a transposed manner
				sum += A[i * N + k] * B[j * N + k];
			}
			res[i * N + j] = sum;
		}
	}

	return res;
}

// Multiply a lower triangular matrix with a normal one
double* multiply_lower_triangular(int N, double *L, double *B, double *res) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0;
			// Only iterate through the elements of L that are not 0
			for (int k = 0; k <= i; k++) {
				sum += L[i * N + k] * B[k * N + j];
			}
			res[i * N + j] = sum;
		}
	}

	return res;
}

// Multiply a normal matrix with an upper triangular one
double* multiply_upper_triangular(int N, double *B, double *U, double *res) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0;
			// Only iterate through the elements of U that are not 0
			for (int k = 0; k <= j; k++) {
				sum += B[i * N + k] * U[k * N + j];
			}
			res[i * N + j] = sum;
		}
	}

	return res;
}

// Transpose an upper triangular matrix
double* transpose_upper_triangular(int N, double *A, double *At) {
	for (int i = 0; i < N; i++) {
		// Get only the elements above the diagonal, the rest are 0
		for (int j = i; j < N; j++) {
			At[j * N + i] = A[i * N + j];
		}
	}

	return At;
}

// Add two matrices
double* add_matrices(int N, double *A, double *B, double *res) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			res[i * N + j] = A[i * N + j] + B[i * N + j];
		}
	}

	return res;
}

double* my_solver(int N, double *A, double* B) {
	double *At = transpose_upper_triangular(N, A, allocate_matrix(N));
	double *AtB = multiply_lower_triangular(N, At, B, allocate_matrix(N));
	double *BA = multiply_upper_triangular(N, B, A, allocate_matrix(N));
	double *paranthesis = add_matrices(N, AtB, BA, allocate_matrix(N));
	double *result = multiply_with_transpose(N, paranthesis, B, allocate_matrix(N));
	
	free(At);
	free(AtB);
	free(BA);
	free(paranthesis);

	return result;
}
