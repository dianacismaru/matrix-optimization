/*
 * Tema 3 ASC
 * 2024 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */

double* allocate_matrix(int N) {
	double *A = (double*) calloc(N * N, sizeof(double));
	if (!A) {
		printf("Memory allocation failed\n");
		exit(1);
	}
	return A;
}

double* multiply_with_transpose(int N, double *A, double *B, double *res) {	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0;
			for (int k = 0; k < N; k++) {
				sum += A[i * N + k] * B[j * N + k];
			}
			res[i * N + j] = sum;
		}
	}

	return res;
}

double* multiply_lower_triangular(int N, double *L, double *B, double *res) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0;
			for (int k = 0; k <= i; k++) {
				sum += L[i * N + k] * B[k * N + j];
			}
			res[i * N + j] = sum;
		}
	}

	return res;
}

double* multiply_upper_triangular(int N, double *B, double *U, double *res) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0;
			for (int k = 0; k <= j; k++) {
				sum += B[i * N + k] * U[k * N + j];
			}
			res[i * N + j] = sum;
		}
	}

	return res;
}

double* transpose_sup_triangular(int N, double *A, double *At) {
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			At[j * N + i] = A[i * N + j];
		}
	}

	return At;
}

double* add_matrices(int N, double *A, double *B, double *res) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			res[i * N + j] = A[i * N + j] + B[i * N + j];
		}
	}

	return res;
}

void print_matrix(int N, double *A) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%f ", A[i * N + j]);
		}
		printf("\n");
	}
}
/*
 * C = (At x B + B x A) x Bt
 */
double* my_solver(int N, double *A, double* B) {
	printf("NEOPT SOLVER\n");

	double *At = transpose_sup_triangular(N, A, allocate_matrix(N));
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
