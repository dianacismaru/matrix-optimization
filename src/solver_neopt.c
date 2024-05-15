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

double* multiply(int N, double *A, double *B, double *res) {	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0;
			for (int k = 0; k < N; k++) {
				sum += A[i * N + k] * B[k * N + j];
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

double* transpose_matrix(int N, double *A) {
	double *At = allocate_matrix(N);
	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			At[i * N + j] = A[j * N + i];
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

	double *At = transpose_matrix(N, A);
	double *AtB = allocate_matrix(N); 
	AtB = multiply_lower_triangular(N, At, B, AtB);
	// AtB = multiply(N, At, B, AtB);
	
	double *BA = allocate_matrix(N);
	// BA = multiply(N, B, A, BA);
	BA = multiply_upper_triangular(N, B, A, BA);

	double *paranthesis = allocate_matrix(N);
	paranthesis = add_matrices(N, AtB, BA, paranthesis);

	double *Bt = transpose_matrix(N, B);
	double *result = allocate_matrix(N);
	result = multiply(N, paranthesis, Bt, result);

	// print_matrix(N, result);

	return result;
}
