/*
 * Tema 3 ASC
 * 2024 Spring
 */
#include "utils.h"
#include "helper.h"

// Multiply a normal matrix with a transposed one
double* multiply_with_transpose(int N, double *A, double *B, double *res) {
	register int index = 0;
	register double sum;
	register double *res_ptr, *orig_pa, *pa, *pb;

	for (int i = 0; i < N; i++, index += N) {
		res_ptr = &res[index];
		orig_pa = &A[index];
		pb = &B[0];

		for (int j = 0; j < N; j++) {
			sum = 0.0;
			pa = orig_pa;

			for (int k = 0; k < N; k++) {
				sum += (*pa) * (*pb);
				pa++;
				pb++;
			}

			res_ptr[j] = sum;
		}
	}

	return res;
}

// Multiply a lower triangular matrix with a normal one
double* multiply_lower_triangular(int N, double *L, double *B, double *res) {
	register int index = 0;
	register double sum;
	register double *res_ptr, *orig_pa, *pa, *pb;

	for (int i = 0; i < N; i++, index += N) {
		orig_pa = &L[index];
		res_ptr = &res[index];

		for (int j = 0; j < N; j++) {
			pa = orig_pa;
			pb = &B[j];
			sum = 0.0;

			for (int k = 0; k <= i; k++) {
				sum += (*pa) * (*pb);
				pa++;
				pb += N;
			}
			res_ptr[j] = sum;
		}
	}

	return res;
}

// Multiply a normal matrix with an upper triangular one
double* multiply_upper_triangular(int N, double *B, double *U, double *res) {
	register int index = 0;
	register double sum;
	register double *res_ptr, *orig_pa, *pa, *pb;

	for (int i = 0; i < N; i++, index += N) {
		orig_pa = &B[index];
		res_ptr = &res[index];

		for (int j = 0; j < N; j++) {
			pa = orig_pa;
			pb = &U[j];
			sum = 0.0;

			for (int k = 0; k <= j; k++) {
				sum += (*pa) * (*pb);
				pa++;
				pb += N;
			}
			res_ptr[j] = sum;
		}
	}

	return res;
}

// Transpose an upper triangular matrix
double* transpose_upper_triangular(int N, double *A, double *At) {
	register int index = 0;
	register double *pa, *pAt;

	for (int i = 0; i < N; i++, index += N) {
		pa = &At[index + i];
		pAt = &A[index];

		for (int j = i; j < N; j++) {
			*pa = pAt[j];
			pa += N;
		}
	}

	return At;
}

// Add two matrices
double* add_matrices(int N, double *A, double *B, double *res) {
	register int index = 0;
	register double *res_ptr, *pAt, *pb;

	for (int i = 0; i < N; i++, index += N) {
		res_ptr = &res[index];
		pAt = &A[index];
		pb = &B[index];

		for (int j = 0; j < N; j++) {
			res_ptr[j] = pAt[j] + pb[j]; 
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
