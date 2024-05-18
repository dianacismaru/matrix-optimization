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

/*
for(i = 0; i < N; i++){
  double *orig_pa = &a[i][0];
  for(j = 0; j < N; j++){
    double *pa = orig_pa;
    double *pb = &b[0][j];
    register double suma = 0;
    for(k = 0; k < N; k++){
      suma += *pa * *pb;
      pa++;
      pb += N;
    }
    c[i][j] = suma;
  }
}
*/

// double* multiply_with_transpose(int N, double *A, double *B, double *res) {	
// 	for (int i = 0; i < N; i++) {
// 		for (int j = 0; j < N; j++) {
// 			register double sum = 0;
// 			for (int k = 0; k < N; k++) {
// 				sum += A[i * N + k] * B[j * N + k];
// 			}
// 			res[i * N + j] = sum;
// 		}
// 	}

// 	return res;
// }

double* multiply_with_transpose(int N, double *A, double *B, double *res) {
	for (int i = 0; i < N; i++) {
		register double *res_ptr = &res[i * N];
		register double *orig_pa = &A[i * N];
		for (int j = 0; j < N; j++) {
			register double sum = 0;
			register double *pa = orig_pa;
			register double *pb = &B[j * N];
			for (int k = 0; k < N; k++) {
				sum += *pa * *pb;
				pa++;
				pb++;
			}
			res_ptr[j] = sum;
		}
	}

	return res;
}

double* multiply_lower_triangular(int N, double *L, double *B, double *res) {
	for (int i = 0; i < N; i++) {
		register double *orig_pa = &L[i * N];
		register double *res_ptr = &res[i * N];

		for (int j = 0; j < N; j++) {
			register double *pa = orig_pa;
			register double *pb = &B[j];
			register double sum = 0.0;
			for (int k = 0; k <= i; k++) {
				sum += *pa * *pb;
				pa++;
				pb += N;
			}
			res_ptr[j] = sum;
		}
	}

	return res;
}

double* multiply_upper_triangular(int N, double *B, double *U, double *res) {
	for (int i = 0; i < N; i++) {
		register double *orig_pa = &B[i * N];
		register double *res_ptr = &res[i * N];
		for (int j = 0; j < N; j++) {
			register double *pa = orig_pa;
			register double *pb = &U[j];
			register double sum = 0.0;
			for (int k = 0; k <= j; k++) {
				sum += *pa * *pb;
				pa++;
				pb += N;
			}
			res_ptr[j] = sum;
		}
	}

	return res;
}

double* transpose_sup_triangular(int N, double *A, double *At) {
	for (int i = 0; i < N; i++) {
		register int index = i * N;
		register double* At_ptr = &At[index + i];
		register double* A_ptr = &A[index];
		for (int j = i; j < N; j++) {
			// At[j * N + i] = A[i * N + j];
			*At_ptr = A_ptr[j];
			At_ptr += N;
		}
	}

	return At;
}


double* add_matrices(int N, double *A, double *B, double *res) {
	for (int i = 0; i < N; i++) {
		register int index = i * N;
		double *res_ptr = &res[index];
		double *A_ptr = &A[index];
		double *B_ptr = &B[index];
		for (int j = 0; j < N; j++) {
			res_ptr[j] = A_ptr[j] + B_ptr[j]; 
		}
	}

	return res;
}

/*
 * C = (At x B + B x A) x Bt
 */
double* my_solver(int N, double *A, double* B) {
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
