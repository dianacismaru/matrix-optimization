#ifndef HELPER_H
#define HELPER_H

// Allocate memory for a matrix of size N x N
double* allocate_matrix(int N) {
	double *A = (double*) calloc(N * N, sizeof(double));

	// Safety check
	if (!A) {
		printf("Memory allocation failed\n");
		exit(1);
	}

	return A;
}

#endif // HELPER_H