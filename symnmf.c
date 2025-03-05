#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

#define MAX_LINE_LENGTH 10000
#define EPSILON 1e-4
#define MAX_ITER 300
#define BETA 0.5

/* Function to calculate Euclidean distance squared */
double euclidean_dist_squared(double* point1, double* point2, int dim) {
    double dist = 0.0;
    int i;
    for (i = 0; i < dim; i++) {
        double diff = point1[i] - point2[i];
        dist += diff * diff;
    }
    return dist;
}

/* Similarity Matrix Calculation */
double** calculate_similarity_matrix(double** data, int n, int dim) {
    double** A;
    int i, j;

    /* Allocate memory for similarity matrix */
    A = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        A[i] = (double*)malloc(n * sizeof(double));
    }

    /* Calculate similarity */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                A[i][j] = 0.0;
            } else {
                double dist = euclidean_dist_squared(data[i], data[j], dim);
                A[i][j] = exp(-dist);
            }
        }
    }

    return A;
}

/* Diagonal Degree Matrix Calculation */
double** calculate_diagonal_degree_matrix(double** similarity, int n) {
    double** D;
    int i, j;

    /* Allocate memory for diagonal degree matrix */
    D = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        D[i] = (double*)calloc(n, sizeof(double));
        D[i][i] = 0.0;
        for (j = 0; j < n; j++) {
            D[i][i] += similarity[i][j];
        }
    }

    return D;
}

/* Normalized Similarity Matrix Calculation */
double** calculate_normalized_similarity(double** A, double** D, int n) {
    double** W;
    int i, j;

    /* Allocate memory for normalized matrix */
    W = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        W[i] = (double*)malloc(n * sizeof(double));
    }

    /* Calculate D^(-1/2) * A * D^(-1/2) */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            W[i][j] = A[i][j] / sqrt(D[i][i] * D[j][j]);
        }
    }

    return W;
}

/* SymNMF H Optimization */
double** symnmf_optimization(double** W, int n, int k, double** initial_H) {
    double** H = initial_H;
    double** H_prev;
    int iter, i, j;
    double frobenius_norm;

    /* Allocate memory for previous H */
    H_prev = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        H_prev[i] = (double*)malloc(k * sizeof(double));
        memcpy(H_prev[i], H[i], k * sizeof(double));
    }

    /* Iterative update */
    for (iter = 0; iter < MAX_ITER; iter++) {
        /* Update H */
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                double numerator = 0.0, denominator = 0.0;
                int l;

                /* Calculate components for update */
                for (l = 0; l < n; l++) {
                    double temp = 0.0;
                    int m;
                    for (m = 0; m < k; m++) {
                        temp += H[l][m] * H[l][m];
                    }
                    numerator += W[i][l] * H[l][j];
                    denominator += temp * H[i][j];
                }

                H[i][j] *= (1.0 - BETA + BETA * numerator / denominator);
            }
        }

        /* Check convergence */
        frobenius_norm = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                double diff = H[i][j] - H_prev[i][j];
                frobenius_norm += diff * diff;
            }
        }

        if (sqrt(frobenius_norm) < EPSILON) {
            break;
        }

        /* Update H_prev */
        for (i = 0; i < n; i++) {
            memcpy(H_prev[i], H[i], k * sizeof(double));
        }
    }

    /* Free H_prev */
    for (i = 0; i < n; i++) {
        free(H_prev[i]);
    }
    free(H_prev);

    return H;
}

int main(int argc, char* argv[]) {
    /* Your implementation for command-line argument processing */
    return 0;
}
