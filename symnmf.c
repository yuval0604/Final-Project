#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

#define EPSILON 1e-4
#define MAX_ITER 300
#define BETA 0.5

/* Memory management: free 2D matrix */
void free_matrix(double** matrix, int rows) {
    int i;
    if (matrix) {
        for (i = 0; i < rows; i++) {
            free(matrix[i]);
        }
        free(matrix);
    }
}

/* Safe memory allocation with error handling */
double** allocate_matrix(int rows, int cols) {
    double** matrix;
    int i;

    matrix = (double**)malloc(rows * sizeof(double*));
    if (!matrix) {
        fprintf(stderr, "An Error Has Occurred\n");
        exit(1);
    }

    for (i = 0; i < rows; i++) {
        matrix[i] = (double*)calloc(cols, sizeof(double));
        if (!matrix[i]) {
            /* Free previously allocated rows */
            free_matrix(matrix, i);
            fprintf(stderr, "An Error Has Occurred\n");
            exit(1);
        }
    }
    return matrix;
}

/* Safe value handling to prevent division by zero */
double safe_divide(double numerator, double denominator) {
    if (fabs(denominator) < EPSILON) {
        return numerator / EPSILON;
    }
    return numerator / denominator;
}

/* Euclidean distance calculation */
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
    double** A = allocate_matrix(n, n);
    int i, j;

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
    double** D = allocate_matrix(n, n);
    int i, j;

    for (i = 0; i < n; i++) {
        D[i][i] = 0.0;
        for (j = 0; j < n; j++) {
            D[i][i] += similarity[i][j];
        }
    }

    return D;
}

/* Normalized Similarity Matrix Calculation */
double** calculate_normalized_similarity(double** A, double** D, int n) {
    double** W = allocate_matrix(n, n);
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double sqrt_deg_i = sqrt(D[i][i]);
            double sqrt_deg_j = sqrt(D[j][j]);
            
            W[i][j] = (sqrt_deg_i < EPSILON || sqrt_deg_j < EPSILON) 
                ? 0.0 
                : A[i][j] / (sqrt_deg_i * sqrt_deg_j);
        }
    }

    return W;
}

/* SymNMF H Optimization */
double** symnmf_optimization(double** W, int n, int k, double** initial_H) {
    double** H = initial_H;
    double** H_prev = allocate_matrix(n, k);
    int iter, i, j;
    double frobenius_norm;

    /* Copy initial H */
    for (i = 0; i < n; i++) {
        memcpy(H_prev[i], H[i], k * sizeof(double));
    }

    /* Iterative update */
    for (iter = 0; iter < MAX_ITER; iter++) {
        /* Update H */
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                double numerator = 0.0, denominator = 0.0;
                int l, m;

                /* Calculate components for update */
                for (l = 0; l < n; l++) {
                    double temp = 0.0;
                    for (m = 0; m < k; m++) {
                        temp += H[l][m] * H[l][m];
                    }
                    
                    numerator += W[i][l] * H[l][j];
                    denominator += temp * H[i][j];
                }

                H[i][j] *= (1.0 - BETA + BETA * safe_divide(numerator, denominator));
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

    /* Free temporary matrix */
    free_matrix(H_prev, n);

    return H;
}
