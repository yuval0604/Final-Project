#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

#define EPSILON 1e-4
#define MAX_ITER 300
#define BETA 0.5

void print_error_and_exit(void) {
    printf("An Error Has Occurred\n");
    exit(1);
}

// free matrix memory
void free_matrix(double** matrix, int rows) {
    int i;
    if (matrix) {
        for (i = 0; i < rows; i++) {
            free(matrix[i]);
        }
        free(matrix);
    }
}

// Allocate matrix memory
double** allocate_matrix(int rows, int cols) {
    double** matrix;
    int i;

    matrix = (double**)malloc(rows * sizeof(double*));
    if (!matrix) {
        print_error_and_exit(void);
    }
    for (i = 0; i < rows; i++) {
        matrix[i] = (double*)calloc(cols, sizeof(double));
        if (!matrix[i]) {
            free_matrix(matrix, i);
            print_error_and_exit(void);
        }
    }
    return matrix;
}

static double** marix_mul(double** mat1, int n, double** mat2, int k) {
    double** res = allocate_matrix(n, k);
    double sum_res = 0.0;
    int row, col, l;
    for (row = 0; row < n; row++) {
        for (col = 0; col < k; col++) {
            for (l = 0; l < n; l++) {
                sum_res += mat1[row][l] * mat2[l][col];
                }
        res[row][col] = sum_res;
        }
    }
    return res;
}

static double** compute_HtH(double** H, int n, int k) {
    double** res = allocate_matrix(k, k);
    double sum_res = 0.0;
    int row, col, rowH;
    for (row = 0; row < k; row++) {
        for (col = 0; col < k; col++) {
            for (rowH = 0; rowH < n; rowH++) {
                sum_res += H[rowH][r] * H[rowH][c];
            }
            res[row][col] = sum_res;
        }
    }
    return res;
}


/*
calculate squared euclidean distance 
input: two points, dimention of the points
output: squared euclidean distance between point1, point2
*/
double euclidean_dist_squared(double* point1, double* point2, int dim) {
    double dist = 0.0;
    int i;
    for (i = 0; i < dim; i++) {
        double diff = point1[i] - point2[i];
        dist += diff * diff;
    }
    return dist;
}

/*
calculate similarity matrix
input: set of points, size of matrix, dim of points
output: similarity matrix
*/
double** calculate_similarity_matrix(double** data, int n, int dim) {
    double** A = allocate_matrix(n, n);
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                A[i][j] = 0.0;
            }
            else {
                double dist = euclidean_dist_squared(data[i], data[j], dim);
                A[i][j] = exp(-dist);
            }
        }
    }
    return A;
}

/*
calculate diagonal matrix
input: similarity matrix, size of matrix
output: diagonal matrix
*/
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

/*
calculate normalized similarity matrix
input: similarity matrix, diagonal matrix, size of matrix
output: normalized similarity matrix
*/
double** calculate_normalized_similarity(double** A, double** D, int n) {
    double** W = allocate_matrix(n, n);
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double sqrt_deg_i = sqrt(D[i][i]);
            double sqrt_deg_j = sqrt(D[j][j]);
            if (sqrt_deg_i < EPSILON || sqrt_deg_j < EPSILON) {
                W[i][j] = 0.0;
            } else {
                W[i][j] = A[i][j] / (sqrt_deg_i * sqrt_deg_j);
            }
        }
    }
    return W;
}

/*
calculate H
input: normalized similarity matrix, size of matrix W, size of H, initial H
output: H
*/
double** symnmf_optimization(double** W, int n, int k, double** initial_H) {
    double** H = initial_H;
    double** H_prev = allocate_matrix(n, k);
    int iter, i, j;
    double norm;

    // Copy initial H to H_prev
    for (i = 0; i < n; i++) {
        memcpy(H_prev[i], H[i], k * sizeof(double));
    }

    for (iter = 0; iter < MAX_ITER; iter++) {
        
        double** WH = marix_mul(W, n, H, k)
        double** HtH = compute_HtH(H, n, k)
        double** HHth = marix_mul(HHt, n, H, k)

        // Update H
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                double hij = H[i][j];
                double numerator   = WH[i][j];
                double denominator = HHtH[i][j];

                // avoid division by zero
                if (abs(denominator) < EPSILON) {
                    denominator = EPSILON;  
                }

                {
                    double frac = numerator / denominator;
                    double updated = hij * ((1.0 - BETA) + BETA * frac);

                    if (updated < 0.0) {
                        updated = 0.0;
                    }

                    H[i][j] = updated;
                }
            }
        }

        free_matrix(WH,   n);
        free_matrix(HtH,  k);
        free_matrix(HHtH, n);

        //Check convergence
        norm = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                double diff = H[i][j] - H_prev[i][j];
                norm += diff * diff;
            }
        }
        if (sqrt(norm) < EPSILON) {
            break;
        }

        // copy H to H_prev
        for (i = 0; i < n; i++) {
            memcpy(H_prev[i], H[i], k * sizeof(double));
        }
    }

    free_matrix(H_prev, n);
    return H;
}




