#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

#define EPSILON 1e-4
#define MAX_ITER 300
#define BETA 0.5

/* --- פונקציות ניהול זיכרון --- */
double** allocate_matrix(int rows, int cols);
void free_matrix(double** mat, int rows);
void print_matrix(double** mat, int rows, int cols);

/* --- חישובי מטריצות --- */
double euclidean_dist_squared(double* point1, double* point2, int dim);
double** calculate_similarity_matrix(double** data, int n, int dim);
double** calculate_diagonal_degree_matrix(double** similarity, int n);
double** calculate_normalized_similarity(double** A, double** D, int n);
double** symnmf_optimization(double** W, int n, int k, double** initial_H);

/* --- ניהול זיכרון --- */
double** allocate_matrix(int rows, int cols) {
    int i;
    double** mat = (double**)malloc(rows * sizeof(double*));
    if (!mat) {
        fprintf(stderr, "An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < rows; i++) {
        mat[i] = (double*)calloc(cols, sizeof(double));
        if (!mat[i]) {
            free_matrix(mat, i);
            fprintf(stderr, "An Error Has Occurred\n");
            exit(1);
        }
    }
    return mat;
}

void free_matrix(double** mat, int rows) {
    int i;
    if (mat) {
        for (i = 0; i < rows; i++) {
            if (mat[i]) {
                free(mat[i]);
            }
        }
        free(mat);
    }
}

void print_matrix(double** mat, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f ", mat[i][j]);
        }
        printf("\n");
    }
}

/* --- חישובי מטריצות --- */
double euclidean_dist_squared(double* point1, double* point2, int dim) {
    double dist = 0.0;
    int i;
    for (i = 0; i < dim; i++) {
        double diff = point1[i] - point2[i];
        dist += diff * diff;
    }
    return dist;
}

double** calculate_similarity_matrix(double** data, int n, int dim) {
    double** A = allocate_matrix(n, n);
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = (i == j) ? 0.0 : exp(-euclidean_dist_squared(data[i], data[j], dim));
        }
    }
    return A;
}

double** calculate_diagonal_degree_matrix(double** similarity, int n) {
    double** D = allocate_matrix(n, n);
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            D[i][i] += similarity[i][j];
        }
    }
    return D;
}

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

/* --- SymNMF אופטימיזציה --- */
double** symnmf_optimization(double** W, int n, int k, double** initial_H) {
    double** H = allocate_matrix(n, k);
    double** H_prev = allocate_matrix(n, k);
    int iter, i, j, l, m;
    double frobenius_norm;

    /* אתחול H */
    for (i = 0; i < n; i++) {
        memcpy(H[i], initial_H[i], k * sizeof(double));
        memcpy(H_prev[i], H[i], k * sizeof(double));  // מעתיק את H ולא את initial_H
    }

    for (iter = 0; iter < MAX_ITER; iter++) {
        /* עדכון H */
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                double numerator = 0.0, denominator = 0.0;

                for (l = 0; l < n; l++) {
                    numerator += W[i][l] * H[l][j];  // חישוב המונה (W * H)

                    // חישוב המכנה (H * H^T * H)
                    double hht_h = 0.0;
                    for (m = 0; m < k; m++) {
                        hht_h += H[l][m] * H[l][m];  // נוסחה מתוקנת
                    }
                    denominator += H[i][l] * hht_h;
                }

                H[i][j] *= (1.0 - BETA + BETA * (numerator / (denominator + EPSILON)));
            }
        }

        /* בדיקת התכנסות */
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

        /* עדכון H_prev */
        for (i = 0; i < n; i++) {
            memcpy(H_prev[i], H[i], k * sizeof(double));
        }
    }

    /* שחרור זיכרון */
    free_matrix(H_prev, n);
    return H;
}
