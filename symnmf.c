#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

#define EPSILON 1e-4
#define MAX_ITER 300
#define BETA 0.5
#define MAX_LINE_LEN 1024


void print_error_and_exit(void) {
    printf("An Error Has Occurred\n");
    exit(1);
}

void free_matrix(double** matrix, int rows) {
    int i;
    if (matrix) {
        for (i = 0; i < rows; i++) {
            free(matrix[i]);
        }
        free(matrix);
    }
}

double** allocate_matrix(int rows, int cols) {
    double** matrix;
    int i;

    matrix = (double**)malloc(rows * sizeof(double*));
    if (!matrix) {
        print_error_and_exit();
    }
    for (i = 0; i < rows; i++) {
        matrix[i] = (double*)calloc(cols, sizeof(double));
        if (!matrix[i]) {
            free_matrix(matrix, i);
            print_error_and_exit();
        }
    }
    return matrix;
}

static double** matrix_mul(double** A, int rows, int inner, double** B, int cols) {
    double** res = allocate_matrix(rows, cols);
    int i, j, k;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            double sum = 0.0;
            for (k = 0; k < inner; k++) {
                sum += A[i][k] * B[k][j];
            }
            res[i][j] = sum;
        }
    }
    return res;
}




static double** compute_HtH(double** H, int n, int k) {
    double** res = allocate_matrix(k, k);
    int i, j, l;
    for (i = 0; i < k; i++) {
        for (j = 0; j < k; j++) {
            double sum = 0.0;
            for (l = 0; l < n; l++) {
                sum += H[l][i] * H[l][j];
            }
            res[i][j] = sum;
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

    for (i = 0; i < n; i++) {
        memcpy(H_prev[i], H[i], k * sizeof(double));
    }

    for (iter = 0; iter < MAX_ITER; iter++) {
        
        double** WH = matrix_mul(W, n, n, H, k);
        double** HtH = compute_HtH(H, n, k);
        double** HHtH = matrix_mul(H, n, k, HtH, k);


        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                double hij = H[i][j];
                double numerator   = WH[i][j];
                double denominator = HHtH[i][j];
                double frac;
                double updated;

                if (fabs(denominator) < EPSILON) {
                    denominator = EPSILON;  
                }

                {
                    frac = numerator / denominator;
                    updated = hij * ((1.0 - BETA) + BETA * frac);

                    if (updated < 0.0) {
                        updated = 0.0;
                    }

                    H[i][j] = updated;
                }
            }
        }

        free_matrix(WH, n);
        free_matrix(HtH, k);
        free_matrix(HHtH, n);

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

        for (i = 0; i < n; i++) {
            memcpy(H_prev[i], H[i], k * sizeof(double));
        }
    }

    free_matrix(H_prev, n);
    return H;
}


int main(int argc, char *argv[]) {
    char *goal, *filename, *token;
    FILE *fp;
    int n=0;
    int d=0;
    int i=0;
    int j = 0;
    int r, c;
    char line[MAX_LINE_LEN];
    double **data, **A, **D; 
    double **result_matrix=NULL;
    
    if (argc != 3) {
        printf("An Error Has Occurred\n");
        return 1;
    }
    
    goal = argv[1];
    filename = argv[2];

    if (strcmp(goal, "sym") != 0 && strcmp(goal, "ddg") != 0 && strcmp(goal, "norm") != 0) {
        printf("An Error Has Occurred\n");
        return 1;
    }

    fp = fopen(filename, "r");
    if (!fp) {
        printf("An Error Has Occurred\n");
        return 1;
    }

   
    if (fgets(line, sizeof(line), fp) == NULL) {
        printf("An Error Has Occurred\n");
        fclose(fp);
        return 1;
    }
   
    n = 1; 
    token = strtok(line, ",");
    while (token != NULL) {
        d++;
        token = strtok(NULL, ",");
    }
    while (fgets(line, sizeof(line), fp) != NULL) {
        n++;
    }
    fclose(fp);

    data = allocate_matrix(n, d);
    if (!data) {
        printf("An Error Has Occurred\n");
        return 1;
    }

    fp = fopen(filename, "r");
    if (!fp) {
        printf("An Error Has Occurred\n");
        free_matrix(data, n);
        return 1;
    }
    
    while (fgets(line, sizeof(line), fp) != NULL && i < n) {
        line[strcspn(line, "\n")] = '\0';
        token = strtok(line, ",");
        while (token != NULL && j < d) {
            data[i][j] = atof(token);
            j++;
            token = strtok(NULL, ",");
        }
        i++;
    }
    fclose(fp);

    if (strcmp(goal, "sym") == 0) {
        result_matrix = calculate_similarity_matrix(data, n, d);
    }
    else if (strcmp(goal, "ddg") == 0) {
        A = calculate_similarity_matrix(data, n, d);
        result_matrix = calculate_diagonal_degree_matrix(A, n);
        free_matrix(A, n);
    }
    else if (strcmp(goal, "norm") == 0) {
        A = calculate_similarity_matrix(data, n, d);
        D = calculate_diagonal_degree_matrix(A, n);
        result_matrix = calculate_normalized_similarity(A, D, n);
        free_matrix(A, n);
        free_matrix(D, n);
    }
    free_matrix(data, n);

    if (!result_matrix) {
        printf("An Error Has Occurred\n");
        return 1;
    }

    for (r = 0; r < n; r++) {
        for (c = 0; c < n; c++) {
            printf("%.4f", result_matrix[r][c]);
            if (c < n - 1)
                printf(",");
        }
        printf("\n");
    }
    free_matrix(result_matrix, n);
    return 0;
}
