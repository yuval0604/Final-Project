#ifndef SYMNMF_H
#define SYMNMF_H

void print_error_and_exit(void);

// free matrix memory
void free_matrix(double** matrix, int rows);


// Allocate matrix memory
double** allocate_matrix(int rows, int cols);

static double** marix_mul(double** mat1, int n, double** mat2, int k);


static double** compute_HtH(double** H, int n, int k);



/*
calculate squared euclidean distance 
input: two points, dimention of the points
output: squared euclidean distance between point1, point2
*/
double euclidean_dist_squared(double* point1, double* point2, int dim);


/*
calculate similarity matrix
input: set of points, size of matrix, dim of points
output: similarity matrix
*/
double** calculate_similarity_matrix(double** data, int n, int dim) ;


/*
calculate diagonal matrix
input: similarity matrix, size of matrix
output: diagonal matrix
*/
double** calculate_diagonal_degree_matrix(double** similarity, int n);


/*
calculate normalized similarity matrix
input: similarity matrix, diagonal matrix, size of matrix
output: normalized similarity matrix
*/
double** calculate_normalized_similarity(double** A, double** D, int n);


/*
calculate H
input: normalized similarity matrix, size of matrix W, size of H, initial H
output: H
*/
double** symnmf_optimization(double** W, int n, int k, double** initial_H);


#endif

