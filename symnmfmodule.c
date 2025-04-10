# define PY_SSIZE_T_CLEAN
#include "py_compat.h"
#include "symnmf.h"

/*
Convert a NumPy 2D-array to a C matrix
input: Numpy 2D array, 2 pointers to store the number of rows and columns.
output: Matrix (or NULL if memory allocation failed)
*/
static double** numpy2D_to_C_matrix(PyArrayObject* array, int* rows, int* cols){
    int i, j;
    double** matrix;

    *rows = (int)PyArray_DIM(array, 0);
    *cols = (int)PyArray_DIM(array, 1);

    matrix = allocate_matrix(*rows, *cols);
    if (!matrix){return NULL;}

    for(i = 0; i < *rows; i++){
        for(j = 0; j < *cols; j++){
            matrix[i][j] = *(double*)PyArray_GETPTR2(array, i, j);
        }
    }

    return matrix;
}

/*
Convert a C matrix back to a numpy 2D-array 
input: C Matrix, num of rows, num of columns.
output: A NumPy array (or NULL if creation failed)
*/
static PyObject* C_matrix_to_numpy2D(double** matrix, int rows, int cols){
    int i, j;
    npy_intp dims[2] = {rows, cols};
    PyObject* np_array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    
    if (!np_array) {return NULL;}
    
    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            *(double*)PyArray_GETPTR2((PyArrayObject*)np_array, i, j) = matrix[i][j];
        }
    }
    
    free_matrix(matrix, rows);
    return np_array;
}

/*
Computes the similarity matrix A
input: Numpy 2D array
output: Similarity matrix A (as a 2D numpy array)
*/
static PyObject* py_sym(PyObject* self, PyObject* args){
    PyArrayObject* data_array;
    double **data, **similarity;
    int rows, cols;
    PyObject* result;
    (void)self;
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &data_array)) {
        return NULL;
    }
    if (PyArray_NDIM(data_array) != 2) {return NULL;}

    data = numpy2D_to_C_matrix(data_array, &rows, &cols); 
    if (!data) {return NULL;}

    similarity = calculate_similarity_matrix(data, rows, cols);
    free_matrix(data, rows);

    if (!similarity) {return NULL;}
    result = C_matrix_to_numpy2D(similarity, rows, rows);
    return result;
}

/*
Computes the diagonal degree matrix D
input: Numpy 2D array
output: Diagonal degree matrix D (as a 2D numpy array)
*/
static PyObject* py_ddg(PyObject* self, PyObject* args) {
    PyArrayObject* data_array;
    double** data,**similarity, **diagonal;
    int rows, cols;
    PyObject* result;
    (void)self;

    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &data_array)) return NULL;

    data = numpy2D_to_C_matrix(data_array, &rows, &cols);
    if (!data) return NULL;

    similarity = calculate_similarity_matrix(data, rows, cols);
    free_matrix(data, rows);
    if (!similarity) return NULL;
    
    diagonal = calculate_diagonal_degree_matrix(similarity, rows);
    free_matrix(similarity, rows);
    if (!diagonal) return NULL;

    result = C_matrix_to_numpy2D(diagonal, rows, rows);
    return result;
}

/*
Computes the normalized similarity matrix W
input: Numpy 2D array
output: Normalized similarity matrix W (as a 2D numpy array)
*/
static PyObject* py_norm(PyObject* self, PyObject* args) {
    PyArrayObject* data_array;
    double** data, **similarity, **diagonal, **normalized;
    int rows, cols;
    PyObject* result;
    (void)self;

    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &data_array)){return NULL;}

    data = numpy2D_to_C_matrix(data_array, &rows, &cols);
    if (!data){return NULL;}

    similarity = calculate_similarity_matrix(data, rows, cols);
    diagonal = calculate_diagonal_degree_matrix(similarity, rows);
    normalized = calculate_normalized_similarity(similarity, diagonal, rows);

    free_matrix(data, rows);
    free_matrix(similarity, rows);
    free_matrix(diagonal, rows);

    if (!normalized) return NULL;

    result = C_matrix_to_numpy2D(normalized, rows, rows);
    return result;
}

/*
Performs full SymNMF on the normalized similarity matrix W, given an initial H
input: Numpy 2D array (representing W), Numpy 2D array (representing H) 
output: Optimized matrix H (as a 2D numpy array)
*/
static PyObject* py_symnmf(PyObject *self, PyObject *args) {
    PyArrayObject *w_array, *h_array;
    double **W, **H, **result_H;
    int n, k;
    PyObject *result;
    (void)self;

    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &w_array, &PyArray_Type, &h_array)) {return NULL;}

    W = numpy2D_to_C_matrix(w_array, &n, &n);
    if (!W) {return NULL;}

    H = numpy2D_to_C_matrix(h_array, &n, &k);
    if (!H) {return NULL;}

    result_H = symnmf_optimization(W, n, k, H);
    free_matrix(W, n);

    result = C_matrix_to_numpy2D(result_H, n, k);
    return result;
}

/*
Method definitions
*/
static PyMethodDef SymnmfMethods[] = {
    {"sym", py_sym, METH_VARARGS, "Calculate similarity matrix"},
    {"ddg", py_ddg, METH_VARARGS, "Calculate diagonal degree matrix"},
    {"norm", py_norm, METH_VARARGS, "Calculate normalized similarity matrix"},
    {"symnmf", py_symnmf, METH_VARARGS, "Perform SymNMF optimization"},
    {NULL, NULL, 0, NULL} 
};

/*
Module definition
*/
static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf",       
    "Symmetric Non-negative Matrix Factorization (SymNMF) implementation",
    -1,             
    SymnmfMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

/*
Initialization function
*/
PyMODINIT_FUNC PyInit_symnmf(void) {
    PyObject *m;
    
    import_array();
    
    m = PyModule_Create(&symnmfmodule);
    if (m == NULL) {
        return NULL;
    }
    
    return m;
}
