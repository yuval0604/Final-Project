#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include "symnmf.h"

/*
Convert a pyhton list of lists to a C matrix
input: Python list, 2 pointers to store the number of rows and columns.
output: Matrix (or NULL if memory allocation failed)
*/
double** pylist_to_cmatrix(PyObject* py_mat, int* n_ptr, int* d_ptr) {
    int n, d, i, j;
    PyObject* row;
    double** matrix;

    if (!PyList_Check(py_mat)) {return NULL;}

    n = PyList_Size(py_mat);
    if (n == 0) {return NULL;}

    row = PyList_GetItem(py_mat, 0);
    if (!PyList_Check(row)) {return NULL;}
    d = PyList_Size(row);

    matrix = allocate_matrix(n, d);
    if (matrix == NULL) {return NULL;}

    for (i = 0; i < n; i++) {
        row = PyList_GetItem(py_mat, i);
        if (!PyList_Check(row) || PyList_Size(row) != d) {
            free_matrix(matrix, n);
            return NULL;
        }
        for (j = 0; j < d; j++) {
            matrix[i][j] = PyFloat_AsDouble(PyList_GetItem(row, j));
        }
    }

    *n_ptr = n;
    *d_ptr = d;
    return matrix;
}

/*
Convert a C matrix back to a python list of lists.
input: C Matrix, num of rows, num of columns.
output: A python list
*/
PyObject* cmatrix_to_pylist(double** matrix, int n, int m) {
    int i, j;
    PyObject* py_mat = PyList_New(n);
    if (py_mat == NULL) {return NULL;}

    for (i = 0; i < n; i++) {
        PyObject* row = PyList_New(m);
        if (row == NULL) {
            Py_DECREF(py_mat);
            return NULL;
        }
        for (j = 0; j < m; j++) {
            PyObject* item = PyFloat_FromDouble(matrix[i][j]);
            if (item == NULL) {
                Py_DECREF(row);
                Py_DECREF(py_mat);
                return NULL;
            }
            PyList_SetItem(row, j, item); 
        }
        PyList_SetItem(py_mat, i, row);
    }
    return py_mat;
}

/*
Computes the similarity matrix A
input: Python object - list of lists
output: Python object - Similarity matrix
*/
static PyObject* py_sym(PyObject* self, PyObject* args) {
    PyObject* py_data;
    double** data;
    double** result;
    int n, d;
    PyObject* output;
    (void)self;

    if (!PyArg_ParseTuple(args, "O", &py_data)) {return NULL;}

    data = pylist_to_cmatrix(py_data, &n, &d);
    if (data == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    result = calculate_similarity_matrix(data, n, d);
    free_matrix(data, n);

    if (result == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    output = cmatrix_to_pylist(result, n, n);
    free_matrix(result, n);

    return output;
}

/*
Computes the diagonal degree matrix D
input: Python object - list of lists
output: Python object - diagonal matrix
*/
static PyObject* py_ddg(PyObject* self, PyObject* args) {
    PyObject* py_data;
    double** data;
    double** similarity;
    double** result;
    int n, d;
    PyObject* output;
    (void)self;

    if (!PyArg_ParseTuple(args, "O", &py_data)) {return NULL;}

    data = pylist_to_cmatrix(py_data, &n, &d);
    if (data == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    similarity = calculate_similarity_matrix(data, n, d);
    free_matrix(data, n);
    if (similarity == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    result = calculate_diagonal_degree_matrix(similarity, n);
    free_matrix(similarity, n);

    if (result == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    output = cmatrix_to_pylist(result, n, n);
    free_matrix(result, n);

    return output;
}

/*
Computes the normalized similarity matrix W
input: Python object - list of lists
output: Python object - normalized similarity matrix
*/
static PyObject* py_norm(PyObject* self, PyObject* args) {
    PyObject* py_data;
    double** data;
    double** similarity;
    double** result;
    int n, d;
    PyObject* output;
    (void)self;

    if (!PyArg_ParseTuple(args, "O", &py_data)) {return NULL;}

    data = pylist_to_cmatrix(py_data, &n, &d);
    if (data == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    similarity = calculate_similarity_matrix(data, n, d);
    free_matrix(data, n);
    if (similarity == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    result = calculate_normalized_similarity(similarity, n);
    free_matrix(similarity, n);

    if (result == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }

    output = cmatrix_to_pylist(result, n, n);
    free_matrix(result, n);

    return output;
}

/* Performs full SymNMF on the normalized similarity matrix W, given an initial H
input: Python object - list of lists(H), Python object - list of lists(W) 
output: Python object - Optimized H */
static PyObject* py_symnmf(PyObject* self, PyObject* args) {
    PyObject* py_H;
    PyObject* py_W;
    double** H;
    double** W;
    int n, k;
    int save_dim;
    PyObject* output;
    (void)self;

    if (!PyArg_ParseTuple(args, "OO", &py_H, &py_W)) {return NULL;}
    H = pylist_to_cmatrix(py_H, &n, &k);
    if (H == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    W = pylist_to_cmatrix(py_W, &save_dim, &save_dim);
    if (W == NULL) {
        free_matrix(H, n);
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    if (save_dim != n) {
        free_matrix(H, n);
        free_matrix(W, save_dim);
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    if (symnmf(H, W, n, k) == NULL) {
        free_matrix(H, n);
        free_matrix(W, n);
        PyErr_SetString(PyExc_RuntimeError, "An Error Has Occurred");
        return NULL;
    }
    output = cmatrix_to_pylist(H, n, k);
    free_matrix(H, n);
    free_matrix(W, n);
    return output;
}

static PyMethodDef symnmfMethods[] = {
    {"sym", py_sym, METH_VARARGS, "Calculate similarity matrix"},
    {"ddg", py_ddg, METH_VARARGS, "Calculate diagonal degree matrix"},
    {"norm", py_norm, METH_VARARGS, "Calculate normalized similarity matrix"},
    {"symnmf", py_symnmf, METH_VARARGS, "Run symnmf algorithm"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf",
    NULL,
    -1,
    symnmfMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_symnmf(void) {
    return PyModule_Create(&symnmfmodule);
}
