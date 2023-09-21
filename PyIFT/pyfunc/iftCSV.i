#include "ift.h"
#include <Python.h>

//! swig(newobject)
iftCSV *PyListToCSV(PyObject *list) {
    int nrows = PyList_Size(list);

    if (nrows <= 0) {
        char error[512];
        sprintf(error, "Number of Rows in Python CSV List for creating a CSV %d is <= 0", nrows);
        SWIG_Error(0, error);
        return NULL;
    }


    PyObject *row0 = PyList_GetItem(list, 0);
    int ncols = PyList_Size(row0);
    
    if (ncols <= 0) {
        char error[512];
        sprintf(error, "Number of Columns in Python CSV List for creating a CSV %d is <= 0", ncols);
        SWIG_Error(0, error);
        return NULL;
    }


    iftCSV *csv = iftCreateCSV(nrows, ncols);


    for (int r = 0; r < nrows; r++) {
        PyObject *row = PyList_GetItem(list, r);
        
        int ncols_of_row = PyList_Size(row);
        if (ncols_of_row != ncols) {
            char error[512];
            sprintf(error, "Number of Columns from Row %d is different from the Number of Columns from the CSV: %d != %d", r, ncols_of_row, ncols);
            SWIG_Error(0, error);
            return NULL;
        }


        for (int c = 0; c < ncols; c++) {
            PyObject *col = PyUnicode_AsEncodedString(PyList_GetItem(row, c), "UTF-8", "strict");
            strcpy(csv->data[r][c], PyString_AsString(col));
        }
    }
    
    return csv;
}