//
// Created by Samuel Martins on 2018-12-16.
//

#ifndef IFT_MATRIX_H
#define IFT_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Struct of a float Matrix.
 */
//! swig(extend = iftMatrixExt.i, destroyer = iftDestroyMatrix, name = iftMatrix)
struct ift_matrix {
    /** Number of rows */
    int nrows;
    /** Number of columns */
    int ncols;
    /** Total number of matrix elements: nrows * ncols */
    long n;
    /** 1D array with the matrix values */
    float *val;
    /** Look-up table to speed up index access */
    long *tbrow;
    bool allocated;
};



/**
 * @brief Creates a (float) Matrix with <nrows> rows and <ncols> cols.
 */
//! swig(newobject)
iftMatrix *iftCreateMatrix(int ncols, int nrows);

iftMatrix *iftCreateMatrix_omp(int ncols, int nrows);

/**
 * @brief Destroys a (float) Matrix.
 */
void iftDestroyMatrix(iftMatrix **M);


/**
 * @brief Reads a (float) Matrix from a 2-D Numpy array.
 *
 * @param npy_path Pathname (*.npy) from the numpy array.
 * @return 2-D numpy array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftMatrix *iftReadMatrix(const char *npy_path, ...);


/**
 * @brief Writes an Integer Matrix as a 1-D Numpy array.
 * @param arr Array.
 * @param npy_path Pathname (*.npy) from the Numpy array file.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteMatrix(const iftMatrix *mat, const char *npy_path, ...);


/**
 * @brief Prints a (float) Matrix.
 */
//! swig()
void iftPrintMatrix(const iftMatrix *M);


#ifdef __cplusplus
}
#endif

#endif //IFT_MATRIX_H
