//
// Created by Samuel Martins on 2018-12-16.
//

#ifndef IFT_INT_MATRIX_H
#define IFT_INT_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Matrix of integer values.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @note Its type definition in iftBasicDataType.h
 */
//! swig(extend = IntMatrixExt.i, destroyer = iftDestroyIntMatrix, name = iftIntMatrix)
struct ift_int_matrix {
    /** Matrix of elements. */
    int *val;
    /** Number of Columns, Rows, and look-up table to speed up index access. */
    int ncols, nrows, *tbrow;
    /** Total number of elements = ncols*nrows. */
    int n;
};



/**
 * @brief Creates a Matrix of Integer values.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 *
 * @param ncols Number of Columns.
 * @param nrows Number of Rows.
 * @return The created Integer Matrix.
 */
//! swig(newobject)
iftIntMatrix *iftCreateIntMatrix(int ncols, int nrows);


/**
 * @brief Destroys a Matrix of Integer values.
 *
 * @author Samuel Martins
 * @date October 15, 2015
 */
void iftDestroyIntMatrix(iftIntMatrix **iM);


/**
 * @brief Reads an Integer Matrix from a 2-D Numpy array.
 *
 * @param npy_path Pathname from the numpy array.
 * @return 2-D numpy array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftIntMatrix *iftReadIntMatrix(const char *npy_path, ...);


/**
 * @brief Writes an Integer Matrix as a 1-D Numpy array.
 * @param arr Array.
 * @param npy_path Pathname from the Numpy array file.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteIntMatrix(const iftIntMatrix *mat, const char *npy_path, ...);


/**
 * @brief Copies an integer Matrix
 * @author Samuka Martins
 * @date Sep 17, 2017
 */
iftIntMatrix *iftCopyIntMatrix(int *val, int nrows, int ncols);


/**
 * @brief Reallocates memory for a matrix and copies the original data. The new size could be higher or lower
 * @warning If the original matrix is larger then some data will be lost, i.e. only ncols columns and nrows rows will be copied
 * @param M Matrix pointer
 * @param ncols New number of columns
 * @param nrows New number of rows
 * @author Cesar Castelo
 * @date Jul 20, 2018
 */
void iftResizeIntMatrix(iftIntMatrix **M, int ncols, int nrows);


/**
 * @brief Rearrange a matrix to a string for printing.
 * @author Samuka
 * @date Sep 16, 2017
 */
char *iftIntMatrixAsString(const iftIntMatrix *mat);


/**
 * @brief Prints an Integer Matrix.
 * @author Samuka
 * @date May 7, 2019
 */
//! swig()
void iftPrintIntMatrix(const iftIntMatrix* M);


/**
 * @brief Merges an array of iftIntMatrix into a new iftIntMatrix (horizontally)
 * @param matArray Array of int matrices
 * @param n Size of the array
 * @author Cesar Castelo
 * @date Feb 16, 2018
 */
iftIntMatrix* iftMergeIntMatrixArrayHorizontal(iftIntMatrix **matArray, int n);

#ifdef __cplusplus
}
#endif

#endif //IFT_INT_MATRIX_H
