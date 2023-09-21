//
// Created by Samuel Martins on 12/12/18.
//

#ifndef IFT_DBL_ARRAY_H
#define IFT_DBL_ARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Array of double values.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
//! swig(extend = DblArrayExt.i, destroyer = iftDestroyDblArray, name = iftDblArray)
struct ift_dbl_array {
    /** Number of elements. */
    long n;
    /** Array of double values. */
    double *val;
};



/**
 * @brief Creates an iftDblArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
iftDblArray *iftCreateDblArray(long n);

/**
 * Copies a Double array, from data.
 * @author Peixinho
 * @date April, 2016
 */
iftDblArray *iftCopyDblArray(const double* array, long n);


/**
 * @brief Destroys an iftDblArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
void iftDestroyDblArray(iftDblArray **darr);


/**
 * @brief Reads a Double Array from a 1-D Numpy array.
 *
 * @param npy_path Pathname from the numpy array.
 * @return 1-D numpy array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftDblArray *iftReadDblArray(const char *npy_path, ...);


/**
 * @brief Writes a Double Array as a 1-D Numpy array.
 * @param  arr Array.
 * @param  npy_path Pathname from the Numpy array file.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteDblArray(const iftDblArray *arr, const char *npy_path, ...);


/**
 * @brief Reallocates memory for an iftDblArray and copies the original data. The new size could be higher or lower
 * @warning If the original array is larger then some data will be lost, i.e. only n elements will be copied
 * @author Cesar Castelo
 * @date Jul 18, 2018
 * @ingroup Memory
 */
void iftResizeDblArray(iftDblArray **iarr, long n);


/**
 * @brief Rearrange an array to a string for printing.
 * @author Samuka;
 * @date Sep 16, 2017
 */
char *iftDblArrayAsString(const iftDblArray *darr);


/**
 * @brief Convert an array (data) of the datatype <dtype> with <length> elements.
 *
 * @warning If the array is not a Double Array, the values are converted to.
 *
 * @param data Pointer to the array.
 * @param length Number of elements from the array.
 * @param cdtype Datatype of the array.
 * @return Converted array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftDblArray *iftConvertToDblArray(void *data, size_t length, iftCDataType cdtype);

#ifdef __cplusplus
}
#endif

#endif //IFT_DBL_ARRAY_H
