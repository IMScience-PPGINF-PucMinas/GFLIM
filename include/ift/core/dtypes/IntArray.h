//
// Created by Samuel Martins on 10/12/18.
//

#ifndef IFT_INTARRAY_H
#define IFT_INTARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"

/**
 * @brief Array of integer values.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 *
 * @note Its type definition in iftBasicDataType.h
 */
//! swig(extend = IntArrayExt.i, destroyer = iftDestroyIntArray, name = iftIntArray)
struct ift_int_array {
    /** Number of elements. */
    long n;
    /** Array of integer values. */
    int *val;
};



/**
 * @brief Creates an iftIntArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
//! swig(newobject)
iftIntArray *iftCreateIntArray(long n);


/**
 * @brief Destroys an iftIntArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
void iftDestroyIntArray(iftIntArray **iarr);


/**
 * @brief Reads an Integer Array from a 1-D Numpy array.
 *
 * @param npy_path Pathname from the numpy array.
 * @return 1-D numpy array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftIntArray *iftReadIntArray(const char *npy_path, ...);


/**
 * @brief Writes an Integer Array as a 1-D Numpy array.
 * @param  arr Array.
 * @param  npy_path Pathname from the Numpy array file.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteIntArray(const iftIntArray *arr, const char *npy_path, ...);


/**
 * @brief Reallocates memory for an iftIntArray and copies the original data. The new size could be higher or lower
 * @warning If the original array is larger then some data will be lost, i.e. only n elements will be copied
 * @author Cesar Castelo
 * @date Jul 18, 2018
 * @ingroup Memory
 */
void iftResizeIntArray(iftIntArray **iarr, long n);


/**
 * @brief Shuffles an integer array by the modern version Fisher-Yates shuffle algorithm.
 * @details See https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 * @param array Integer array
 * @param n Size of the array
 */
void iftShuffleIntArray(int* array, int n);


/**
 * @brief Converts a voxel v to an Integer Array of size 3.
 * [0] = x, [1] = y, [2] = z
 *
 * @author Samuka
 * @date Nov 16, 2017
 */
iftIntArray *iftVoxelToIntArray(iftVoxel v);


/**
 * @brief Converts an integer array of size 3 to a voxel
 * [0] = x, [1] = y, [2] = z
 *
 * @author Samuka
 * @date Nov 16, 2017
 */
iftVoxel iftIntArrayToVoxel(const iftIntArray *arr);


/**
 * @brief Rearrange an array to a string for printing.
 * @author Samuka;
 * @date Sep 16, 2017
 */
char *iftIntArrayAsString(const iftIntArray *iarr);


/**
 * @brief Convert an array (data) of the datatype <dtype> with <length> elements.
 *
 * @warning If the array is not an Integer Array, the values are converted to.
 *
 * @param data Pointer to the array.
 * @param length Number of elements from the array.
 * @param cdtype Datatype of the array.
 * @return Converted array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftIntArray *iftConvertToIntArray(void *data, size_t length, iftCDataType cdtype);

#ifdef __cplusplus
}
#endif

#endif //IFT_INTARRAY_H
