//
// Created by Samuel Martins on 12/12/18.
//

#ifndef IFT_LONG_ARRAY_H
#define IFT_LONG_ARRAY_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Array of long values.
 * @author Thiago Vallin spina
 * @date Mar 1, 2016
 * @ingroup Memory
 */
typedef struct ift_long_array {
    /** Number of elements. */
    long n;
    /** Array of integer values. */
    long *val;
} iftLongArray;




/**
 * @brief Creates an iftLongArray.
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016
 * @ingroup Memory
 */
iftLongArray *iftCreateLongArray(long n);


/**
 * @brief Destroys an iftLongArray.
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016
 * @ingroup Memory
 */
void iftDestroyLongArray(iftLongArray **iarr);


/**
 * @brief Reads a Long Integer Array from a 1-D Numpy array.
 *
 * @param npy_path Pathname from the numpy array.
 * @return 1-D numpy array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftLongArray *iftReadLongArray(const char *npy_path, ...);


/**
 * @brief Writes a Long Integer Array as a 1-D Numpy array.
 * @param  arr Array.
 * @param  npy_path Pathname from the Numpy array file.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteLongArray(const iftLongArray *arr, const char *npy_path, ...);


/**
 * @brief Reallocates memory for an iftLongArray and copies the original data. The new size could be higher or lower
 * @warning If the original array is larger then some data will be lost, i.e. only n elements will be copied
 * @author Cesar Castelo
 * @date Jul 18, 2018
 * @ingroup Memory
 */
void iftResizeLongArray(iftLongArray **iarr, long n);


/**
 * @brief Convert an array (data) of the datatype <dtype> with <length> elements.
 *
 * @warning If the array is not a Long Integer Array, the values are converted to.
 *
 * @param data Pointer to the array.
 * @param length Number of elements from the array.
 * @param cdtype Datatype of the array.
 * @return Converted array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftLongArray *iftConvertToLongArray(void *data, size_t length, iftCDataType cdtype);

#ifdef __cplusplus
}
#endif

#endif //IFT_LONG_ARRAY_H
