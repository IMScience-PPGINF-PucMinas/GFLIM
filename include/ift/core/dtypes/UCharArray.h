//
// Created by Samuel Martins on 11/12/18.
//

#ifndef IFT_UCHAR_ARRAY_H
#define IFT_UCHAR_ARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Array of uchar values (can be used to store uint8 values).
 * @author Thiago Vallin Spina
 * @date Mar 13, 2016
 * @ingroup Memory
 */
typedef struct ift_uchar_array {
    /** Number of elements. */
    long n;
    /** Array of uchar values. */
    uchar *val;
} iftUCharArray;



/**
 * @brief Creates an iftUCharArray.
 * @author Thiago Vallin Spina
 * @date Mar 13, 2016
 * @ingroup Memory
 */
iftUCharArray *iftCreateUCharArray(long n);


/**
 * @brief Destroys an iftUCharArray.
 * @author Thiago Vallin Spina
 * @date Mar 13, 2016
 * @ingroup Memory
 */
void iftDestroyUCharArray(iftUCharArray **iarr);


/**
 * @brief Reads an Unsigned Char Array from a 1-D Numpy array.
 *
 * @param npy_path Pathname from the numpy array.
 * @return 1-D numpy array.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
iftUCharArray *iftReadUCharArray(const char *npy_path, ...);


/**
 * @brief Writes an Unsigned Char Array as a 1-D Numpy array.
 * @param  arr Array.
 * @param  npy_path Pathname from the Numpy array file.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteUCharArray(const iftUCharArray *arr, const char *npy_path, ...);


/**
 * @brief Reallocates memory for an iftUCharArray and copies the original data. The new size could be higher or lower
 * @warning If the original array is larger then some data will be lost, i.e. only n elements will be copied
 * @author Cesar Castelo
 * @date Jul 18, 2018
 * @ingroup Memory
 */
void iftResizeUCharArray(iftUCharArray **iarr, long n);


/**
 * @brief Convert an array (data) of the datatype <dtype> with <length> elements.
 *
 * @warning If the array is not an Unsigned Char Array, the values are converted to.
 *
 * @param data Pointer to the array.
 * @param length Number of elements from the array.
 * @param cdtype Datatype of the array.
 * @return Converted array.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
iftUCharArray *iftConvertToUCharArray(void *data, size_t length, iftCDataType cdtype);

#ifdef __cplusplus
}
#endif

#endif //IFT_UCHAR_ARRAY_H
