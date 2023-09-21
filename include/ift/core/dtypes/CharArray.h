#ifndef IFT_CHAR_ARRAY_H
#define IFT_CHAR_ARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Array of char values (can be used to store int8 values).
 * @author Samuel Martins
 * @date Dec 15, 2018
 * @ingroup Memory
 *
 * @note Its type definition in iftBasicDataType.h
 */
//! swig(destroyer = iftDestroyCharArray, name = iftCharArray)
struct ift_char_array {
    /** Number of elements. */
    long n;
    /** Array of char values. */
    char *val;
};


/**
 * @brief Creates an iftCharArray.
 * @author Samuel Martins
 * @date Dec 15, 2018
 * @ingroup Memory
 */
iftCharArray *iftCreateCharArray(long n);


/**
 * @brief Destroys an iftCharArray.
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftDestroyCharArray(iftCharArray **carr);


/**
 * @brief Reads a Char Array from a 1-D Numpy array.
 *
 * @param npy_path Pathname from the numpy array.
 * @return 1-D numpy array.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
iftCharArray *iftReadCharArray(const char *npy_path, ...);


/**
 * @brief Writes a Char Array as a 1-D Numpy array.
 * @param  arr Array.
 * @param  npy_path Pathname from the Numpy array file.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteCharArray(const iftCharArray *arr, const char *npy_path, ...);


/**
 * @brief Convert an array (data) of the datatype <dtype> with <length> elements.
 *
 * @warning If the array is not a Char Array, the values are converted to.
 *
 * @param data Pointer to the array.
 * @param length Number of elements from the array.
 * @param cdtype Datatype of the array.
 * @return Converted array.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
iftCharArray *iftConvertToCharArray(void *data, size_t length, iftCDataType cdtype);

#ifdef __cplusplus
}
#endif

#endif //IFT_CHAR_ARRAY_H
