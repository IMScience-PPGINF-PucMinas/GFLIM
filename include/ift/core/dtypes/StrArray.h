//
// Created by Samuel Martins on 12/12/18.
//

#ifndef IFT_STR_ARRAY_H
#define IFT_STR_ARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Array of strings.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 *
 * @note Its type definition in iftBasicDataType.h
 */
struct ift_str_array {
    /** Number of elements. */
    long n;
    /** Array of strings. */
    char **val;
};


/**
 * @brief Creates an Array of Strings with at most 2048 characters.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
iftStrArray *iftCreateStrArray(long n);


/**
 * @brief Destroys an iftStrArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
void iftDestroyStrArray(iftStrArray **sarr);


/**
 * @brief Reallocates memory for an iftStrArray and copies the original data. The new size could be higher or lower
 * @warning If the original array is larger then some data will be lost, i.e. only n elements will be copied
 * @author Cesar Castelo
 * @date Jul 18, 2018
 * @ingroup Memory
 */
void iftResizeStrArray(iftStrArray **iarr, long n);


/**
 * @brief Copies an iftStrArray.
 * @author Samuel Martins
 * @date Sep 17, 2017
 * @ingroup Memory
 */
iftStrArray *iftCopyStrArray(char **str_arr, long n);


/**
 * @brief Rearrange an array to a string for printing.
 * @author Samuka;
 * @date Sep 16, 2017
 */
char *iftStrArrayAsString(const iftStrArray *sarr);


/**
 * @brief Creates a string array with specific values.
 * @author Peixinho
 * @param n Number of values to search.
 * @param ... Sequence of values to search.
 */
iftStrArray *iftStrValues(int n, ...);

#ifdef __cplusplus
}
#endif

#endif //IFT_STR_ARRAY_H
