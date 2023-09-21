//
// Created by Samuel Martins on 11/12/18.
//

#ifndef IFT_ULONG_ARRAY_H
#define IFT_ULONG_ARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Array of unsigned long values.
 * @author Thiago Vallin spina
 * @date Mar 1, 2016
 * @ingroup Memory
 */
typedef struct ift_ulong_array {
    /** Number of elements. */
    long n;
    /** Array of integer values. */
    ulong *val;
} iftULongArray;


/**
 * @brief Creates an iftULongArray.
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016
 * @ingroup Memory
 */
iftULongArray *iftCreateULongArray(long n);


/**
 * @brief Destroys an iftULongArray.
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016
 * @ingroup Memory
 */
void iftDestroyULongArray(iftULongArray **iarr);


/**
 * @brief Reallocates memory for an iftULongArray and copies the original data. The new size could be higher or lower
 * @warning If the original array is larger then some data will be lost, i.e. only n elements will be copied
 * @author Cesar Castelo
 * @date Jul 18, 2018
 * @ingroup Memory
 */
void iftResizeULongArray(iftULongArray **iarr, long n);

#ifdef __cplusplus
}
#endif


#endif //IFT_ULONG_ARRAY_H
