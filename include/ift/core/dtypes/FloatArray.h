//
// Created by Samuel Martins on 12/12/18.
//

#ifndef IFT_FLOAT_ARRAY_H
#define IFT_FLOAT_ARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/DblArray.h"

/**
 * @brief Array of float values.
 * @author Thiago Vallin Spina
 * @date Mar 4, 2016
 * @ingroup Memory
 */
//! swig(extend = FloatArrayExt.i, destroyer = iftDestroyFloatArray)
typedef struct ift_flt_array {
    /** Number of elements. */
    long n;
    /** Array of float values. */
    float *val;
} iftFloatArray;



/**
 * @brief Creates an iftFloatArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
iftFloatArray *iftCreateFloatArray(long n);


/**
 * @brief Destroys an iftDblArray.
 * @author Samuel Martins
 * @date Oct 15, 2015
 * @ingroup Memory
 */
void iftDestroyFloatArray(iftFloatArray **darr);


/**
 * @brief Reads a Float Array from a 1-D Numpy array.
 *
 * @param npy_path Pathname from the numpy array.
 * @return 1-D numpy array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftFloatArray *iftReadFloatArray(const char *npy_path, ...);


/**
 * @brief Writes a Float Array as a 1-D Numpy array.
 * @param  arr Array.
 * @param  npy_path Pathname from the Numpy array file.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteFloatArray(const iftFloatArray *arr, const char *npy_path, ...);


/**
 * @brief Reallocates memory for an iftFloatArray and copies the original data. The new size could be higher or lower
 * @warning If the original array is larger then some data will be lost, i.e. only n elements will be copied
 * @author Cesar Castelo
 * @date Jul 18, 2018
 * @ingroup Memory
 */
void iftResizeFloatArray(iftFloatArray **iarr, long n);


/**
 * @brief Rearrange an array to a string for printing.
 * @author Samuka;
 * @date Dec 28, 2017
 */
//! swig(newobject, stable)
char *iftFloatArrayAsString(const iftFloatArray *farr);


/**
 * @brief Convert an array (data) of the datatype <dtype> with <length> elements.
 *
 * @warning If the array is not a Float Array, the values are converted to.
 *
 * @param data Pointer to the array.
 * @param length Number of elements from the array.
 * @param cdtype Datatype of the array.
 * @return Converted array.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftFloatArray *iftConvertToFloatArray(void *data, size_t length, iftCDataType cdtype);

/**
 * @brief Creates a iftDblArray from a iftFloatArray and returns it
 *
 * @param array_src The source float array.
 * @ingroup Memory
 * 
 * @author Cesar Castelo
 * @date Feb 8, 2019
 */
iftDblArray *iftFloatArrayToDblArray(iftFloatArray *array_src);

/**
 * @brief Creates a iftFloatArray from a iftDblArray and returns it
 *
 * @param array_src The source double array.
 * @ingroup Memory
 * 
 * @author Cesar Castelo
 * @date Feb 8, 2019
 */
iftFloatArray *iftDblArrayToFloatArray(iftDblArray *array_src);



#ifdef __cplusplus
}
#endif

#endif //IFT_FLOAT_ARRAY_H
