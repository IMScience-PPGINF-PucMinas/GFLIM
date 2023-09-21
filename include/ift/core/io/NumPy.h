//
// Created by Samuel Martins on 2018-12-14.
//

#ifndef IFT_NUMPY_H
#define IFT_NUMPY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/FloatArray.h"
#include "ift/core/io/Stream.h"


/**
 * @brief Struct to deal with the reading and writing of numpy array's headers.
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
typedef struct ift_numpy_header {
    /** Numpy datatype: e.g.: i1, u4, f8, ... */
    char *dtype;
    /** Flag to indicate if the array is stored in memory in the Fortran Ordem (column-major) */
    bool fortran_order;
    /** Number of dimensions from the array. */
    size_t n_dims;
    /** Shape of the array. */
    long *shape;
    /** Size of the array: multiplication of shape elements. */
    size_t size;
} iftNumPyHeader;




/**
 * @brief Creates a numpy header of a numpy array.
 *
 * @param cdtype CDataType from the numpy array.
 * @param shape Array's shape.
 * @param n_dims Number of dimensions of the array.
 * @return Created numpy header.
 *
 * @author Samuel Martins
 * @date Dec 16, 2018
 */
iftNumPyHeader *iftCreateNumPyHeader(iftCDataType cdtype, const long *shape, size_t n_dims);


/**
 * @brief Destroys a numpy header.
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftDestroyNumPyHeader(iftNumPyHeader **header);


/**
 * @brief Reads the a numpy array file. It returns the data and the numpy header.
 *
 * @note Based on http://pyopengl.sourceforge.net/pydoc/numpy.lib.format.html
 *
 * @param npy_path Numpy array pathname (*.npy).
 * @param header_out Reference to return the numpy header.
 * @return Numpy array data stream.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void *iftReadNumPy(const char *npy_path, iftNumPyHeader **header_out, ...);


/**
 * @brief Writes a NumPy array (header and data).
 *
 * @param header Numpy header.
 * @param data Array data stream.
 * @param npy_path Pathname to write the data (*.npy).
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
void iftWriteNumPy(const iftNumPyHeader *header, const void *data, const char *npy_path, ...);


/**
 * @brief Gets the corresponding CDataType from the numpy dtype.
 * @note Supported numpy dtypes: i1, i2, i4, i8, u1, u2, u4, u8, f4, f8
 *
 * @param dtype Numpy dtype
 * @return Corresponding CDataType.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
iftCDataType iftNumPyDTypeAsCDataType(const char *dtype);


/**
 * @brief Gets the corresponding NumPy dtype from a given CDataType.
 *
 * @param cdtype CDataType.
 * @return Corresponding numpy dtype.
 *
 * @author Samuel Martins
 * @date Dec 15, 2018
 */
char *iftCDataTypeAsNumPyDType(iftCDataType cdtype);

#ifdef __cplusplus
}
#endif


#endif //IFT_NUMPY_H
