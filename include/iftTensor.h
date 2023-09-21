#ifndef _IFTTENSOR_H_
#define _IFTTENSOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/DblArray.h"
#include "ift/core/dtypes/Matrix.h"
#include "ift/core/dtypes/FloatArray.h"
#include "iftCommon.h"

#include "iftMatrix.h"


/**
 * @brief Tensor multidimensional matrix
 * @author Peixinho
 * @date Nov, 2016
 */
typedef struct ift_tensor {
    float* val;/*Tensor values*/
    int* dimension;/*Tensor dimension*/
    int* accDimension;

    bool allocated;

    int ndimension;
    long n;/*Total number of elements*/
} iftTensor;

/**
 * @brief creates a Tensor with the dimensions specified
 * @param ndim
 * @param dim
 * @author Sep, 2017
 * @return tensor
 * @author Peixinho
 */
iftTensor* iftCreateTensorWithDimArray(int ndim, int* dim);
iftTensor* iftCreateTensorPointerWithDimArray(float *val, int ndim, int *dim);

/**
 * @brief creates a Tensor
 * @param ndimension Number of dimension, followed by the specified dimensions
 * @date Nov, 2016
 * @author Peixinho
 * @return a tensor object
 */
iftTensor* iftCreateTensor(int ndimension, ...);
iftTensor* iftCreateTensorPointer(float *val, int ndimension, ...);

/**
 * @brief Destroy a tensor object
 * @author Peixinho
 * @date Nov, 2016
 * @param tensor
 */
void iftDestroyTensor(iftTensor** tensor);

/**
 * @brief Destroy a tensor pointer object
 * @author Felipe Galvao
 * @date Nov, 2019
 * @param t Tensor
 */
void iftDestroyTensorPointer(iftTensor** t);

/**
 * @brief Create Tensor Copy
 * @author Deangeli
 * @date may, 2017
 * @param tensor
 */
iftTensor* iftCreateTensorCopy(iftTensor* src);

#define iftTensorElem(m, ...) (m)->val[iftGetTensorIndex((m),((m)->ndimension), __VA_ARGS__)]

#define iftTensorElemWithDimArray(m, a) (m)->val[iftGetTensorIndexWithDimArray((m), (a))]

void iftReshapeTensorWithDimensionArray(iftTensor *T, int ndimension, int *dimension);

void iftReshapeTensor(iftTensor *T, int ndimension, ...);

int iftGetTensorIndexDebug(iftTensor* tensor, int ndimension, ...);

int iftGetTensorIndex(iftTensor* tensor, int ndimension, ...);

int iftGetTensorIndexWithDimArray(iftTensor* tensor, int *dimension);

iftTensor* iftReadTensorCsv(char* tensorFilename, char* dimArrayFilename, char separator);

void iftWriteTensorCsv(iftTensor *tensor, char *tensorFilename, char *dimArrayFilename, char separator);

iftTensor *iftMultTensors(iftTensor *A, iftTensor *B);

iftTensor *iftSliceTensor(iftTensor *A, int **sliceIdx);

void iftPrintTensorDim(iftTensor* T);

char *iftTensorDimStr(iftTensor* T);

#ifdef __cplusplus
}
#endif

#endif
