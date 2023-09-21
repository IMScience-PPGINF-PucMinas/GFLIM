#ifndef IFT_MATRIX_CUH
#define IFT_MATRIX_CUH

#ifdef __cplusplus
extern "C" {
#endif

#include "iftMatrix.h"
#include "iftTensor.h"



/**
 * @brief Multiplies matrices A*B in GPU device.
 * @author Peixinho
 * @date Jan 2017
 * @warning This function is intended to be used two multiply large matrices (Prefereably with all the dimensions above 600). If iftInitGPU() is not called, the initialization of GPU will be called here, what may cause some overhead.
 */
void iftMultMatricesInPlaceGPU(const iftMatrix *A, const iftMatrix *B, bool transposeA, bool transposeB, iftMatrix **C);

void iftCopyTensorToGPU(iftTensor* dst, iftTensor* src, int n);

iftTensor* iftAllocTensorGPU(int n);

#ifdef __cplusplus
}
#endif

#endif//IFT_MATRIX_CUH
