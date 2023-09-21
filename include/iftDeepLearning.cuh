#ifndef IFT_DEEPLEARNING_CUH
#define IFT_DEEPLEARNING_CUH

#ifdef __cplusplus
extern "C" {
#endif

#include "iftMatrix.cuh"

/**
 * @brief Init the GPU devices.
 * @warning If not called, the GPU devices will be started in the first call for a GPU function.
 * @author Peixinho
 */

iftTensor* iftCreateTensorCopy(iftTensor* src);

iftTensor* iftCreateTensorCopyGPU(iftTensor* src);

//void iftBatchToImage(iftTensor* dst, iftTensor* src, int id_image);

void iftImageToBatch(iftTensor* dst, iftTensor* src, int id_image);

void iftComputeConvolutionGPU(iftTensor* image, iftTensor* kernel, iftTensor* out);
void iftComputeCorrelationGPU(iftTensor* image, iftTensor* kernel, iftTensor* out);

void iftNeuralCorrelationGPU(iftTensor* data, iftTensor* kernel, iftTensor* out);
void iftNeuralConvolutionGPU(iftTensor* data, iftTensor* kernel, iftTensor* out);
void iftNeuralConvolutionUpdateGPU(iftTensor* data, iftTensor* kernel, iftTensor* out);


#ifdef __cplusplus
}
#endif

#endif//IFT_DEEPLEARNING_CUH
