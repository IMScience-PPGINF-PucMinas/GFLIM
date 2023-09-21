#include "iftDeepLearning.cuh"

#include "ift/core/io/Stream.h"
#include "iftMemory.cuh"

#define iftTensorElemGPU3(m, ...) (m)->val[iftGetTensorIndexGPU3((m),((m)->ndimension), __VA_ARGS__)]
#define iftTensorElemGPU4(m, ...) (m)->val[iftGetTensorIndexGPU4((m),((m)->ndimension), __VA_ARGS__)]


__host__ __device__ int iftGetTensorIndexDebugGPU4(iftTensor* tensor, int ndimension, int d0, int d1, int d2, int d3) {

	if(ndimension == 3){
		printf("Invalid size to GPU @iftTensorElemGPU4.\n");
		return 0;
	}

	int dim;
    int index = 0;
    for (int i = 0; i < ndimension; ++i) {
		if(i == 0) dim = d0;
		if(i == 1) dim = d1;
		if(i == 2) dim = d2;
		if(i == 3) dim = d3;

        if(dim<0 || dim>=tensor->dimension[i]) {
            printf("Invalid Index. %dth dimension = %d @iftTensorElemGPU4.\n", i, dim);
        }
        index += dim * tensor->accDimension[i];
    }

    return index;
}

__host__ __device__ int iftGetTensorIndexDebugGPU3(iftTensor* tensor, int ndimension, int d0, int d1, int d2) {

	if(ndimension == 4){
		printf("Invalid size to GPU @iftTensorElemGPU3.\n");
		return 0;
	}

	int dim;
    int index = 0;
    for (int i = 0; i < ndimension; ++i) {
        if(i == 0) dim = d0;
		if(i == 1) dim = d1;
		if(i == 2) dim = d2;

        if(dim<0 || dim>=tensor->dimension[i]) {
            printf("Invalid Index. %dth dimension = %d @iftTensorElemGPU3.\n", i, dim);
        }
        index += dim * tensor->accDimension[i];
    }

    return index;
}

__host__ __device__ int iftGetTensorIndexGPU4(iftTensor* tensor, int ndimension, int d0, int d1, int d2, int d3) {

	if(ndimension == 3){
		printf("Invalid size to GPU @iftTensorElemGPU4.\n");
		return 0;
	}

	int dim;
    int index = 0;
    for (int i = 0; i < ndimension; ++i) {
		if(i == 0) dim = d0;
		if(i == 1) dim = d1;
		if(i == 2) dim = d2;
		if(i == 3) dim = d3;
        index += dim * tensor->accDimension[i];
    }

    return index;
}

__host__ __device__ int iftGetTensorIndexGPU3(iftTensor* tensor, int ndimension, int d0, int d1, int d2) {

	if(ndimension == 4){
		printf("Invalid size to GPU @iftTensorElemGPU3.\n");
		return 0;
	}

    int dim;
    int index = 0;
    for (int i = 0; i < ndimension; ++i) {
		if(i == 0) dim = d0;
		if(i == 1) dim = d1;
		if(i == 2) dim = d2;
        index += dim * tensor->accDimension[i];
    }

    return index;
}

__global__ void KernelCorrelationGPU(iftTensor* image, iftTensor* kernel, iftTensor* out) { 

	int col = blockIdx.x * blockDim.x +threadIdx.x;
	int row = blockIdx.y * blockDim.y +threadIdx.y;
	
	int batchSize = image->dimension[0];
	int width = image->dimension[2];
	int height = image->dimension[3];

	int nkernels = kernel->dimension[0];
	int nbands = kernel->dimension[1];
	int kwidth = kernel->dimension[2];
	int kheight = kernel->dimension[3];
	
	int kwidthHalf = (kwidth - 1)/2;
	int kheightHalf = (kheight - 1)/2;

	if((col >= kwidthHalf) && (col < (width-kwidthHalf)) && (row >= kheightHalf) && (row < height-kheightHalf)){
		for (int m = 0; m < batchSize; ++m) {		
		//kernels
			for (int k = 0; k < nkernels; ++k) {
				//bands
				double sum = 0.0;
				for (int b = 0; b < nbands; ++b) {
					//convolution
					for (int i = - kwidthHalf; i <= kwidthHalf; ++i) {
						for (int j = - kheightHalf; j <= kheightHalf ; ++j) {
							sum += iftTensorElemGPU4(image, m, b, col + i, row + j) * iftTensorElemGPU4(kernel, k, b, i+kwidthHalf, j+kheightHalf);
						}
					}
					iftTensorElemGPU4(out, m, k, col-kwidthHalf, row-kheightHalf) = sum;
				}
			}
		}
	}
}

__global__ void KernelConvolutionGPU(iftTensor* error, iftTensor* kernel, iftTensor* out) { 

	int x = blockIdx.x * blockDim.x +threadIdx.x;
	int y = blockIdx.y * blockDim.y +threadIdx.y;

  int batchSize = error->dimension[0];

//  int nbands = kernel->dimension[1];
  int width = error->dimension[2];
  int height = error->dimension[3];

  int nkernels = kernel->dimension[0];
  int nbands = kernel->dimension[1];
  int kwidth = kernel->dimension[2];
  int kheight = kernel->dimension[3];

  int kwidthHalf = (int)(kwidth - 1)/2;
  int kheightHalf = (int)(kheight - 1)/2;

	if((x >= 0) && (x < width) && (y >= 0) && (y < height)){

	  //image in batch
		for (int m = 0; m < batchSize; ++m) {
			//kernels
			for (int b = 0; b < nbands; ++b) {
						double sum = 0.0;
						//bands
						for (int i = 0; i < kwidth; ++i) {
							for (int j = 0; j < kheight ; ++j) {
								for (int k = 0; k < nkernels; ++k) {
									int row = x - (i - kwidthHalf);
									int col = y - (j - kheightHalf);
																
									if((row >= 0) && (col >= 0) && (row < width) && (col < height)){
										sum += iftTensorElemGPU4(kernel, k, b, i, j) * iftTensorElemGPU4(error, m, b, row, col);
									}
								}
							}
							iftTensorElemGPU4(out, m, b, x, y) = sum;						
		      			}
			}
		}
	}

}

__global__ void KernelConvolutionUpdateGPU(iftTensor* last, iftTensor* cur, iftTensor* out){

	int x = blockIdx.x * blockDim.x +threadIdx.x;
	int y = blockIdx.y * blockDim.y +threadIdx.y;

  int batchSize = last->dimension[0];

  int nbands = cur->dimension[1];
  int width = cur->dimension[2];
  int height = cur->dimension[3];

  int lwidth = last->dimension[2];
  int lheight = last->dimension[3];

  int nkernels = out->dimension[0];
  // int kbands = out->dimension[1];
  int kwidth = out->dimension[2];
  int kheight = out->dimension[3];

  int kwidthHalf = (int)(kwidth - 1)/2;
  int kheightHalf = (int)(kheight - 1)/2;


	if((x >= 0) && (x < width) && (y >= 0) && (y < height)){
	  //image in batch
		for (int k = 0; k < nkernels; ++k) {
			for (int m = 0; m < batchSize; ++m) {
				//kernels
				for (int b = 0; b < nbands; ++b) {
							//bands
							for (int i = 0; i < kwidth; ++i) {
								for (int j = 0; j < kheight ; ++j) {
									int row = x - (i - kwidthHalf);
									int col = y - (j - kheightHalf);
									//printf("i: %d j: %d\n", i, j);
									if((row >= 0) && (col >= 0) && (row < lwidth) && (col < lheight)){
										iftTensorElemGPU4(out, k, b, i, j) += iftTensorElemGPU4(cur, m, b, x, y) * iftTensorElemGPU4(last, m, k, row, col);
									}
								}
							
							}						
				}
			}
		}
	}


}

iftTensor* iftCreateTensorCopyGPU(iftTensor* src){

	iftTensor *dst = (iftTensor*) iftAlloc(1, sizeof(iftTensor));
	iftTensor *d_dst = (iftTensor*) iftAllocGPU(1, sizeof(iftTensor));

	dst->ndimension = src->ndimension;
    dst->allocated = true;
    dst->n = src->n;

	dst->dimension = iftAllocIntArrayGPU(src->ndimension);
	dst->accDimension = iftAllocIntArrayGPU(src->ndimension);
	dst->val = iftAllocFloatArrayGPU(src->n);

	iftCopyIntArrayToGPU(dst->dimension, src->dimension, src->ndimension);	
	iftCopyIntArrayToGPU(dst->accDimension, src->accDimension, src->ndimension);
	iftCopyFloatArrayToGPU(dst->val, src->val, src->n);
    
	iftCopyToGPU((void*)d_dst,(void*)dst, sizeof(iftTensor));

	iftFree(dst);

	return d_dst;
}

void iftDestroyTensorGPU(iftTensor* gpuTensor) {
	iftTensor* host = (iftTensor*)iftAlloc(1, sizeof(iftTensor));

	iftCopyFromGPU((void*)host, (void*)gpuTensor, sizeof(iftTensor));

	iftFreeGPU(host->dimension);
	iftFreeGPU(host->accDimension);
	iftFreeGPU(host->val);
	iftFreeGPU(gpuTensor);
	iftFree(host);
}

iftTensor* iftCreateTensorCopyFromGPU(iftTensor* src){

	iftTensor *dst = (iftTensor*) iftAlloc(1, sizeof(iftTensor));
	iftCopyFromGPU((void*)dst,(void*)src, sizeof(iftTensor));
	
	int *dimension = dst->dimension;
	int *accDimension = dst->accDimension;
	float *val = dst->val;
	
	dst->dimension = iftAllocIntArray(dst->ndimension);
	dst->accDimension = iftAllocIntArray(dst->ndimension);
	dst->val = iftAllocFloatArray(dst->n);
	
	iftCopyIntArrayFromGPU(dst->dimension, dimension, dst->ndimension);	
	iftCopyIntArrayFromGPU(dst->accDimension, accDimension, dst->ndimension);
	iftCopyFloatArrayFromGPU(dst->val, val, dst->n);

	iftFreeGPU(src);
	iftFreeGPU(dimension);
	iftFreeGPU(accDimension);
	iftFreeGPU(val);
	
	return dst;
}


void iftBatchToImageGPU(iftTensor* dst, iftTensor* src, int id_image){
	for (int i = 0; i < dst->ndimension; ++i)
		dst->dimension[i] = src->dimension[i+1];

	for (int b = 0; b < src->dimension[1]; ++b)
		for (int i = 0; i < src->dimension[2]; ++i)
		    for (int j = 0; j < src->dimension[3]; ++j)
		        iftTensorElem(dst, b, i, j) = iftTensorElem(src, id_image, b, i, j);
}


void iftImageToBatchGPU(iftTensor* dst, iftTensor* src, int id_image) {
	
	for (int b = 0; b < dst->dimension[1]; ++b)
		for (int i = 0; i < dst->dimension[2]; ++i)
		    for (int j = 0; j < dst->dimension[3]; ++j)
		         iftTensorElem(dst, id_image, b, i, j) = iftTensorElem(src, b, i, j);
}

void iftNeuralConvolutionGPU(iftTensor* error, iftTensor* kernel, iftTensor* out) {

	int tile = 32;
	//cudaError_t status; 
	//Cria cópia host
	iftTensor *h_image = iftCreateTensorCopy(error);
	iftTensor *h_kernel = iftCreateTensorCopy(kernel);
	iftTensor *h_out = iftCreateTensorCopy(out);

//	//Aloca Memória na GPU
//	iftTensor *d_image = (iftTensor*)iftAllocGPU(1, sizeof(iftTensor));
//	iftTensor *d_kernel = (iftTensor*)iftAllocGPU(1, sizeof(iftTensor));
//	iftTensor *d_out = (iftTensor*)iftAllocGPU(1, sizeof(iftTensor));
	
	//Transfere dado pra GPU
	iftTensor* d_image = iftCreateTensorCopyGPU(h_image);
	iftTensor* d_kernel = iftCreateTensorCopyGPU(h_kernel);
	iftTensor* d_out = iftCreateTensorCopyGPU(h_out);
	
	//Kernel
	dim3 dimGrid(ceil((float)error->dimension[2]/tile), ceil((float)error->dimension[3]/tile), 1);
	dim3 dimBlock(tile, tile, 1);
	KernelConvolutionGPU<<<dimGrid, dimBlock>>>(d_image, d_kernel, d_out);
//	cudaDeviceSynchronize();

	/*status = cudaGetLastError();

	if(status==cudaErrorMemoryAllocation)//remove if slow down the function
	{
		iftError("Could not allocate enough memory.", "iftAlloc");
	}
    else if(status!=cudaSuccess) {
		iftError("Check Cuda documentation for error %d.\n", "iftAllocGPU", status);
    }*/


	//Transfere dado para CPU
	iftTensor *temp;
	temp = iftCreateTensorCopyFromGPU(d_out);

	for(int i = 0; i < temp->n; i++){
		out->val[i] = temp->val[i];
		//printf("temp: %f out: %f\n", temp->val[i], out->val[i]);
	}

	iftDestroyTensor(&h_image);
	iftDestroyTensor(&h_kernel);
	iftDestroyTensor(&h_out);
	iftDestroyTensor(&temp);
	iftDestroyTensorGPU(d_image);
	iftDestroyTensorGPU(d_kernel);
//	iftDestroyTensorGPU(d_out);
}

void iftNeuralConvolutionUpdateGPU(iftTensor* error, iftTensor* kernel, iftTensor* out) {

	int tile = 32;
	//cudaError_t status; 
	//Cria cópia host
	iftTensor *h_image = iftCreateTensorCopy(error);
	iftTensor *h_kernel = iftCreateTensorCopy(kernel);
	iftTensor *h_out = iftCreateTensorCopy(out);

	//Aloca Memória na GPU
//	iftTensor *d_image = (iftTensor*)iftAllocGPU(1, sizeof(iftTensor));
//	iftTensor *d_kernel = (iftTensor*)iftAllocGPU(1, sizeof(iftTensor));
//	iftTensor *d_out = (iftTensor*)iftAllocGPU(1, sizeof(iftTensor));
	
	//Transfere dado pra GPU
	iftTensor *	d_image = iftCreateTensorCopyGPU(h_image);
	iftTensor *	d_kernel = iftCreateTensorCopyGPU(h_kernel);
	iftTensor *	d_out = iftCreateTensorCopyGPU(h_out);
	
	//Kernel
	dim3 dimGrid(ceil((float)error->dimension[2]/tile), ceil((float)error->dimension[3]/tile), 1);
	dim3 dimBlock(tile, tile, 1);
	KernelConvolutionUpdateGPU<<<dimGrid, dimBlock>>>(d_image, d_kernel, d_out);
//	cudaDeviceSynchronize();

	/*status = cudaGetLastError();

	if(status==cudaErrorMemoryAllocation)//remove if slow down the function
	{
		iftError("Could not allocate enough memory.", "iftAlloc");
	}
    else if(status!=cudaSuccess) {
		iftError("Check Cuda documentation for error %d.\n", "iftAllocGPU", status);
    }*/


	//Transfere dado para CPU
	iftTensor *temp;
	temp = iftCreateTensorCopyFromGPU(d_out);

	for(int i = 0; i < temp->n; i++){
		out->val[i] = temp->val[i];
		//printf("temp: %f out: %f\n", temp->val[i], out->val[i]);
	}

	iftDestroyTensor(&h_image);
	iftDestroyTensor(&h_kernel);
	iftDestroyTensor(&h_out);
	iftDestroyTensor(&temp);
	iftDestroyTensorGPU(d_image);
	iftDestroyTensorGPU(d_kernel);
//	iftDestroyTensorGPU(d_out);
}

void iftComputeCorrelationGPU(iftTensor* image, iftTensor* kernel, iftTensor* out) {

	int tile = 32;

	//Cria cópia host
	iftTensor *h_image = iftCreateTensorCopy(image);
	iftTensor *h_kernel = iftCreateTensorCopy(kernel);
	iftTensor *h_out = iftCreateTensorCopy(out);

//	//Aloca Memória na GPU
//	iftTensor *d_image = (iftTensor*)iftAllocGPU(1, sizeof(iftTensor));
//	iftTensor *d_kernel = (iftTensor*)iftAllocGPU(1, sizeof(iftTensor));
//	iftTensor *d_out = (iftTensor*)iftAllocGPU(1, sizeof(iftTensor));
	
	//Transfere dado pra GPU
	iftTensor* d_image = iftCreateTensorCopyGPU(h_image);
	iftTensor* d_kernel = iftCreateTensorCopyGPU(h_kernel);
	iftTensor* d_out = iftCreateTensorCopyGPU(h_out);
	
	//Kernel
	dim3 dimGrid(ceil((float)image->dimension[2]/tile), ceil((float)image->dimension[3]/tile), 1);
	dim3 dimBlock(tile, tile, 1);
	KernelCorrelationGPU<<<dimGrid, dimBlock>>>(d_image, d_kernel, d_out);
//	cudaDeviceSynchronize();

	/*error = cudaGetLastError();

	if(error==cudaErrorMemoryAllocation)//remove if slow down the function
	{
		iftError("Could not allocate enough memory.", "iftAlloc");
	}
    else if(error!=cudaSuccess) {
		iftError("Check Cuda documentation for error %d.\n", "iftAllocGPU", error);
    }*/


	//Transfere dado para CPU
	iftTensor *temp;
	temp = iftCreateTensorCopyFromGPU(d_out);

	for(int i = 0; i < temp->n; i++){
		out->val[i] = temp->val[i];
		//printf("temp: %f out: %f\n", temp->val[i], out->val[i]);
	}

	iftDestroyTensor(&h_image);
	iftDestroyTensor(&h_kernel);
	iftDestroyTensor(&h_out);
	iftDestroyTensor(&temp);
	iftDestroyTensorGPU(d_image);
	iftDestroyTensorGPU(d_kernel);
//	iftDestroyTensorGPU(d_out);
}

void iftNeuralCorrelationGPU(iftTensor* data, iftTensor* kernel, iftTensor* out) {

	//int batchSize = data->dimension[0];

	//for (int m = 0; m < batchSize; ++m) {
		//iftTensor *result;
		//iftTensor *image = iftCreateTensor(3, data->dimension[1], data->dimension[2], data->dimension[3]);
		//iftTensor *image_out = iftCreateTensor(3, out->dimension[1], out->dimension[2], out->dimension[3]);
		//iftBatchToImageGPU(image, data, m); // m
		
		//result = iftComputeConvolutionGPU(image, kernel, image_out);
		iftComputeCorrelationGPU(data, kernel, out);

		//iftImageToBatchGPU(out, result, m);//m
		//iftDestroyTensor(&image);
		//iftDestroyTensor(&result);
		//iftDestroyTensor(&image_out);
  	//}
}




//int main() {
//
//    iftMatrix* M = iftRandomMatrix(5,5,0.0,1.0);
//
//    iftMatrix* N = iftRandomMatrix(5,5,0.0,1.0);
//
//    iftPrintMatrix(iftMultMatrices(M, N));
//
//    return 0;
//}

