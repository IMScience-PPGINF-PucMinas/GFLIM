#include "iftMemory.cuh"

#include "ift/core/tools/Dialog.h"
// these includes CANNOT be inside an extern "C" {
#include <cuda.h>
#include <cublas.h>


/* device=0, 1, etc */

void iftStartGPU(int device) {
    int nDevices;
    cudaError_t err = cudaGetDeviceCount(&nDevices);
    if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
    if ((device > nDevices-1) || (device < 0))
        iftError("Invalid device number","iftStartGPU");	

    err = cudaSetDevice(device);
    if (err != cudaSuccess) iftError("Check Cuda documentation for error: %s.\n", "iftStartGPU", err);
    
    cublasInit();
    struct cudaDeviceProp props;
    err = cudaGetDeviceProperties(&props, device);
    if (err != cudaSuccess) iftError("Check Cuda documentation for error: %s.\n", "iftStartGPU", err);

    printf("GPU-%d: %s @ %.0fMHz. %.0fGB/%.0fKB Global/Shared memory.\n\n", device, props.name, (float)props.clockRate/1000.0, (float)props.totalGlobalMem/pow(1024.0,3), props.sharedMemPerBlock/1024.0);

    size_t free, total;

    err = cudaMemGetInfo(&free, &total);
    if (err != cudaSuccess) {
        printf("Error cudaMemGetInfo\n");
    } else {
        printf("Free=%lu, total=%lu\n", free, total);
    }
}

void iftStopGPU() {
    cudaDeviceReset();
}

void iftCopyToGPU(void* dst, void* src, size_t size) {
    cudaError_t status = cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice);
    if(status != cudaSuccess) {
        iftError("Check Cuda documentation for error %d.\n", "iftCopyToGPU", status);
    }
}

void iftCopyFromGPU(void* dst, void* src, size_t size) {
    cudaError_t status = cudaMemcpy(dst, src, size, cudaMemcpyDeviceToHost);
    if(status != cudaSuccess) {
        iftError("Check Cuda documentation for error %d.\n", "iftCopyFromGPU", status);
    }
}

void iftCopyIntArrayFromGPU(int* dst, int* src, int n) {
    iftCopyFromGPU((void*)dst, (void*)src, n*sizeof(int));
}

void iftCopyIntArrayToGPU(int* dst, int* src, int n) {
    iftCopyToGPU((void*)dst, (void*)src, n*sizeof(int));
}

int* iftAllocIntArrayGPU(int n) {
    return (int*) iftAllocGPU(n, sizeof(int));
}

void iftCopyFloatArrayFromGPU(float* dst, float* src, int n) {
    iftCopyFromGPU((void*)dst, (void*)src, n*sizeof(float));
}

void iftCopyFloatArrayToGPU(float* dst, float* src, int n) {
    iftCopyToGPU((void*)dst, (void*)src, n*sizeof(float));
}

float* iftAllocFloatArrayGPU(int n) {
    return (float*) iftAllocGPU(n, sizeof(float));
}

void *iftAllocGPU(int n, size_t size) {
    void *mem = NULL;
    
    cudaError_t status = cudaMalloc(&mem, n * size);
    if(status==cudaErrorMemoryAllocation)//remove if slow down the function
    {
        iftError("Could not allocate enough memory.", "iftAllocGPU");
    }
    else if(status!=cudaSuccess) {
        iftError("Check Cuda documentation for error %d.\n", "iftAllocGPU", status);
    }

    return mem;
}

void iftFreeGPU(void *mem) {
    cudaFree(mem);
}

float iftGetFreeMemoryGPU(int device)
{
    size_t freeMem, totalMem;
    cudaError_t err = cudaMemGetInfo(&freeMem, &totalMem);
    if (err != cudaSuccess)
        printf("%s\n", cudaGetErrorString(err));

    return freeMem;
}


