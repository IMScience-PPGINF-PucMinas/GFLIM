#include "iftFLIM.cuh"

#include <cuda.h>
#include <cublas.h>


int iftNumberOfDevices()
{
  int nDevices;
  cudaError_t err = cudaGetDeviceCount(&nDevices);
  if (err != cudaSuccess){
    iftWarning("%s","NumberOfDevices",cudaGetErrorString(err));
    return(0);
  }
  return(nDevices);
}

bool iftStartDevice(int device)
{
  cudaError_t err = cudaSetDevice(device);
  if (err != cudaSuccess){
    iftWarning("%s","StartDevice",cudaGetErrorString(err));
    return(false);
  }
  cublasInit();
  return(true);
}

bool iftStopDevice(int device)
{
  cudaError_t err = cudaDeviceReset();
  if (err != cudaSuccess){
      iftWarning("%s","StopDevice",cudaGetErrorString(err));
    return(false);
  }
  return(true);
}

int iftFLIMBatchSizeGPU(iftFLIMArch *arch, int input_image_nvoxels, int input_image_nchannels, int device)
{
  float freeMemoryMb, percMemory=0.95; 

  if (device < 0) device = 0;
  freeMemoryMb = (float)iftGetFreeMemoryGPU(device)/1024.0/1024.0;
  struct cudaDeviceProp props;
  cudaError_t err = cudaGetDeviceProperties(&props, device);
  if (err == cudaSuccess){
    printf("GPU-%d: %s @ %.0fMHz. %.0fGB of GPU memory.\n\n", device, props.name, (float)props.clockRate/1000.0, (float)props.totalGlobalMem/pow(1024.0,3));
  }
  float nMbytesPerDouble   = sizeof(double)/1024.0/1024.0;
  float ninput_channels   = input_image_nchannels;
  float nMbytesInputImage = input_image_nvoxels*ninput_channels*nMbytesPerDouble;
  float max_requiredMemory = 0.0;
  for (int l=0; l < arch->nlayers; l++){
    int KS2 = 1;
    if (arch->layer[l].kernel_size[2] != 0)
       KS2 = arch->layer[l].kernel_size[2];
    float kernelsize =    arch->layer[l].kernel_size[0]*arch->layer[l].kernel_size[1]*KS2;
    float nMbytesKernels      = arch->layer[l].noutput_channels*kernelsize*ninput_channels*nMbytesPerDouble;
    float nMbytesOutputImage  = input_image_nvoxels*arch->layer[l].noutput_channels*nMbytesPerDouble;
    float requiredMemory  = nMbytesKernels + nMbytesOutputImage + nMbytesInputImage;
     if (requiredMemory > max_requiredMemory){ 
       max_requiredMemory  = requiredMemory; 
     } 
    ninput_channels        = arch->layer[l].noutput_channels;
    nMbytesInputImage      = nMbytesOutputImage/arch->layer[l].pool_stride;
    input_image_nvoxels    = input_image_nvoxels/arch->layer[l].pool_stride;
  }
  
  int batchsize = (int)floor(percMemory*freeMemoryMb/(1.5*max_requiredMemory));

  return(batchsize);
}
