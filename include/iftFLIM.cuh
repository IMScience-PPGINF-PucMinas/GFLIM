/**
 * @file 
* @brief FLIM: Feature Learning From Image Markers -- auxiliary functions
* that depend on cuda.
 * 
 * @note <b>Programs:</b>
 * * @ref iftFLIM-LearnModel.c = Learns the CNN model
 * * @ref iftFLIM-ExtractFeatures.c = Extracts image features.
 */

#ifndef IFT_FLIM_CUH
#define IFT_FLIM_CUH

#ifdef __cplusplus
extern "C" {
#endif

#include "iftFLIM.h"
#include "iftMemory.cuh"

int  iftNumberOfDevices();
bool iftStartDevice(int device);
int  iftFLIMBatchSizeGPU(iftFLIMArch *arch, int input_image_nvoxels, int input_image_nchannels, int device);
bool iftStopDevice(int device);
    
#ifdef __cplusplus
}
#endif


#endif

