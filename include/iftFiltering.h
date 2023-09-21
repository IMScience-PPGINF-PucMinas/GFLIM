#ifndef IFT_FILTERING_H_
#define IFT_FILTERING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftFImage.h"
#include "iftMImage.h"
#include "iftAdjacency.h"
#include "iftSort.h"
#include "iftMatrix.h"
#include "iftKernel.h"
#include "ift/imgproc/dtypes/Roi.h"

//! swig(newobject, stable)
iftImage  *iftLinearFilter(const iftImage *img, iftKernel *K);

//! swig(newobject, stable)
iftImage  *iftLinearFilterInRegion(iftImage *img, iftImage *mask, iftKernel *K);
//! swig(newobject)
iftImage  *iftFastLinearFilter(iftImage *img, iftKernel *K, char crop);
iftImage  *iftCroppedFastLinearFilterByMatrixMult(iftImage *img, iftKernel *K);
//! swig(newobject)
iftFImage *iftFLinearFilter(iftFImage *img, iftKernel *K);
//! swig(newobject)
iftFImage  *iftFLinearFilterInRegion(iftFImage *img, iftImage *mask, iftKernel *K);
//! swig(newobject)
iftFImage *iftFastFLinearFilter(iftFImage *img, iftKernel *K, char crop);
//! swig(newobject, stable)
iftMImage *iftMLinearFilter(iftMImage *img, iftMKernel *K);
//! swig(newobject, stable)
iftMImage *iftMMLinearFilter(iftMImage *img, iftMMKernel *k_bank);

/**
 * @brief Apply the Median Filter on an Input Image.
 * @param  img Input Image to be filtered.
 * @param  Ain Input adjacent used in filtering. If NULL, a 4-neighborhood (2D) or
 *             6-neighborhood (3D) is considered.
 * @return     Filtered Image.
 */
//! swig(newobject, stable)
iftImage *iftMedianFilter(const iftImage *img, iftAdjRel *Ain);


/**
 * @brief Applies a Median filter on an Input Image by considering all neighbors from an adjacency
 * relation in the computation.
 * 
 * @param img Image to be filtered.
 * @param Ain Adjacency Relation. If NULL, a 4-neighborhood (2D) or
 *            neighborhood (3D) is considered.
 * 
 * @return Filtered image.
 * 
 * @author Samuel Martins
 * @date Sep 14, 2018
 */
//! swig(newobject, stable)
iftImage *iftMeanFilter(const iftImage *img, const iftAdjRel *Ain);
//! swig(newobject, stable)
iftMImage *iftMMedianFilter(iftMImage *img, iftAdjRel *A);
//! swig(newobject, stable)
iftImage  *iftModaFilter(iftImage *img, iftAdjRel *A);
//! swig(newobject, stable)
iftImage  *iftSobelGradientMagnitude(iftImage *img);
//! swig(newobject, stable)
iftMatrix *iftImageToMatrix(iftImage *img, iftFastAdjRel *F, char crop);
//! swig(newobject, stable)
iftMatrix *iftImageAdjToMatrix(iftImage *img, iftAdjRel *A);
//! swig(newobject, stable)
iftImage  *iftMatrixToImage(iftMatrix *M, int xsize, int ysize, int zsize);
//! swig(newobject)
iftMatrix *iftFImageToMatrix(iftFImage *img, iftFastAdjRel *F, char crop);
//! swig(newobject)
iftFImage  *iftMatrixToFImage(iftMatrix *M, int xsize, int ysize, int zsize);
//! swig(newobject, stable)
iftMatrix *iftMImageToMatrix(iftMImage *img, iftFastAdjRel *F);
//! swig(newobject, stable)
iftMImage *iftMatrixToMImage(iftMatrix *M, int xsize, int ysize, int zsize, int nbands, char band_orientation);
//! swig(newobject, stable)
iftMImageArray *iftMatrixToMImageArray(iftMatrix *M, int nImgs, int xsize, int ysize, int zsize, int nbands, char band_orientation);
//! swig(newobject, stable)
iftMatrix *iftMImageToFeatureMatrix(iftMImage *mimg, iftAdjRel *A, float *fill_band_with);
//! swig(newobject, stable)
iftMatrix *iftMImageArrayToFeatureMatrix(iftMImageArray *mimgArray, iftAdjRel *A);
//! swig(newobject, stable)
iftMatrix *iftKernelToMatrix(iftKernel *K);
//! swig(newobject, stable)
iftMatrix *iftMKernelToMatrix(iftMKernel *K);
//! swig(newobject, stable)
iftMatrix *iftMMKernelToMatrix(iftMMKernel *k_bank);
//! swig(newobject, stable)
iftImage  *iftSmoothImage(iftImage *img, iftAdjRel *A, float sigma);
//! swig(newobject, stable)
iftImage  *iftSmoothImageInRegion(iftImage *img, iftImage *mask, iftAdjRel *A, float sigma);

void iftFastBilateralFilter2DAux(const float * __restrict__ _input,
                                 float       * __restrict__ _output,
                                 int   width,
                                 int   height,
                                 int 	channels,
                                 int   s_sigma,
                                 float r_sigma);

//! swig(newobject, stable)
iftImage *iftFastBilateralFilter2D(iftImage *img, int s_sigma, float r_sigma);
//! swig(newobject, stable)
iftMImage *iftFastMBilateralFilter2D(iftMImage *img, int s_sigma, float r_sigma);

//! swig(newobject, stable)
iftImage *iftNormalizeImage(iftImage *img, iftAdjRel *A, int Imax);
//! swig(newobject, stable)
iftImage *iftAlphaPooling(iftImage *img, iftAdjRel *A, int stride, float alpha);
//! swig(newobject, stable)
iftMatrix *iftRandomKernelBankAsMatrix(int size, int nbands, int nKernels);

/**
 * @brief Applies local batch centralization to a set of images and returns the fileset with the output images paths
 * 
 * @param fs Fileset with the input images
 * @param kernelSize Kernel size
 * @param outputDir Output directory
 * @param colSpace Color space
 * 
 * @return Fileset with the output images paths
 * 
 * @author Cesar Castelo
 * @date Mar 25, 2019
 */
iftFileSet *iftLocalBatchCentralization(iftFileSet *fs, int kernelSize, char *outputDir, char colSpace);

/**
 * @brief Computes the convolution between a mimage and a kernel bank and then applies a given stride
 * 
 * @param mimg Input multiband image
 * @param kernels Kernel bank (matrix)
 * @param stride Stride
 * 
 * @return Mimage with the resulting convolutions
 * 
 * @author Cesar Castelo
 * @date Apr 22, 2019
 */
//! swig(newobject, stable)
iftMImage *iftMConvolution(iftMImage *mimg, iftMatrix *kernels, int stride);

/**
 * @brief Computes the convolution between an array of mimages and a kernel bank and then applies a given stride
 * 
 * @param mimg Input multiband image array
 * @param kernels Kernel bank (matrix)
 * @param stride Stride
 * 
 * @return Mimage array with the resulting convolutions
 * 
 * @author Cesar Castelo
 * @date Apr 26, 2019
 */
//! swig(newobject, stable)
iftMImageArray *iftMConvolutionArray(iftMImageArray *mimgArray, iftMatrix *kernels, int stride);

/**
 * @brief Applies max pooling to a mimage given a patch size and a stride
 * 
 * @param mimg Input multiband image
 * @param radius of the pooling adjacency
 * @param stride Stride
 * 
 * @return Mimage with max-pooling applied
 * 
 * @author Cesar Castelo
 * @date Mar 25, 2019
 */
//! swig(newobject, stable)
iftMImage *iftMMaxPooling(iftMImage *mimg, float radius, int stride);

/**
 * @brief Applies max pooling to an array of mimages given a patch size and a stride
 * 
 * @param mimg Input multiband image array
 * @param radius of the pooling adjacency
 * @param stride Stride
 * 
 * @return Mimage array with max-pooling applied
 * 
 * @author Cesar Castelo
 * @date Apr 26, 2019
 */
//! swig(newobject, stable)
iftMImageArray *iftMMaxPoolingArray(iftMImageArray *mimgArray, float radius, int stride);

/**
 * @brief Applies max pooling to a mimage given an array of ROIs. The max value is computed inside each ROI
 * and all the pixels in that ROI get the computed max value. Finally a stride is applied to create the new image
 * 
 * @param mimg Input multiband image
 * @param roiArray Array of ROIs to apply the max-pooling operation
 * @param stride Stride
 * 
 * @return Mimage with max-pooling applied
 * 
 * @author Cesar Castelo
 * @date Mar 25, 2019
 */
//! swig(newobject, stable)
iftMImage *iftMMaxPoolingRoi(iftMImage *mimg, iftRoiArray *roiArray, int stride);

/**
 * @brief Applies max pooling to an array of mimages given an array of ROIs. The max value is computed inside each ROI
 * and all the pixels in that ROI get the computed max value. Finally a stride is applied to create the new image
 * 
 * @param mimg Input multiband image array
 * @param roiArray Array of ROIs to apply the max-pooling operation
 * @param stride Stride
 * 
 * @return Mimage array with max-pooling applied
 * 
 * @author Cesar Castelo
 * @date Apr 26, 2019
 */
//! swig(newobject, stable)
iftMImageArray *iftMMaxPoolingRoiArray(iftMImageArray *mimgArray, iftRoiArray *roiArray, int stride);

/**
 * @brief Applies min pooling to a mimage given a patch size and a stride
 * 
 * @param mimg Input multiband image
 * @param radius of the pooling adjacency
 * @param stride Stride
 * 
 * @return Mimage with min-pooling applied
 * 
 * @author Cesar Castelo
 * @date May 15, 2019
 */
//! swig(newobject, stable)
iftMImage *iftMMinPooling(iftMImage *mimg, float radius, int stride);

/**
 * @brief Applies min pooling to an array of mimages given a patch size and a stride
 * 
 * @param mimg Input multiband image array
 * @param radius of the pooling adjacency
 * @param stride Stride
 * 
 * @return Mimage array with min-pooling applied
 * 
 * @author Cesar Castelo
 * @date Apr 26, 2019
 */
//! swig(newobject, stable)
iftMImageArray *iftMMinPoolingArray(iftMImageArray *mimgArray, float radius, int stride);

/**
 * @brief Applies ReLU to a mimage.
 * 
 * @param mimg Input multiband image
 * 
 * @return Mimage with ReLU applied
 * 
 * @author Cesar Castelo
 * @date Apr 22, 2019
 */
//! swig(newobject, stable)
iftMImage *iftMReLU(iftMImage *mimg);

/**
 * @brief Applies ReLU to an array of mimages.
 * 
 * @param mimg Input multiband image array
 * 
 * @return Mimage array with ReLU applied
 * 
 * @author Cesar Castelo
 * @date Apr 26, 2019
 */
//! swig(newobject, stable)
iftMImageArray *iftMReLUArray(iftMImageArray *mimg);

#ifdef __cplusplus
}
#endif

#endif
