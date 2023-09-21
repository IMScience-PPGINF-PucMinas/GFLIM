#ifndef IFT_RADIOMETRIC_H_
#define IFT_RADIOMETRIC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/Color.h"
#include "iftCommon.h"
#include "iftImage.h"
#include "iftFImage.h"

/* Radiometric Transformations */

//! swig(newobject, stable)
iftImage *iftLinearStretch(iftImage *img, double f1, double f2, double g1, double g2);
iftFImage *iftFLinearStretch(iftFImage *img, double f1, double f2, double g1, double g2);


/**
 * @brief Normalize an Image by MinMax normalization only inside a Region of Interest to the
 * range [minval, maxval].
 * 
 * @param  img    Input image to be normalized.
 * @param  region Region of Interest. If NULL, the entire image is considered.
 * @param  minval Inferior bound of the normalization.
 * @param  maxval nferior bound of the normalization.
 * @return        Normalized Image.
 *
 * @author Samuka Martins
 * @date Dec 21, 2017
 */
iftImage *iftNormalizeInRegion(const iftImage *img, const iftImage *region, double minval, double maxval);


/**
 * @brief Standardize the luminance histogram of gray and color images (i.e., it preserves matiz and saturation) to 
 * be a Gaussian with desired mean and standard deviation. 
 * 
 * @param  orig   input image to be standardized.
 * @param  mean   desired mean.
 * @param  stdev  desired standard deviation.
 * @return        standardized image
 *
 * @author Alexandre Falcao
 * @date March, 9th 2022
 */
iftImage *iftGaussianStandardization(iftImage *orig, float mean, float stdev);
  
/**
 * @brief Normalize an Image by MinMax normalization to the range [minval, maxval].
 * 
 * @param  img    Input image to be normalized.
 * @param  minval Inferior bound of the normalization.
 * @param  maxval nferior bound of the normalization.
 * @return        Normalized Image.
 *
 * @author Samuka Martins
 * @date Dec 21, 2017
 */
//! swig(newobject, stable);
iftImage *iftNormalize(const iftImage *img, double minval, double maxval);


/**
 * @brief Apply divisive normalization to an image. For each pixel it divides its value by the sum of all the
 * squared values of the adjacent pixels
 * 
 * @param img           Input image to be normalized.
 * @param A             Adjacency to be used to apply the normalization
 * @param minValNorm    Minimum value for the normalized image
 * @param maxValNorm    Maximum value for the normalized image
 * @return              Normalized Image.
 *
 * @author Cesar Castelo
 * @date Set 10, 2018
 */
iftImage *iftImageDivisiveNormalization(iftImage *img, iftAdjRel *A, int minValNorm, int maxValNorm);


/**
 * @brief Normalize an Image only inside a Region of Interest to the range [minval, maxval]
 * ignoring voxel outliers with high brightness.
 *
 * It computes the normalized accumulated histogram to figure out the percentage of voxels from 0 until 
 * each voxel v only inside the Region of Interest. When attaining a given percentage perc,
 * it will ignore the brightness of the remaining voxels (outliers)
 * 
 * @param  img    Input Image
 * @param  region Region of Interest. If NULL, the entire image is considered.
 * @param  minval Inferior bound of the normalization.
 * @param  maxval Inferior bound of the normalization.
 * @param  perc   Percentage of true voxels to be considered in normalization (Suggestion: 0.98)
 * @return        Normalized Image.
 *
 * @author Samuka Martins
 * @date Jun 20, 2017
 */
//! swig(newobject)
iftImage *iftNormalizeWithNoOutliersInRegion(const iftImage *img, const iftImage *region,
                                            int minval, int maxval, float perc);


/**
 * @brief Normalize an Image to the range [minval, maxval] ignoring voxel outliers with high brightness.
 *
 * It computes the normalized accumulated histogram to figure out the percentage of voxels from 0 until 
 * each voxel v. When attaining a given percentage perc, it will ignore the brightness of the 
 * remaining voxels (outliers)
 * 
 * @param  img    Input Image
 * @param  minval Inferior bound of the normalization.
 * @param  maxval Inferior bound of the normalization.
 * @param  perc   Percentage of true voxels to be considered in normalization (Suggestion: 0.99)
 * @return        Normalized Image.
 *
 * @author Samuka Martins
 * @date Jun 20, 2017
 */
//! swig(newobject)
iftImage *iftNormalizeWithNoOutliers(const iftImage *img, int minval, int maxval, float perc);


iftImage *iftWindowAndLevel(iftImage *img, int width, int level, int maxval);
iftImage *iftGaussianStretch(iftImage *img, double mean, double stdev, int maxval);
iftImage *iftExponenStretch(iftImage *img, double f1, double f2, double g1, double g2);

/**
 * @brief Equalizes (uniformly) the histogram of an image by the sorting method.
 *
 * Input image will be equalized uniformly in the range of [0, max_value].
 * 
 * 1) Sort all pixels from <img> in ascending order of brightness; \n
 * 2) Split the sorted pixels in (max_val + 1) buckets, k = 0, 1, ..., max_val, with the same number of voxels; \n
 * 3) Equalized Image <out>, such that out[p] = the bucket k where the voxel p has its brightness img[p] mapped to. \n\n
 *
 * POSSIBLE PROBLEM:
 * If the passed max_val <= 0, it will be considered the maximum value of the image.
 * Noises with high brightness values can result in a wrong equalization, stretching the brightness range incorrectly.
 * 
 * @param  img     Image to be equalized.
 * @param  max_val Maximum value of the output image. If <= 0, it is considered the maximum value of the image. 
 * @return         Equalized Image.
 *
 * @author Samuka
 * @date Jan 3, 2017
 */
iftImage *iftEqualize(const iftImage *img, int max_val);


/**
 * @brief Matches the histogram of an image <img> with a reference image <ref> according to
 * the Traditional Algorithm.
 *
 * This traditional algorithm works as follows:
 * 1) builds the normalized accumulated histograms for the images (himg and href)
 * 2) For each pixel p in img
 *    tries to find the bucket of href with value exactly equal to himg[img->val[p]]
 *    this will be value of p in the matched output image. 
 *
 * A mask for the image and another one for the reference image can be passed to limitate
 * the histogram computation for them.
 * 
 * OBS: Equal images (img == ref) does not necessarily return an exact copy of ref.
 * A set of voxels can have values slightly different from the original ones, because it gets
 * the first ref img's hist bucket with the same value of the img's hist.
 * Because it is an accumulated histogram, different sequential buckets can have the same value.
 * 
 * @param  img Image whose histogram will be matched.
 * @param  img_mask Mask that defines the histogram computation for the Image. Use NULL to consider the entire image.
 * @param  ref Reference image.
 * @param  ref_mask Mask that defines the histogram computation for the Reference Image. Use NULL to consider the entire image.
 * @return     Resulting Image with Matched Histogram.
 *
 * @author Samuka
 * @date Jan 3, 2017
 */
//! swig(newobject)
iftImage *iftMatchHistogram(const iftImage *img, const iftImage *img_mask,
                            const iftImage *ref, const iftImage *ref_mask);


/* Radiometric Resolution */
int 	  iftRadiometricResolution(iftImage *img);


#ifdef __cplusplus
}
#endif

#endif
