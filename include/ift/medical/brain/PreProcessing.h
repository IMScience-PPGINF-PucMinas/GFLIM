#ifndef IFT_BRAIN_PREPROCESSING_H_
#define IFT_BRAIN_PREPROCESSING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftImage.h"
#include "iftPlane.h"



/**
 * @brief Set the environment to run ANTs programs (available in ift_dir/externals/ANTs)
 * @note Supported by Linux and Mac OS
 *
 * @author Samuka Martins
 * @date Jul 18, 2017
 */
void iftSetANTsEnvironment(void);



/**
 * @brief Apply the N4 algorithm for Inhomogeneity Correct in a MRI Image.
 *
 * It is used the default parameters, suggested by 3D Slicer tool [2], for N4 [1].\n
 * [1] Tustison, Nicholas J., et al. \"N4ITK: improved N3 bias correction.\" IEEE transactions on medical imaging 29.6 (2010): 1310-1320.\n
 *
 * @param  img           Image to be corrected.
 * @param  mask          Binary mask that defines the structure of your interest. If NULL, N4 considers the entire image.
 * @param  shrink_factor Image Resampling/Shrinking Factor (1, 2, 3, 4, ...) to decrease the computation time.
 *                       Shrink factors <= 4 are commonly used.
 * @param  out_bias      The resulting Bias Field Image by the N4 correction. If NULL, it is not considered.
 * @return               Corrected image by N4.
 *
 * @author Samuka Martins
 * @date Mar 13, 2017
 */
iftImage *iftN4BiasFieldCorrection(const iftImage *img, const iftImage *mask, int shrink_factor,
                                   iftImage **out_bias);


/**
 * @brief A good scheme to apply the N4 over the input image in a more effective way.
 *
 * Firtly, N4 is applied with shrink factor 4 and the default arguments.
 * Then, another N4 is applied with shrink factor 6 and the default arguments, in order to make
 * a fine tunning in the correction.
 *
 * @param  img      Image to be corrected.
 * @param  mask     Binary mask that defines the structure of your interest. If NULL, N4 considers the entire image.
 * @param  out_bias The resulting Bias Field Image by the N4 correction. If NULL, it is not considered.
 * @return          Corrected image by N4.
 *
 * @author Samuka Martins
 * @date Mar 13, 2017
 */
iftImage *iftEffectiveN4BiasFieldCorrection(const iftImage *img, const iftImage *mask, iftImage **out_bias);


/**
 * @brief Apply a set of Pre-Processing operation on to a Brain MR Image.
 *
 * The pre-processing consists of the following ordered steps:
 * 1 - Bias Field Correction by N4;
 * 2 - Median Filter using a Spherical Adjacency of Radius 1;
 * 3 - Normalization within [0, 2ˆnbits - 1] (nbits is passed);
 * 4 - Extraction and image's alignment with the Mid-Sagittal Plain (MSP). A MSP previously extracted can be passed;
 * 5 - Histogram Matching with a given Reference Image (a binary mask can be passed to define the regions for Hist. Matching)
 *
 * Obs: If a valid MSP (!= NULL) is passed, the function uses it to align the image. Otherwise, the function extracts the MSP.
 *
 * @param mri Brain MR Image
 * @param nbits Number of Bits for normalization within [0, 2ˆnbits - 1]. If nbits <= 0, no normalization is done.
 * @param msp_in Mid-Sagittal Plane (previously extracted) used to align the image's alignment.
 * @param mask Binary mask defined for the input mri image used to restrict the Histogram Matching.
 * @param ref_mri Reference Brain MR Image used for histogram matching.
 * @param ref_mask Binary mask defined for the reference image used to restrict the Histogram Matching.
 * @param skip_n4 Skip the N4 Bias Field Correction.
 * @param skip_median_filter Skip the Median Filter.
 * @param skip_msp_alignment Skip the MSP alignment.
 * @param skip_hist_matching Skip the Histogram Matching.
 * @param msp_out Reference to return the MSP for the alignment. If NULL, nothing is returned.
 * @return Pre-processed MR Brain Image.
 *
 * @author Samuka
 * @date Dec 20, 2017
 */
//! swig(newobject)
iftImage *iftBrainMRIPreProcessing(const iftImage *mri, int nbits, const iftPlane *msp_in, const iftImage *mask,
                                   const iftImage *ref_mri, const iftImage *ref_mask, bool skip_n4,
                                   bool skip_median_filter, bool skip_msp_alignment, bool skip_hist_matching,
                                   iftPlane **msp_out);


#ifdef __cplusplus
}
#endif

#endif
