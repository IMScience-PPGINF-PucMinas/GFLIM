#ifndef IFT_BRAINMSP_H_
#define IFT_BRAINMSP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftImage.h"
#include "iftInterpolation.h"
#include "iftPlane.h"


/**
 * @brief Finds the Mid-Sagittal Plane (MSP) from a brain image and then aligns it based on the MSP.
 *
 * Algorithm from
 * [1] Ruppert, Guilherme CS, et al. "A new symmetry-based method for mid-sagittal plane extraction in neuroimages."
 * Biomedical Imaging: From Nano to Macro, 2011 IEEE International Symposium on. IEEE, 2011.
 *
 * @param img Brain image.
 * @param msp_out Return by reference the found msp. If NULL, nothing is returned.
 * @return Aligned brain image based on its found MSP. The MSP becomes the central sagittal slice from the image.
 *
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
//! swig(newobject)
iftImage *iftAlignBrainByMSP(const iftImage *img, iftPlane **msp_out);


/**
 * @brief Rotates a brain image based on a Mid-Sagittal Plane.
 *
 * @note It forces that the returned rotated image has the same domain and voxel sizes from the original ones.
 *
 * @param img Brain Image to be rotated.
 * @param msp Mid-Sagittal Plane.
 * @param interp_type Type of Interpolation. Use IFT_NEAREST_NEIGHBOR_INTERP to rotate label images, otherwise use IFT_LINEAR_INTERP.
 * @return Rotated Image.
 *
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
//! swig(newobject)
iftImage *iftRotateImageToMSP(const iftImage *img, const iftPlane *msp, iftInterpolationType interp_type);


/**
 * @brief Converts a MSP into a Labeled Image with the same resolution of a given reference image.
 *
 * @param msp Mid-Sagittal Plane.
 * @param ref_img Reference Image.
 * @return Labeled image with the voxels on the MSP set to 1.
 *
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
//! swig(newobject)
iftImage *iftMSPAsImage(const iftPlane *msp, const iftImage *ref_img);

#ifdef __cplusplus
}
#endif

#endif
