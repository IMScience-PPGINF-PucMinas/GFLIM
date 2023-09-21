/**
 * @file iftRegistration.h
 * @brief Image Registration module.
 * @author Samuel Martins
 * @date Nov 5th, 2015
 * @ingroup Registration
 *
 * @note Examples of Image Registration can be found in ift/demo/Registration/iftRegisterImageByElastix.c and
 * ift/demo/Registration/iftRegisterImageSetByElastix.c.
 * note Examples of Image Registration given deformation fields can be found in ift/demo/Registration/iftTransformImageByTransformix.c and
 * ift/demo/Registration/iftTransformImageSetByTransformix.c
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef IFT_REGISTRATION_H
#define IFT_REGISTRATION_H

#include "iftCommon.h"
#include "iftImage.h"
#include "iftMatrix.h"
#include "iftRadiometric.h"
#include "iftRepresentation.h"
#include "iftSeeds.h"
#include "iftGeometric.h"


/**
 * @brief Computes a rotation matrix to align images based on the principal axes of a given object. The object is given in a binary image. 
 * @author Alexandre Falcao
 * @date Jan 19th, 2019
 * @ingroup Registration
 *
 * Computes a rotation matrix to align images by using
 * iftTransformImageByMatrix. The rotation matrix is computed based on
 * the principal axes of a given object. The object is given in a
 * binary image. As a collateral effect, images may be flipped due to
 * simmetry of the object points in two or more directions.
 * 
 *
 * @param binary_image: A binary image containing the reference object. 
 * @return A rotation matrix. 
 */
  iftMatrix *iftRotationMatrixToAlignByPrincipalAxes(iftImage *bin);

  
/**
 * @brief Centralizes a set of Label Images in a same domain. The new images are written in an output directory.
 * @author Samuel Martins
 * @date Oct 7, 2015
 * @ingroup Registration
 *
 * Centralizes a set of Label Images into a same domain. The new images are saved in the <output_dir>.\n
 * The resulting images have the same filenames from their original ones.\n\n
 * 
 * In order to centralize and put all label images into a same domain, the following steps are executed for each image:\n
 * 1. Gets the Minimum Bounding Box (MBB) and its Image Geometric Center, considereing all objects inside it;
 * 2. Centralize all MBB by their Geometric Centers
 * 3. The maximum bounding box, after the centralization, correspond to the resulting image domain of the
 * centralized label images
 *
 * @param label_imgs_paths A File Array with the pathnames from the label images to be centralized.
 * @param output_dir Pathname from the output directory where the new image will be saved/stored.
 *                   If NULL, a temporary directory is automatically created.
 * @return A file set with the paths from the centralized image domains. 
 */
iftFileSet *iftCentralizeLabelImages(const iftFileSet *label_imgs_paths, const char *output_dir);


void iftAlignImagesForSimilarity(char *dirIn, char *dirOut);


/**
 * @brief Register a cloud point in order to maximize the given score image.
 *
 * Given an array of points this function tries to apply affine transformations in order to better fit these points maximizing the score image.
 * The score image can be computed as the Inverted Euclidean Distance Transform as in iftEuclideanScoreImage(), or any other maximization score.
 *
 * @author Peixinho
 *
 * @date 12 Aug 2015
 *
 */
iftMatrix* iftShapeRegister(iftPoint* orig, int norig, iftImage* score);


/**
 * @brief Computes the inverse of Euclidean Distance Transform in a border Image.
 *
 * @param img The input binary border image.
 * @param decay The decay factor for the score, larger decay factors creates a score image that penalizes points too far from the border.
 *
 * @author Peixinho
 * @date 12 Aug 2015
 */
iftImage* iftEuclideanScoreImage(iftImage* img, float decay);
iftMatrix *iftShapeTransform(iftPoint *orig, int norig, float rx, float rz, iftImage* score);

/**
 * @brief Computes the Mean Square Error (MSE) of the fixed image and moved image.
 *
 * @param img1 Input registered image.
 * @param img2 The decay factor for the score, larger decay factors creates a score image that penalizes points too far from the border.
 *
 * @author Azael Sousa
 * @date 12 Aug 2015
 */
float iftRegistrationRMSE(iftImage *fixed_img, iftImage *moved_img);

#endif

#ifdef __cplusplus
}
#endif

