//
// Created by Cesar Castelo on Mar 21, 2019.
//

#ifndef IFT_ROI_H
#define IFT_ROI_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/imgproc/dtypes/VoxelArray.h"
#include "ift/imgproc/dtypes/BoundingBox.h"
#include "ift/core/io/Stream.h"
#include "iftImage.h"
#include "iftMImage.h"
#include "iftImageMath.h"

/**
 * @brief Region of interest.
 * @author Cesar Castelo
 * @date Mar 21, 2019
 * @ingroup Roi
 */
typedef struct ift_roi {
    /** Number of Elements */
    long n;
    /** Array of Voxels */
    iftVoxel *val;
    /** Label */
    int label;
} iftRoi;

/**
 * @brief Array of ROIs.
 * @author Cesar Castelo
 * @date Mar 21, 2019
 * @ingroup Roi
 */
typedef struct ift_roi_array {
    /** Number of Elements */
    long n;
    /** Array of ROIs */
    iftRoi **val;
} iftRoiArray;

/**
 * @brief Creates a region of interest
 * 
 * @param n Number of elements.
 * @return Region of interest
 * 
 * @author Cesar Castelo
 * @date Mar 21, 2019
 */
iftRoi *iftCreateRoi(long n);

/**
 * @brief Destroys a region of interest.
 * 
 * @param roi Region of interest
 * 
 * @author Cesar Castelo
 * @date Mar 21, 2019
 */
void iftDestroyRoi(iftRoi **roi);

/**
 * @brief Creates a copy of a region of interest
 * 
 * @param roi Region of interest
 * 
 * @author Cesar Castelo
 * @date Mar 22, 2019
 */
iftRoi *iftCopyRoi(iftRoi *roi);

/**
 * @brief Creates an array of ROIs
 * 
 * @param n Number of elements.
 * @return Array of ROIs
 * 
 * @author Cesar Castelo
 * @date Mar 21, 2019
 */
iftRoiArray *iftCreateRoiArray(long n);

/**
 * @brief Destroys an array of ROIs
 * 
 * @param arr Array of ROIs
 * 
 * @author Cesar Castelo
 * @date Mar 21, 2019
 */
void iftDestroyRoiArray(iftRoiArray **arr);

/**
 * @brief Creates a copy of an array of ROIs
 * 
 * @param arr Array of ROIs
 * 
 * @author Cesar Castelo
 * @date Mar 22, 2019
 */
iftRoiArray *iftCopyRoiArray(iftRoiArray *arr);

/**
 * @brief Writes an array of ROIs to the disk
 * 
 * @param arr Array of ROIs
 * @param filename Filename to be created
 * 
 * @author Cesar Castelo
 * @date Nov 29, 2019
 */
void iftWriteRoiArray(iftRoiArray *arr, char *filename);

/**
 * @brief Read an array of ROIs from the disk
 * 
 * @param filename Filename to be read
 * 
 * @author Cesar Castelo
 * @date Nov 29, 2019
 */
iftRoiArray *iftReadRoiArray(char *filename);

/**
 * @brief Returns the index of the biggest ROI in an array of ROIs
 * 
 * @param arr Array of ROIs
 * 
 * @author Cesar Castelo
 * @date Mar 22, 2019
 */
int iftBiggestRoiInArray(iftRoiArray *arr);

/**
 * @brief Creates a ROI from a bounding box
 * 
 * @param bb Bounding box
 * 
 * @author Cesar Castelo
 * @date Mar 22, 2019
 */
iftRoi *iftRoiFromBoundingBox(iftBoundingBox bb);

/**
 * @brief Creates an image that contains the minimum bounding box of the given ROI extracted from the image.
 * The minimum bounding box contains the original pixel values of the ROI and the background is black.
 * 
 * @param img Image
 * @param roi Region of interest
 * 
 * @author Cesar Castelo
 * @date Mar 22, 2019
 */
iftImage* iftExtractRoiNoBkgd(iftImage *img, iftRoi *roi);

/**
 * @brief Creates a mimage that contains the minimum bounding box of the given ROI extracted from the mimage.
 * The minimum bounding box contains the original pixel values of the ROI and the background is black.
 * 
 * @param mimg Mimage
 * @param roi Region of interest
 * 
 * @author Cesar Castelo
 * @date Jun 13, 2019
 */
iftMImage* iftMExtractRoiNoBkgd(iftMImage *mimg, iftRoi *roi);

/**
 * @brief Creates a binary mask from a ROI in an image
 * 
 * @param roi Region of interest
 * @param maskSize Size of the resulting mask image
 * 
 * @author Cesar Castelo
 * @date Mar 22, 2019
 */
iftImage *iftMaskFromRoi(iftRoi *roi, iftSize maskSize);

/**
 * @brief Extracts the center voxel of a ROI (center of gravity)
 * 
 * @param img Reference image
 * @param roi Region of interest
 * 
 * @author Cesar Castelo
 * @date Mar 22, 2019
 */
iftVoxel iftRoiCenterVoxel(iftImage *img, iftRoi *roi);

/**
 * @brief Defines a ROI with the pixels that belong to a given superpixel inside a superpixel labels image
 * 
 * @param spixLabels Image containing the tuperpixel labels
 * @param label Chosen superpixel label
 * 
 * @author Cesar Castelo
 * @date Mar 24, 2019
 */
iftRoi *iftRoiFromSuperpixel(iftImage *spixLabels, int label);

#ifdef __cplusplus
}
#endif

#endif //IFT_ROI_H
