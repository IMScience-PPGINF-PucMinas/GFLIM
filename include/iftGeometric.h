#ifndef IFT_GEOMETRIC_H_
#define IFT_GEOMETRIC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftMatrix.h"
#include "iftInterpolation.h"
#include "iftDataSet.h"



/**
 * @brief Transforms/Maps an Image from a Transformation Matrix.
 * @author afalcao
 * @author Samuel Martins
 * @date Dec 9th, 2015
 *
 * @param img Image to be Transformed/Mapped.
 * @param M Transformation Matrix.
 * @return The Transformated/Mapped Image.
 */
iftImage *iftTransformImageByMatrix(const iftImage *img, const iftMatrix *M);
  
/**
 * @brief Rotates a Gray 2D image around the z-axis (tilt [0, 180])
 *
 *
 * @param img Image to be rotated.
 * @param theta Rotation angle on z-axis (tilt [0, 180)).
 * @return Rotated Image.
 *
 * @author Alexandre Falcao
 * @date Jan, 2015
 */
//! swig(newobject, stable)
iftImage *iftRotateImage2D(iftImage *img, float theta);

/**
 * @brief Align binary image by PCA
 *
 *
 * @param bin Binary image for alignment
 * @return Rotated Image.
 *
 * @author Alexandre Falcao
 * @date Feb, 2022
 */
//! swig(newobject, stable)
iftImage *iftAlignObject(iftImage *bin);

/**
 * @brief Scales a 2D image by the factors sx and sy for the x-axis and y-axis, respectively.
 * @author Samuka Martins
 * @date Aug 17, 2017
 */
iftImage *iftScaleImage2D(const iftImage *img, float sx, float sy);


/**
 * @brief Rotates a Gray 3D image on the x-axis (tilt [0, 180]) and y-axis (spin [0, 180]).
 *
 * @note The resulting rotated Image is cubic whose side has the size of the main diagonal of the input image.
 *
 * @param img Image to be rotated.
 * @param theta_x Rotation angle on x-axis (tilt [0, 180)).
 * @param theta_y Rotation angle on y-axis (spin [0, 180)).
 * @return Rotated Image.
 *
 * @author Samuel Martins
 * @date Nov 23, 2018
 */
//! swig(newobject, stable)
iftImage *iftRotateImage(const iftImage *img, float theta_x, float theta_y);


/**
 * @brief Scale an image by the factors sx, sy, and sz for the x-axis, y-axis, and z-axis, respectively.
 * @author Samuka Martins
 * @date Aug 17, 2017
 */
//! swig(newobject, stable)
iftImage *iftScaleImage(const iftImage *img, float sx, float sy, float sz);



/**
 * @brief Definition for a function that flips a given voxel <u> on some coordinate.
 * @param u Voxel to be flipped.
 * @param xsize Image's width used for flipping on X-Axis
 * @param ysize Image's height used for flipping on Y-Axis
 * @param zsize Image's depth used for flipping on Z-Axis
 * 
 * @author Samuka Martins
 * @date Sep 2, 2018
 */
typedef iftVoxel (*iftFlipVoxelFunc)(iftVoxel u, int xsize, int ysize, int zsize);


/**
 * @brief Flips a voxel <u> on X-Axis.
 * @note It has the same signature of the function definition iftFlipVoxelFunc
 * @author Samuka Martins
 * @date Sep 2, 2018
 */
iftVoxel iftFlipVoxelOnXAxis(iftVoxel u, int xsize, int ysize, int zsize);


/**
 * @brief Flips a voxel <u> on Y-Axis.
 * @note It has the same signature of the function definition iftFlipVoxelFunc
 * @author Samuka Martins
 * @date Sep 2, 2018
 */
iftVoxel iftFlipVoxelOnYAxis(iftVoxel u, int xsize, int ysize, int zsize);


/**
 * @brief Flips a voxel <u> on Z-Axis.
 * @note It has the same signature of the function definition iftFlipVoxelFunc
 * @author Samuka Martins
 * @date Sep 2, 2018
 */
iftVoxel iftFlipVoxelOnZAxis(iftVoxel u, int xsize, int ysize, int zsize);


/**
 * @brief Flip an Image on a given axis (IFT_AXIS_X, IFT_AXIS_Y, or AXIS_Z)
 * @author Samuel Martins
 * @date Mar 12, 2018
 */
//! swig(newobject, stable)
iftImage *iftFlipImage(const iftImage *img, iftAxis axis);


/**
 * @brief Flip a Multi-Band Image on a given axis (IFT_AXIS_X, IFT_AXIS_Y, or AXIS_Z)
 * @author Samuel Martins
 * @date Sep 3, 2018
 */
iftMImage *iftFlipMImage(const iftMImage *mimg, iftAxis axis);



iftImage  *iftTransformImageClipping(iftImage *img, iftMatrix *InvE, int xsize, int ysize, int zsize);
iftFImage *iftTransformFImageClipping(iftFImage *img, iftMatrix *InvE, int xsize, int ysize, int zsize);
iftImage  *iftTransformImageClip(iftImage *img, iftMatrix *M);

#ifdef __cplusplus
}
#endif

#endif
