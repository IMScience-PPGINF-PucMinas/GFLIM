#ifndef IFT_INTERPOLATION_H_
#define IFT_INTERPOLATION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "ift/core/dtypes/Color.h"
#include "iftImage.h"
#include "iftFImage.h"
#include "iftMImage.h"
#include "iftDataSet.h"
#include "iftAdjacency.h"
#include "iftPlane.h"
#include "iftMatrix.h"
#include "iftSegmentation.h"
#include "iftRepresentation.h"
#include "iftMSPS.h"


/**
 * @brief Types of Interpolations. The flag IFT_LINEAR_INTERP must be used for BILINEAR or TRILINEAR interpolation type.
 * @author Samuel Martins
 * @date Nov 24, 2018
 */
//! swig()
typedef enum {
  IFT_NEAREST_NEIGHBOR_INTERPOLATION, IFT_LINEAR_INTERPOLATION
} iftInterpolationType;


/**
 * @brief Definition/template for a function that interpolates a point P inside the image img.
 *
 * One can use the iftInterpolationType to assign a given interpolation function to this function pointer (template).
 *
 * @param img Image used for interpolation.
 * @param P point to be interpolated.
 * @return Interpolated intensity.
 *
 * @author Samuka Martins
 * @date Nov 25, 2018
 */
typedef int (*iftInterpPointFunc)(const iftImage *img, iftPoint P);



// Returns the coordinate on the border of the image for the given
// dimension that is closest to the original coordinate
static inline int iftBorderInterpolateCoord(int dim_size, int p) {
  return (p < 0) ? 0 : ((p >= dim_size)? dim_size - 1 : p);
}

// Interpolates a voxel outside the float image's domain to a border
// voxel
static inline iftVoxel iftBorderInterpolateVoxel(iftImage *img, iftVoxel v) {
  v.x = iftBorderInterpolateCoord(img->xsize, v.x);
  v.y = iftBorderInterpolateCoord(img->ysize, v.y);
  v.z = iftBorderInterpolateCoord(img->zsize, v.z);

  return v;
}

// Interpolates a voxel outside the float image's domain to a border
// voxel
static inline iftVoxel iftFBorderInterpolateVoxel(iftFImage *img, iftVoxel v) {
  v.x = iftBorderInterpolateCoord(img->xsize, v.x);
  v.y = iftBorderInterpolateCoord(img->ysize, v.y);
  v.z = iftBorderInterpolateCoord(img->zsize, v.z);

  return v;
}


// Interpolates a voxel outside the multi-band image's domain to a border
// voxel
static inline iftVoxel iftMBorderInterpolateVoxel(iftMImage *img, iftVoxel v) {
  v.x = iftBorderInterpolateCoord(img->xsize, v.x);
  v.y = iftBorderInterpolateCoord(img->ysize, v.y);
  v.z = iftBorderInterpolateCoord(img->zsize, v.z);

  return v;
}

iftPlane *iftFindBestCutPlane(iftImage *weight, iftPoint pos, int xviewsize, int yviewsize);
iftPlane *iftFindBestObjectCutPlane(iftImage *obj, iftImage *weight);

/**
 * @brief Reslice a 3D image by following a specified order among axes
 * x, y, and z, and with increments dx, dy, dz from the corresponding
 * octant, which is found from dx, dy, and dz values. The input 3D
 * image starts from octant (0,0,0) and its voxels are accessed from
 * that octant by following the order XYZ with increments
 * dx=dy=dz=1. In order to reslice it from octant
 * (xsize-1,ysize-1,zsize-1), for instance, the increments must be
 * specified as dx=dy=dz=-1. If we also wish to change the voxel
 * access order to z, y, and x, then the axis order must be ZYX.
 * @author Alexandre Falcao 
 * @date Feb, 23rd 2016
 * @ingroup Interpolation
 * @param  img: A 3D input image.
 * @param  dx : increment in {-1,1} along the axis x.
 * @param  dy : increment in {-1,1} along the axis y.
 * @param  dz : increment in {-1,1} along the axis z.
 * @param  axis_order : new order among axes from {XYZ, XZY, YXZ, YZX, ZXY, ZYX}
 * @return Image that results from the reslicing of the input image
 * with increments dx, dy, dz from the corresponding octant and with a
 * given axis order.
 */

iftImage  *iftResliceImage(iftImage *img, iftPlane *pl, int xsize, int ysize, int zsize);
iftImage *iftResliceImageSimple(const iftImage *img, int dx, int dy, int dz, iftAxisOrder axis_order);
iftImage  *iftGetSlice(iftImage *img, iftPlane *pl, int xsize, int ysize);

/**
 * @brief Interpolate an image (color or gray) by using the Nearest Neighbor Algorithm.
 * @param  img Image to be interpolated.
 * @param  sx  Scale Factor in x-axis.
 * @param  sy  Scale Factor in y-axis.
 * @param  sz  Scale Factor in z-axis.
 * @return     Interpolated Image
 */
//! swig(newobject, stable)
iftImage *iftInterpByNearestNeighbor(const iftImage *img, float sx, float sy, float sz);

/**
 * @brief Interpolate a 2D image (color or gray) by using the Nearest Neighbor Algorithm.
 * @param  img Image to be interpolated.
 * @param  sx  Scale Factor in x-axis.
 * @param  sy  Scale Factor in y-axis.
 * @return     Interpolated Image
 */
//! swig(newobject, stable)
iftImage *iftInterp2DByNearestNeighbor(const iftImage *img, float sx, float sy);


/**
 * @brief Resize a 3D image (color or gray) by using linear interpolation 
 * @param  img Image to be interpolated.
 * @param  sx  Scale Factor in x-axis.
 * @param  sy  Scale Factor in y-axis.
 * @param  sz  Scale Factor in z-axis.
 * @return     Interpolated Image
 */
//! swig(newobject, stable)
iftImage  *iftInterp(const iftImage *img, float sx, float sy, float sz);

/**
 * @brief Resize a 2D image (color or gray) by using linear interpolation 
 * @param  img Image to be interpolated.
 * @param  sx  Scale Factor in x-axis.
 * @param  sy  Scale Factor in y-axis.
 * @return     Interpolated Image
 */
//! swig(newobject, stable)
iftImage *iftInterp2D(const iftImage *img, float sx, float sy);

/**
 * @brief Resize a 2D/3D image (color or gray) by using linear interpolation 
 * @param  img Image to be interpolated.
 * @param  xsize  desired size in x-axis.
 * @param  ysize  desired size in y-axis.
 * @param  zsize  desired size in z-axis.
 * @return     Interpolated Image
 */
//! swig(newobject, stable)
iftImage* iftResizeImage(const iftImage *img, int xsize, int ysize, int zsize);

iftFImage  *iftFGetSlice(iftFImage *img, iftPlane *pl, int xsize, int ysize);
iftFImage *iftFInterp2D(iftFImage *img, float sx, float sy);
iftMImage *iftMInterp2D(iftMImage *mimg, float sx, float sy);
iftMImage *iftMInterp(iftMImage *img, float sx, float sy, float sz);
iftMImage* iftResizeMImage(iftMImage* img, int xsize, int ysize, int zsize);
iftFImage *iftFInterp(iftFImage *img, float sx, float sy, float sz);
iftImage  *iftShapeBasedInterp2D(iftImage *label, float sx, float sy);
iftImage  *iftShapeBasedInterp(iftImage *label, float sx, float sy, float sz);
iftImage  *iftResliceOnPrincipalAxis(iftImage *img, iftImage *bin);

/**
 * @brief Reslice an image by using a 4X4 Transformation Matrix R on homogeneous coordinates.
 *
 * R must only be a Rotation Matrix on homogeneous coordinates with only rotations of 0, 90, 270,
 * ie, only values -1, 0, +1 and only one value != 0 (-1, +1) per column and per row.
 * Otherwise, the reslice will not work correctly.
 *  
 * @param  img Image to be resliced.
 * @param  R   4x4 Rotation Matrix.
 * @return     Resliced Image.
 *
 * @author Samuka Martins
 * @date Apr 22, 2018
 */
iftImage *iftResliceImageByTransMatrix(const iftImage *img, const iftMatrix *R);



/**
 * @brief Gets the interpolated voxel value in a given (float) point by trilinear interpolation.
 *
 * @attention If the point is out of the image domain, it returns 0.
 *
 * @param img Target image.
 * @param P Point to be interpolated.
 * @return Interpolated intensity.
 *
 * @author Samuel Martins
 * @date Nov 24, 2018
 */
int iftImageValueAtPoint(const iftImage *img, iftPoint P);


/**
 * @brief Gets the 3D image interpolated pixel value in a given position according to its nearest neighbor.
 *
 * If the point is out of the image domain, it returns 0.
 *
 * @attention It considers that the point P is inside the image's domain. There is no validation for this.
 *
 * @param img Target image.
 * @param P Point to be interpolated.
 * @return Interpolated intensity.
 *
 * @author Samuel Martins
 * @date Nov 24, 2018
 */
int iftImageValueAtPointNearestNeighbor(const iftImage *img, iftPoint P);



/**
 * @brief Gets the 2D image interpolated pixel value in a given position.
 *
 * @warning This function is for 2D images.
 *
 * @param img Target image.
 * @param p position to be interpolated.
 * @return Interpolated intensity.
 *
 * @author Samuel Martins
 * @date Jan 21, 2019
 */
int iftImageValueAtPoint2D(const iftImage *img, iftPoint P);


/**
 * @brief Gets the 2D image interpolated pixel value in a given position according to its nearest neighbor.
 *
 * @warning This function is for 2D images.
 *
 * @param img Target image.
 * @param p position to be interpolated.
 * @return Interpolated intensity.
 */
int iftImageValueAtPoint2DNearestNeighbor(const iftImage *img, iftPoint P);



#ifdef __cplusplus
}
#endif

#endif

