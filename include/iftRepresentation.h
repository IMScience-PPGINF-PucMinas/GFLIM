#ifndef IFT_REPRESENTATION_H_
#define IFT_REPRESENTATION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftAdjacency.h"
#include "iftImage.h"
#include "iftRadiometric.h"
#include "iftFImage.h"
#include "iftSeeds.h"
#include "iftSegmentation.h"


/**
 * @brief It defines the API of functions to initialize the cost of a pixel with value <val> in EDT algorithm.
 * Check an example in function iftIntializeDistTransCostInterior in iftRepresentation.h
 *
 * @param  val Value to be compared.
 * @return     A given return.
 *
 * @author Samuka Martins
 * @date Dec 7, 2017
 */
typedef int (*iftIntializeDistTransCost)(int val);


/**
 * @brief Initialize the cost of a given pixel with value <val>, in order
 * to compute the EDT only in Object Interior.
 * If the val > 0 (object pixel), return IFT_INFINITY, otherwise return 0.
 *
 * @author Samuka Martins
 * @date Dec 1, 2017
 */
int iftIntializeDistTransCostInterior(int val);


/**
 * @brief Initialize the cost of a given pixel with value <val>, in order
 * to compute the EDT only in Background (Exterior of the Objects).
 * If the val == 0 (background pixel), return IFT_INFINITY, otherwise return 0.
 *
 * @author Samuka Martins
 * @date Dec 1, 2017
 */
int iftIntializeDistTransCostExterior(int val);

/**
 * @brief Initialize the cost of a given pixel to compute the EDT in both Background and Object pixels.
 * It returns IFT_INFINITY regardless of the input pixel value.
 *
 * @author Samuka Martins
 * @date Dec 1, 2017
 */
int iftIntializeDistTransCostBoth(int val);


/**
 * @brief Gets the Initialization EDT Cost Function for a given side (IFT_INTERIOR, IFT_EXTERIOR, IFT_BOTH)
 * @author Samuka Martins
 * @date Dec 1, 2017
 */
iftIntializeDistTransCost iftGetIntializeDistTransCost(iftSide side);

//! swig(newobject, stable)
iftFImage *iftGeodesicDistTrans(const iftSet *S, const iftImage *mask, const iftAdjRel *A);

//! swig(newobject, stable)
iftImage  *iftBorderDistTrans(const iftImage *label, iftAdjRel *A);
//! swig(newobject, stable)
  iftImage  *iftSetDistTrans(iftSet *S, int xsize, int ysize, int zsize);
  
iftImage       *iftShellDistTrans(iftImage *bin, iftAdjRel *A, char side, float max_dist);
void            iftDistTransRootMap(iftImage *bin, iftAdjRel *A, char side, iftImage **dist, iftImage **root);
iftFImage      *iftSignedDistTrans(iftImage *bin, iftAdjRel *A);
iftFImage      *iftShellSignedDistTrans(iftImage *bin, iftAdjRel *A,float max_dist);

//! swig(newobject, stable)
iftFImage      *iftMSSkel(iftImage *bin); /* multiscale skeletons by geodesic length */

  /* compute surface skeleton by thresholding the geodesic skeleton
     and then select a given number of components */

  iftImage *iftSurfaceSkeleton(iftFImage *skel, float thres, int number_of_components);

  // This function computes the multiscale skeleton and returns the distance transformed
  // used for computing it. Both maps can be used to compute a Medial Axis Transform.
iftFImage *iftMSSkel2DDistMap(iftImage *label_img, iftAdjRel *A, iftSide side, iftImage **dist_out,
                              iftImage **relabel_img_out);

  // This function should be used when the euclidean distance to the skeleton is irrelevant.
  // It is a wrapper for function iftMSSkel2DDistMap
//! swig(newobject, stable)
iftFImage *iftMSSkel2D(iftImage *label_img, iftAdjRel *A, iftSide side, iftImage **dist_out,
                       iftImage **relabel_img_out);

iftImage       *iftIntMSSkel2D(iftImage *bin, iftAdjRel *A, iftSide side);
  // This function computes the border of all labels in the label image and then uses the
  // given set of border pixels to compute a distance transform from them.
  // These functions are a superset of the binary case.
  void  	  iftMultiLabelDistTransFromBorders(iftImage *label, iftAdjRel *A, char side, iftImage **dist, iftImage **root);

  /**
   * @brief Computes the distance transform for a multi-label image from a given set S.
   * S is usually composed of the boundary pixels for all labels. This function generalizes
   * iftDistTrans and iftShellDistTrans since it accepts a maximum distance parameter.
   *
   * @author Thiago Vallin Spina
   *
   * @param S An iftSet representing the object's border or any set of points from which the EDT will be computed.
   * @param label A label image used to determine the EDT's size.
   * @param A An adjacency relation used to compute the EDT.
   * @param side Determines whether the EDT will be computed inside, outside, or on both sides of the union of labels.
   * @param max_dist The maximum distance until which the EDT will be computed.
   * @param dist Returns the computed EDT.
   * @param root Returns the root value in S for each spel.
   */
  void            iftMultiLabelShellDistTransFromSet(iftSet *S, iftImage *label, iftAdjRel *A, char side, double max_dist, iftImage **dist, iftImage **root);


  void 		  iftMultiLabelDistTransFromSet(iftSet *S, iftImage *label, iftAdjRel *A,	char side, iftImage **dist, iftImage **root);

/**
 * @brief Computes the Euclidean Distance Transform from a Labeled Image (binary or multi-label)
 * with root, label, and path propagation.
 * The Labeled Image is relabed, since there could have different objects (connected components) with the same label.
 * They are returned if the passed pointers are different from NULL.
 *
 * @param  label_img     Input Labeled Image (binary or multi-label).
 * @param  Ain           Input adjacent used to find the borders. If NULL, a 4-neighborhood (2D) or
 *                       6-neighborhood (3D) is considered.
 * @param  side          Side for EDT computation (IFT_INTERIOR, IFT_EXTERIOR, IFT_BOTH).
 * @param  root_out      If != NULL, returns the propagated root.
 * @param  edt_label_out If != NULL, returns the propagated labeled image.
 * @param  pred_out      If != NULL, returns the propagated predecessor map.
 * @return               The distance map after EDT computation.
 */
//! swig(newobject, stable)
iftImage *iftEuclDistTrans(const iftImage *label_img, const iftAdjRel *Ain, iftSide side, iftImage **root_out,
                           iftImage **edt_label_out, iftImage **pred_out);

/**
 *
 * @param mask          mask delimiting region
 * @param set           set to compute de edt from
 * @param root_out      If != NULL, returns the propagated root.
 * @return
 */
//! swig(newobject, stable)
iftImage *iftGeodesicDistTransFromSet(const iftImage *mask, const iftSet *set, iftImage **root_out);

  /**
   * @brief Computes the number of descendants each node of the forest
   * has in a target set, as represented by a binary image.
   *
   * @author Alexandre Falcao
   *
   * @param pred:  input spanning forest.
   * @param target_bin: binary image that represents the target set.
   */

  iftImage *iftNumberOfDescendantsInTarget(iftImage *pred, iftImage *target_bin);
  /**
   * @brief Computes the number of ascendants each node has in the forest.
   *
   * @author Alexandre Falcao
   *
   * @param pred:  input spanning forest.
   */

  iftImage *iftNumberOfAscendants(iftImage *pred);

  /**
   * @brief Computes the number of descendants each node has in the forest.
   *
   * @author Alexandre Falcao
   *
   * @param pred:  input spanning forest.
   */

  iftImage *iftNumberOfDescendants(iftImage *pred);

/**
 * @brief Signed version of function iftMultiLabelDistTransFromBorders.
 *
 * @author Thiago Vallin Spina
 * @date Feb 15, 2016
 *
 * @param label Label image.
 * @param A Adjacency relation to establish local neighborhood.
 * @param dist Signed euclidean distance transform result.
 * @param root The root map. (May be NULL if not desired).
 *
 */
  iftFImage     *iftMultiLabelSignedDistTransFromBorders(iftImage *label, iftAdjRel *A, char side, iftImage **root);
  void           iftLabelRootPropagation(iftImage *bin, iftAdjRel *A, char side, iftImage **root, iftImage **label, iftImage **dist);
  iftImage       *iftRootPropagation(iftImage *bin, iftAdjRel *A, char side, float max_dist);

/**
 * @brief Compute the Geodesic length for each pixel in the objects' contour.
 * It respects the label from the input Image to find the objects' borders,
 * so that two adjacent objects with different labels in the input image will also be
 * different objects in the labeled objects' contours.
 *
 * @param  label_img Labeled Image (binary or multi-label)
 * @return           Geodesic length Map.
 */
iftFImage *iftLabelContourPixelByGeoLen(iftImage *label_img);

iftImage       *iftLifmssketImage(iftImage *img);
iftImage       *iftDropImage(iftImage *bin);
iftFImage      *iftIntegralImage(iftImage *img);
float           iftGetIntegralValueInRegion(const iftFImage *integ, iftVoxel *v, int npts);
iftImage       *iftMarkGeometricCenters(iftImage *bin);

iftImage       *iftComponentSizes(iftImage *bin, iftAdjRel *A);
iftImage       *iftBorderSize(iftImage *bin);

void iftObjectAreaFromPixel(const iftImage* comp, int index, iftImage *area_img);

/**
 * @brief Runs the medial axis transform in 2D, such that the output
 * image is a 2D skeleton image with the squared Euclidean distance
 * from the skeleton points to the objects' boundaries.
 * @author Alexandre Falcao
 * @date   January 14th, 2015
 * @param  binary image with one or multiple objects.
 * @param  scale threshold for binarization of the multiscale skeleton.
 * @param  object side in which the skeleton will be computed: IFT_INTERIOR, IFT_EXTERIOR, or IFT_BOTH.
 * @return The skeleton image whose values are the squared EDT values for shape reconstruction.
 *
 */

//! swig(newobject, stable)
iftImage       *iftMedialAxisTrans2D(iftImage *bin, float scale_thres, iftSide side);

/**
 * @brief Reconstructs the shape from the medial axis image. It is a
 * shape smoothing filter.
 * @author Alexandre Falcao
 * @date   January 14th, 2015
 * @param  medial axis image
 * @param  reconstruction value
 * @return Binary image with the reconstructed shapes.
 *
 */

//! swig(newobject, stable)
iftImage       *iftShapeReconstruction(iftImage *medial_axis, int value);


/**
 * @brief Terminal points of a binary skeleton.
 * @author Alexandre Falcao
 * @date   January 14th, 2015
 * @param  Binary skeleton.
 * @return Binary image with the terminal points of the input skeleton.
 *
 */

//! swig(newobject, stable)
iftImage       *iftTerminalPoints2D(iftImage *skel);


/**
 * @brief  Branch points of a binary skeleton.
 * @author Alexandre Falcao
 * @date   Feb, 2022
 * @param  Binary skeleton.
 * @return Binary image with the branch points of the input skeleton.
 *
 */

//! swig(newobject, stable)
iftImage *iftBranchPoints2D(iftImage *skel);
  
/**
 * @brief  Extract the set of terminal points of a binary skeleton.
 * @author Alexandre Falcao
 * @date   January 14th, 2015
 * @param  Binary skeleton.
 * @return Set of terminal points of the input skeleton.
 *
 */

//! swig(newobject, stable)
iftSet         *iftTerminalPointSet2D(iftImage *skel);

float iftWeightedIntegralValue(const iftFImage *int_img, iftVoxel center, int width, int height, float weight);

iftMImage *iftExtractRegionOfInterest(iftMImage *morig, iftMImage *mlabel, float radius);

//! swig(newobject)
iftVoxelArray *iftContourToArray(const iftImage *contour);

//! swig(newobject)
iftVoxelArray *iftApproxContour(const iftImage *contour, double epsilon);

//! swig(newobject)
iftVoxelArray *iftApproxVoxelArray(const iftVoxelArray *array, double epsilon);

//! swig(newobject)
iftSet *iftVoxelArrayToSet(const iftImage *image, const iftVoxelArray *array);

//! swig(newobject)
iftSet *iftNearestInContour(const iftImage *contour, const iftSet *set, const iftImage *mask);


#ifdef __cplusplus
}
#endif

#endif
