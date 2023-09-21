/**
 * @file iftSegmentation.h
 * @brief Definitions and functions about Segmentation Methods
 * @author
 * @date
 * @ingroup Segmentation
 */

#ifndef IFT_SEGMENTATION_H_
#define IFT_SEGMENTATION_H_

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BMap.h"
#include "ift/imgproc/dtypes/BoundingBox.h"

#include "iftAdjacency.h"
#include "iftClassification.h"
#include "iftClustering.h"
#include "iftCommon.h"
#include "iftDataSet.h"
#include "iftFiltering.h"
#include "iftFImage.h"
#include "iftImage.h"
#include "iftImageForest.h"
#include "iftImageMath.h"
#include "iftMathMorph.h"
#include "iftMetrics.h"
#include "iftMImage.h"
#include "iftRadiometric.h"
#include "iftRepresentation.h"
#include "iftSeeds.h"

//struct iftDataSet; // declaration added to avoid headers cross-reference


/**
 * @brief Enumeration to indicate different kinds of Gradient Algorithms.
 * @author Samuel Martins
 * @date Aug 22, 2016
 * @ingroup Segmentation
 */
typedef enum {
    IFT_NONE_GRAD, IFT_IMAGE_BASINS, IFT_BRAIN_GRAD, IFT_IMAGE_GRAD_MAGNITUDE
} iftGradientAlg;

/**
 * @brief Enumeration to indicate different kinds of mask join operations
 * @author Cesar Castelo
 * @date Set 17, 2018
 * @ingroup Segmentation
 */
typedef enum {
    IFT_MASK_JOIN_OPERATION_UNION, IFT_MASK_JOIN_OPERATION_INTERSECTION
} iftMaskJoinOperation;

/**
 * @brief Smooth Frontier Data Structure to be used in the Relaxation Methods
 * @author
 * @date
 * @ingroup Segmentation
 */
typedef struct ift_smooth_frontier {
  iftFImage *border_weight;
  iftFImage *norm_factor;
  iftFImage *prev_weight;
  iftFImage *next_weight;
  iftImage  *prev_label;
  iftImage  *next_label;
  iftImage  *prev_marker;
  iftImage  *next_marker;
  float      smooth_factor;
  int        smooth_iterations;
} iftSmoothBorder;

/**
 * @brief Converts a iftMaskJoinOperation name (str) into its corresponding iftMaskJoinOperation value
 * @author Cesar Castelo
 * @date Mar 19, 2019
 * @ingroup Segmentation
 */
iftMaskJoinOperation iftMaskJoinOperStrToMaskJoinOper(char *maskJoinOperStr);

/**
 * @brief Allocate a memory region to store the smooth frontier data structure
 * @author
 * @date
 * @ingroup
 *
 * @param basins Input Image.
 * @param A Adjacence Relation.
 * @param smooth_iterations Number of smooth interations.
 * @param smooth_factor Smooth's factor.
 * @return Pointer to the memory region storing the smooth frontier data structure.
 */
iftSmoothBorder *iftCreateSmoothBorder(iftImage *basins, iftAdjRel *A, int smooth_iterations, float smooth_factor);

/**
 * @brief Free the memory region used to store the smooth frontier data structure
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param smooth Reference to the smooth frontier data structure to be destroyed.
 */
void iftDestroySmoothBorder (iftSmoothBorder **smooth);

/**
 * @brief Filter used to enhance image's object
 * @author
 * @date
 * @ingroup Segmentation
 *
 * Filter used to enhance image's object by computing object map based on the distances between each voxel and
 * its closest roots in each set: object root set and background root set. This method is used in
 * pre-processment for image segmentation by watershed transforms.
 *
 * @param img Input Image.
 * @param seed Labeled Seeds.
 * @param obj Object's Label.
 * @return Image with enhanced object regarding its label.
 */
iftImage      *iftEnhanceObject(iftImage *img, iftLabeledSet *seed, int obj);

/**
 * @brief Filter used tp enhance image's edges
 * @author
 * @date
 * @ingroup Segmentation
 *
 * Filter all objects in the image and then combine them resulting in the
 * enhancement of all edges.
 *
 * @param img Input Image.
 * @param A Adjacence Relation.
 * @param seed Labeled Seeds.
 * @param alpha Contrast in the edge area is enhanced
 * @return Image with enhanced edges.
 */
//! swig(newobject, stable)
iftImage      *iftEnhanceEdges(iftImage *img, iftAdjRel *A, iftLabeledSet *seed, float alpha);

/**
 * @brief Filter used to enhance a MImage's object
 * @author
 * @date
 * @ingroup Segmentation
 *
 * Filter used to enhance a MImage's object by computing object map based on the distances between each voxel and
 * its closest roots in each set: object root set and background root set. This method is used in
 * pre-processment for image segmentation by watershed transforms.
 *
 * @param img Input MImage.
 * @param seed Labeled Seeds.
 * @param obj Object's Label.
 * @return Image with enhanced object regarding its label.
 */
iftImage      *iftMEnhanceObject(iftMImage *img, iftLabeledSet *seed, int obj);

/**
 * @brief Filter used tp enhance MImage's edges
 * @author
 * @date
 * @ingroup Segmentation
 *
 * Filter all objects in the MImage and then combine them resulting in the
 * enhancement of all edges.
 *
 * @param img Input MImage.
 * @param A Adjacence Relation.
 * @param seed Labeled Seeds.
 * @param alpha Contrast in the edge area is enhanced
 * @return Image with enhanced edges.
 */
iftImage      *iftMEnhanceEdges(iftMImage *img, iftAdjRel *A, iftLabeledSet *seed, float alpha);



/**
 * @brief Computes the Specialized Gradient for a Brain Image.
 * @author Samuel Martins
 * @date Jun 30, 2016
 * 
 * @param  brain_img Brain Image.
 * @return           Brain Gradient.
 */
iftImage *iftBrainGrad(const iftImage *brain_img);

/**
 * @brief Applies a SShape on the Input Image Image returning an enhanced image.
 * @author Samuel Martins
 * @date Jun 30, 2016
 */
iftImage *iftApplySShape(const iftImage *img, int a, int b, int c);

/**
 * @brief Computes the gradient vector field at voxels inside a given
 * mask. If the mask is not provided, it computes the gradient vector
 * field of the entire image.  
 * @author Alexandre Falcao
 * @date Nov 6, 2016
 *
 */
iftFImage **iftGradientVectorField(iftImage *img, iftImage *mask);

/**
 * @brief Gets the Gradient of Texture from an Image.
 * @author Samuel Martins
 * @date Jun 30, 2016
 *
 * @note This is an adaptation of the function TextGradient3 from old ift.
 */
iftImage *iftTextGradient(const iftImage *img);


/**
 * @brief      Compute the gradient of the input image according to a given algorithm:
 *
 * @param[in]  img       Input image.
 * @param[in]  grad_alg  Considered gradient algorithm.
 *
 * @return     The gradient of the input image.
 */
iftImage *iftComputeGradient(const iftImage *img, iftGradientAlg grad_alg);


/**
 * @brief Combine two gradient images linearly: (alpha * grad_img1) + ((1-alpha) * grad_img2)
 * 
 * @param  grad_img1     First Gradient Image
 * @param  grad_img2     Second Gradient Image
 * @param  alpha         Factor for combination.
 * @param  max_img_range Maximum image range used to normalize, by linear stretch, the input gradient
 *                       images (eg: 255, 4095, ...)
 * @return               Linear Combined Gradient Image.
 */
iftImage *iftCombineGradImages(iftImage *grad_img1, iftImage *grad_img2, float alpha, int max_img_range);

/**
 * @brief Computes the Watershed Transform from a set of labeled markers computing the predecessor map
 * @author Azael Sousa
 * @date
 * @ingroup Segmentation
 *
 * @param basins    Input Gradient Image.
 * @param Ain       The adjacency relation used to consider neighboring voxels.
 *                  If NULL, the function considers 8-neighborhood (if 2D image) or 26-neighborhood (if 3D image)
 * @param seeds     Set of labeled markers
 * @param forbidden Set of voxel indices from the image that cannot be achieved/segmented.
 *                  If NULL, all image's voxels can be segmented.
 * @param pred      The predecessor map
 * @return          Image with objects segmented by the watershed transform.
 */
iftImage *iftWatershedWithPredMap(const iftImage *basins, iftAdjRel *Ain, iftLabeledSet *seeds, iftSet *forbidden, iftImage **pred);

/**
 * @brief Computes the Watershed Transform from a set of labeled markers
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param basins    Input Gradient Image.
 * @param Ain       The adjacency relation used to consider neighboring voxels.
 *                  If NULL, the function considers 4-neighborhood (if 2D image) or 6-neighborhood (if 3D image)
 * @param seeds     Set of labeled markers
 * @param forbidden Set of voxel indices from the image that cannot be achieved/segmented.
 *                  If NULL, all image's voxels can be segmented.
 * @return          Image with objects segmented by the watershed transform.
 */
//! swig(newobject, stable)
iftImage *iftWatershed(const iftImage *basins, iftAdjRel *Ain, iftLabeledSet *seeds, iftSet *forbidden);


/**
 * @brief Computes GCmax, equivalent to the Watershed Transform, arc-weight = |I(p) - I(q)|
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param mimg      Multi-band Image
 * @param Ain       The adjacency relation used to consider neighboring voxels.
 *                  If NULL, the function considers 4-neighborhood (if 2D image) or 6-neighborhood (if 3D image)
 * @param seeds     Set of labeled markers
 * @param forbidden Set of voxel indices from the image that cannot be achieved/segmented.
 *                  If NULL, all image's voxels can be segmented.
 * @return          Image with objects segmented by the watershed transform.
 */
 //! swig(newobject, stable)
iftImage *iftWaterCut(iftMImage *mimg, iftAdjRel *Ain, iftLabeledSet *seeds, iftSet *forbidden);


/**
 * @brief Computes the Watershed Transform from a set of labeled markers, and also returns the
 * mean gradient value along the boundaries of the objects.
 * @author Samuel Martins
 * @date Jul 4, 2016
 * @ingroup Segmentation
 *
 * @param basins    Input Gradient Image.
 * @param Ain       The adjacency relation used to consider neighboring voxels.
 *                  If NULL, the function considers 8-neighborhood (if 2D image) or 26-neighborhood (if 3D image)
 * @param seeds     Set of labeled markers
 * @param forbidden Set of voxel indices from the image that cannot be achieved/segmented.
 *                  If NULL, all image's voxels can be segmented.
 * @param result    Return by Reference the mean gradient value along the boundaries of the objects.
 * @return          Image with objects segmented by the watershed transform.
 */
iftImage *iftWatershedMeanCut(const iftImage *basins, iftAdjRel *Ain, const iftLabeledSet *seeds,
							  const iftSet *forbidden, float *result);


/**
 * @brief Computes the watershed transform in a forest structure.
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param  fst  iftImageForest.    Forest structure created with a gradient image.
 * @param  seed iftLabeledSet.     List of spels with image index and seed label.
 * @param  removal_markers iftSet. List of spels marked for removal. NULL if empty

 * @return void. All forest maps are updated by reference.
 */
void iftDiffWatershed(iftImageForest *fst, iftLabeledSet *seed, iftSet *trees_for_removal);

/**
 * @brief Applies a relaxation in the label map.
 * @author
 * @date
 * @ingroup Segmentation
 *
 * This function relies on the bitmap 'processed' inside the iftImageForest structure. The bitmaps is automatically
 * created with the forest but it is only updated in iftDiffWatershed function. The object's border is computed only
 * in the region marked with 1 in the bitmap. To apply the relaxation to the entire label map, fill the processed
 * bitmap with 1 before invoking this function.
 *
 * Besides the object relaxation, this function also corrects any possible inconsistencies in the forest that could be
 * caused during the relaxation. More details about forest inconsistency during post-processing filters are described
 * in Nikolas Moya's masters dissertation.
 *
 * @param  fst  iftImageForest updated using iftDiffWatershed.
 * @param  smooth iftSmoothBorder. Structure with the relaxation parameters.

 * @return void. Consistent forest with relaxation roots.
 */
void iftRelaxObjects(iftImageForest *fst, iftSmoothBorder *smooth);

/**
 * @brief Computes a single shot watershed transform and applies relaxation.
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param  basins  iftImage    Gradient image created in the client.
 * @param  A iftAdjRel.     Adjacency Relation structure.
 * @param  seed iftLabeledSet. Set of labeled seeds from the client.
 * @param  num_smooth_iterations. Integer with the number of relaxation iterations.
 * @param  smooth_factor. Float number within [0, 1] with the relaxation factor.
 *
 * @return void. All forest maps are updated by reference.
 */
iftImage  *iftRelaxedWatershed(iftImage *basins, iftAdjRel *A, iftLabeledSet *seed, int num_smooth_iterations, float smooth_factor);

/**
 * @brief Enchance image when voxel's Cr is greater than Cb
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param img Colored Input Image
 *
 * @return Enhanced Image.
 */
iftImage  *iftEnhanceWhenCrIsGreaterThanCb(iftImage *img);

/**
 * @brief Computes the Watershed Transform Based on the path weight of a given dataset
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param dataset Set of data samples
 * @param A Adjacence Relation
 * @param seed Set of Seeds
 *
 * @return Image with the watershed transform applied.
 */
iftImage  *iftWatershedOnVoxelDist(iftDataSet *dataset, iftAdjRel *A,iftLabeledSet *seed);

/**
 * @brief Compute the normalize factor to a voxel and its neighbors
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param weight Weight map of each voxel
 * @param A Adjacence Relation (voxel's neighbors)
 *
 * @return Image with voxel's normalize factor.
 */
iftFImage *iftWeightNormFactor(const iftFImage *weight, iftAdjRel *A);

/**
 * @brief Computes the Watergray Transform from a marker image
 * @author
 * @date
 * @ingroup Segmentation
 *
 * Watershed transforms from grayscale markers: Note that, in order to obtain a connected operator for some graph
 * topology, the adjacency relation must be the same used to compute the input image.
 *
 * @param img Input Image.
 * @param marker Image Marker
 * @param seed Set of labeled markers
 *
 * @return Image with objects segmented by the watershed transform.
 */

//! swig(newobject)  
iftImage  *iftWaterGray(iftImage *basins, iftImage *marker, iftAdjRel *A);

/**
 * @brief Computes the Watergray Transform from a marker image in a forest structure
 * @author
 * @date
 * @ingroup Segmentation
 *
 * Watershed transforms from grayscale markers: Note that, in order to obtain a connected operator for some graph
 * topology, the adjacency relation must be the same used to compute the input image.
 *
 * @param fst Image Forest
 * @param marker Image Marker
 *
 * @return Void.
 */
 //! swig()
void iftWaterGrayForest(iftImageForest *fst, iftImage *marker);

/**
 * @brief Watershed Transform in a Distance Transformed Image
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param dist Image with distance transform
 * @param label Labeled Image
 * @param H Distance Weight
 * @param A Adjacency Relation
 *
 * @return Image with Watershed Transform.
 */
iftImage  *iftWaterDist(iftImage *dist, iftImage *label, int H, iftAdjRel *A);

/**
 * @brief Computes the dual WaterGray segmentation algorithm
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param domes The image domes.
 * @param marker The marker image.
 * @param A The adjacency relation of the graph.
 *
 * @return Void.
 */
iftImage  *iftDualWaterGray(iftImage *domes, iftImage *marker, iftAdjRel *A);

/**
 * @brief Computes the dual WaterGray segmentation algorithm inside a given mask.
 *
 * @param domes The image domes.
 * @param marker The marker image.
 * @param mask The mask
 * @param A The adjacency relation of the graph.
 * @param roots If not NULL, the set of roots is returned.
 * @return The segmentation label.
 *
 * @note The marker image is modified inside this function.
 */
iftImage *iftDualWaterGrayOnMask(iftImage *domes, iftImage *marker, iftImage *mask, iftAdjRel  *A, iftSet **roots);

/**
 * @brief Computes the Watershed Transform based on a given orientation
 * @author Alexandre Falcao
 * @date   Jul 5th 2016
 * @ingroup Segmentation
 *
 * @param img: Original image, only brightness band
 * @param grad: Gradient image used to compute the graph arc-weight
 * @param A: Adjacency relation of the image graph (optional, type NULL for default options: spheric for 3D images with radius 1.0 or circular for 2D images with radius 1.5)
 * @param seed: Set of seeds (internal and external markers)
 * @param alpha: array of values between (-1, 1), when alpha > 0 it penalizes the path cost from bright to dark nodes, otherwise dark to bright.
 * @param forbidden: Image region excluded from the watershed computation (optional, type NULL if there is none)
 *
 * @return Labeled image that results from the watershed computation.
 */
//! swig(newobject, stable)
iftImage *iftOrientedWatershed(const iftImage *img, const iftImage *grad, iftAdjRel *Ain,
                               iftFloatArray *alpha, iftLabeledSet *seeds, iftSet *forbidden);

/**
 * @brief Threshold in the Image Histogram
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param img Input Image
 * @param perc Percentage Threshold
 *
 * @return Position in the histogram whose value surpass the Threshold.
 */
int iftCumHistogramThres(iftImage *img, float perc);

/**
 * @brief Define a Threshold value through Otsu Method
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param img Input Image
 *
 * @return Threshold Value.
 */
//! swig(stable)
int iftOtsu(const iftImage *img);

/**
 * @brief Applies the Otsu Method in Regions defined by a Mask Image.
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param img Input Image
 * @param mask Image Mask
 *
 * @return Threshold value.
 */
int iftOtsuInRegion(iftImage *img, iftImage *mask);

/**
 * @brief Computes the Threshold in an image
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param img Input Image
 * @param lowest Threshold Base
 * @param highest Threshold Ceil
 * @param value Number to be assign to the voxels
 *
 * In case the image's voxel has intensity between the base and ceil Threshold, then, the image's voxel receives
 * the parameter value, otherwise, it receives 0.
 *
 * @param img Input Image.
 * @param lower Base value for the Threshold.
 * @param highest Ceil value fot the Threshold.
 * @param value Value to be assign to the image's voxel.
 * @return Image with the Threshold applied.
 */
//! swig(newobject)
iftImage *iftThreshold(const iftImage *img, int lowest, int highest, int value);

/**
 * @brief Computes the Threshold in a Fimage
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param img Input FImage
 * @param lowest Threshold Base
 * @param highest Threshold Ceil
 * @param value Number to be assign to the voxels
 *
 * In case the Fimage's voxel has intensity between the base and ceil Threshold, then, the Fimage's voxel receives
 * the parameter value, otherwise, it receives 0.
 *
 * @param img Input FImage.
 * @param lower Base value for the Threshold.
 * @param highest Ceil value fot the Threshold.
 * @param value Value to be assign to the Fimage's voxel.
 * @return Image with the Threshold applied.
 */
iftImage  *iftFThreshold(const iftFImage *img, float lowest, float highest, int value);

/**
 * @brief Use of OPF method to binarize an image
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param orig Original Image
 * @param enha Enhanced Image
 * @param init_thres Initial Threshold
 * @param train_perc Considered percentage of the Enhanced image
 *
 * @return Binarized Image.
 */
iftImage  *iftBinarizeByOPF(iftImage *orig, iftImage *enha, int init_thres, float train_perc);

/**
 * @brief Computes the mean value within the adjacency region A and selects voxels that are above perc*mean, by
 * assigning value (positive number) to them. Repeats this process for a certain number of iterations (niters),
 * excluding voxels already selected.
 * @author
 * @date
 * @ingroup Segmentation
 *
 * Repeat the following process during a certain number of iterations (niters), excluding previously selected voxels:
 * For each voxel p, compute the mean value among its adjacent voxels. If the intensity of p is above a percentile of
 * the mean value, then p is selected (i.e., assign value to the output intensity of p).
 *
 * @param img Input Image
 * @param A Adjacency Relation
 * @param mask An optional mask or NULL to constrain the thresholding
 * @param perc Percentage
 * @param niters Number of iterations
 * @param value Value to be assign
 *
 * @return Threshold Image.
 */
//! swig(newobject)  
  iftImage  *iftAboveAdaptiveThreshold(iftImage *img, iftImage *mask, iftAdjRel *A, float perc, int niters, int value);

/**
 * @brief Computes the mean value within the adjacency region A and selects voxels that are below perc*mean, by
 * assigning value (positive number) to them. Repeats this process for a certain number of iterations (niters),
 * excluding voxels already selected.
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param img Input Image
 * @param A Adjacency Relation
 * @param mask An optional mask or NULL to constrain the thresholding
 * @param perc Percentage
 * @param niters Number of iterations
 * @param value Value to be assign
 *
 * @return Threshold Image.
 */
//! swig(newobject)    
  iftImage  *iftBelowAdaptiveThreshold(iftImage *img, iftImage *mask, iftAdjRel *A, float perc, int niters, int value);

/**
 * @brief Computes a smooth weight image
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param basins Input Image Basins
 * @param beta Contrast Factor
 *
 * @return Weighted FImage.
 */
iftFImage *iftSmoothWeightImage(const iftImage *basins, float beta);

/**
 * @brief Perform a fast smooth on the objects
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param labelIn Input Label Image
 * @param weight Weight Map
 * @param ninter Number of Iterations
 *
 * @return Void.
 */
iftImage *iftFastSmoothObjects(const iftImage *labelIn, const iftFImage *weight, int n_iters);

/**
 * @brief It transforms an image of borders into an image of regions. The regions are the interior of the closed
 * contours and the lines disappear.
 * @author Alexandre Falcao.
 * @date 4/12/2015
 * @ingroup Segmentation
 *
 * @param  border Binary image with value 0 for background and value >= 1 for border voxels.
 *
 * @return An image of regions labeled from 1 to N.
 */
iftImage* iftBorderImageToLabelImage(iftImage* border);

/**
 * @brief Preserves the input label of border voxels inside and outside the object and sets the label of the remaining
 * voxels to zero.
 * @author Alexandre Falcao
 * @date Nov 30, 2015
 * @ingroup Segmentation
 *
 * @param label: it can be an image with labeled objects or labeled supervoxels. If the label image contains a single
 * object, the function will preserve only the internal border voxels.
 * @param get_margins: set if the margins of the original image are needed. 1 for yes and 0 for not. Parameter added by Adan Echemendia on May, 2017.
 *
 * @return border: A labeled image of border voxels. Each border voxel preserves the orginal label of its supervoxel/object
 * in the input image.
 */
//! swig(newobject)
iftImage *iftBorderImage(const iftImage *label, bool get_margins);


/* Similarity measures for binary segmentation */
/**
 * @brief Compute the Beta-Score for a given binary image, a groundtrut and a beta factor
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param bin Imput Binary Image
 * @param  gt Correct segmentation. 0 for background and 1 for object.
 * @param beta Score Factor
 *
 * @return Void.
 */
float    iftFBetaScore(iftImage *bin, iftImage *gt, float beta);

/**
 * @brief F1-Score for binary segmentation.\n Given a groundtruth and binary label image, compute the F1-Score.\n
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param  bin iftImage.    Binary label image. 0 for background and 1 for object.
 * @param  gt  iftImage.    Correct segmentation. 0 for background and 1 for object.
 *
 * @return A float with the FScore [0, 1].
 */
float iftFScoreError(iftImage *bin, iftImage *gt);

/**
 * @brief F1-Score extension for multilabels.\n Given a multilabel image and a multilabel ground truth, it computes the
 * Fscore measure indivudally. The average is stored at index 0. The other i indexes store the fscore error for the ith
 * object.The perfect segmentation produces error 1.\n
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param  mlabel                iftImage.    Multi label image resultant from a segmentation.
 * @param  mlgt                  iftImage.    Multi label groundtruth read from disk.
 * @param  number_of_objects     Integer.     Number of objets. For binary segmentation, this value must be 1.
 *
 * @return A float array of size (number_of_objects + 1) with the fscore of each object. Index 0 is the average
 * between all objects.
 */
float *iftFScoreMultiLabel (iftImage *mlabel, iftImage *mlgt, int number_of_objects);

/* Error statistics */
/**
 * @brief Computes the Statistics Error on the Segmentation
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param gt_image Groundtruth image
 * @param cl_image Segmentation image
 *
 * @note This Function only works on Groundtruth images and Segmentation Images
 *
 * @return Segmentation Statistic Error.
 */
iftErrorClassification iftSegmentationErrors(iftImage* gt_image, iftImage* cl_image);
 
/**
 * @brief Computes the recall between two border images, the ground truth and the resulting image from segmentation
 * @author Alexandre Falcao
 * @date Nov 30, 2015
 * @ingroup Segmentation
 *
 * @param  gt: the ground truth image of borders.
 * @param  border: the border image resulting from segmentation.
 * @param  tolerance_dist: a distance of tolerance from the ground truth, since the borders of ground-truth masks.
 * created by manual tracing do not follow the border of the object in the image.
 *
 * @return a number within [0,1] that represents the boundary recall.
 */
float iftBoundaryRecall(iftImage *gt, iftImage *border, float tolerance_dist);

/**
 * @brief Computes the under segmentation error --- total of leaking voxels across the boundaries of a ground truth image.
 * Implementation according to the SLIC paper
 * @author Alexandre Falcao and John Vargas
 * @date Dec, 7th 2015
 * @ingroup Segmentation
 *
 * @param  gt: the ground truth image of regions@param  label: the label image resulting from segmentation.
 * @param  perc_of_intersection: Percentage of intersection between region and ground truth segment to consider the
 * region as part of the object (true intersection).
 *
 * @return A number that represents the under segmentation error.
 */
float iftUnderSegmentationSLIC(iftImage *gt, iftImage *label, float perc_of_intersection);


/**
 * @brief Computes the under segmentation error --- total of leaking voxels across the boundaries of a ground truth image.
 * Implementation according to the paper "Superpixel Benchmark and Comparison"
 * @author Alexandre Falcao and John Vargas
 * @date Dec, 7th 2015
 * @ingroup Segmentation
 *
 * @param  gt: the ground truth image of regions
 * @param  label: the label image resulting from segmentation.
 * @param  perc_of_intersection: Percentage of intersection between region and ground truth segment to consider the
 * region as part of the object (true intersection).
 *
 * @return A number that represents the under segmentation error.
 */
float iftUnderSegmentation(iftImage *gt_image, iftImage *label_image);

/**
 * @brief Compute the Topology Measure
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param label Input Image Label
 *
 * @return Void.
 */
float iftTopologyMeasure(iftImage *label);

/**
 * @brief Compress 2D Image
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param label Input Label Image
 *
 * @return float.
 */
float iftCompactness2D(iftImage *label);

/**
 * @brief Select and Propagate Regions Above a given area
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param label Input Label Image
 * @param area Area Threshold
 *
 * @return Image with all regions above the area propagated.
 */
  
iftImage *iftSelectAndPropagateRegionsAboveArea(iftImage *label, int area);

/**
 * @brief Join superpixels whose size is less than <br>area</br> with their most color similar adjacent superpixel
 * @author Adan Echemendia Montero
 * @date 19/06/2017
 * @ingroup Segmentation
 *
 * @param img An image
 * @param label Input Label Image
 * @param area Area Threshold
 *
 * @return Image with all regions above the area propagated.
 */
//! swig(newobject)  
iftImage *iftSelectAndPropagateRegionsAboveAreaByColor(iftImage *img, iftImage *label, int area);

/**
 * @brief Smoothes Regions by Diffusion
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param label Label Image
 * @param orig Original Image
 * @param smooth_factor Smooth Factor
 * @pparam niters Number of Iterations
 *
 * @return Image with Smoothed Regions.
 */
//! swig(newobject)
iftImage *iftSmoothRegionsByDiffusion(const iftImage *label_img, const iftImage *orig_img,
                                      float smooth_factor, int n_iters);

/**
 * @brief Check Segmentation Consistency
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param fst Image Forest Structure
 *
 * @return 1 for consistent segmentation and 0 for inconsistent segmentation.
 */
char iftIsSegmentationConsistent(iftImageForest *fst);

/**
 * @brief Checks if the label root is connected
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param pred Input Image
 * @param label Label Image
 * @param p Voxel position
 *
 * @return 1, if it is connect, 0 otherwise.
 */
int iftIsLabelRootConnected (iftImage *pred, iftImage *label, int p);

/**
 * @brief
 * @author
 * @date
 * @ingroup Segmentation
 *
 * @param
 * @param
 *
 * @return Void.
 */
iftLabeledSet *iftLabelToForestGeodesicRobot(iftImage *gradient,iftImage *label, int seeds_per_iteration, int number_of_objects, int min_distance_border, int max_marker_size, int min_marker_size);


/**
 * @brief Shifts the voxel coordinates from a Labeled Set. Ignores shifted voxels out of the Image Domain.
 * @author Samuel Martins
 * @date Jul 4, 2016
 * 
 * @param  S   Labeled Set with the coordinates of the voxels to be shifted.
 * @param  img Image from which the coordinates belongs to.
 * @param  dx  Displacement in x-axis.
 * @param  dy  Displacement in y-axis.
 * @param  dz  Displacement in z-axis.
 * @return     The shifted Labeled Set.
 */
iftLabeledSet *iftShiftCoordsFromLabeledSet(const iftLabeledSet *S, const iftImage *img, int dx, int dy, int dz);


/**
 * @brief Shifts the voxel coordinates from a Set. Ignores shifted voxels out of the Image Domain.
 * @author Samuel Martins
 * @date Jul 4, 2016
 * 
 * @param  S   Set with the coordinates of the voxels to be shifted.
 * @param  img Image from which the coordinates belongs to.
 * @param  dx  Displacement in x-axis.
 * @param  dy  Displacement in y-axis.
 * @param  dz  Displacement in z-axis.
 * @return     The shifted Set.
 */
iftSet *iftShiftCoordsFromSet(const iftSet *S, const iftImage *img, int dx, int dy, int dz);

/**
 * @brief Extends a binary marker by selecting the connected voxels whose intensities fall in a given interval.
 * @author Alexandre Falcao
 * @date Nov 10th, 2016
 * 
 * @param  img Input image
 * @param  A   Input adjacency relation
 * @param  marker Input binary marker
 * @param  minthres Minimum intensity
 * @param  maxthres Maximum intensity
 * @return extended binary marker
 */
iftImage *iftConnectedThreshold(const iftImage *img, const iftAdjRel *A, iftImage *marker, const int minthres, const int maxthres);


/**
 * @brief Gets the most frequent label in as matrix of occurrences for each voxel from an image by Majority Voting.
 *
 * label_occurs is an integer matrix n_voxels x n_labels. \n
 * Each voxel p has the occurrences for the labels from 0 to (n_labels-1).
 * 
 * @param  img          Image analyzed.
 * @param  label_occurs Matrix (n_voxels x n_labels) of occurrences of each label for each voxel of input image.
 * @return              Image with the most frequent label for each voxel.
 */
iftImage *iftMajorityVoting(const iftImage *img, const iftIntMatrix *label_occurs);

/**
 * @brief Computes the descendant map of an optimum-path forest for Tree Pruning
 * @author Thiago Vallin Spina
 *
 * @param  pred             Predecessor map/optimum-path forest.
 * @param  A                Adjacency relation used during segmentation.
 * @param set_of_interest   (Optional, may be NULL) Binary map indicating the set of voxels from which the descendant map should be computed.
 * @return                  Image with the number of descendants for each node/voxel in the forest.
 */
iftImage* iftDescendantMap(iftImage *pred, iftAdjRel *A, iftBMap *set_of_interest);


/**
 * @brief Label the connected components (objects) of a binary or labeled image.
 * @param  label_img [description]
 * @param  n_objs    [description]
 * @return           [description]
 */
iftImage *iftConnCompLabeling(const iftImage *label_img, int *n_objs);

/**
 * @brief Computes the Oriented Image Foresting Transform based on Oriented Image Foresting Transform Segmentation
by Seed Competition paper by Paulo A. V. Miranda and Lucy A. C. Mansilla.
 * @author Jordão Bragantini
 * @date   Jun 21rst 2018
 * @ingroup Segmentation
 *
 * @param img: Original image, only brightness band
 * @param A: Adjacency relation of the image graph (optional, type NULL for default options: spheric for 3D images with radius 1.0 or circular for 2D images with radius 1.5)
 * @param seed: Set of seeds (internal and external markers)
 * @param alpha: array of values between (-1, 1), when alpha > 0 it penalizes the path cost from bright to dark nodes, otherwise dark to bright.
 * @param forbidden: Image region excluded from the watershed computation (optional, type NULL if there is none)
 *
 * @return Labeled image that results from the watershed computation.
 */
//! swig(newobject, stable)
iftImage *iftOrientedWaterCut(const iftImage *img, iftAdjRel *Ain, iftFloatArray *alpha,
                                          iftLabeledSet *seeds, iftSet *forbidden);

//! swig(newobject)
iftImage *iftOrientedColorWaterCut(iftMImage *mimg, iftImage *orient, iftAdjRel *Ain, float beta, iftLabeledSet *seeds);

//! swig(newobject, stable)
iftImage *iftEnhancedWaterCut(iftMImage *mimg, iftImage *objmap, iftAdjRel *Ain, iftLabeledSet *seeds, float alpha);

//! swig(newobject, stable)
iftImage *iftEnhancedWatershed(iftImage *basins, iftImage *objmap, iftAdjRel *Ain, iftLabeledSet *seeds, float alpha);

//! swig(newobject)
iftImage *iftSuperPixelMajorityVote(iftImage *comp, iftImage *objmap, float threshold);

/**
 * @brief Creates a mask by joining a set of masks in a file dir. The join operation could be union or intersection
 * @author César Castelo
 * @date   Set 17, 2018
 * @ingroup Segmentation
 *
 * @param maskDir: Directory with masks
 * @param operation: Operation to the applied in the masks: Union/Intersection
 * @param maskVal: Value to represent the mask pixels (e.g. 1 or 255)
 *
 * @return Joint mask
 */
iftImage *iftCreateJointMaskFromFileDir(char *maskDir, iftMaskJoinOperation operation, int maskVal);
iftImage *iftCreateJointMaskFromFileset(iftFileSet *maskFileset, iftMaskJoinOperation operation, int maskVal);

/**
 * @brief Creates a mask with the superpixel boundaries
 * @author César Castelo
 * @date   Mar 14, 2019
 * @ingroup Segmentation
 *
 * @param spixLabels: Image containing the superpixel labels
 * @param maskVal: Value to be set in the mask
 *
 * @return Mask containing the superpixel boundaries
 */
iftImage *iftCreateSuperpixelBoundariesMask(iftImage *spixLabels, int maskVal);


//! swig(newobject)
iftImage *iftBoundingBoxArrayToLabel(const iftImage *img, const iftBoundingBoxArray *bb_ary);

//! swig(newobject)
iftSet *iftFindOptPath(const iftMImage *mimg, int src, int dst, float sigma, iftFImage *pathval, iftImage *pred);

//! swig(newobject)
iftSet *iftFindOptPathWithProb(const iftMImage *mimg, const iftFImage *obj, const iftFImage *bkg, int src, int dst,
                               float sigma, float gamma, iftFImage *pathval, iftImage *pred);

//! swig(newobject)
iftImage *iftSinglePathContour(const iftImage *label);

//! swig(newobject)
iftImage *iftFillContour(const iftImage *contour);

/**
* @brief Obtains a label image by executing the DISF algorithm
* @author Felipe Belem
* @date July 21, 2020
* @ingroup Segmentation
*
* @param mimg - Multiband Image (CIELAB, preferably)
* @param A - Adjacency Relation (8- or 26-adjacency, for 2D and 3D)
* @param num_init_seeds - Initial number of seeds N_0 (N_0 > 0)
* @param num_superpixels - Final number of superpixels N_f (N_0 >> N_f)
* @param mask - Mask indicating spels to be conquered (NULL for any spel within mimg)
* @return Label image with labels l in [1,N_f]
*/

//! swig(newobject)
iftImage *iftDISF(iftMImage *mimg, iftAdjRel *A, int num_init_seeds, int num_superpixels, iftImage *mask);


#ifdef __cplusplus
}
#endif

#endif // IFT_SEGMENTATION_H_

