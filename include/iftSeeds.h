#ifndef IFT_SEEDS_H_
#define IFT_SEEDS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BMap.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "ift/core/dtypes/List.h"
#include "ift/core/dtypes/Set.h"
#include "ift/core/dtypes/DHeap.h"

#include "iftAdjacency.h"
#include "iftCommon.h"
#include "iftDataSet.h"
#include "iftImage.h"
#include "iftMatrix.h"
#include "iftMathMorph.h"
#include "iftMImage.h"
#include "iftRegion.h"

//! swig(newobject)
iftSet 	      *iftExtractRemovalMarkers(iftLabeledSet **s);
//! swig(newobject)
iftLabeledSet *iftLabelObjBorderSet(iftImage *bin, iftAdjRel *A);
//! swig(newobject)
iftLabeledSet *iftImageBorderLabeledSet(iftImage *img);
//! swig(newobject)
iftLabeledSet *iftLabelCompSet(iftImage *bin, iftAdjRel *A);
//! swig(newobject)
iftLabeledSet *iftFuzzyModelToLabeledSet(iftImage *model);
//! swig(newobject)
iftLabeledSet *iftMAdjustSeedCoordinates(iftLabeledSet *S, iftMImage *input, iftMImage *output);
//! swig(newobject)
iftLabeledSet *iftAdjustSeedCoordinates(iftLabeledSet *Sin, iftImage *orig, iftMImage *output);

//! swig(newobject)
iftSet        *iftImageBorderSet(const iftImage *img);
//! swig(newobject)
iftLabeledSet *iftMultiObjectBorderLabeledSet(iftImage *img, iftAdjRel *A);


/**
 * @brief Gets the voxel's index from all object voxels (val > 0), returning them into a set
 * @param  label_img Label Image.
 * @return           Set with the voxel's indices from all objects.
 */
//! swig(newobject)
iftSet *iftObjectToSet(const iftImage *label_img);

/**
 * @brief Gets the voxel's index from all object voxels (val > 0), returning them into a labeledset
 * @param  label_img Label Image.
 * @return LabeledSet with the voxel's indices from all objects.
 */
//! swig(newobject)
iftLabeledSet *iftObjectToLabeledSet(const iftImage *label_img);

/**
 * @brief Gets the voxel's index from all object voxels (val > 0), returning them into a list
 * @param  label_img Label Image.
 * @return           Set with the voxel's indices from all objects.
 *
 * @author Samuka Martins
 * @date Nov 16, 2017
 */
//! swig(newobject)
iftList *iftObjectsToList(const iftImage *label_img);


/**
 * @brief Finds and returns an iftSet with the coordinates of all border spels of a label image for a given adjacency <A>.
 * @param label Image label.
 * @param A Adjacency Relation.
 * @return The iftSet with the coordinates of all border spels of the label image <label>.
 */
//! swig(newobject)
iftSet *iftObjectBorderSet(const iftImage *label_img, iftAdjRel *Ain);
//! swig(newobject)
iftSet *iftBackgroundBorderSet(const iftImage *label_img, iftAdjRel *Ain);


/**
 * @brief Finds and returns an iftBMap with the coordinates of all border spels of a label image for a given adjacency <A>.
 * @param label Image label.
 * @param A Adjacency Relation.
 * @param n_border_spels Number of border pixels/voxels of the label image that will be returned by reference.
 * @return The iftBMap with the coordinates of all border spels of the label image <label>.
 */
//! swig(newobject)
iftBMap *iftObjectBorderBitMap(const iftImage *label, const iftAdjRel *A, int *n_border_spels);


/**
 * @brief Stores all pixels from a Mask (Bin Image, except the background 0) to an input Set.
 * @author Samuel Martins
 * @dataset Jul 12, 2016
 * 
 * @param  mask Input Mask (Binary Image).
 * @return     The resulting Set.
 *
 * @note The input Set S may already have elements.
 */
//! swig()
void iftMaskImageToSet(const iftImage *mask, iftSet **S);

/**
 * @brief Create a binary image with values 1 at voxels from a given set. 
 * @author Alexandre Falcao 
 * @dataset Sep 10, 2021
 * 
 * @param  Seed set.
 * @return Resulting binary mask.
 *
 */
//! swig(newobject)  
  iftImage *iftMaskImageFromSet(iftSet *S, int xsize, int ysize, int zsize);


  
/**
 * @brief Sets the value <obj> in all voxels in the input image (which must have been allocated previously) from the input set.
 * @author Samuka Martins
 * @date Nov 16, 2017
 */
//! swig()  
void iftSetToImage(const iftSet *S, iftImage *img, int obj);


/**
 * @brief Sets the value <obj> in all voxels in the input image (which must have been allocated previously) from the input List.
 * @author Samuka Martins
 * @date Apr 19, 2018
 */
//! swig()
void iftListToImage(const iftList *L, iftImage *img, int obj);

/**
 * @brief Create a list with the indices of all objects voxels from a mask (binary image).
 * @author Samuel Martins
 * @date Apr 24, 2018
 */
//! swig(newobject)
iftList *iftMaskToList(const iftImage *mask);


/**
 * @brief Insert an Labeled Set into image <img>. In order to save background seeds (which has 0-value)
 * to 1, use increment_label = true.
 * 
 * @author Samuka Martins
 * @date Dec 27, 2017
 */
//! swig()
void iftLabeledSetToImage(const iftLabeledSet *S, iftImage *img, bool increment_label);


/**
 * @brief Sets the value <obj> in all voxels in the input image (which must have been allocated previously) from the input integer array.
 * @author Samuel Martins
 * @date Apr 24, 2018
 */
//! swig()
void iftIntArrayToImage(const iftIntArray *iarr, iftImage *img, int obj);


/**
 * @brief Create an integer array with the indices of all objects voxels from a mask (binary image).
 * @author Samuel Martins
 * @date Apr 24, 2018
 */
//! swig(newobject)  
iftIntArray *iftMaskToIntArray(const iftImage *mask);


/**
 * @brief   Copy all labeled region voxels to a Labeled Set. 
 * @author  Alexandre Falcao (Updated by Jorda)
 * @dataset Aug 1st, 2018 (Oct. 30th, 2018)
 * 
 * @param  label: Labeled region image, usually from 1 to n regions.
 * @param  decrement_label : Boolean value to indicate when the labels must be decremented by 1. It is used when label 1 indicates background.
 * @return The resulting labeled seed set.
 *
 */

//! swig(newobject)
iftLabeledSet *iftRegionsToLabeledSet(const iftImage *label, bool decrement_label);
  


//! swig(newobject)
iftSet        *iftEndPoints(iftImage *skel, iftAdjRel *A);
//! swig(newobject)
iftSet        *iftFindPathOnSkeleton(iftImage *skel, iftAdjRel *A, int src, int dst);
//! swig(newobject)
iftSet        *iftSkeletonPoints(iftImage *skel);


/**
 * @brief Return the Label Image's borders for a given Adjacency.
 * 
 * A border voxel is one with some adjacent voxel with different value or a frame voxel.
 * If keep_border_labels is true, it keeps the original label from the input label image.
 * Otherwise, all borders are binarized with value 1.
 * If include_image_frame is true, object pixels on image frame are included.
 * 
 * 
 * @param  label_img          (Multi) Label Image.
 * @param  Ain                Input adjacent used to find the borders. If NULL, a 4-neighborhood (2D) or 6-neighborhood (3D)
 *                            considered.
 * @param keep_border_labels  Keeps the original label from the input label image or binarize them with value 1.
 * @param include_image_frame Includes the object pixels on image frame.
 * @return                    An image with the found borders.
 *
 * @author Samuka Martins, Cesar Castelo
 * @date Nov 30, 2017
 */
//! swig(newobject)
iftImage *iftObjectBorders(const iftImage *label_img, const iftAdjRel *Ain, bool keep_border_labels,
                           bool include_image_frame);


/**
 * @brief Find the object borders of a labeled image and ALWAYS relabel them, even if the image
 * has multiple labels, since more than one object could have the same label.
 * It respects the label from the input Label Image to propagate the new label for the objects,
 * so that two adjacent objects with different labels in the input image will be also
 * different objects in the relabeled image.
 * 
 * @param  label_img Labeled Image (binary or multi-label)
 * @param  Ain       Input adjacent used to label the connected components within the borders,
 *                   If NULL, a 8-neighborhood (2D) or 26-neighborhood (3D) is considered.
 * @return           Image with the labeled object borders.
 * 
 * @author Samuka Martins, Cesar Castelo
 * @date Nov 30, 2017
 */
//! swig(newobject)  
iftImage *iftFindAndLabelObjectBorders(const iftImage *label_img, const iftAdjRel *Ain);
//! swig(newobject)  
iftImage      *iftEasyLabelComp(iftImage *bin, iftAdjRel *A);
//! swig(newobject)  
iftImage      *iftLabelComp(iftImage *bin, iftAdjRel *A);


/**
 * @brief Select the Largest Component of a Binary Image.
 * 
 * @param bin Binary Image (0 = background, any value > 0 to represent the object).
 * @param Ain   Adjacency Relation. If NULL, it uses 4-neighborhood (2D) or 6-neighborhood (3D)
 * @return The largest component of the binary image with the object's label.
 */
//! swig(newobject)
iftImage *iftSelectLargestComp(const iftImage *bin, const iftAdjRel *Ain);


//! swig(newobject)
iftImage      *iftSelectSmallestComp(iftImage *bin, iftAdjRel *A);
//! swig(newobject)
iftImage      *iftSelectKLargestComp(iftImage *bin, iftAdjRel *A, int K);
//! swig(newobject)
iftImage      *iftSelectKSmallestComp(iftImage *bin, iftAdjRel *A, int K);
//! swig(newobject)
iftImage      *iftComponentArea(iftImage *bin, iftAdjRel *A);
//! swig(newobject)
iftImage      *iftSelectCompAboveArea(iftImage *bin, iftAdjRel *A, int thres);
//! swig(newobject)
iftImage      *iftSelectCompBelowArea(iftImage *bin, iftAdjRel *A, int thres);
//! swig(newobject)
iftImage      *iftSelectCompInAreaInterval(iftImage *bin, iftAdjRel *A, int thres_min, int thres_max);

//! swig(newobject)
iftImage      *iftRegionalMaxima(iftImage *img);
//! swig(newobject)
iftImage      *iftRegionalMinima(iftImage *img);
//! swig(newobject)
iftImage      *iftRegionalMaximaInRegion(iftImage *img, iftImage *mask);
//! swig(newobject)
iftImage      *iftRegionalMinimaInRegion(iftImage *img, iftImage *mask);
//! swig(newobject)
iftImage      *iftRootVoxels(iftImage *pred);
//! swig(newobject)
iftImage      *iftLeafVoxels(iftImage *pred, iftAdjRel *A);
//! swig(newobject)
iftImage      *iftLeafVoxelsOnMask(iftImage *pred, iftAdjRel *A, iftImage *mask);

  /* Assume that label contains regions with values in
     {1,2,..,n}. That is, label equal to zero is considered background
     and it does not define any region */

//! swig(newobject)
iftImage      *iftRegionArea(iftImage *label);


/**
 * @brief Select the Largest Component of a Labeled Image whose objects has label from 1 to n.
 * 
 * @param bin Binary Image
 * @param A   Adjacency Relation. If NULL, it uses 4-neighborhood (2D) or 6-neighborhood (3D)
 * @return The largest component of the binary image. Its label is the original region's label. 
 */
//! swig(newobject)
iftImage *iftSelectLargestRegion(const iftImage *label_img);

//! swig(newobject)
iftImage      *iftSelectSmallestRegion(iftImage *label);
//! swig(newobject)
iftImage      *iftSelectKLargestRegions(iftImage *label, int K);
//! swig(newobject)
iftImage      *iftSelectKLargestRegionsAndPropagateTheirLabels(iftImage *label, iftAdjRel *A, int K);
//! swig(newobject)
iftImage      *iftSelectKSmallestRegions(iftImage *label, int K);
//! swig(newobject)
iftImage      *iftSelectRegionsAboveArea(iftImage *label, int thres);
//! swig(newobject)
iftImage      *iftSelectRegionsAboveAreaAndPropagateTheirLabels(iftImage *label, int thres);
//! swig(newobject)
iftImage      *iftSelectRegionsBelowArea(iftImage *label, int thres);
/**
 * @brief Selects the regions with area size within the given interval.
 *
 * The regions are relabeled automatically in the output image.
 *
 * @param in_label The input label.
 * @param min_thres The minimum area size.
 * @param max_thres The maximum area size.
 * @return The relabeled image with removed regions.
 *
 */

//! swig(newobject)
iftImage      *iftSelectRegionsInAreaInterval(iftImage *label, int min_area, int max_area);
/**
 * @brief Selects the regions with area size within the given interval.
 *
 * The regions are relabeled automatically in-place, as opposed to iftSelectRegionsInAreaInterval.
 *
 * @param in_label The input label.
 * @param min_thres The minimum area size.
 * @param max_thres The maximum area size.
 * @sa iftSelectRegionsInAreaInterval
 *
 */
//! swig(newobject)
void           iftSelectRegionsInAreaIntervalInplace(iftImage *label, int min_area, int max_area);

  /* ---------------------------- */

//! swig(newobject)
iftImage  *iftLabelContPixel(iftImage *bin);
char       iftValidArc(iftBMap *bndr, iftImage *pred, iftImage *bin, iftAdjRel *A, iftAdjRel *L, iftAdjRel *R, iftVoxel u, int i, int q);
  char       iftValidStartingPixel(iftBMap *bndr, iftImage *pred, iftImage *bin, iftAdjRel *A, iftAdjRel *L, iftAdjRel *R, int p, int *nvalidarcs);

//! swig(newobject, stable)
iftLabeledSet *iftReadSeeds(const iftImage *img, const char *filename, ...);

//! swig(newobject, stable)
iftLabeledSet *iftMReadSeeds(iftMImage *img, char *filename);



iftLabeledSet *iftReadSeedsGraph(char *filename);
void iftWriteSeedsGraph(iftLabeledSet* seed, const char *_filename, ...);



//! swig(stable)
void iftMWriteSeeds(iftLabeledSet *seed, iftMImage *image, char *filename);

//! swig(stable)
void iftWriteSeeds(iftLabeledSet *seed, const iftImage *image, const char *filename, ...);

//! swig(newobject)
iftImage	  *iftSeedImageFromLabeledSet(iftLabeledSet* labeled_set, iftImage *image);

//! swig(newobject, stable)
iftLabeledSet *iftLabeledSetFromSeedImage(iftImage* seed_image, bool decrement);
//! swig(newobject)
iftLabeledSet *iftLabeledSetFromSeedImageMarkersAndHandicap(iftImage* seed_image, iftImage *marker, iftImage *handicap);

/**
 * @brief This function computes the segmentation error components and relabels them to ensure that each one is given a
 * unique id (refactored from iftBorderMarkersForPixelSegmentation)
 *
 * @author Thiago Vallin Spina
 *
 * @param gt_image Ground truth image
 * @param label Segmentation result. It may be NULL, in which case we return an error component image with one component
 * per label
 * @param adj_relabeling Adjacency relation considered for relabeling the error components
 */
iftImage *iftRelabelSegmentationErrorComponents(iftImage *gt_image, iftImage *label, iftAdjRel *adj_relabeling);
iftLabeledSet *iftBorderMarkersForPixelSegmentation(iftImage *grad_image, iftImage *gt_image, float border_distance);
iftLabeledSet *iftGeodesicMarkersForSegmentation(iftImage *gt_image, iftImage *label);
iftLabeledSet* iftBorderMarkersForSuperpixelSegmentation(iftImage* label_image,iftImage* gt_image, iftDataSet* dataset);

//Pops the first "nelem" elements from "lset" with label "label"
//! swig(newobject)
iftLabeledSet* iftGetSeeds(iftLabeledSet* S, int nelem, int label);
//! swig(newobject)
iftLabeledSet* iftGetMisclassifiedSeeds(iftLabeledSet* S, int nelem, int label, iftImage* gt_image, iftImage* cl_image);

// Binary Segmentation only
int iftCheckNewSeeds(int *nelem, int length);

int iftSelectCircularRobotSeeds(iftImage *seed_image, iftBMap *used_seeds, iftImage *gt, double dist_border,
                                double max_marker_radius, double min_marker_radius, iftAdjRel *distance_border,
                                int center_seed);
//int iftMarkersFromMisclassifiedSeeds(iftImage* seed_image, iftLabeledSet* all_seeds, iftBMap* used_seeds, int nseeds,iftImage* gt_image, iftImage* cl_image, int dist_border, int max_marker_radius, int min_marker_radius);
int iftMarkersFromMisclassifiedSeeds(iftImage *seed_image, iftLabeledSet *all_seeds, iftBMap *used_seeds, int nseeds,
                                     int number_of_labels, iftImage *gt, iftImage *label, int dist_border,
                                     int max_marker_radius, int min_marker_radius);

void iftWriteSeedsOnImage(iftImage* image, iftLabeledSet* seed);

int iftRootVoxel(iftImage *pred, int p);

//! swig(newobject)
iftImage *iftFastLabelComp(const iftImage *bin, const iftAdjRel *Ain);

//! swig(newobject)
iftSet *iftBinaryMaskToSet(iftImage *mask);

//! swig(newobject)
iftImage *iftHBasins(iftImage *img, int H);
//! swig(newobject)
iftImage *iftHDomes(iftImage *img, int H);

/**
 * @brief Applies a grid sampling on a binary mask image according to a radius, selecting randomly n_samples.
 *
 * If the number of required samples is greater than the total number of found samples in the grid, 
 * all samples are returned.
 * If n_samples is <= 0, it also returns all samples from the grid.
 * 
 * @warning This function only works with a binary mask with a single connect component.
 * @warning Your program should have the statement 'iftRandomSeed(time(NULL));' in order to guarantee
 * a true random selection.
 * 
 * @note There is a program to compute this function: ift/demo/Miscellaneous/Sampling/iftGridSamplingOnMask.c
 * @note See also: ift/demo/Miscellaneous/Sampling/iftExtractPatchesByGeodesicGridSamplingOnMask.c
 *  
 * @param  bin_mask  Binary mask where the topological grid sampling is applied.
 * @param  radius    Radius/stride between the samples of the grid.
 * @param  initial_obj_voxel_idx Index of the initial object voxel to start the grid sampling.
 *                               If < 0, the first object voxel found in the image is considered.
 * @param  n_samples Number of required samples randomly extracted from the grid.
 *                   If it is greater than the total number of samples from the grid, all are returned.
 *                   If it <= 0, all samples is also returned.
 * @return           An array with the n_samples randomly chosen from the grid.
 * 
 * @author Samuka Martins
 * @date Apr 13, 2018
 */
//! swig(newobject)
iftIntArray *iftGridSamplingOnMask(const iftImage *bin_mask, float radius,
                                   int initial_obj_voxel_idx, long n_samples);

float iftEstimateGridOnMaskSamplingRadius
(const iftImage *binMask, int initialObjVoxelIdx, int nSamples);

/**
* @brief Performs a grid sampling on an multilabel image by assigning a proportional
* amount of seeds to each components area. The seeds are equally distributed within
* the object
* @author Felipe Belem
* @date June 28,2018
*
* @param label - Labeled image
* @param mask - Mask indicating reachable elements
* @param nSeeds - Number of desired seeds
*
* @return iftSet of seeds
*/
iftSet *iftMultiLabelGridSamplingOnMaskByArea
(const iftImage *label, iftImage *mask, int nSeeds );

/**
* @brief Given a labeled image, for each component, the method samples the centroid
*        if, and only if, the component area percentage is larger than the threshold
* @details In order to accept all kinds of labeled images, the method relabels the image
*          and computes the non-background area (thus, label > 0). Then, for every
*          connected component, the algorithm computes it proportional area (compared
*          with the total object area) and, if it is sufficiently large, samples the 
*          components' centroid.
* @author Felipe Belem
* @date Dec 4, 2018
* @ingroup Seeds
* 
* @param[in] label - Labeled image
* @param[in] mask - Mask indicating reachable elements (x > 0)
* @param[in] thresh - Area percentage threshold (0 <= x <= 1)
*
* @return Set of relevant (by area) centroid indexes
* 
* @see iftGeometricCentersFromLabelImage
*/
iftSet *iftMultiLabelCentroidSamplingOnMaskByArea
(const iftImage *label, iftImage *mask, float thresh);

/**
* @brief Performs a geodesic grid sampling on the thresholded grayscale object
*        saliency map, with respect to the percentage of seeds within the object.
* @details The method, through the thresholding, delimits the components belonging
*          to the background, and those belonging to the object(s). Thus, it performs
*          a geodesic grid sampling on each component, with respect to the percentage
*          of seeds within the object (also, indirectly, within the background).
* @author Felipe Belem
* @date Dec 4, 2018
* @ingroup Seeds
*
* @param[in] objsm - Grayscale object probability map
* @param[in] mask - Mask indicating reachable elements (x > 0)
* @param[in] k - Number of desired seeds
* @param[in] thresh - Threshold value for the map (0 <= x <= 1)
* @param[in] obj_perc - Percentage of seeds within the object (0 <= x <= 1)
*
* @return A binary image, where elements whose values are > 0,
*         indicates seeds
* @see iftMultiLabelGridSamplingOnMaskByArea for more detais
*/
iftImage *iftGrayObjMapGridSamplOnMaskByArea
(iftImage *objsm, iftImage *mask, int k, float thresh, float obj_perc );

/**
* @brief It oversegments the background, while for every object component, samples
*        its centroid if, and only if, its area is sufficiently large (in comparison
*        with other object components).
* @details The method, through the thresholding, delimits the components belonging
*          to the background, and those belonging to the object(s). Thus, for the
*          background, it samples, using a geodesic grid sampling approach, all the
*          desired seeds. For the object, it samples the component's centroid if its
*          proportional area is higher than the threshold established in parameter.
* @author Felipe Belem
* @date Dec 4, 2018
* @ingroup Seeds
* 
* @attention It may take too long if has too many noises! Filtering is advised!
*
* @param[in] objsm - Grayscale object probability map
* @param[in] mask - Mask indicating reachable elements (x > 0)
* @param[in] k - Number of desired seeds
* @param[in] map_thr - Threshold value for the map (0 <= x <= 1)
* @param[in] seed_thr - Area percentage threshold for the seed selection (0 <= x <= 1)
*
* @return A binary image, where elements whose values are > 0,
*         indicates seeds
* @see iftMultiLabelGridSamplingOnMaskByArea
* @see iftMultiLabelCentroidSamplingOnMaskByArea 
*/
iftImage *iftGrayObjMapCentrSamplOnMaskByArea
(iftImage *objsm, iftImage *mask, int k, float map_thr, float seed_thr );

/**
 * @brief Performs a grid sampling by choosing the voxels with maximum path values when running several
 * differential IFTs.
 * 
 * Initially, the central voxel of the image (or the geometric center of the object, if bin_mask != NULL)
 * is chosen to be the first voxel of the grid.
 * Then, it performs a differential IFT (DIFT) only with this voxel as seed.
 * At the end of the DIFT execution, the voxel with maximum cost in the cost map is added to the grid
 * and to the seed set and a new DIFT is executed by keeping the computed cost map.
 * This process is repeated until choosing the <n_samples> points of the grid.
 * 
 * @attention The final number of grid samples can be less than the input <n_samples>, since
 * the radius of the first grid sampling is automatically defined and there is no a strong guarantee that
 * there will be n_samples/2 samples with this radius inside the mask.
 * 
 * @param img Image where the grid sampling is performed.
 * @param bin_mask (Optional) Binary mask used to define the grid sampling only inside the object.
 *                 If NULL, the entire image is considered.
 * @param n_samples Number of required samples for the grid.
 * @return An array with the <n_samples> grid voxels chosen.
 * 
 * @author Samuel Martins
 * @date Sep 20, 2018
 */
iftIntArray *iftGridSamplingByMaxPathValues(const iftImage *img, const iftImage *bin_mask, long n_samples);


/**
 * @brief Performs a hybrid grid sampling on a binary mask.
 * 
 * This method selects half of the required samples by using the common topological grid sampling on mask (@see iftGridSamplingOnMask),
 * and from these points, it selects the remaining ones by using the max path values (@see iftGridSamplingByMaxPathValues).
 * For example, suppose we want n_samples = 100 samples.
 * 50 samples are initially selected using the common grid sampling. Then, we run a differential IFT (DIFT)
 * with these points as the seeds. The voxel of maximum path value of each tree is selected for the final grid.
 * Therefore, we have 50 + 50 = 100 selected samples.
 * 
 * The first grid is more conservative in the selection, and the second is a kind of "fine-tuning" to
 * select points with very different texture from the first ones.
 * 
 * @param img Image where the grid sampling is performed.
 * @param bin_mask Binary mask where the topological grid sampling is applied.
 * @param n_samples n_samples Number of required samples for the grid.
 * @return An array with the <n_samples> grid voxels chosen.
 * 
 * @author Samuel Martins
 * @date Sep 20, 2018
 */
iftIntArray *iftHybridGridSamplingOnMask(const iftImage *img, const iftImage *bin_mask, int n_samples);


/**
 * @brief Performs a hybrid grid sampling on a binary mask usind the cost function of ISF algorithm.
 * 
 * This method is similar to @see iftHybridGridSamplingOnMask, but with the cost function of the ISF algortihm:
 * Given a pixel p, with root r, trying to conquer a neighbor q, the cost function of the path extension is:
 * f(path_r-->p, <p, q>) = ((alpha * |I(q) - I(p)|) ^ beta) + dist(p, q)
 * 
 * where I(p) is the value of the pixel p, and dist(p, q) is the euclidean distance between p and q
 * 
 * @param img Image where the grid sampling is performed.
 * @param bin_mask Binary mask where the topological grid sampling is applied.
 * @param n_samples n_samples Number of required samples for the grid.
 * @param alpha Alpha factor of cost function.
 * @param alpha Beta factor of cost function.
 * @return An array with the <n_samples> grid voxels chosen.
 * 
 * @author Samuel Martins
 * @date Sep 28, 2018
 */
iftIntArray *iftHybridGridSamplingOnMaskISFCostFunction(const iftImage *img, const iftImage *bin_mask, int n_samples,
                                                        float alpha, float beta);


/**
 * @brief Apply a Grid Sampling for Patch Extraction where each grid point will be the central point
 * from its patch.
 * 
 * Each point of the grid is spaced each other by (stride_x, stride_y, stride_z).
 * The first grid point is the first one that fits the patch of size (patch_xsize, patch_ysize, patch_zsize).
 * It is guaranteed that the all patches centralized into the grid points will be entirely inside the image domain.
 * 
 * @param  img_dom  Domain of a given Image for grid sampling. Use the function iftGetImageDomain to get the domain of a given image.
 * @param  patch_xsize Patch's xsize.
 * @param  patch_ysize Patch's ysize.
 * @param  patch_zsize Patch's zsize.
 * @param  stride_x Stride (spacing) between the grid points on the x-axis.
 * @param  stride_y Stride (spacing) between the grid points on the y-axis.
 * @param  stride_z Stride (spacing) between the grid points on the z-axis.
 * @return          Integer array with the grid points.
 *
 * @author Samuel Martins
 * @date Jun 20, 2018.
 */
//! swig(newobject)
iftIntArray *iftGridSamplingForPatchExtraction(iftImageDomain img_dom, int patch_xsize, int patch_ysize,
                                               int patch_zsize, int stride_x, int stride_y, int stride_z);


/**
 * @brief Get a set of bounding boxes along an input image from an array of voxels.
 *
 * Each voxel, whose indices are in the integer array <voxel_indices>, corresponds to the central voxel
 * of a cubic bounding box of side <size>.
 * If a given bounding box is out of the image domain, its size and position are fit to it.
 * 
 * @param  img           Image used to get the voxel coordinate of the voxel indices.
 * @param  voxel_indices Indices of the central voxels on the image domain.
 * @param  size          Size of the cubic bounding boxes.
 * @return               Array with the found cubic bounding boxes along the image.
 *
 * @author Samuka Martins
 * @date Apr 25, 2018.
 */
//! swig(newobject)
iftBoundingBoxArray *iftBoundingBoxesAroundVoxels(const iftImage *img, const iftIntArray *voxel_indices, int size);

//! swig(newobject)
iftImage *iftLabelComponentsBySeeds(iftImage *comp, iftLabeledSet *seeds, bool incr);


/**
 * @brief Convert a voxel/region dataset to labeled seeds (voxels) by
 * using their true labels as seed label. In the case of region
 * (component) datasets, it extends the seed set to include all voxels
 * of regions marked with a true label in the dataset.
 *
 * @param  Z             the training dataset set
 * @param  comp          Image with labeled regions (samples in the dataset)
 * @return               Labeled seed set
 *
 * @author Alexandre Falcao
 * @date Aug 10, 2018.
 */

//! swig(newobject, stable)
iftLabeledSet *iftDataSetToLabeledSeeds(iftDataSet *Z, iftImage *comp);

/**
 * @brief Connect seed through the optimum path between them.
 *
 * @param mimg          Multi-band image
 * @param seeds         Input seed markers
 * @return              New labeled set with additional seeds
 *
 * @author Jordao Bragantini
 * @date Aug. 13, 2018
 */
//! swig(newobject, stable)
iftLabeledSet *iftConnectObjectSeeds(iftImage *img, iftLabeledSet *seeds);

/**
 * @brief Propagate true labeled samples to samples within the same cluster, given that the cluster is labeled with
 * truelabel above a given percentage
 *
 * @param Z             DataSet with group and true labels
 * @param truelabel     Label to be propagated
 * @param purity        Minimum percentage o labeling in each group
 * @return              New labeled set with additional, must extend to rest of the image
 *
 * @author Jordao Bragantini
 * @date Aug 13, 2018.
 */

//! swig(newobject)
iftLabeledSet *iftPropagateSeedsToCluster(const iftDataSet *Z, int truelabel, float purity);


/**
 * @brief Select seeds within a cluster/group which mixture percentage are above a given threshold
 * @param seeds         original markers
 * @param groups        groups mapping
 * @param label         label to be enhanced
 * @param threshold     threshold
 * @return              a subset of markers given restrictions
 *
 * @author Jordao Bragantini
 * @date Aug 20, 2018.
 */

//! swig(newobject, stable)
iftLabeledSet *iftSelectSeedsForEnhancement(iftLabeledSet *seeds, iftImage *groups, int label, float threshold);

  
/**
* @brief    Performs a sampling on the object saliency map, through the
*           selection of the extremities of the ordered values.
* @details  This method samples 'num_seeds' seeds by considering the values of the
*           grayscale object saliency map, in such a way that 'obj_perc'% of those
*           consists on the highest values (in the map). Therefore, the number of
*           seeds in the background consists on the remaining quantity of seeds 
*           in order to obtain the exact number of samples.
*           If a mask is given, it will not select the spels outside the mask. However,
*           it will still try to maintain the amount of seeds, by selecting those within it.
* @author   Felipe Belem
* @author   Leonardo Melo
* @author   Alexandre Falcao
* @date     Jan 18, 2019
*
* @param    objsm - Grayscale object saliency map.
* @param    mask - Mask inidicating the reachable spels (or NULL --- every spel is reachable).
* @param    num_seeds - Number of overall seeds.
* @param    obj_perc - Percentage of seeds in the object (i.e. higher values).
* @param    dist_penalty - Penalization factor greater than 0 (higher values --> higher concentrations).
* @return   Image in which the seeds have value greater than 0.
*
* @warning  This function was not tested for 3D images
* @warning  If the number of samples is extremely high (close to the number of
*           spels), some seeds might be selected twice (for bkg, and for obj)
* @see      iftObjSalMapSamplByHighestValue
*/
iftImage *iftSamplingByOSMOX(iftImage* objsm, iftImage *mask, int num_seeds, float obj_perc, float dist_penalty);

/**
* @brief    Performs a sampling on the object saliency map, taking into
*           account the highest values in it.
* @details  This method samples seeds by considering the values within
*           the object saliency map. First, it orders the values (decrease order),
*           and, for each spel selected, a gaussian penalization is performed in the 
*           adjacents (estimated by a pre-computed patch size --- of equal dimensions)
*           in order to spread the seeds (and avoid the selection of two extremelu close
*           seeds).
*           If a mask is given, it will not select the spels outside the mask. However,
*           it will still try to maintain the amount of seeds, by selecting those within it.
* @author   Felipe Belem
* @author   Leonardo Melo
* @author   Alexandre Falcao
* @date     Jan 18, 2019
*
* @param    objsm - Grayscale object saliency map.
* @param    mask - Mask inidicating the reachable spels (or NULL -- every spel is reachable).
* @param    num_seeds - Number of overall seeds.
* @param    dist_penalty - Penalization factor greater than 0 (higher values --> higher concentrations).
* @return   Set of seed indexes.
*
* @warning  This function was not tested for 3D images
*/
iftSet *iftObjSalMapSamplByHighestValue(iftImage *objsm, iftImage *mask, int num_seeds, float dist_penalty);

//! swig(newobject)
iftDoubleMatrix *iftSeedsFeatures(iftLabeledSet *seeds, const iftMImage *mimg, int label);

//! swig()
void iftSeedsSuperpixelBins(iftLabeledSet *seeds, const iftMImage *mimg, const iftImage *superpixel,
                            iftDoubleMatrix *data, iftIntArray *label);

//! swig(newobject)
iftImage *iftDrawDilatedSeeds(const iftImage *image, const iftLabeledSet *seeds,
                              const iftAdjRel *A, const iftColorTable *ctb_rgb);


/*
* @author   Alexandre Falcao
* @date     Feb, 2022
*
* @param    A binary image with components connected to its frame
* @return   new binary image without those components
*
*/
//! swig(newobject)
  iftImage *iftRemoveFrameComp(iftImage *bin);

void  iftRemoveFrameCompInPlace(iftImage *bin);



#ifdef __cplusplus
} 
#endif

#endif
