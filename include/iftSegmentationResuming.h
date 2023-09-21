//
// Created by tvspina on 1/18/16.
//

#ifndef IFT_IFTSEGMENTATIONRESUMING_H
#define IFT_IFTSEGMENTATIONRESUMING_H

#include "ift.h"

///**
// * @brief This function converts an input label image into an optimum-path forest by computing markers whose result
// * using DIFT produce the original label as close as possible. The markers are selected using a robot user that add
// * seeds close to the objects' boundaries in increasing order of gradient value. Hence, seeds are added on regions
// * with low boundary information first, which are problematic for DIFT. Candidate pixels are selected from a given
// * distance to the boundary, and circular markers are drawn centered on them with a radius that at most touches the
// * objects boundaries.
// *
// * @author Thiago Vallin Spina
// *
// * @param gradient Gradient image values
// * @param label Original input label image
// * @param seeds_per_iteration Number of seeds that are added for each label at each iteration (e.g., 2 seeds per object per iteration)
// * @param border_erosion_radius_for_candidate_seeds Radius used to erode the objects' masks in order to select seeds.
// * @param min_safe_distance_to_border This is a padding value that constrains the maximum marker radius. If greater than
// * 0 then the markers will not touch the objects' boundaries
// * @param max_marker_radius Maximum marker radius
// * @param min_marker_radius Minimum markers radius (it may be 0, in which case the markers may be a single pixel)
// * @param stopping_threshold Accuracy/Error threshold selected to stop the algorithm
// * @param iftStoppingCriterion This function evaluates the segmentation accuracy of all objects and returns true if
// * the resuming algorithm should stop because all objects have reached the given accuracy/error threshold
// */
//iftLabeledSet *iftLabelToForestPixelRobot(iftImage *gradient, iftImage *label, iftAdjRel *A, int seeds_per_iteration,
//                                          float border_dilation_radius_for_candidate_seeds,
//                                          int min_safe_distance_to_border, int min_marker_radius, int max_marker_radius,
//                                          iftLabeledSet *optional_input_seeds, double stopping_threshold,
//                                          int secondary_stopping_threshold,
//                                          uchar (*iftStoppingCriterion)(iftImage *, iftImage *, double),
//                                          uchar (*iftSecondaryStoppingCriterion)(int, int, int));

/**
 * @brief This function computes the maximum Dice accuracy among all objects in GT and returns true if it is more than or
 * equal to the threshold. It may be used with iftLabelToForestPixelRobot to stop the label reconstruction.
 *
 * @param label Label image being computed
 * @param gt Original label image
 * @param threshold Minimum accepted accuracy for all objects
 */
uchar iftStopLabelReconstructionByDice(iftImage *label, iftImage *gt, double threshold);
/**
 * @brief This function computes the minimum ASSD error among all objects in GT and returns true if it is less than or
 * equal to the threshold. It may be used with iftLabelToForestPixelRobot to stop the label reconstruction.
 *
 * @param label Label image being computed
 * @param gt Original label image
 * @param threshold Maximum accepted error for all objects
 */
uchar iftStopLabelReconstructionByASSD(iftImage *label, iftImage *gt, double threshold);


/**
 * @brief Stops label reconstruction when a certain number of iterations is performed.
 *
 * @author Thiago Vallin Spina
 * @date Jan 27, 2016
 *
 * @param cur_iteration The current iteration number.
 * @param cur_num_seeds The current number of added seeds.
 * @param threshold The maximum number of allwed iterations.
 *
 * @return True or false when the criterion is met.
 */
uchar iftStopLabelReconstructionByNumInterations(int cur_iteration, int cur_num_seeds, int threshold);

/**
 * @brief Stops label reconstruction when a certain number of seeds is added.
 *
 * @author Thiago Vallin Spina
 * @date Jan 27, 2016
 *
 * @param cur_iteration The current iteration number.
 * @param cur_num_seeds The current number of added seeds.
 * @param threshold The maximum number of seeds that is allowed to be added.
 *
 * @return True or false when the criterion is met.
 */
uchar iftStopLabelReconstructionByNumSeeds(int cur_iteration, int cur_num_seeds, int threshold);

/* EDT-based robot */

/**
 * @brief Finds and returns an iftLabeledSet with the coordinates of all border spels of a label image for a given adjacency <A>,
 * which include spels that belong to *any* label including the background. That is, background spels adjancet to object
 * spels will also be considered as border spels. If a spel is on the image's border but belongs to the background
 * then they are not considered as border spels, while the opposite is true for object spels.
 *
 * @author Thiago Vallin Spina
 *
 * @param label Label image.
 * @param A Adjacency relation for considering the frontier between labels
 * @return The iftSet with the coordinates of all border spels of the label image <label>
 */
iftLabeledSet * iftLabelBorderSet(iftImage *label, iftAdjRel *A);


iftLabeledSet *iftBorderMarkersDistAndRootForPixelSegmentation(iftImage *grad_image, iftImage *gt_image,
                                                               double border_distance, iftImage **dist,
                                                               iftImage **root);

/**
 * @brief This function selects seed voxels on the border between foreground and background according to the euclidean
 * distance that each pixel has to it. It takes as input a voxel on the border of the object or the background, which
 * serves as the "center" of the marker. Then, it selects as seed voxels those whose root on the border is within a
 * given radius to the center voxel and belong to the same error component. Hence, the resulting marker follows the
 * contour instead of being a simple circle/sphere.
 *
 * @author Thiago Vallin Spina
 *
 * @param center_marker_voxel The voxel on the object's or background's border that has been selected as the marker's "center".
 * @param gt_image The ground truth image.
 * @param dist The EDT distance image. NOTE: This function selects ALL voxels with finite distance as seeds, hence, <dist>
 * should be created only until a given radius from the border to ensure that the markers will be small. Ideally, this
 * distance should also be <marker_radius>.
 * @param root The EDT's root image.
 * @param seed_image The seed image where the seeds will be placed on.
 * @param error_components The image with labeled error components.
 * @param marker_radius The radius used to select the roots around the center voxel.
 * @param min_border_distance A minimum distance to the border that must be respected in order to select seed voxels.
 *
 */
ulong iftSelectSeedsBasedOnEDT(iftVoxel center_marker_voxel, iftImage *gt_image, iftImage *dist, iftImage *root,
                               iftImage *error_components, iftImage *seed_image, iftImage *seed_image_mk_id,
                               double marker_length, double min_border_distance, int new_seed_id,
                               bool select_seeds_on_axes_only);
/**
 * @brief This function recursively removes seeds from a given set if they are within a given radius of a given voxel.
 *
 *
 * @param Thiago Vallin Spina
 *
 * @param S Seed set.
 * @param center The center voxel.
 * @param label The center voxel's label
 * @param img The image in which the seed voxels have been selected.
 * @param radius The radius for removal.
 */
void iftRemoveSeedsWithinRadius(iftLabeledSet **S, iftVoxel center, int label, iftImage *img, double radius);
iftSet **iftErrorComponentsPerSegmentationLabelSortedByArea(iftImage *gt, iftImage *error_components);
int iftMarkersFromMisclassifiedComponentsSortedByAreaAndSeedEDT(iftImage *seed_image, iftImage *seed_image_mk_id,
                                                                iftLabeledSet **all_seeds, iftImage *gt,
                                                                iftImage *label, iftImage *dist, iftImage *root,
                                                                int nseeds_per_object, double marker_length,
                                                                double min_border_distance, iftAdjRel *A,
                                                                bool select_seeds_on_axes_only);

void iftDrawSeeds(iftImage *img, iftImage *seed_image, iftColorTable *cmap);
iftLabeledSet *iftLabelToForestPixelRobotEDT(iftImage *gradient, iftImage *label, iftAdjRel *A, int nseeds_per_label_per_iteration,
                                             double min_safe_distance_to_border, double max_marker_width, double max_marker_length,
                                             iftLabeledSet *optional_input_seeds, double stopping_threshold,
                                             int secondary_stopping_threshold,
                                             uchar (*iftStoppingCriterion)(iftImage *, iftImage *, double),
                                             uchar (*iftSecondaryStoppingCriterion)(int, int, int), bool select_seeds_on_axes_only);

iftLabeledSet *iftLabelToForestPixelRobotEDTOrientedWatershed(iftImageForest *fst, iftImage *gt_image,
                                                              int nseeds_per_label_per_iteration,
                                                              double min_safe_distance_to_border,
                                                              double max_marker_width, double max_marker_length,
                                                              double gamma, iftLabeledSet *optional_input_seeds,
                                                              double stopping_threshold,
                                                              int secondary_stopping_threshold,
                                                              uchar (*iftStoppingCriterion)(iftImage *, iftImage *,
                                                                                            double),
                                                              uchar (*iftSecondaryStoppingCriterion)(int, int, int),
                                                              bool select_seeds_on_axes_only);

iftLabeledSet *iftLabelToForestPixelRobotEDTIGraphOrientedWatershed(iftIGraph *igraph, iftImage *gt_image,
                                                                    int nseeds_per_label_per_iteration,
                                                                    double min_safe_distance_to_border,
                                                                    double max_marker_width, double max_marker_length,
                                                                    double gamma, iftLabeledSet *optional_input_seeds,
                                                                    double stopping_threshold,
                                                                    int secondary_stopping_threshold,
                                                                    uchar (*iftStoppingCriterion)(iftImage *,
                                                                                                  iftImage *, double),
                                                                    uchar (*iftSecondaryStoppingCriterion)(int, int,
                                                                                                           int),
                                                                    bool select_seeds_on_axes_only);

//iftLabeledSet *iftLabelToForestGeodesicRobotByACC(iftImage *gradient, iftImage *label, iftAdjRel *A,
//                                                  int seeds_per_iteration, int min_distance_border, int max_marker_size,
//                                                  int min_marker_size, iftLabeledSet *optional_input_seeds,
//                                                  double stopping_threshold, int secondary_stopping_threshold,
//                                                  uchar (*iftStoppingCriterion)(iftImage *, iftImage *, double),
//                                                  uchar (*iftSecondaryStoppingCriterion)(int, int, int));

iftLabeledSet *
        iftLabelToForestISF_Root(const iftImage *img, const iftImage *input_label, double alpha, double beta,
                                         double gamma, int niters, int *nseeds, iftIGraph **igraph_out, iftDHeap **Q_out);

iftImage *iftRespSystemGradient(iftImage *img, iftAdjRel *A);

void iftResumeImageSegmentation(iftImage *orig, iftImage *input_label, iftImage *basins, iftAdjRel *A, iftDict *json,
                                iftLabeledSet *optional_input_seeds, const char *out_dir, char *out_img_basename);

void iftIGraphDiffISF_Resuming_Root(iftIGraph *igraph, const iftImage *input_label, iftLabeledSet *new_seeds,
                                    iftSet **trees_for_removal, iftDHeap *Q, double alpha, double beta, double gamma,
                                    bool allow_label_fixing);

void iftDiffOrientedWatershedResuming(iftImageForest *fst, iftImage *input_label, iftLabeledSet *seed,
                                      iftSet * removal_markers, double gamma);

void iftIGraphResumingDiffOrientedWatershed(iftIGraph *igraph, iftImage *input_label, iftLabeledSet *seeds,
                                            iftSet *trees_for_removal, iftDHeap *Q, double gamma);

//! swig()
int iftLocalMaximum(const iftImage *weights, int p, const iftAdjRel *disk);

//! swig()
int iftFurthestInError(const iftImage *source, const iftImage *target, const iftImage *mask, int p);

//! swig(newobject)
iftSet *iftQueryForAnchorsPosition(const iftImage *gt_contour, const iftImage *source_mask,
                                   const iftImage *gt_mask, const iftSet *anchors);

//! swig(newobject)
iftSet *iftFurtherThanThreshold(const iftVoxelArray *anchors,
                                const iftImage *mask,
                                float threshold);


#endif //IFT_IFTSEGMENTATIONRESUMING_H
