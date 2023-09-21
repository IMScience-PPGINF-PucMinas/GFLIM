#include "iftSaliency.h"


typedef struct ift_tensor_scale {
    iftFloatArray *anisotropy;
    iftFloatArray *orientation;
    iftFloatArray *minor_axis;
    iftFloatArray *major_axis;
    iftFloatArray *area;
    iftVoxelArray *pos_focus, *neg_focus;
    int m_pairs;
    long n;
} iftTensorScale;

iftSaliencyArray *iftFocusPrior(iftImage *superpixel_img, iftMImage *mimg, float variance);
iftSaliencyArray *iftFocusXRayPrior(iftImage *superpixel_img, iftMImage *mimg, float variance, float border_threshold);
iftSaliencyArray *iftImageCenterSaliencyPrior(iftImage *superpixel_img, float variance);
iftSaliencyArray *iftForegroundDistanceSaliencyPrior(iftImage *superpixel_img, iftImage *saliency, float variance);
iftSaliencyArray *iftForegroundDistanceSaliencyPriorNew(iftImage *superpixel_img, iftImage *saliency, float variance);
iftSaliencyArray *iftForegroundFeatureDistanceSaliencyPrior(iftImage *superpixel_img, iftImage *saliency, iftMImage *normalized_feats, float variance);
iftSaliencyArray *iftEllipseMatchingSaliencyPrior(iftImage *superpixel_img, float variance, int min_size, int max_size);
iftSaliencyArray *iftImageCenterSaliencyPrior(iftImage *superpixel_img, float variance);
iftSaliencyArray *iftImageScribbleDistanceSaliencyPrior(iftImage *superpixel_img, iftLabeledSet *scribbles, float variance);
void addAllPriorsToGraph(iftSaliencyGraph **saliency_graph_ref, iftImage *superpixel_img, iftMImage *features, iftImage *saliency, iftITSELFParameters *params);


/**
 * @brief Add a saliency array as a prior in the saliency graph.
 * @date March 02, 2020
 * @author Leo Melo
 * @note You can add multiple priors to the graph and combine them using iftCuboidPriorIntegrationImplicit.
 *
 * @param  ref_saliency_graph The reference saliency graph to add the prior to.
 * @param  new_prior The prior array containing a saliency value to each superpixel.
 * @return      void.
 *
 */
void iftAddPriorToGraph(iftSaliencyGraph **ref_saliency_graph, iftSaliencyArray *new_prior);

/**
 * @brief Auxiliary function to iftAddPriorToGraph. Add a saliency array as a prior in the saliency graph.
 * @date March 02, 2020
 * @author Leo Melo
 * @note This function is called by iftAddPriorToGraph.
 *
 * @param  ref_saliency_graph The reference saliency graph to add the prior to.
 * @param  new_prior The prior array containing a saliency value to each superpixel.
 * @param  npriors The current number of priors inside the graph.
 * @return      void.
 *
 */
void iftAppendPriorToList(iftSaliencyArray ***prior_list, iftSaliencyArray *new_prior, int npriors);

/**
 * @brief Remove all priors from a saliency graph.
 * @date March 02, 2020
 * @author Leo Melo
 *
 * @param  ref_saliency_graph The reference saliency graph to remove the priors from.
 * @return      void.
 *
 */
void iftResetPriorList(iftSaliencyGraph **ref_saliency_graph);

void iftDestroyPriorList(iftSaliencyGraph **ref_saliency_graph);

iftSet *iftBoundingBoxToSeedSet(iftBoundingBox bb, iftImage *image);
iftImage *iftFillEllipse(iftImage *label_img, int label, iftVoxel *new_centroid, float *new_area);
iftTensorScale *iftSuperpixelToTensorScale(iftImage *label_img, int m_pairs, int min_size, int max_size);
void iftDestroyTensorScale(iftTensorScale **tensor_scale_ref);