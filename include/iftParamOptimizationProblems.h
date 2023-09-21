/**
 * @file iftParamOptimizationProblems.h
 * @brief Definitions and functions for building problems for parameter optimization.
 * @author Samuel Martins
 * @date Jan 01, 2017
 */

#ifndef IFT_PARAM_OPTIMIZATION_PROBLEMS_H
#define IFT_PARAM_OPTIMIZATION_PROBLEMS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/LabeledSet.h"
#include "ift/imgproc/dtypes/BoundingBox.h"
#include "ift/medical/segm/SOSM-S.h"
    
#include "iftCommon.h"
#include "iftImage.h"
#include "iftMSPS.h"
#include "iftSimilarity.h"



/**
 * @brief Problem struct for Gradient Matching optimization.
 * @author Samuka
 * @date Jan 2, 2017
 */
typedef struct ift_grad_match_problem {
    /** Fixed gradient image */
    const iftImage *fixed_grad_img;
    /** Moving Gradient Image that will try to match with the fixed one */
    const iftImage *mov_grad_img;
    /** Array with the coordinates of all voxels with value > 0 (for the Moving Gradient Image) */
    iftIntArray *mov_grad_elems;
    iftVector max_disp;
} iftGradMatchProb;


/**
 * @brief Problem struct for optimizations with Shape Models.
 * @author Samuka
 * @date Jan 2, 2017
 */
typedef struct ift_obj_model_mean_arc_weight {
    const iftImage *grad_test_img;
    int label;
    const iftLabeledSet *seeds;
    const iftSet *certain_obj_region;
    const iftSet *forbidden;
    iftVector max_disp;
} iftObjModelMeanArcWeight;


typedef struct ift_patch_localization_nmi_prob {
    iftImage *test_img;
    iftImage *template_img;
    iftBoundingBox true_template_patch;
    iftBoundingBox test_patch;
    iftBoundingBox search_region;
} iftPatchLocalizationNMIProb;


/**
 * @brief Build the MSPS Delta Matrix with regular stride for each scale, according to: 1 + scale_index*stride.
 *
 * Ex: n_params = 3, n_scales = 4, stride = 2 \n
 * msps->delta = [[1, 1, 1],          \n
 *                [3, 3, 3],          \n
 *                [5, 5, 5],          \n
 *                [7, 7, 7]]          \n
 * 
 * @param  n_params Number of Parameter to be optimized by MSPS.
 * @param  n_scales Number of Scales used by MSPS.
 * @param  stride   Stride factor used to compute the Regular Delta Matrix.
 * @return          Regular MSPS Delta Matrix.
 *
 * @author Samuka
 * @date Jan 2, 2017
 */
iftMatrix *iftRegularMSPSDeltaMatrix(int n_params, int n_scales, int stride);


/**
 * @brief Allocates a struct for Gradient Matching Problem.
 * Fixed and Moving images are only assigned to struct (not copied). They must have the same domain and voxel size.
 * 
 * @param  fixed_grad_img Fixed Gradient Image.
 * @param  mov_grad_img   Moving Gradient Image that will try to match with the fixed one.
 * @return                Created Gradient Matching Problem Struct.
 *
 * @author Samuka
 * @date Jan 2, 2017
 */
iftGradMatchProb *iftCreateGradMatchProb(const iftImage *fixed_grad_img, const iftImage *mov_grad_img, iftVector max_disp);


/**
 * @brief Destroys a Grad Matching Problem Struct.
 * @author Samuka
 * @date Jan 2, 2017
 */
void iftDestroyGradMatchProb(iftGradMatchProb **prob);


iftPatchLocalizationNMIProb *iftCreatePatchLocalizationNMIProb(const iftImage *test_img,
                                    const iftImage *ref_img, iftBoundingBox true_ref_patch,
                                    iftBoundingBox patch, iftBoundingBox search_region);


void iftDestroyPatchLocalizationNMIProb(iftPatchLocalizationNMIProb **prob);


float iftPatchNMI(void *prob, float *theta);


iftBoundingBox iftMSPSMaxPatchNMI(const iftImage *img, const iftImage *template_img, iftBoundingBox true_template_patch,
                                        iftBoundingBox patch, iftBoundingBox search_region);

/**
 * @brief Fitness function to find the best matching between the gradients whose score is the
 * mean value of the intersection between the gradients.
 *
 * For all translated voxels from moving gradient image with values > 0, it takes the mean value of the Fixed Image. \n
 * Such function must be maximized to find the optimum matching.\n
 * The array <b>theta</b> contains the displacement for the coord x, y, and z (only for 3D images). \n
 * Each voxels from moving image is translated before matching.
 * 
 * @param  problem Gradient Matching Problem.
 * @param  theta   Array with the displacement for the coordinates x, y, and z (only for 3D images) which
 *                 is used to translate the moving gradient's voxels before matching.
 * @return         Resulting score.
 *
 * @author Samuka
 * @date Jan 2, 2017
 */
float iftMatchGradProbMean(void *problem, float *theta);


/**
 * @brief Find the best gradient matching by translating (in x, y, z) the fixed grad image over the
 * moving one and maximizing the mean value of their intersections.
 * Images must have the same domain and voxel size.
 *
 * Stride is used to compute a regular translation for each scales, according to: 1 + scale_index*stride. \n
 * For example, n_params = 3 (x, y, z), n_scales = 4, stride = 2 \n\n
 * msps->theta = [0, 0, 0] // initial translation vector \n\n
 * 
 * // perturbation of each parameter (coordinate/column) in each scale 
 * msps->delta = [[1, 1, 1],          \n
 *                [3, 3, 3],          \n
 *                [5, 5, 5],          \n
 *                [7, 7, 7]]          \n
 * 
 * @param  fixed_grad_img Fixed Gradient Image.
 * @param  mov_grad_img   Moving Gradient Image that will try to match with the fixed one.
 * @param  n_scales       Number of scales for optimization.
 * @param  stride         Stride used to compute a regular translation (displacement) of the moving
 *                        image's coordinates for each scale.
 * @return                Resulting displacement vector of the best matching.
 *
 * @author Samuka
 * @date Jan 2, 2017
 */
iftVector iftMSPSMatchGradients(const iftImage *fixed_grad_img, const iftImage *mov_grad_img, int n_scales, int stride, iftVector max_disp);



iftObjModelMeanArcWeight *iftCreateObjModelMeanArcWeightProb(const iftImage *grad_test_img, int label,
                                                             const iftLabeledSet *seeds, const iftSet *certain_obj_region,
                                                             const iftSet *forbidden, iftVector max_disp);


/**
 * @brief Destroys an Optimization Problem with Object Model.
 * author Samuka
 * @date Jan 2, 2017
 */
void iftDestroyObjModelMeanArcWeightProb(iftObjModelMeanArcWeight **prob);



iftVector iftMSPSObjModel(const iftImage *grad_test_img, int label, const iftLabeledSet *seeds,
                          const iftSet *certain_obj_region, const iftSet *forbidden, 
                          int n_scales, int stride, iftVector max_disp);


float iftComputeObjModelMeanArcWeight(void *problem, float *theta);


// iftBoundingBox iftMSPSBestNMIBoundingBox







#ifdef __cplusplus
}
#endif

#endif
