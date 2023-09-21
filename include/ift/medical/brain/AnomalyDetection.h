
#ifndef IFT_ANOMALYDETECTION_H
#define IFT_ANOMALYDETECTION_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/DblArray.h"
#include "ift/core/dtypes/IntArray.h"
#include "iftDataSet.h"
#include "iftImage.h"
#include "iftMImage.h"



/**
 * @brief Computes the magniture of the registration error --- voxel-wise absolute difference -- between a
 * registered image and its template.
 *
 * After computing the error, the function can SUBTRACT a bias image
 * (e.g. the normal registration error magnitude map) to attenuate brain differences specially on the cortex.
 * Voxels with negative values after this attenuation are set to 0.
 *
 * @attention The function expects that the image is registered on the template.
 *
 * @param img Registered Image.
 * @param template_img Template Image (Reference Image).
 * @param bias (Optional) Bias image used to attenuate the registration error Magnitude. Use NULL to ignore it.
 * @return Resulting Registration Error Magnitude.
 *
 * @author Samuel Martins
 * @date Sep 16, 2017
 */
//! swig(newobject)
iftImage *iftRegErrorMagnitude(const iftImage *img, const iftImage *template_img, const iftImage *bias);



/**
 * @brief Computes the Mean Registration Error Magnitude for a given registered image set.
 *
 * For each registered image, it computes the registration error magnitude between the image and its template.
 * The resulting map is the mean registration error magnitude from the image set.
 * Optionally, it can add the standard deviation from the registration error magnitudes to this map.
 *
 * @warning All images must be registered in the template.
 *
 * @param img_set Set of registered images
 * @param template_img Template Image (Reference Image).
 * @param add_stdev_error Add the standard deviation reg error magnitudes to the resulting mean map.
 * @return Mean Registration Error Magnitude.
 *
 * @author Samuel Martins
 * @date Sep 16, 2017
 */
//! swig(newobject)
iftImage *iftMeanRegErrorMagnitude(const iftFileSet *img_set, const iftImage *template_img, bool add_stdev_error);


/**
 * @brief Creates a 2-band image by stacking an input registered image and its template.
 * @author Samuel Martins
 * @date Sep 16, 2017
 */
iftMImage *iftBuildImageTemplateMImage(const iftImage *img, const iftImage *template_img);

/**
 * Perform ISF supervoxel segmentation on each object of a label image whose seeds are generated from registration
 * errors between a registered image and its template.
 *
 * For each object in the label image, the function performs ISF on it. Its initial seeds are obtained from
 * the maxima of the registration error magnitudes and from a set of seeds on the remaining region.
 * ISF is performed on a 2-band image obtained by stacking the registered image and its template.
 * All supervoxel maps are combined into a single image so that the supervoxels are relabeled.
 *
 * @param img Input Registered image.
 * @param reg_error_mag Registration Error Magnitudes from the input Registered image.
 * @param template_img Template image.
 * @param label_img Label image with the target objects.
 * @param alphas Array with the ISF alpha factors, one for each object. Index [i] results the value for object i.
 * @param betas Array with the ISF beta factors, one for each object. Index [i] results the value for object i.
 * @param thres_factors Array with the factors, one for each object, used to increase/decrease the otsu threshold.
 *                      Index [i] results the value for object i.
 * @param min_dist_to_borders Array with the min. distance to the object's border that the corresponding initial seeds
 *                            must have. Index [i] results the value for object i.
 * @param n_seeds_on_correct_regions Array with the number of seeds on the corrected region of each target object.
 *                                   Index [i] results the value for object i.
 * @return Supervoxel Image.
 *
 * @author Samuka Martins
 * @date Sep 18, 2019
 */
//! swig(newobject)
iftImage *iftISFOnRegErrors(const iftImage *img, const iftImage *reg_error_mag, const iftImage *template_img,
                            const iftImage *label_img, const iftDblArray *alphas, const iftDblArray *betas,
                            const iftDblArray *thres_factors, const iftDblArray *min_dist_to_borders,
                            const iftIntArray *n_seeds_on_correct_region);

//! swig(newobject)
iftImage *iftISFOnRegErrorsFast(const iftImage *img, const iftImage *reg_error_mag, const iftImage *template_img,
                                const iftImage *label_img, const iftDblArray *alphas, const iftDblArray *betas,
                                const iftDblArray *thres_factors, const iftDblArray *min_dist_to_borders,
                                const iftIntArray *n_seeds_on_correct_region);

    /**
 * @brief Applies a Grid Sampling by using the maxima on domes of an image..
 *
 * Given an image, the function first threshold it by using a passed threshold factor, resulting in a binary image
 * with the domes of the image.
 * Such domes are then eroded (removing small ones) to decrease the number of grid points.
 * The final grid points consists of: (i) the maxima of each dome;  (ii) a set of n points on the flat region (i.e. the complement
 * of the binary image) is uniformly distributted.
 *
 * @param img Image.
 * @param bin_mask Binary mask that delimits where the grid points will be selected.
 * @param n_samples_on_flat_region Number of grid points selected on flat regions.
 * @param thres_factor Threshold factor used to binarize the image and obtain the domes.
 * @return Grid points.
 *
 * @author Samuka Martins
 * @date Sep 18, 2019
 */
    //! swig(newobject)
    iftIntArray *iftGridSamplingOnDomes(const iftImage *img, const iftImage *bin_mask, int n_samples_on_flat_region,
                                        float thres_factor, iftImage **domes_out);

void iftSupervoxelHistFeats(const iftImage *img, const iftImage *svoxels_img, iftDataSet **Zarr, int s, int n_bins,
                            int n_supervoxels);


//! swig(newobject)
iftDataSet **iftExtractSupervoxelHistRegErrorsFeats(const iftImage *test_reg_error_mag, const iftImage *test_svoxels_img,
                                                    const iftFileSet *train_set, int n_bins, int *n_svoxels_out);

//! swig(newobject)
iftDataSet **iftExtractSupervoxelHistRegErrorsFeatsDilation(const iftImage *test_img, const iftImage *test_svoxels_img,
                                                            const iftFileSet *train_set, const iftImage *template_img,
                                                            int n_bins, float radius, const iftImage *bias, int *n_svoxels_out);

//! swig(newobject)
iftDataSet **iftExtractSupervoxelBandHistFeats(const iftMImage *test_filt_img, const iftImage *test_svoxels_img,
                                               const iftFileSet *train_set, int n_bins, int *n_svoxels_out);

//! swig(newobject)
void iftWriteSupervoxelDataSets(const iftDataSet **Zarr, int n_supervoxels, const char *out_dir);





//! swig(newobject)
iftFImage *iftComputeLinearAttenuationWeightsByEDT(const iftImage *label_img, float max_attenuation_factor);

//! swig(newobject)
iftFImage *iftComputeExponentialAttenuationWeightsByEDT(const iftImage *label_img, float max_attenuation_factor, float exponent);

//! swig(newobject)
iftImage *iftWeightedRegErrorMagnitude(const iftImage *img, const iftImage *template_img, const iftFImage *weights);

//! swig(newobject)
iftImage *iftRemoveSVoxelsByVolAndMeanRegError(const iftImage *svoxels_img, const iftImage *reg_error_mag,
                                               int min_vol, float min_mean_reg_error_on_svoxel);

//! swig(newobject)
iftImage *iftISFOnAttentionMap(const iftImage *img, const iftImage *attention_map, const iftImage *target_img,
                               const iftImage *label_img, const iftDblArray *alphas, const iftDblArray *betas,
                               const iftDblArray *thres_factors, const iftDblArray *min_dist_to_borders,
                               const iftIntArray *n_seeds_on_correct_region);


//! swig(newobject)
iftDataSet **iftExtractSupervoxelBICFeats(const iftImage *test_img, const iftImage *test_svoxels_img,
                                          const iftFileSet *train_set, int n_bins_per_channel,
                                          int *n_svoxels_out);

//! swig(newobject)
iftDataSet **iftExtractSupervoxelAttentionMapFeats(const iftImage *test_attention_map, const iftImage *test_svoxels_img,
                                                   const iftFileSet *train_set, int n_bins, int *n_svoxels_out);


//! swig(newobject)
iftDataSet **iftExtractSupervoxelLBPFeats(const iftImage *test_img, const iftImage *test_svoxels_img,
                                          const iftFileSet *train_set, int n_bins, int *n_svoxels_out);

//! swig(newobject)
iftDataSet **iftExtractSupervoxelVLBPFeats(const iftImage *test_img, const iftImage *test_svoxels_img,
                                           const iftFileSet *train_set, int n_bins, int *n_svoxels_out);


//! swig(newobject)
iftDataSet **iftExtractSupervoxelTextureFeats(const iftImage *test_img, const iftImage *test_attention_map,
                                              const iftImage *test_svoxels_img, const iftFileSet *train_set,
                                              const iftFileSet *train_attention_map_set,
                                              int n_bins, int *n_svoxels_out);


#ifdef __cplusplus
}
#endif


#endif //IFT_ANOMALYDETECTION_H
