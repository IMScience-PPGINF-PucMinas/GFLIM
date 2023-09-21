/**
 * @file AsymmetryAnalysis.h
 * @brief Definitions and functions for Brain Asymmetry Analysis.
 * @author Samuel Martins
 * @date Apr 28, 2018
 */

#ifndef IFT_ASYMMETRY_ANALYSIS_H
#define IFT_ASYMMETRY_ANALYSIS_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"
#include "iftDataSet.h"
#include "iftImage.h"
#include "iftMImage.h"


/**
 * @brief Enum with flags for Right and Left Brain Side
 * @author Samuel Martins
 * @date Dec 11, 2018
 */
//! swig()
typedef enum {
    IFT_RIGHT_BRAIN_SIDE, IFT_LEFT_BRAIN_SIDE
} iftBrainSide;


/**
 * @brief Computes the voxel-wise absolute asymmetries between the brain right and left sides according to its Mid-Sagittal Plane (MSP).
 *
 * After computing the brain absolute asymmetries, it can SUBTRACT a bias image
 * (e.g. the normal asymmetry map) to attenuate brain asymmetries specially on the cortex.
 * Voxels with negative values after this attenuation are set to 0.
 *
 * @attention The function assumes that the image's MSP is its central sagittal slice (xslice).
 *
 * @param img Brain image to compute its asymmetries.
 * @param bias (Optional) Bias image used to attenuate brain asymmetries. Use NULL to ignore it.
 * @return Resulting brain absolute asymmetries.
 *
 * @author Samuel Martins
 * @date Dec 12, 2018
 */
//! swig(newobject)
iftImage *iftBrainAsymMap(const iftImage *img, const iftImage *bias);



/**
 * @brief Computes the Mean Brain Absolute Asymmetries for a given image set.
 * 
 * For each image, it computes the absolute asymmetries between the brain sides based on its Mid-Sagittal Plane (MSP).
 * The resulting map is the mean brain absolute asymmetries from the image set.
 * Optionally, it can add the standard deviation from the brain asymmetries to this map.
 *
 * @warning All images must have the same domain, voxel sizes, orientation, and their MSPs must be the central sagittal slice.
 * 
 * @param img_set Set of images;
 * @param add_stdev_asymmetries Add the standard deviation asymmetries to the resulting mean asymmetry map. 
 * @return Mean Asymmetry Map.
 * 
 * @author Samuel Martins
 * @date Dec 4, 2018
 */
//! swig(newobject)
iftImage *iftMeanBrainAsymMap(const iftFileSet *img_set, bool add_stdev_asymmetries);


/**
 * @brief Computes the Mean Brain Difference Asymmetries for a given image set and a template.
 *
 * For each image, it computes the absolute difference between it and a template.
 * The resulting map is the mean brain absolute differences from the image set.
 * Optionally, it can add the standard deviation from the brain asymmetries to this map.
 *
 * @warning All images must be registered in the template.
 *
 * @param img_set Set of images;
 * @param template_img Template where all images are registered;
 * @param add_stdev_asymmetries Add the standard deviation asymmetries to the resulting mean asymmetry map.
 * @return Mean Asymmetry Map.
 *
 * @author Samuel Martins
 * @date Jul 14, 2019
 */
//! swig(newobject)
iftImage *iftMeanBrainDiffMap(const iftFileSet *img_set, const iftImage *template_img, bool add_stdev_asymmetries);


/**
 * @brief Applies a Grid Sampling by using the brain asymmetries.
 * 
 * Given an asymmetry map of a brain, it is thresholded in order to get the asymmetries as binary objects.
 * Such objects are eroded (removing small objects) to decrease the number of grid points on the asymmetric
 * regions. The geometric centers of the n_points_on_asymmetries resulting objects are selected as grid points.
 * 
 * The remaining points are selected on the symmetric region. For that, the functions tries to sample
 * n_points_on_symmetries = (n_samples - n_points_on_asymmetries) points equally distanced on this regions.
 * If n_points_on_symmetries < min_samples_on_symmetries, then n_points_on_symmetries = min_samples_on_symmetries.
 * 
 * The final grid points are the union of the points on the symmetric and asymmetric regions.
 * 
 * @param asym_map Asymmetry Map of a brain.
 * @param bin_mask Mask of a side of the Hemisphere.
 * @param n_samples Number of required samples/points. Note that resulting number of samples
 *                  can be less or greater than this number. 
 * @param min_samples_on_symmetries Minimum number of super samples/points on the symmetric region.
 * @param thres_factor Factor used to multiply the otsu threshold for binarizing the asymmetry map in the seed initialization.
 * @return Resulting grid points.
 * 
 * @author Samuel Martins
 * @date Oct 1, 2018
 */
//! swig(newobject)
iftIntArray * iftGridSamplingByBrainAsymmetries(const iftImage *asym_map, const iftImage *bin_mask,
                                                int min_samples_on_symmetries, float thres_factor);


/**
 * @brief Builds a Multi-Band Image by stacking the right and the flipped left brain sides from an input image.
 * 
 * @attention The function assumes that the image's Mid-Sagittal Plane is its central xslice.
 * 
 * The resulting Multi-Band Image has 2 bands: the first is the input image, and the second is the flipped input image.
 * By assuming that the mid-sagittal plane is the central xslice of the input image, then each voxel
 * will have two values: one from the right brain side and the other from the flipped left brain side.
 * 
 * @param img Brain image.
 * @return Multi-Band Image with the stacked brain hemispheres.
 * 
 * @author Samuel Martins
 * @date Sep 1, 2018
 */
 //! swig(newobject)
iftMImage *iftBuildBrainHemisphereMImage(const iftImage *img);


/**
 * @brief Builds a Multi-Band Image by stacking the right and the flipped left brain sides from an input image, and their
 * absolute asymmetries (differences).
 *
 * @attention The function assumes that the image's Mid-Sagittal Plane is its central xslice.
 *
 * The resulting Multi-Band Image has 3 bands: the first is the input image, the second is the flipped input image,
 * and the third is the absolute asymmetry map between the first 2 bands.
 * By assuming that the mid-sagittal plane is the central xslice of the input image, then each voxel
 * will have 3 values.
 *
 * @param img Brain image.
 * @return Multi-Band Image with the stacked brain hemispheres.
 *
 * @author Samuel Martins
 * @date Dec 14, 2018
 */
//! swig(newobject)
iftMImage *iftBuildBrainHemisphereMImageAsym(const iftImage *img, const iftImage *asym_map);


/**
 * @brief Extracts the HAA feats (Histogram of Absolute Asymmetries) from an asymmetry map.
 *
 * This functions extracts HAA feats from an asymmetry map for each of its supervoxels.
 * The extracted supervoxels' feature vectors are stored in a set of datasets, one per supervoxel.
 * By assuming a set of <n_supervoxels> supervoxels, an array <Zarr> of <n_supervoxels> + 1 (for convenience), and
 * that the sample index from the asymmetry map is <s>, the extracted feature vector for each supervoxel i is stored
 * in the dataset Zarr[i] at the sample <s>.
 *
 * @warning The number of features of all supervoxel datasets must be equal to the number of histogram bins.
 *
 * @param asym_map Map with the absolute asymmetries between the brain sides.
 * @param svoxels_brain_side Label image with the symmetric supervoxels (from 1 to n_supervoxels) only from one brain side.
 * @param Zarr Array of datasets, one per supervoxel.
 * @param s Sample index in the datasets for the input asymmetry map.
 * @param n_bins Number of bins of the histograms.
 * @param max_val Maximum value used for histogram computing.
 * @param n_supervoxels Number of supervoxels from the Supervoxel label image.
 *
 * @author Samuel Martins
 * @date Dec 11, 2018
 */
void iftSupervoxelHAAFeatsInAsymMap(const iftImage *asym_map, const iftImage *svoxels_brain_side,
                                    iftDataSet **Zarr, int s, int n_bins, int max_val, int n_supervoxels);


/**
 * @brief Extract symmetrical supervoxels on the Brain.
 * 
 * @attention The function assumes that the image's Mid-Sagittal Plane is its central xslice.
 * 
 * It expects a binary mask with the target object from a brain side (according to its
 * Mid-Sagittal Plane - MSP), e.g. a mask with right hemisphere. Then, it extracts supervoxels, so that the
 * same supervoxels inside the binary mask are mirrored in the other brain side.
 * 
 * To find/determine corresponding regions in both brain sides (based on the target binary mask),
 * it firstly builds a multi-band image with 2 bands by stacking the right brain side with the mirrored left one,
 * so that both sides have the same orientation.
 * Such MImage has voxels with 2 values, one for each brain side.
 * 
 * Then, it extracts supervoxels from this MImage by using ISF only inside the binary mask.
 * The resulting supervoxels/regions that associates both brain sides are mirrored so that
 * the resulting label image has the same supervoxels in both brain sides.
 * 
 * The SymmISF's seed initialization uses the absolute asymmetries between the brain sides (based on its MSP).
 * If a normal asymmetry map is passed, such asymmetries are attenuated by it.
 * 
 * @param img Input Brain Image.
 * @param bin_mask Binary Mask with the target object from a given brain side.
 * @param alpha Alpha factor of SymmISF.
 * @param beta Beta factor of SymmISF.
 * @param thres_factor Factor used to multiply the otsu threshold for binarizing the asymmetry map in the seed initialization.
 * @param min_dist_to_border Minimum euclidean distance to the binary mask borders that the initial seeds must have.
 * @param n_seeds_on_symmetric_regions Number of Seeds on Symmetryc Regions.
 * @param normal_asymmap (Optional) Normal asymmetry map used to attenuate the brain asymmetries for
 *                                  the seed initialization for SymmISF. Pass NULL to ignore it.
 * @return Label image with the resulting supervoxels.
 * 
 * @author Samuel Martins
 * @date Sep 1, 2018
 */
//! swig(newobject)
iftImage *iftSymmISF(const iftImage *img, const iftImage *bin_mask, float alpha, float beta, float thres_factor,
                     float min_dist_to_border, int n_seeds_on_symmetric_regions, const iftImage *normal_asymmap);


/**
 * @brief Extract symmetrical supervoxels on the Brain by using OISF.
 *
 * @attention The function assumes that the image's Mid-Sagittal Plane is its central xslice.
 *
 * This function is very similar to @see iftSymmISF, but using the OISF supervoxel method instead ISF.
 * Note that OISF requires an extra argument: gamma.
 *
 * @param img Input Brain Image.
 * @param bin_mask Binary Mask with the target object from a given brain side.
 * @param n_supervoxels Required number of supervoxels.
 * @param alpha Alpha factor of OISF.
 * @param beta Beta factor of OISF.
 * @param beta Gamma factor of OISF.
 * @param thres_factor Factor used to multiply the otsu threshold for binarizing the asymmetry map in the seed initialization.
 * @param normal_asymmap (Optional) Normal asymmetry map used to attenuate the brain asymmetries for
 *                                  the seed initialization for SymmISF. Pass NULL to ignore it.
 * @return Label image with the resulting supervoxels.
 *
 * @author Samuel Martins
 * @date May 6, 2019
 */
//! swig(newobject)
iftImage *iftSymmOISF(const iftImage *img, const iftImage *bin_mask, int n_supervoxels, float alpha, float beta,
                      float gamma, const iftImage *normal_asymmap, float thres_factor);



/**
 * @brief Extracts a brain side (right or left based on its mid-sagittal plane - MSP) from a brain image.
 *
 * The resulting image has the same domain and voxel sizes from input one.
 *
 * @attention It assumes that the MSP is the central sagittal slice (xslice) from the image.
 *
 * @param img Input Brain Image.
 * @param side Chosen brain side.
 * @return Image with the extracted brain side
 *
 * @author Samuel Martins
 * @date Dec 11, 2018
 */
//! swig(newobject)
iftImage *iftExtractBrainSide(const iftImage *img, iftBrainSide side);

/**
 * @brief Builds a dataset for each pair of symmetric supervoxels (extracted from a test image) by computing its
 * Histogram of the Absolute Asymmetries (HAA) as feature vector.
 *
 * For each pair of symmetric supervoxels in the label image <test_sym_supervoxels> extracted from the test image <test_img>,
 * a dataset is built by computing the Histogram of Absolute Asymmetries (HAA feats) with <n_bins> bins between the pair of
 * symmetric supervoxes for all training images and the test image.
 * Then, each supervoxel dataset has n_train_samples + 1 samples, where n_train_samples is the number of training images
 * in <train_asym_set>.
 * The testing sample is the last sample in all supervoxel datasets.
 *
 * The asymmetry maps are computed by the absolute difference between the brain sides from the images according to
 * their mid-sagittal plane (MSP).
 * A normal asymmetry map can be passed to only attenuate the TEST image asymmetries.
 *
 * @note For programming convenience, the resulting array has n_supervoxels + 1 datasets, where n_supervoxels
 * is the number of pair of symmetric supervoxels, so that the index [i] contains the dataset for the symmetric
 * supervoxels with label i.
 * @note The testing sample is the last sample in all supervoxel datasets.
 * @attention The function assumes that the test image's MSP is its central sagittal slice (xslice).
 *
 * @param test_img Test Image.
 * @param test_sym_svoxels Label image with pairs of symmetric supervoxels extracted from the test image.
 * @param train_set Training image set.
 * @param n_bins Number of the histogram bins.
 * @param normal_asym_map (Optional) Normal asymmetry map used to only attenuate the TEST image asymmetries. Pass NULL to ignore it.
 * @param n_svoxels_out Reference to return the number of supervoxels in the test label image.
 * @return Array of datasets, one for each pair of symmetric supervoxels.
 * 
 * @author Samuel Martins
 * @date Sep 1, 2018
 */
//! swig(newobject)
iftDataSet **iftExtractSupervoxelHAAFeats(const iftImage *test_img, const iftImage *test_sym_svoxels,
                                          const iftFileSet *train_set, int n_bins, const iftImage *normal_asym_map,
                                          int *n_svoxels_out);

/**
 * @brief Builds a dataset for each pair of symmetric supervoxels (extracted from a test image) by computing its
 * BIC descriptor as feature vector.
 *
 * For each pair of symmetric supervoxels in the label image <test_sym_supervoxels> extracted from the test image <test_img>,
 * a dataset is built by computing the BIC desciptor with <n_bins_per_channel> number of bins per image channel.
 * Then, each supervoxel dataset has n_train_samples + 1 samples, where n_train_samples is the number of training images
 * in <train_asym_set>.
 * The testing sample is the last sample in all supervoxel datasets.
 *
 * The asymmetry maps are computed by the absolute difference between the brain sides from the images according to
 * their mid-sagittal plane (MSP).
 * A normal asymmetry map can be passed to only attenuate the TEST image asymmetries.
 *
 * @note For programming convenience, the resulting array has n_supervoxels + 1 datasets, where n_supervoxels
 * is the number of pair of symmetric supervoxels, so that the index [i] contains the dataset for the symmetric
 * supervoxels with label i.
 * @note The testing sample is the last sample in all supervoxel datasets.
 * @attention The function assumes that the test image's MSP is its central sagittal slice (xslice).
 *
 * @param test_img Test Image.
 * @param test_sym_svoxels Label image with pairs of symmetric supervoxels extracted from the test image.
 * @param train_set Training image set.
 * @param n_bins_per_channel Number of bins per channel.
 * @param normal_asym_map (Optional) Normal asymmetry map used to only attenuate the TEST image asymmetries. Pass NULL to ignore it.
 * @param n_svoxels_out Reference to return the number of supervoxels in the test label image.
 * @return Array of datasets, one for each pair of symmetric supervoxels.
 * 
 * @author Samuel Martins
 * @date Oct 2, 2019
 */
//! swig(newobject)
iftDataSet **iftExtractSupervoxelBICAsymmFeats(const iftImage *test_img, const iftImage *test_sym_svoxels,
                                               const iftFileSet *train_set, int n_bins_per_channel, const iftImage *normal_asym_map,
                                               int *n_svoxels_out);

/**
 * @brief Given a pre-computed array of supervoxel datasets built from a training set and a single test image, it builds and
 * sets a CSV reference data to each dataset.
 *
 * The function assumes that each supervoxel dataset were built (feature extraction) with exactly the same training
 * image set, test image, and label image with the pair of symmetric supervoxels from the test image.
 * Then, it creates a CSV for each supervoxel dataset that contains the following columns:
 * "image_path,symmetric_supervoxel_path,supervoxel_label,status"
 *
 * Therefore, each sample [s] in the supervoxel datasets contains its original image path, the considered symmetric
 * supervoxel image and its corresponding label used for feature extraction.
 *
 * @note For programming convenience, the input array has n_supervoxels + 1 datasets, where n_supervoxels
 * is the number of pair of symmetric supervoxels, so that the index [i] contains the dataset for the symmetric
 * supervoxels with label i.
 *
 * @param Zarr Array of datasets, one for each pair of symmetric supervoxels.
 * @param n_supervoxels Number of supervoxels considered.
 * @param test_img_path Test image's pathname.
 * @param test_supervoxels_path Pathname from the label image with the pair of symmetric supervoxels.
 * @param train_set Training image set used for feature extraction.
 */
//! swig(newobject)
void iftBuildRefDataSupervoxelDataSets(iftDataSet **Zarr, int n_supervoxels, const char *test_img_path,
                                       const char *test_supervoxels_path, const iftFileSet *train_set);

#ifdef __cplusplus
}
#endif

#endif
