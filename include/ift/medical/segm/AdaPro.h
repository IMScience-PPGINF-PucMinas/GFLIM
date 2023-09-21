//
// Created by Samuel Martins on 06/12/18.
//

#ifndef IFT_ADAPRO_H
#define IFT_ADAPRO_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftFImage.h"
#include "iftImage.h"
#include "iftSVM.h"



/**
 * @brief Object Model of Adaptive Probabilistic Atlas AdaPro.
 * @author Samuka
 * @date Dec 17, 2017
 * @ingroup ObjectModel
 */
typedef struct ift_obj_model {
    /** Label of the target Object */
    int label;
    /** (Cropped) Prior Probability Atlas */
    iftFImage *prob_atlas;
    /** Shape (domain) of the Template Image, where the object model is built */
    iftImageDomain template_shape;
    /** Begin voxel of the cropped prob atlas on template image coordinate space */
    iftVoxel begin;
    /** Radius used to erode the prob. atlas to estimate the inner seeds */
    float erosion_radius;
    /** Radius used to dilate the prob. atlas to estimate the outer seeds */
    float dilation_radius;
} iftObjModel;


/**
 * @brief Adaptive Probabilist Atlas.
 * @author Samuka Martins
 * @date Mar 13, 2018
 */
typedef struct ift_adapro {
    /** Array of labels of each Object Model */
    iftIntArray *labels;
    /** Array of Object Models. */
    iftObjModel **obj_models;
    /** Rough Segmentation Mask resulting from the union of all prob. atlases after binarization and morphological closing */
    iftImage *rough_mask;
    /** Template (Reference Image) where the object models are trained (if required). */
    iftImage *template_img;
    /** Image with the training voxels (labeled voxels) for the linear SVM training. */
    iftImage *train_voxels_img;
    /** Linear SVM for texture Classification. */
    iftSVM *svm;
    /** Elastix Registration Files */
    iftStrArray *elastix_files;
} iftAdaPro;




/**
 * @brief Put the (cropped) prob atlas of an object model on its template image's domain from its begining voxel.
 *
 * @author Samuka
 * @date Nov 16, 2017
 * @ingroup ObjectModels
 */
iftFImage *iftProbAtlasOnTemplateImageDomain(const iftFImage *prob_atlas, iftImageDomain template_shape,
                                             iftVoxel begin);


/**
 * @brief Destroy an AdaPro Model.
 * @author Samuka
 * @date Mar 13, 2018
 */
void iftDestroyAdaPro(iftAdaPro **adapro);


/**
 * @brief Read an AdaPro Model.
 * @author Samuka
 * @date Mar 13, 2018
 */
iftAdaPro *iftReadAdaPro(const char *format, ...);


/**
 * @brief Write an AdaPro Model.
 * @author Samuka
 * @date Mar 13, 2018
 */
void iftWriteAdaPro(const iftAdaPro *adapro, const char *path, ...);


/**
 * @brief Train an Adaptive Probabilistic Atlases (AdaPro).
 *
 * Each object is indexed at a position i, so that each one of is parameter (erosion radius, ...)
 * is in the same position of the corresponding array.
 * If the label image with the training voxels <train_voxels_img> is NULL, the model won't be adaptive,
 * only relying on the shape constraints of the probabilistic atlases during segmentation
 *
 * @param atlas_set        Set of Registered Atlases (Label Images) on the template <template_img>.
 * @param template_img     Template Image (Reference Image Space) where the atlases are registered.
 * @param train_voxels_img Image with the markers chosen on template for linear SVM Training. BG voxels has label 0 and objects 1.
 *                         If NULL, the model won't be adaptive, only relying on the shape constraints of the probabilistic atlases.
 * @param labels           Array with the labels of the target objects for training the AdaPro.
 * @param e_radius_arr     (Optional) Array with the erosion radius to generate the inner seeds for each object.
 *                         If NULL, nothing is stored.
 * @param d_radius_arr     (Optional) Array with the dilation radius to generate the outer seeds for each object.
 *                         If NULL, nothing is stored.
 * @param C                Parameter C for linear SVM training.
 * @param elastix_fset     Elastix files used to register images to the adapro's template.
 * @return                 The trained AdaPro.
 *
 * @author Samuka Martins
 * @date Mar 13, 2018
 */
iftAdaPro *iftTrainAdaPro(const iftFileSet *atlas_set, const iftImage *template_img, const iftImage *train_voxels_img,
                          const iftIntArray *labels, const iftFloatArray *e_radius_arr,
                          const iftFloatArray *d_radius_arr, double C, const iftFileSet *elastix_fset);


/**
 * @brief Segment a test image by AdaPro. The test image and the adapro's template must be in the same
 * coordinate space.
 * @param  img                          Image to be segmented.
 * @param  adapro                       Adaptive Probabilistic Atlas.
 * @param  skip_texture_classification  If true, the a texture classification is ignored.
 * @param  aux_basename                 Auxiliar Basename used to save the intermediate steps (seed estimation, gradient,
 *                                      classification mask, etc) of the segmentation. If NULL, nothing is saved.
 *                                      Otherwise, only the shape-based segmentation is performed.
 * @return                              Segmented Image.
 *
 * @author Samuka
 * @date Mar 13, 2018
 */
iftImage *iftSegmentByAdaPro(const iftImage *img, const iftAdaPro *adapro, const char *aux_basename);


/**
 * @brief Writes the Elastix Parameter Files of an AdaPro model into the disk.
 *
 * @param adapro AdaPro model with the Elastix Parameter Files.
 * @param basename Optional basename for the elastix parameter files.
 * @return File set with the pathnames of the Elastix Parameter files.
 *
 * @author Samuel Martins
 * @date Jan 3, 2019
 */
iftFileSet *iftWriteAdaProElastixParamFiles(const iftAdaPro *adapro, const char *basename);


/**
 * @brief Register (by Elastix) a AdaPro Model on a test image's space. Resulting registrations
 * and mapping are assigned in the own input model.
 *
 * Reference image is registerd with the test one. This will be the new reference image. \n
 * Then, all object models are mapped to new space using the deformation fields. \n
 * Deformation fields are not saved.
 *
 * @param adapro       AdaPro Model.
 * @param test_img     Testing Image.
 * @param elastix_files Elastix Parameter Files.
 * @param def_fields_out Optional reference to return the resulting deformation fields. Pass NULL to ignore it.
 *
 * @author Samuka
 * @date Mar 18, 2018
 * @ingroup ObjectModels
 */
void iftRegisterAdaProOnTestImageByElastix(iftAdaPro *adapro, const iftImage *test_img, const iftFileSet *elastix_files,
                                           iftFileSet **def_fields_out);


/**
 * @brief Finds the Object Model's Seeds for delineation.
 * It assumes that the test image is registered on the template image (or vice-versa), used to
 * build the object model
 *
 * The borders of the dilated and eroded certainty region of the prob. atlas form the background and object seeds, respectively.
 * A texture classification mask/membership map can be passed to identify the regions to be adapted into the model.
 * Thus, voxels inside the (dilated) uncertainty and certainty regions (rough segmentation) and classified as background in the classification mask
 * are considered "forbidden for delineation".
 *
 * @param  obj_model          Object Model.
 * @param  clf_mask           Classification mask (Membership Map) of the test image used to filter the seeds. If NULL, no filtering is applied.
 * @param  rough_mask         Rough segmentation (Volume of Interest) used to define forbidden voxels inside it
 *                            classified as background by the classification mask. If NULL or clf_mask is NULL, nothing is done.
 * @param  forbidden          Reference to save the set with the forbidden voxels identified by the texture classification mask.
 *
 * @author Samuel Martins
 * @date Mar 18, 2018
 * @ingroup ObjectModels
 */
iftLabeledSet *iftFindObjModelSeeds(const iftObjModel *obj_model, const iftImage *clf_mask, iftSet **forbidden);


/**
 * @brief Sets the elastix files to an adapro model.
 *
 * @param adapro AdaPro model where the elastix files are set.
 * @param elastix_fset Elastix files to be set.
 *
 * @author Samuel Martins
 * @date Nov 3, 2018
 * @ingroup ObjectModels
 */
void iftSetElastixFilesToAdaPro(iftAdaPro *adapro, const iftFileSet *elastix_fset);


/**
 * @brief Dilates the AdaPro's rough segmentation by its largest dilation radius.
 *
 * @param adapro AdaPro Model.
 * @return Dilated Rough Segmentation.
 *
 * @author Samuel Martins
 * @date Nov 10, 2018
 * @ingroup ObjectModels
 */
iftImage *iftDilateAdaProRoughSegmentation(const iftAdaPro *adapro);

#ifdef __cplusplus
}
#endif

#endif //IFT_ADAPRO_H
