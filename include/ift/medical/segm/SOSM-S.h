//
// Created by Samuel Martins on 06/12/18.
//

#ifndef IFT_SOSM_S_H
#define IFT_SOSM_S_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftFImage.h"
#include "iftImage.h"
#include "ift/core/dtypes/LabeledSet.h"



/**
 * @brief SOSM-S Model of a single Object.
 * @author Samuka
 * @date Dec 13, 2016
 * @ingroup ObjectModel
 */
typedef struct ift_sosms_obj_model {
    /** Label of the target Object */
    int label;
    /** (Cropped) Prior Probability Atlas */
    iftFImage *prob_atlas;
    /** Shape (domain) of the Template Image, where the object model is built */
    iftImageDomain template_shape;
    /** Begin voxel of the cropped prob atlas on template image coordinate space */
    iftVoxel begin;
    // search region for the object localization from the reference voxel
    iftBoundingBox search_region;
} iftSOSMSObjModel;


/**
 * @brief SOSM-S Model.
 * @author Samuka
 * @date Dec 13, 2016
 * @ingroup ObjectModel
 */
typedef struct ift_sosm_s {
    /** Array of labels of each SOSM-S Object Model */
    iftIntArray *labels;
    /** Array of Object Models. */
    iftSOSMSObjModel **obj_models;
    /** Template (Reference Image) where the object models are trained (if required). */
    iftImage *template_img;
} iftSOSMS;



/**
 * @brief Destroys a SOSM-S Model.
 * @author Samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
void iftDestroySOSMS(iftSOSMS **sosm_s);


/**
 * @brief Stores on disk a SOSM-S Model.
 * @author Samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
void iftWriteSOSMS(const iftSOSMS *sosm_s, const char *path);


/**
 * @brief Reads from disk a SOSM-S Model.
 * @author Samuka
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 */
iftSOSMS *iftReadSOSMS(const char *path);


/**
 * @brief Trains a SOSM-S Model [1] from a set of registered segmentation masks.
 *
 * [1] Phellan, 2016 - Medical physics - Medical image segmentation via atlases and fuzzy object models
 *
 * @param  template_img Template Image (Reference Image).
 * @param  labels       Array with the object labels
 * @return              SOSM-S Model.
 *
 * @author sAMUKA
 * @date Dec 20, 2016
 * @ingroup ObjectModels
 *
 */
iftSOSMS *iftTrainSOSMS(const iftFileSet *atlas_set, const iftImage *template_img,
                        const iftIntArray *labels);


/**
 * @brief Segments a Testing Image using the Statistical Multi-Object Shape Model SOSM-S (Phellan, 2016).
 *
 * It applies an Object Location by MSPS translating the seed models over the test image gradient. \n
 * Phellan, 2016 - Medical physics - Medical image segmentation via atlases and fuzzy object models
 *
 * @param  test_img Test Image to be segmented.
 * @param  sosm_s   Statistical Object Shape Model.
 * @return          Segmented Image.
 *
 * @author Samuka
 * @date Jan 3, 2016
 * @ingroup ObjectModels
 */
iftImage *iftSegmentBySOSMS(const iftImage *img, const iftImage *grad_img_in, const iftSOSMS *sosm_s);


/**
 * @brief Finds the SOSM-S Object Model's Seeds for delineation.
 * It assumes that the test image is registered on the template image (or vice-versa), used to
 * build the object model
 *
 * Given a Model (probabilistic map), background's seeds are those with prob 0 and object's seeds
 * with prob 1. \n
 *
 * @param  test_img           Testing Image.
 * @param  obj_model          SOSM-S Object Model.
 * @return                    Labeled seeds for delineation.
 *
 * @author Samuka
 * @date Jan 3, 2017
 * @ingroup ObjectModels
 */
iftLabeledSet *iftFindSOSMSSeeds(const iftImage *test_img, const iftSOSMSObjModel *obj_model);

#ifdef __cplusplus
}
#endif


#endif //IFT_SOSM_S_H
