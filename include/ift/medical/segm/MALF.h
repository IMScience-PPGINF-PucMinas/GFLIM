//
// Created by Samuel Martins on 06/12/18.
//

#ifndef IFT_MALF_H
#define IFT_MALF_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/io/FileSet.h"
#include "iftImage.h"


/**
 * @brief Selects the best atlases to a test image by Normalized Mutual Information [1].
 * Train set and test image must be already registered.
 *
 * @warning All images should be registered in the same coordinate space.
 *
 * If masks are passed for the training images and testing image, they are masked using their corresponding mask
 * before normalized mutual information. \n *
 * The higher the NMI, more similar (better) the train image is. \n
 * [1] Aljabar, 2009 - Neuroimage - Multi-atlas based segmentation of brain images: atlas selection and its effect on accuracy
 *
 * @param test_img        Testing Image.
 * @param train_img_set   Training Image Set.
 * @param n_atlases       Number of atlases to be selected.
 * @param test_img        (Optional) Testing Mask.
 * @param train_mask_set  (Optional) Set of masks for the training set. The pathnames must be in the same order of the training images.
 * @param selected_idxs_out (Optional) Reference to return the indices from the selected atlases.
 * @return     A file set with the pathnames from the selected images.
 *
 * @author Samuka
 * @date Jan 27, 2017
 * @ingroup ObjectModels
 */
iftFileSet *iftAtlasSelectionByNMI(const iftImage *test_img, const iftFileSet *train_img_set, int n_atlases,
                                   const iftImage *test_mask, const iftFileSet *train_mask_set,
                                   iftIntArray **selected_idxs_out);


/**
 * @brief Segments an Image by Classical MALF. This function considers the training label set is already registered on testing image.
 *
 * Classical MALF registers all atlases (img + label img) on test image's space. \n
 * The label of each voxel from segmented image is the most frequent label (Majority Voting) \n
 *
 * @param  test_img        Testing Image to be segmented.
 * @param  train_label_set Training label set already registered on testing image's space.
 * @return                 Segmented image.
 *
 * @author Samuka
 * @date Jan 27, 2017
 * @ingroup ObjectModels
 */
iftImage *iftSegmentByClassicalMALF(const iftImage *test_img, const iftFileSet *train_label_set);


/**
 * @brief Segments an Image by MALF using STAPLE as Label Fusion.
 * This function considers that the training label set is already registered on testing image.
 *
 * @warning The CRKIT, which has the STAPLE binary program, must be installed before.
 *
 * @param  test_img        Testing Image to be segmented.
 * @param  train_label_set Training label set already registered on testing image's space.
 * @return                 Segmented image.
 *
 * @author Samuka
 * @date Jan 27, 2017
 * @ingroup ObjectModels
 */
iftImage *iftSegmentByMALFSTAPLE(const iftImage *test_img, const iftFileSet *train_atlas_set);

#ifdef __cplusplus
}
#endif


#endif //IFT_MALF_H
