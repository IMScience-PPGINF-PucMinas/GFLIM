/**
 * @file iftSimilarity.h
 * @brief Functions to compute Similarities or Distances between images.
 * @author Samuel Martins
 * @ingroup ImageScore
 */

#ifndef IFT_SIMILARITY_H
#define IFT_SIMILARITY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftMatrix.h"
#include "iftMetrics.h"
#include "iftSeeds.h"


/**
 * @brief Computes the Similarity between two Binary Images using Dice.
 * @author Samuel Martins
 * @date Sep 1st, 2015
 * @ingroup ImageMetrics
 *
 * Computes the Dice Similarity between two Binary Images using Dice.\n
 * The binary image must have only the value 0 and another value (not necessarily 0 and 1).\n
 *
 * @param bin_source First Binary Image.
 * @param bin_target Second Binary Image. Commonly the ground truth.
 * @return The dice similarity between bin_source and bin_target.
 *
 * @note The dice similarity goes from 0 (no similarity) and 1 (perfect similarity).
 * 
 * @warning The function does not check if the input images are really Binary Images.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
 //! swig()
double iftDiceSimilarity(const iftImage *bin_source, const iftImage *bin_target);


/**
 * @brief Computes the similarity between two Label Images using Dice.
 * @author Samuel martins
 * @date Feb 18, 2016
 * @ingroup ImageMetrics
 *
 * The (Binary) Dice Similarity is computed for each object (1..N).\n
 * The result is stored in the ith position of a double array, that indeed is its the object index (label value).\n
 * Index 0 holds the average for all objects.\n\n
 * 
 * E.g:\n
 * output[0] = Average Dice; \n
 * output[1] = Dice of the object 1 (its voxel value is 1)\n
 * output[N] = Dice of the object N (its voxel value is N)\n
 * 
 * @param label_img_source Source Label Image.
 * @param label_img_target Target Label Image. Commonly the ground truth.
 * @param n_objects Number of objects of the input label images.
 * @return A double array with the dices. Index 0 has the average dice and ith index has the dice of
 * the ith object.
 *
 * @note The dice similarity goes from 0 (no similarity) and 1 (perfect similarity).
 * 
 * @warning The function does not check if the input images are really Label Images with <b>n_objects</b> objects,
 * whose labels goes from 0 to n_objects.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
iftDblArray *iftDiceSimilarityMultiLabel(const iftImage *label_source,
                                         const iftImage *label_target, int n_objects);


/**
 * @brief Computes the Average Surface Distance between two Binary Images.
 * @author Samuel Martins
 * @date Feb 21, 2016
 * @ingroup ImageMetrics
 *
 * Computes the Mean Euclidean Distance between the boundary spels from the Source Binary Image to the boundary
 * spels of the Target Binary Image.\n
 * For each boundary spel from the Source Binary Image, it computes the minimum euclidean distance to a boundary spel
 * of the Target Binary Image. Then, the average minimum distance is computed and returned.
 * EDT (implemented with IFT) is used to compute the euclidean distances.\n
 * 
 * @param bin_source Source Binary Image.
 * @param bin_target Target Binary Image.
 * @return The Mean Euclidean Distance.
 * 
 * @note The perfect segmentation produces ASD <b>0</b> (no error - perfect similarity).
 * @note It is not symmetric because it only computes the distance <b>one-way</b> (from the sourceto the target binary image).
 * @note Reference: http://mbi.dkfz-heidelberg.de/grand-challenge2007/sites/eval.htm
 * 
 * @warning The function does not check if the input images are really Binary Images.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
double iftASD(const iftImage *bin_source, const iftImage *bin_target);


/**
 * @brief Computes the Average Surface Distance between two Label Images with label from 0 to N.
 * @author Samuel Martins
 * @date Feb 21, 2016
 * @ingroup ImageMetrics
 *
 * The (Binary) ASD is computed for each object (1..N).\n
 * The result is stored in the ith position of a double array, that indeed is its the object index (label value).\n
 * Index 0 holds the average for all objects.\n\n
 * 
 * E.g:\n
 * output[0] = Average ASD; \n
 * output[1] = ASD of the object 1 (its voxel value is 1)\n
 * output[N] = ASD of the object N (its voxel value is N)\n
 *
 * @param label_source Source Label Image.
 * @param label_target Target Label Image.
 * @param n_objects Number of Objects of the two Label Images.
 * @return A double array with all ASD results. Index 0 has the average ASD and ith index has the ASD of
 * the ith object.
 *
 * @note The perfect segmentation produces ASD <b>0</b> (no error - perfect similarity).\n
 * @note It is not symmetric because it only computes the distance <b>one-way</b> (from the sourceto the target binary image).
 * @note Reference: http://mbi.dkfz-heidelberg.de/grand-challenge2007/sites/eval.htm
 * 
 * @warning The function does not check if the input images are really Label Images with <b>n_objects</b> objects,
 * whose labels goes from 0 to n_objects.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
iftDblArray *iftASDMultiLabel(const iftImage *label_source, const iftImage *label_target, int n_objects);


/**
 * @brief Computes the Average Symmetric Surface Distance (ASSD) between two Binary Images.
 * @author Samuel Martins
 * @date Sep 14th, 2015
 * @ingroup ImageMetrics
 *
 * Computes the Average Distance from Source Binary Image to the Target Binary Image and vice versa (Symmetric).\n
 * The binary image must have only the value 0 and another value.
 * 
 * @param bin_source Source Binary Image.
 * @param bin_target Target Binary Image.
 * @return The resulting ASSD.
 * 
 * @note The perfect segmentation produces ASD <b>0</b> (no error - perfect similarity).
 * @note The binary image must have only the value 0 and another value.
 * @note Reference: http://mbi.dkfz-heidelberg.de/grand-challenge2007/sites/eval.htm
 * 
 * @warning The function does not check if the input images are really Binary Images.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
 //! swig(newobject)
double iftASSD(const iftImage *bin_source, const iftImage *bin_target);


/**
 * @brief Computes the Average Symmetric Surface Distance (ASSD) between two Label Images.
 * @author Samuel Martins
 * @date Feb 20, 2016
 * @ingroup ImageMetrics
 *
 * The (Binary) ASSD is computed for each object (1..N).\n
 * The result is stored in the ith position of a double array, that indeed is its the object index (label value).\n
 * Index 0 holds the average for all objects.\n\n
 * 
 * E.g:\n
 * output[0] = Average ASSD; \n
 * output[1] = ASSD of the object 1 (its voxel value is 1)\n
 * output[N] = ASSD of the object N (its voxel value is N)\n
 *
 * @param label_source Source Label Image.
 * @param label_target Target Label Image.
 * @param n_objects Number of Objects of the two Label Images.
 * @return A double array with the ASD. Index 0 has the average ASD and ith index has the ASD of
 * the ith object.
 * 
 * @note The perfect segmentation produces ASD <b>0</b> (no error - perfect similarity).
 * @note Reference: http://mbi.dkfz-heidelberg.de/grand-challenge2007/sites/eval.htm
 * 
 * @warning The function does not check if the input images are really Label Images.
 * To do that, use the program demo/Miscellaneous/iftCheckLabelImages.c
 */
iftDblArray *iftASSDMultiLabel(const iftImage *label_source, const iftImage *label_target, int n_objects);


/**
 * // TODO: Refactor this function and create a function to Multilabel
 * @brief Returns the Sum of the Gradient of a given border points.
 * @author Samuel Martins
 * @date Mar 8, 2016 
 */
double iftBorderGradMatchingScore(iftSet *borders, iftImage *grad);

/**
 * // TODO: Refactor this function and create a function to Multilabel
 * @brief Returns the sum of the maximum gradient values on the region of influence of a given set of border points.
 *
 * The region of influence is given by the euclidean distance transform using the set of border voxels as roots.
 *
 * @date Mar 11, 2016
 * @author Thiago Vallin Spina
 *
 * @param borders The set of border voxels.
 * @param grad The original image gradient.
 * @param A The adjacency relation to be used by the EDT.
 * @param max_dist The maximum distance of voxels to the border set to be considered for computing the region of influence.
 */
double iftBorderRegionOfInfluenceGradMatchingScore(iftSet *borders, iftImage *grad, iftAdjRel *A, double max_dist);


/**
 * @brief Computes the Score Matrix from a set of Label Images using the Dice Similarity.
 * 
 * It considers that the input label images have <n_objs> objects with labels from [1, n_objs]. \n
 * 
 * The matrix element [i, j] is the score between the images i and j, which is the average dice
 * between all dice similarities from the objects. \n
 * 
 * @param  label_img_files Pathnames from the label images to be computed.
 * @param  n_objs          Number of Objects from the label images.
 * @return                 The Image Score Matrix.
 *
 * @author Samuel Martins
 * @date Nov 7, 2017
 * @ingroup ImageMetrics
 */
iftMatrix *iftComputeScoreMatrixByDice(const iftFileSet *label_img_files, int n_objs);


/**
 * @brief Computes the Score Matrix from a set of Label Images using the ASSD Distance.
 * 
 * It considers that the input label images have <n_objs> objects with labels from [1, n_objs]. \n
 * 
 * The matrix element [i, j] is the score between the images i and j, which is the average dice
 * between all ASSD Distance from the objects. \n
 * 
 * @param  label_img_files Pathnames from the label images to be computed.
 * @param  n_objs          Number of Objects from the label images.
 * @return                 The Image Score Matrix.
 *
 * @author Samuel Martins
 * @date Nov 7, 2017
 * @ingroup ImageMetrics
 */
iftMatrix *iftComputeScoreMatrixByASSD(const iftFileSet *label_img_files, int n_objs);


/**
 * @brief Computes the IMage Euclidean Distance (IMED) between two images as proposed in [1].
 * 
 * [1] Wang, Liwei, Yan Zhang, and Jufu Feng. "On the Euclidean distance of images."
 * IEEE transactions on pattern analysis and machine intelligence 27.8 (2005): 1
 * 
 * @warning This function is extremely slow: O(n^2)
 * 
 * @param img1 First Image.
 * @param img2 Second Image.
 * @return IMage Euclidean Distance.
 * 
 * @author Samuel Martins
 * @date Aug 29, 2018
 * @ingroup ImageMetrics
 */
double iftImageEuclideanDistance(const iftImage *img1, const iftImage *img2);


/**
 * @brief Computes the IMage Euclidean Distance (IMED) between two images for each object defined
 * by a label image.
 * 
 * It returns a double array with the IMED for each object, so that the array has size of [0, n_labels],
 * where n_labels is the label with highest value in the label image.
 * 
 * IMED was proposed in:
 * [1] Wang, Liwei, Yan Zhang, and Jufu Feng. "On the Euclidean distance of images."
 * IEEE transactions on pattern analysis and machine intelligence 27.8 (2005): 1
 * 
 * @param img1 First Image.
 * @param img2 Second Image.
 * @param label_img Label Image that defines the objects.
 * @return Double array with IMED for each object, so that [i] is equal to the imed for the object with
 *         label i. 
 * 
 * @author Samuel Martins
 * @date Aug 30, 2018
 * @ingroup ImageMetrics
 */
iftDblArray *iftImageEuclideanDistanceForLabels(const iftImage *img1, const iftImage *img2,
                                                const iftImage *label_img);



/**
 * @brief Computes thes Achievable Segmentation Accuracy (ASA) between a superpixel map and a ground truth.
 *
 * Achievable segmentation accuracy (ASA) is a performance upperbound measure.
 * It gives the highest accuracy achievable for object segmentation that utilizes superpixels as units.
 *
 * To compute ASA we label each superpixel with the label of the ground truth segment that has the largest overlap,
 * including the background (label 0).
 * The fraction of correctly labeled pixels is the achievable accuracy.
 *
 * @param super_map Superpixel Map. It must only have values >= 1.
 * @param gt Ground-Truth Image. It must only have values >= 0.
 * @return Computed Achievable Segmentation Accuracy.
 *
 * @author Samuel Martins
 * @date May 7, 2019
 * @ingroup ImageMetrics
 */
double iftASA(const iftImage* super_map, const iftImage* gt);


/* @brief (O - FP)(O - FN) / O^2
 * O: object size
 * FP: false positives
 * FN: false negatives
 *  @InProceedings{Cappabianco:2019:GBRAC,
 *  author = {Cappabianco, F. A. M. and Ribeiro, Pedro F. O. and Miranda, P. A. V. and Udupa, J. K.},
 *  title = {A GENERAL AND BALANCED REGION-BASED METRIC FOR EVALUATING MEDICAL IMAGE SEGMENTATION ALGORITHMS},
 *  booktitle = {IEEE International Conference on Image Processing},
 *  year = {2019},
 *  note = {To appear}
 *  }
 * @author Jordão Bragantini
 * @date July 31, 2019
 * @ingroup ImageMetrics
 */
//! swig(newobject)
double iftGeneralBalancedCoeff(const iftImage *bin_source, const iftImage *bin_target);

//! swig(newobject)
double          iftShannonEntropy(const iftImage *image);
double          iftJointEntropy(const iftImage *image1, const iftImage *image2);
double          iftNormalizedMutualInformation(const iftImage *image1, const iftImage *image2);
float           iftMeanDistBetweenBoundaries(iftImage *bin1, iftImage * bin2);

#ifdef __cplusplus
}
#endif

#endif


