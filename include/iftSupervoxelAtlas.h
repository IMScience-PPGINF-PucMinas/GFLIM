//
// Created by azaelmsousa on 12/11/20.
//

#ifndef IFT_IFTSUPERVOXELATLAS_H
#define IFT_IFTSUPERVOXELATLAS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftMatrix.h"
#include "iftMImage.h"
#include "iftSegmentation.h"
#include "ift/core/io/FileSet.h"
#include "ift/core/tools/String.h"
#include "ift/core/dtypes/Dict.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/imgproc/basic/Histogram.h"
  
  
typedef struct ift_supervoxel_atlas {
    iftImage *svoxel; //Supervoxel image
    iftMImage *mimg; //Image containing the features used for the gaussian distribution
    iftMatrix **covariance; //Covariance matrix for each supervoxel
    iftMatrix **mean; //Mean vector for each supervoxel
    iftMatrix **stdev; //Standard deviation vector for each supervoxel
    int n_svoxels; //Number of supervoxels
} iftSupervoxelAtlas;



iftMatrix *iftStandardizeTrainMatrix(iftMatrix *M, char axis, iftMatrix **mean, iftMatrix **stdev);

/**
 * @brief Aligns the input mimg features of each supervoxel with the origin. For instance, suppose we have these
 * two histograms:
 *
 *          |                                        |
 *          |          *                             |
 *          |         * *                            |                        *
 *          |        *   *                           |                       * *
 *          |       *     *                          |                      *   *
 *          |     *        *                         |                     *      *
 *          |___*___________*_______________         |____________________*_________*__
 *
 * The output will be:
 *
 *          |                                        |
 *          |*                                       |
 *          * *                                     *|
 *         *|  *                                   * *
 *        * |   *                                 *  |*
 *      *   |    *                               *   |  *
 *    *     |_____*_________________________    *    |____*_____________________________
 *
 * This way, the histograms can be analyzed at the same range.
 *
 *
 * @param svoxel Supervoxel model.
 * @param mimg Multiband image to be aligned.
 * @param mimg_label Label of the Multiband image.
 *
 * @return Aligned multiband image
 *
 * @author Azael M Sousa
 * @date Nov 13, 2020
 */
iftMImage *iftAlignSupervoxelsHistogram(iftImage *svoxel, iftMImage *mimg, iftImage *mimg_label);

/**
 * @brief Allocates memory and creates the Supervoxel Atlas data structure
 *
 * @param mimg Multiband image used to compute the supervoxel model.
 * @param svoxel Supervoxel model.
 *
 * @return Supervoxel Atlas data structure
 *
 * @author Azael M Sousa
 * @date Nov 13, 2020
 */

iftSupervoxelAtlas *iftCreateSupervoxelAtlas(iftMImage *mimg, iftImage *svoxel);

/**
 * @brief Computes the gray histogram of a single band from an specific supervoxel
 *
 * @param atlas Supervoxel atlas structure
 * @param band Band to which compute the histogram
 * @param svoxel_id Supervoxel to compute the histogram
 * @return
 */
iftHist *iftComputeSupervoxelHistogram(iftSupervoxelAtlas *atlas, int band, int svoxel_id);

/**
 * @brief Frees a Supervoxel Atlas instance
 *
 * @param atlas A supervoxel atlas pointer.
 *
 * @author Azael M Sousa
 * @date Nov 13, 2020
 */
void iftDestroySupervoxelAtlas(iftSupervoxelAtlas **atlas);

/**
 * @brief Writes the supervoxel atlas into a zip file
 *
 * @param atlas Supervoxel atlas to be written
 * @param filename Name of the output file with extension .zip
 *
 * @author Azael M Sousa
 * @date Nov 13, 2020
 */
void iftWriteSupervoxelAtlas(const iftSupervoxelAtlas *atlas, const char *filename, ...);

/**
 * @brief Reads a supervoxel atlas model
 *
 * @param filename Path to the saved model.
 *
 * @return Supervoxel Atlas data structure
 *
 * @author Azael M Sousa
 * @date Nov 13, 2020
 */
iftSupervoxelAtlas *iftReadSupervoxelAtlas(const char *filename);

/**
 * @brief Computes the multivariate gaussian pdf. The multivariate gaussian distribution is comuted as follows:
 *
 *                                            1
 *            p(x | phi, sigma) = ------------------------- * e^( -1/2 * (x-phi)^{T}*sigma^(-1)*(x-phi) ),
 *                                (2pi)^(n/2)*sqrt(|sigma|)
 *
 * where phi is the nD mean, sigma is the covariance matrix, n is the number of features and (.)^{T} refers to the
 * transpose of a matrix.
 *
 * @param sample Matrix contianing the voxels of the supervoxels to be segmented
 * @param covariance Covariance matrix of the corresponding supervoxel.
 * @param mean Mean matrix of the corresponding supervoxel.
 *
 * @return Probability for the given sample.
 *
 * @author Azael M Sousa
 * @date Nov 13, 2020
 */
double iftMultivariateGaussPDF(iftMatrix *sample, iftMatrix *covariance, iftMatrix *mean);

/**
 * @brief Segments a test MImage according to a previous trained Multivariate Gaussian Distribution.
 *
 * @param test Multiband image to be segmented.
 * @param atlas Trained Supervoxel Atlas
 * @param bottom_threshold Bottom threshold to define if the supervoxel must be segmented
 * @param top_threshold Top threshold to define if the supervoxel must be segmented
 *
 * @return Image with 0s for the supervoxels that are outside the threshold factor and 1s for the
 * supervoxels inside the threshold.
 *
 * @author Azael M Sousa
 * @date Nov 13, 2020
 */
iftImage *iftSegmentationBySupervoxelAtlas(iftMImage *test, iftSupervoxelAtlas *atlas, double bottom_threshold, double top_threshold);

/**
 * @brief Computes the probabilistic map of a test MImage according to a previous trained Multivariate Gaussian Distribution.
 *
 * @param test Multiband image to be segmented.
 * @param altas Trained Supervoxel Atlas
 *
 * @return Image with the probability of each voxel to belong in the gaussian distribution.
 *
 * @author Azael M Sousa
 * @date Nov 13, 2020
 */
iftImage *iftProbMapBySupervoxelAtlas(iftMImage *test, iftImage *test_label, iftSupervoxelAtlas *atlas);

#ifdef __cplusplus
}
#endif
  
#endif //IFT_IFTDISTRIBUTION_H

