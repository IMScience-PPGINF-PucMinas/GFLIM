/**
 * @file
 * @brief Bag of Features module
 * @author Peixinho
 * @date  April, 15
 */

#ifndef IFT_IFTBOW_H
#define IFT_IFTBOW_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftDataSet.h"
#include "iftKmeans.h"
#include "iftMImage.h"
#include "iftMetrics.h"

//#include "iftGenericVector.h"
//#include "iftGenericMatrix.h"
//#include "iftGenericVectorUtil.h"
//#include "iftGenericMatrixUtil.h"
//#include "iftClassification.h"

/** @addtogroup Descriptor
 * @{ */

/** @addtogroup BagOfFeatures
 * @brief Bag of Visual Words descriptor
 * @{ */


/**
 * @brief Features Extractor Prototype
 * Interface for patch feature extractor
 * @param img Image to extract patch features
 * @param begin Patch starting voxel
 * @param patchSize Size of patch to extract features
 * @param[out] featsOut Extracted features
 * @return Number of extracted features
 * @date Jun, 2016
 * @author Peixinho
 **/
typedef int (*iftPatchFeatsExtractor)(iftImage* img, iftVoxel begin, iftSize patchSize, float* featsOut);

/** @brief Patches Extractor Prototype
 * Interface for patch sampling
 * @param img Image to select samples
 * @param patchSize Size of extracted patches
 * @param nSamples Number of extracted patches
 * @param[out] patchesOut Starting voxel of selected patches in the image
 * @return Number of selected patches
 **/
typedef int (*iftPatchesSampler)(iftImage* img, iftSize patchSize, int nSamples, int *patchesOut);

/**
 * @brief Bag of Words (kernels) estimator Prototype.
 * Interface to learn dictionary kernels
 * @param Z Dataset containing the extracted features of selected patches
 * @param nkernels Dictionary size, number of kernels to be learned
 * @return Dataset containing the dictionary kernels
 **/
typedef iftDataSet*(*iftBowKernelsEstimator)(iftDataSet* Z, int nkernels);

/** @brief Coding and Pooling over iftMImage Prototype
 *  Interface to coding and pooling strategies on Bag of Features.
 *  Takes as input a multiband image, with the patches responses to each dictionary kernel (The higher the response, the closer to that kernel).
 *  @param mimg Multiband image of responses per kernel
 *  @param[out] featsOut Features after pooling and coding
 *  @param n Number of features
 *  @return The number of features (should be the same as the parameter n, only for asserting)
 **/
typedef int (*iftBowCodingPooling)(iftMImage* mImg, float* featsOut, int n);

/**
 * @brief Bag of Words Feature Extractor (Bow)
 * Contains the Bow architecture information along with visual words learned.
 * @author Peixinho
 * @date Jun, 2016
 */
typedef struct ift_bag_of_features
{
    /** True if the image has colors, False otherwise. */
    bool isColored;

    /** Number of patches sampled per image. */
    int nPatchfeats;
    /** Number of samples per image. */
    int samplesPerImage;
    /** Dictionary size. Number of visual words.*/
    int dictionarySize;

    /** Input image size.*/
    iftSize imgSize;
    /** Size for sampled image patches. */
    iftSize patchSize;
    /** Stride size along. */
    iftSize stride;

    /** Number of features extracted in each image sampled path. */
    int nfeats;

    /** Dictionary. Visual words. */
    iftDataSet *dictionary;

    /** Function to extract features from sampled patches.  */
    iftPatchFeatsExtractor featsExtractor;
    /** Function to sample patches in image. */
    iftPatchesSampler sampler;
    /** Function to select the best K visual words. */
    iftBowKernelsEstimator dictionaryEstimator;
    /** Function to apply coding and pooling in the final Bow multiband image representation. */
    iftBowCodingPooling codingPooling;

} iftBagOfFeatures;

/**
 * @brief Extract the Raw pixels brightness.
 * @implements iftPatchFeatsExtractor()
 */
int iftBowPatchRawFeatsExtractor(iftImage* img, iftVoxel begin, iftVoxel patchSize, float* featsOut);

/**
 * @brief Extract the Raw pixels YCbCr.
 * @implements iftPatchFeatsExtractor()
 */
int iftBowPatchColorRawFeatsExtractor(iftImage* img, iftVoxel begin, iftVoxel patchSize, float* featsOut);

/**
 * @brief Extract BIC features from image patch.
 * @implements iftPatchFeatsExtractor()
 */
int iftBowPatchBICFeatsExtractor(iftImage* img, iftVoxel begin, iftVoxel patchSize, float* featsOut);

/**
 * @brief Extract random samples from image
 * @implements iftPatchesSampler()
 */
int iftBowPatchesRandomSampler(iftImage* img, iftVoxel size, int nsamples, int *patchesOut);

/**
 * Extracts samples from image in a dense way
 * @implements iftPatchSampler()
 */
int iftBowPatchesDenseSampler(iftImage* img, iftVoxel size, int nsamples, int *patchesOut);


/**
 * @brief KMeans Bag of Words estimator. Implements iftBowKernelsEstimator()
 */
iftDataSet* iftBowKMeansKernels(iftDataSet* z, int nkernels);

/**
 * @brief Supervised KMeans Bag of Words estimator. Implements iftBowKernelsEstimator()
 *
 * Computes a KMeans Bag of Words estimator for samples in each class.
 */
iftDataSet* iftBowSupKMeansKernels(iftDataSet* z, int nkernels);

/**
 * @brief Soft coding followed of histogram pooling
 * @implements iftBowCodingPooling()
 */
int iftBowSoftCoding(iftMImage* mImg, float *featsOut, int n);

/**
 * @brief Hard coding followed of histogram pooling
 * @implements iftBowCodingPooling()
 */
int iftBowHardCoding(iftMImage* mImg, float *featsOut, int n);

/**
 * @brief Don't realize coding nor histogram pooling, just grab the raw iftMImages
 * @implements iftBowCoding()
 */
int iftNoCoding(iftMImage* mImg, float *featsOut, int n);

/**
 * @brief Horizontal and vertical sum of intensities. Preserves spatial information.
 * @implements iftBowCoding()
 */
int iftBowHVSumPooling(iftMImage *mImg, float *featsOut, int n);

/**
 * @brief Create a Bag of Features extractor.
 * @param color True if color image inputs, False otherwise
 * @param xImgSize Input image x size
 * @param yImgSize Input image y size
 * @param zImgSize Input image z size
 * @param xPatchSize Image patch x size
 * @param yPatchSize Image patch y size
 * @param zPatchSize Image patch z size
 * @param xStrideSize Stride in x
 * @param yStrideSize Stride in y
 * @param zStrideSize Stride in z
 * @param sampler Patches sampling strategy (implementing iftPatchesSampler())
 * @param featsExtractor Patch feature extractor (implementing iftPatchFeatsExtractor())
 * @param kernelsEstimator Dictionary Learning strategy (implementing iftBowKernelsEstimator())
 * @param coding Coding and pooling strategy (iftBowCodingPooling)
 * @param samplesPerImage Number of patches per training image
 * @param dictSize Number of dictionary kernels
 * @return The Bag of Features Feature Extractor
 */
iftBagOfFeatures *iftCreateBow(bool color, int xImgSize, int yImgSize, int zImgSize, int xPatchSize, int yPatchSize, int zPatchSize,
                     int xStrideSize, int yStrideSize, int zStrideSize, iftPatchesSampler sampler,
                     iftPatchFeatsExtractor featsExtractor, iftBowKernelsEstimator kernelsEstimator,
                     iftBowCodingPooling coding, int samplesPerImage, int dictSize);

/**
 * @brief Validate Bow params, to check if the architecture is feasible
 * @param bow The feature extraction architecture
 * @return True if valid, False otherwise
 */
bool iftBowSetup(iftBagOfFeatures *bow);

/**
 * @brief Destroy a Bow feature extractor.
 * @param bow Bag of Features extractor to be destroyed
 * @author Peixinho
 * @date Jun, 2016
 */
void iftDestroyBow(iftBagOfFeatures **bow);

/**
 * @brief Save Bow feature extractor. (Not implemented)
 */
void iftBowWrite(iftBagOfFeatures *bow);

/**
 * @brief Load Bow feature extractor. (Not implemented)
 */
iftBagOfFeatures *iftBowRead(const char *filename);

/**
 * @brief Learn Visual Words over image dataset.
 * @param bow Bag of Features Extractor
 * @param images Set of images to learn
 * @author Peixinho
 * @date Jun, 2016
 */
void iftBowLearn(iftBagOfFeatures *bow, iftFileSet *images);

/**
 * @brief Compute Bow multiband image from image.
 * Computes the multiband image as the response of all image patches for each dictionary kernel
 * @param bow Bag of Features Extractor
 * @param img Input image
 * @author Peixinho
 * @date Jun, 2016
 */
iftMImage* iftBowTransform(iftBagOfFeatures *bow, iftImage *img);

/**
 * @brief Extract Bow features from one single image
 * @param bow Bag of Features architecture
 * @param img Input image
 * @date Jun, 2016
 * @author Peixinho
 */
iftFeatures* iftBowExtractFeatures(iftBagOfFeatures * bow, iftImage* img);

/**
 * @brief Compute Bow features from a batch of images.
 * @param bow Bag of Features extractor
 * @param imgs Image files list
 * @return The dataset of extracted features
 * @author Peixinho
 * @date Jun, 2016
 */
iftDataSet *iftBowExtractFeaturesBatch(iftBagOfFeatures *bow, iftFileSet *imgs);

/** @} */

/** @} */


//Generic Bow
typedef struct _bagOfVisualWordsManager BagOfVisualWordsManager;



//typedef iftGenericVector* (*iftImageSamplerFunction)(iftImage* image, BagOfVisualWordsManager* bagOfVisualWordsManager);
//typedef iftGenericMatrix* (*iftFeatureExtractorFunction)(iftGenericVector* outputSampler, BagOfVisualWordsManager* bagOfVisualWordsManager);
//typedef iftGenericMatrix* (*iftClusteringFunction)(iftGenericMatrix* outputFeatureExtractor_allSamples, BagOfVisualWordsManager* bagOfVisualWordsManager);
//typedef iftGenericVector* (*iftMountHistogramFunction) (iftGenericMatrix* outputFeatureExtractor_singleSample,BagOfVisualWordsManager* bagOfVisualWordsManager);
//
//typedef struct _bagOfVisualWordsManager {
//    //paths
//    iftGenericVector* pathsToImages_dictionary;
//    iftGenericVector* pathsToImages_train;
//    iftGenericVector* pathsToImages_test;
//
//    //dictionery stuffs
//    iftGenericMatrix* dataVisualWords;
//    iftGenericVector* labelsVisualWords;
//    iftGenericMatrix* dictionary;
//
//    //machine learning stuffs
//    iftGenericMatrix* histogramsTraining;
//    iftGenericVector* labelsTraining;
//    iftGenericMatrix* histogramsPredictSamples;
//    iftGenericVector* labelsPredicted;
//
//    bool storeVisualWordsData;
//    bool storeTrainData;
//    bool storePredictedData;
//
//    void* classifier;
//
//    FreeFunction freeFunction2SamplerOutput;
//    FreeFunction freeFunctionClassifier;
//
//    iftArgumentList* argumentListOfSampler;
//    iftArgumentList* argumentListOfFeatureExtractor;
//    iftArgumentList* argumentListOfClustering;
//    iftArgumentList* argumentListOfDistanceFunction;
//    iftArgumentList* argumentListOfHistogramMounter;
//
//    iftFeatureExtractorFunction featureExtractorFunction;
//    iftImageSamplerFunction imageSamplerFunction;
//    iftDistanceFunction distanceFunction;
//    iftClusteringFunction clusteringFunction;
//    iftMountHistogramFunction mountHistogramFunction;
//    iftFitFunction fitFunction;
//    iftPredictFunction predictFunction;
//
//    iftImage* currentImage;
//}BagOfVisualWordsManager;

#ifdef __cplusplus
}
#endif

#endif //IFT_IFTBOW_H
