/**
 * @file iftBagOfVisualWords.h
 * @brief Bag of Visual Words module
 * @author Cesar Castelo
 * @date  January 4, 2018
 */

#ifndef IFT_IFTBAGOFVISUALWORDS_H
#define IFT_IFTBAGOFVISUALWORDS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "stdarg.h"
#include "ift/imgproc/dtypes/Roi.h"
#include "iftImage.h"
#include "iftDataSet.h"
#include "iftMemory.cuh"
#include "iftClustering.h"
#include "iftIGraph.h"

typedef enum {
    BOVW_INT_POINT_DETECTOR_RANDOM, BOVW_INT_POINT_DETECTOR_GRID,
    BOVW_INT_POINT_DETECTOR_UNSUP_SPIX_ISF, BOVW_INT_POINT_DETECTOR_SUP_SPIX_ISF
} iftBovwIntPointDetectorId;

typedef enum {
    BOVW_LOCAL_FEAT_EXTRACTOR_RAW, BOVW_LOCAL_FEAT_EXTRACTOR_BIC,
    BOVW_LOCAL_FEAT_EXTRACTOR_LBP, BOVW_LOCAL_FEAT_EXTRACTOR_BRIEF,
    BOVW_LOCAL_FEAT_EXTRACTOR_CONV, BOVW_LOCAL_FEAT_EXTRACTOR_CONV_MULTI_LAYER,
    BOVW_LOCAL_FEAT_EXTRACTOR_DEEP_FEATS_MIMG
} iftBovwLocalFeatExtractorId;

typedef enum {
    BOVW_DICT_ESTIMATOR_UNSUP_KMEANS, BOVW_DICT_ESTIMATOR_SUP_KMEANS_PIX_LABEL,
    BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMG_CLASS, BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMAGE,
    BOVW_DICT_ESTIMATOR_SUP_KMEANS_POSITION, BOVW_DICT_ESTIMATOR_SUP_KMEANS_IMG_CLASS_AND_POSITION,
    BOVW_DICT_ESTIMATOR_UNSUP_OPF, BOVW_DICT_ESTIMATOR_SUP_OPF_PIX_LABEL,
    BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS, BOVW_DICT_ESTIMATOR_SUP_OPF_IMAGE,
    BOVW_DICT_ESTIMATOR_SUP_OPF_POSITION, BOVW_DICT_ESTIMATOR_SUP_OPF_IMG_CLASS_AND_POSITION,
    BOVW_DICT_ESTIMATOR_SUP_MANUAL_IMG_CLASS, BOVW_DICT_ESTIMATOR_SUP_MANUAL_POSITION,
    BOVW_DICT_ESTIMATOR_SUP_MANUAL_IMG_CLASS_AND_POSITION
} iftBovwDictEstimatorId;

typedef enum {
    BOVW_COD_FUNC_HARD_ASGMT, BOVW_COD_FUNC_SOFT_ASGMT, BOVW_COD_FUNC_HARD_ASGMT_BATCH,
    BOVW_COD_FUNC_SOFT_ASGMT_BATCH, BOVW_COD_FUNC_2ND_ORDER_STATS_FISHER_VECTORS
} iftBovwCodFuncId;

typedef enum {
    BOVW_DIST_FUNC_1, BOVW_DIST_FUNC_2, BOVW_DIST_FUNC_3, BOVW_DIST_FUNC_4, BOVW_DIST_FUNC_5,
    BOVW_DIST_FUNC_6, BOVW_DIST_FUNC_7, BOVW_DIST_FUNC_8, BOVW_DIST_FUNC_9, BOVW_DIST_FUNC_10,
    BOVW_DIST_FUNC_11, BOVW_DIST_FUNC_12
} iftBovwDistFuncId;

typedef enum {
    BOVW_LOCAL_FEATS_ORDERING_PIX_LABEL, BOVW_LOCAL_FEATS_ORDERING_IMG_CLASS,
    BOVW_LOCAL_FEATS_ORDERING_IMAGE, BOVW_LOCAL_FEATS_ORDERING_POSITION
} iftBovwLocalFeatsOrderingId;

typedef enum {
    BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_INVER, BOVW_COD_FUNC_SOFT_ASGMT_WGT_FUNC_GAUSS
} iftBovwCodFuncSoftAsgmtWgtFuncId;

/* ---------------------------------------------------------------------------------------------------------------------*/
/*                                                Structure Definition                                                  */
/* ---------------------------------------------------------------------------------------------------------------------*/

/** @addtogroup Descriptor
 * @{ */

/** @addtogroup BagOfVisualWords
 * @brief Bag of Visual Words descriptor
 * @{ */

/** @brief Interest Points Detector Prototype
 * Interface for interest points detector
 * @param img Image to select samples
 * @param imgMask pixel-level labels of img (mask)
 * @param imgFilePtr Pointer to the image file
 * @param saveImgIntPts Whether or not to save the image with the interest points
 * @param funcParams Extra parameters to be received (for the instanced function)
 * @return iftRoiArray that represents the regions of interest
 * @author Cesar Castelo
 * @date Jan 04, 2018
 **/
typedef iftRoiArray* (*iftBovwIntPointDetector)(iftImage* img, iftImage* imgMask, iftFile *imgFilePtr, bool saveImgIntPts, iftDict *funcParams);

/**
 * @brief Local Features Extractor Prototype
 * Interface for local feature extractor
 * @param img Image to extract patch features
 * @param roiArray Regions of interest
 * @param imgFilePtr Pointer to the image file
 * @param patchDirPath Directory to save the image patches that will be created (reference data). If NULL the reference data will not be created
 * @param funcParams Extra parameters to be received (for the instanced function)
 * @return iftDataSet that contains the local features
 * @author Cesar Castelo
 * @date Jan 04, 2018
 **/
typedef iftDataSet* (*iftBovwLocalFeatExtractor)(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams);

/**
 * @brief Dictionary Estimator Prototype.
 * Interface to dictionary estimator
 * @param localFeats iftDataset containing the set of local features extracted from a set of images
 * @param funcParams Extra parameters to be received (for the instanced function)
 * @param isFisherVectorsEncoding Whether or not the method will be used for Fisher Vectors encoding 
 * @return void* containing the clustering model that represents the dictionary (opf-clustering, k-means, etc)
 * @author Cesar Castelo
 * @date Jan 04, 2018, updated Oct 9, 2019
 **/
typedef void* (*iftBovwDictEstimator)(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding);

/**
 * @brief Coding Function Prototype
 * Interface to coding functions for a new image
 * @param clustModel void* containing the clustering model that represents the dictionary (opf-clustering, k-means, etc)
 * @param localFeats Set of local features extracted from the image to be coded
 * @param funcParams Extra parameters to be received (for the instanced function)
 * @return iftFeatures that represents the activations for each regionOfInterest against every visual word on the dictionary
 * @author Cesar Castelo
 * @date Jan 04, 2018, updated Oct 9, 2019
 **/
typedef iftFeatures* (*iftBovwCodFunc)(void* dict, iftDataSet* localFeats, iftDict *funcParams);

/**
 * @brief Bag of Visual Words Structure
 * Contains the BoVW architecture information along with the visual dictionary learned.
 * @author Cesar Castelo
 * @date Jan 04, 2018, updated Oct 9, 2019
 */
typedef struct ift_bag_of_visual_words
{
    /* set of local features extracted from all the images */
    iftDataSet *localFeats;

    /* generic pointer that represents the dictionary. If the BoVW uses fisher vectors to encode
    the local feats, the pointer will be a clustering model, otherwise, it will be a dataset
    containing the centroids of the clusters */
    void *dict;

    /* order-hierarchy-related properties */
    iftIntMatrix *localFeatsOrderHrchy; // order hierarchy for the local feats from all images
    iftIntArray *orderHrchies; // IDs of the order hierarchies

    /* images used to create the dictionary */
    iftFileSet *dictLearnFileset;
    iftFileSet *dictLearnMasksFileset;

    /* Function to detect interest points in the image */
    iftBovwIntPointDetector intPointDetector;
    iftBovwIntPointDetectorId intPointDetectorId;

    /* Function to extract local features from the interest point detected  */
    iftBovwLocalFeatExtractor localFeatExtractor;
    iftBovwLocalFeatExtractorId localFeatExtractorId;

    /* Function to estimate the visual dictionary */
    iftBovwDictEstimator dictEstimator;
    iftBovwDictEstimatorId dictEstimatorId;

    /* Function to apply the coding operation */
    iftBovwCodFunc codFunc;
    iftBovwCodFuncId codFuncId;

    /* distance function used to compute the dictionary and to perform the coding operation */
    iftArcWeightFun arcWeightFunc;
    iftBovwDistFuncId distFuncId;

    /* function pointers' params */
    iftDict *funcParams;

} iftBagOfVisualWords;

/**
 * @brief OPF Clustering model for FV codification
 * Contains the OPF clustering model (knn-graph) along with the model params (means, covariances and weights per group)
 * @author Cesar Castelo
 * @date Oct 1, 2019
 */
typedef struct ift_fv_opf_clust_model
{
    /* OPF knn-graph model */
    iftKnnGraph *graph;

    /* set of means for each cluster */
    iftMatrix *means;
    
    /* set of covariance matrices for each cluster */
    iftMatrix **covariances;

    /* set of covariance inverse matrices for each cluster */
    iftMatrix **covariancesInv;

    /* set of covariance matrices determinants for each cluster */
    iftFloatArray *covariancesDet;

    /* weights for the clusters */
    iftFloatArray *weights;

    /* number of features per word */
    int nFeatsPerWord;

} iftBovwOPFClustModelFV;

/* ---------------------------------------------------------------------------------------------------------------------*/
/*                                                Initialization functions                                              */
/* ---------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Creates a BoVW structure.
 * @return The Bag of Visual Words structure
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
iftBagOfVisualWords *iftCreateBovw();

/**
 * @brief Sets the function pointers according to the chosen parameters
 * @param bovw BoVW structure
 * @param intPointDetector Function to detect interest points in the image
 * @param localFeatExtractor Function to extract local features from the interest point detected
 * @param dictEstimator Function to estimate the visual dictionary
 * @param codFunc Function to apply the coding operation
 * @param distFunc Distance function to be used
 * @param funcParams Parameters for the function pointers
 * @author Cesar Castelo
 * @date Jan 30, 2018
 */
void iftBovwSetFunctionPointers(iftBagOfVisualWords *bovw, iftBovwIntPointDetectorId intPointDetector,
                                iftBovwLocalFeatExtractorId localFeatExtractor, iftBovwDictEstimatorId dictEstimator,
                                iftBovwCodFuncId codFunc, iftBovwDistFuncId distFunc, iftDict *funcParams);

/**
 * @brief Destroys a BoVW structure.
 * @param bovw BoVW structure to be destroyed
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
void iftDestroyBovw(iftBagOfVisualWords **bovw);

/**
 * @brief Saves a BoVW structure.
 * @param bovw BoVW structure to be saved
 * @param filename Filename to save the bovw structure
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
void iftBovwWrite(iftBagOfVisualWords *bovw, char *filename);

/**
 * @brief Reads a BoVW structure.
 * @param fileName File that contains the BoVW structure
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
iftBagOfVisualWords *iftBovwRead(const char *filename);

/**
 * @brief Creates a OPF clustering model structure for FV codification.
 * @param graph knn-graph structure containing the OPF graph
 * @param means matrix containing the means per group
 * @param covariances set of matrices containing the covariance matrices per group
 * @param covariancesInv set of matrices containing the inverse covariance matrices per group
 * @param covariancesDet set of matrices containing the covariance matrices determinants per group
 * @param weights float array containing the weights per group
 * @return The OPF clustering model structure
 * @author Cesar Castelo
 * @date Oct 9, 2019
 */
iftBovwOPFClustModelFV *iftCreateBovwOPFClustModelFV(iftKnnGraph *graph, iftMatrix *means, iftMatrix **covariances, iftMatrix **covariancesInv, iftFloatArray *covariancesDet, iftFloatArray *weights);

/**
 * @brief Destroys a OPF clustering model structure for FV codification.
 * @param model The OPF clustering model structure to be destroyed
 * @author Cesar Castelo
 * @date Oct 9, 2019
 */
void iftDestroyBovwOPFClustModelFV(iftBovwOPFClustModelFV **model);

/**
 * @brief Saves a OPF clustering model structure for FV codification.
 * @param model OPF clustering model to be saved
 * @param filename Filename to save the model
 * @author Cesar Castelo
 * @date Oct 9, 2019
 */
void iftWriteBovwOPFClustModelFV(iftBovwOPFClustModelFV *model, char *filename);

/**
 * @brief Reads a OPF clustering model structure for FV codification.
 * @param filename Filename to read the model from
 * @return OPF clustering model
 * @author Cesar Castelo
 * @date Oct 9, 2019
 */
iftBovwOPFClustModelFV *iftReadBovwOPFClustModelFV(char *filename);

/* ---------------------------------------------------------------------------------------------------------------------*/
/*                                                Learning functions                                               */
/* ---------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Learns a Bag of Visual Words structure for image classification
 * @param bovw Bag of Visual Words structure to be trained
 * @param dictLearnFileset image set that will be used to learn the dictionary
 * @param dictLearnMasksFileset image set that contains the pixel labels (mask) of dictLearnFileset
 * @param patchDirPath Directory to save the image patches that will be created (reference data)
 * @param saveImgIntPts Whether or not to save the image with the interest points
 * @return The time to learn the visual dictionary
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
float iftBovwLearnImgClassif(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset, iftFileSet *dictLearnMasksFileset, char *patchDirPath, bool saveImgIntPts);

/**
 * @brief Extract the interest points from a set of images given a Bag of Visual Words structure
 * @param bovw Bag of Visual Words structure
 * @param dictLearnFileset image set that will be used to learn the dictionary
 * @param dictLearnMasksFileset image set that contains the pixel labels (mask) of dictLearnFileset
 * @param patchDirPath Directory to save the image patches that will be created (reference data)
 * @param saveImgIntPts Whether or not to save the image with the interest points
 * @return The time to extract the interest points
 * @author Cesar Castelo
 * @date Dec 17, 2019
 */
float iftBovwLearnIntPointDetec(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset, iftFileSet *dictLearnMasksFileset, char *patchDirPath, bool saveImgIntPts);

/**
 * @brief Computes the feature vector for an image using a previously trained dictionary (it uses the coding and pooling functions)
 * @param bovw Bag of Visual Words structure
 * @param img Image to extract the feature vector from
 * @param imgMask pixel-level labels of img (mask)
 * @param imgFilePtr Pointer to the image file
 * @param saveImgIntPts Whether or not to save the image with the interest points
 * @return The computed feature vector
 * @author Cesar Castelo
 * @date Feb 08, 2018
 */
iftFeatures *iftBovwComputeFeatVectFromImage(iftBagOfVisualWords *bovw, iftImage *img, iftImage*imgMask, iftFile *imgFilePtr, bool saveImgIntPts);

/**
 * @brief Computes the feature vector for an image fileset using a previously trained dictionary
 * It calls the iftBovwComputeFeatVectFromImage function for every image
 * @param bovw Bag of Visual Words structure
 * @param imgFileSet Image fileset to extract the feature vectors from
 * @param imgMasksFileSet fileset containing the pixel-levels of imgFileSet
 * @param saveImgIntPts Whether or not to save the image with the interest points
 * @return A dataset containing the feature vectors
 * @author Cesar Castelo
 * @date Feb 08, 2018
 */
iftDataSet *iftBovwComputeFeatVectsFromImgFileSet(iftBagOfVisualWords *bovw, iftFileSet *imgFileSet, iftFileSet *imgMasksFileSet, bool saveImgIntPts);

/**
 * @brief Computes the feature vector for an image fileset using a previously trained dictionary.
 * It computes the feature vectors dividing the image set in batches to speed up the process
 * @param bovw Bag of Visual Words structure
 * @param imgFileSet Image fileset to extract the feature vectors from
 * @param imgMasksFileSet fileset containing the pixel-levels of imgFileSet
 * @param batchSize Number of images per batch
 * @param saveImgIntPts Whether or not to save the image with the interest points
 * @return A dataset containing the feature vectors
 * @author Cesar Castelo
 * @date Feb 06, 2019
 */
iftDataSet *iftBovwComputeFeatVectsFromImgFileSetBatch(iftBagOfVisualWords *bovw, iftFileSet *imgFileSet, iftFileSet *imgMasksFileset, int batchSize, bool saveImgIntPts);

// iftMImage* iftBowTransform(iftBagOfFeatures *bow, iftImage *img); // Compute Bow multiband image from image.

/* ---------------------------------------------------------------------------------------------------------------------*/
/*                                            Function pointers Implementation                                          */
/* ---------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Detects random interest points from image
 * @implements intPointDetector()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - nPoints: Number of points (int)
 * - patch_size: ize of each patch (int)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
iftRoiArray* iftBovwRandomIntPointDetector(iftImage* img, iftImage* imgMask, iftFile *imgFilePtr, bool saveImgIntPts, iftDict *funcParams);

/**
 * @brief Detects interest points from an image using grid sampling
 * @implements intPointDetector()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - patch_size: Size of each patch (int)
 * - stride: Stride of each patch (int)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
iftRoiArray* iftBovwGridIntPointDetector(iftImage* img, iftImage* imgMask, iftFile *imgFilePtr, bool saveImgIntPts, iftDict *funcParams);

/**
 * @brief Divides the images in superpixels using the ISF algorithm and without using class information. It then takes points from
 * the supervoxels' boundaries using grid sampling
 * @implements intPointDetector()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - n_superpixels: Number of superpixels (int)
 * - patch_size: Size of each patch (int)
 *
 * @author Cesar Castelo
 * @date Mar 15, 2019
 */
iftRoiArray* iftBovwUnsupSuperpixelISFIntPointDetector(iftImage* img, iftImage* imgMask, iftFile *imgFilePtr, bool saveImgIntPts, iftDict *funcParams);

/**
 * @brief Divides the image set by class and executes the ISF algorithm using all the images from each class, obtaining 
 * different sets of superpixels for each class. Then, using the superpixels' boundaries, a (different) number of points
 * are obtained using grid sampling for each class.
 * @implements intPointDetector()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - patch_size: Patch size for grid sampling (int)
 * - stride: Stride for grid sampling, to both create the initial seeds for ISF and then obtain the points from the spix mask (int)
 * - adj_rel_rad: Adjacency relation radius for the ISF algorithm (float)
 * - alpha: Alpha value for the ISF algorithm (float)
 * - beta: Beta value for the ISF algorithm (float)
 * - n_iter: Number of iterations for the ISF algorithm (int)
 * - col_space_in: Input color space of the images (str)
 * - col_space_out: Output color space, to be used in the ISF algorithm (str)
 *
 * @author Cesar Castelo
 * @date Mar 15, 2019
 */
iftRoiArray* iftBovwSupSuperpixelISFIntPointDetector(iftImage* img, iftImage* imgMask, iftFile *imgFilePtr, bool saveImgIntPts, iftDict *funcParams);

/**
 * @brief Extracts the Raw pixels color/brightness information from a set of regions of interest
 * @implements iftBovwLocalFeatExtractor()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - col_space_in: Input color space (str)
 * - col_space_out: Output color space (str)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
iftDataSet* iftBovwRawFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams);

/**
 * @brief Extract BIC features from a set of regions of interest
 * @implements iftBovwLocalFeatExtractor()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - n_bins: Number of bins for image quantization (int)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
iftDataSet* iftBovwBICFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams);

/**
 * @brief Extract LBP features from a set of regions of interest
 * @implements iftBovwLocalFeatExtractor()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - circ_sym_neighb_rad: Radius for the circular symmetric neighborhood (double)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018
 */
iftDataSet* iftBovwLBPFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams);

/**
 * @brief Extracts BRIEF features from a set of regions of interest
 * @implements iftBovwLocalFeatExtractor()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - patch_size: Patch size for the descriptor (int)
 * - descriptor_size: Number of features for each local vector (int)
 *
 * @author Cesar Castelo
 * @date Jan 22, 2019
 */
iftDataSet* iftBovwBriefFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams);

/**
 * @brief Extracts features from a set of regions of interest by using convolutions with a kernel bank
 * @implements iftBovwLocalFeatExtractor()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - col_space: Color space to be used (str)
 * - n_kernels: Number of kernels (int)
 * - kernel_size: Size of each kernel, e.g., 5 for 5x5 kernel (int)
 * - conv_stride: Stride for the convolution operation, e.g., 4 for 2x2 convolution stride (int)
 * - pool_size: Size for the max-pooling operation, e.g., 5 for 5x5 max-pooling (int)
 * - pool_stride: Stride for the max-pooling operation, e.g., 4 for 2x2 max-pooling stride (int)
 * - centralized_imgs_dir: Folder containing the batch-centralized images (str)
 * - kernels_dir: Folder containing the kernel matrix (str)
 * - max_pool_type: Max-pooling type: 'max_pool','max_pool_roi'
 * - resize_image: Whether or not to resize the resulting filtered image to match the original image size (str).
 *   If the image is not resized, the interest points previously detected are no longer used (all the positions from the reduced image are used instead).
 * - vect_constr: Feature vector construction method (str): 'pixel_vals_in_patch','mean_val_in_patch','bic_in_patch'.
 * - n_bins_per_band: Number of bins to extract BIC features. This is used only if vect_constr='bic_in_patch' (int).
 *
 * @author Cesar Castelo
 * @date Mar 11, 2019
 */
iftDataSet* iftBovwConvolutionalFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams);

/**
 * @brief Extracts features from a set of regions of interest by using multi-layer convolutions with a set of kernel banks
 * @implements iftBovwLocalFeatExtractor()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - col_space: Color space to be used (str)
 * - network_params: Dictionary containg the params for each layer (iftDict). The structure must be:
 *     {'n_layers',
 *      'layer_1': {
 *          'conv': {'n_bands', 'n_kernels', 'kernel_size', 'stride'},
 *          'relu': 'yes/no',
 *          'max-pool': {'size','stride'}
 *      },
 *      ...,
 *      'layer_n': {
 *          ...
 *      }}
 * - centralized_imgs_dir: Folder containing the batch-centralized images (str)
 * - kernels_dir: Folder containing the kernel matrices (str)
 * - max_pool_type: Max-pooling type: 'max_pool','max_pool_roi'
 * - resize_image: Whether or not to resize the resulting filtered image to match the original image size (str).
 *   If the image is not resized, the interest points previously detected are no longer used (all the positions from the reduced image are used instead).
 * - vect_constr: Feature vector construction method (str): 'pixel_vals_in_patch','mean_val_in_patch','bic_in_patch'.
 * - n_bins_per_band: Number of bins to extract BIC features. This is used only if vect_constr='bic_in_patch' (int).
 *
 * @author Cesar Castelo
 * @date Apr 22, 2019
 */
iftDataSet* iftBovwConvolutionalMultiLayerFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams);

/**
 * @brief Reads features previously extracted and saved as mimages.
 * @implements iftBovwLocalFeatExtractor()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - feature_mimgs_dir: CSV/Directory containing the feature mimages previously computed.
 * - resize_image: Whether or not to resize the feature image to match the original image size (str).
 *   If the image is not resized, the interest points previously detected are no longer used (all the positions from the reduced image are used instead).
 * - vect_constr: Feature vector construction method (str): 'pixel_vals_in_patch','mean_val_in_patch','bic_feats_from_patch'.
 * - n_bins_per_band: Number of bins to extract BIC features. This is used only if vect_constr='bic_feats_from_patch' (int).
 *
 * @author Cesar Castelo
 * @date Jun 18, 2019
 */
iftDataSet* iftBovwDeepFeaturesMImageFeatExtractor(iftImage* img, iftRoiArray* roiArray, iftFile *imgFilePtr, char *patchDirPath, iftDict *funcParams);

/**
 * @brief Unsupervised KMeans dictionary estimator for BoVW
 * @implements iftBovwDictEstimator()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - ret_real_cent: Whether to use the real centroids (present in the data) or the mean vector of each group (str)
 * - n_groups: Number of groups (int)
 * - max_iter: Maximum number of iterations (int)
 * - min_improv: Minimum improvement before stopping the k-means algorithm (double)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018, updated Oct 9, 2019
 */
void* iftBovwUnsupKMeansDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding);

/**
 * @brief Supervised one-level order-hierarchy KMeans dictionary estimator for BoVW (e.g. pixel, img-class, etc)
 * It computes a KMeans dictionary estimator for each ordering at a certain one-level order-hierarchy, e.g., object and
 * background (for pixel), image category (for img-class), etc.
 * @implements iftBovwDictEstimator()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - ret_real_cent: Used to be passed to iftBovwUnsupKMeansDictEstimator
 * - max_iter: Used to be passed to iftBovwUnsupKMeansDictEstimator
 * - min_improv: Used to be passed to iftBovwUnsupKMeansDictEstimator
 * - n_groups: Number of total groups (int)
 * - local_feats_order_hrchy: Order hierarchy for all the local features (iftIntMatrix*)
 *
 * @author Cesar Castelo
 * @date Feb 02, 2018, updated Oct 9, 2019
 */
void* iftBovwSupOneLevelOrderHrchyKMeansDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding);

void* iftBovwSupTwoLevelOrderHrchyKMeansDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding);

/**
 * @brief Unsupervised OPF dictionary estimator for BoVW
 * @implements iftBovwDictEstimator()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - ret_real_cent: Whether to use the real centroids (present in the data) or the mean vector of each group (str)
 * - clust_samp_perc: Percentage of the dataset to be used to compute the clusters (double)
 * - knn_graph_n_neighb: Number of samples to be used to mount the knn graph.
 *   If it is between [0,1], it is used as a percentage of the training set (after clustSampPerc),otherwise it is
 *   used as the number of neighbors (however, if the value is greater than 30% of the samples after clustSampPerc,
 *   30% is used instead)
 * - known_n_groups: Whether or not we know the number of groups (str)
 * - n_groups: Number of groups (int)
 * - use_prot_as_cent: Whether to use the prototypes of the OPF graph as the centroids (str)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018, updated Oct 9, 2019
 */
void* iftBovwUnsupOPFDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding);

/**
 * @brief Supervised one-level order-hierarchy OPF dictionary estimator for BoVW (e.g. pixel, img-class, position, etc)
 * It computes an OPF dictionary estimator for each ordering at a certain one-level order-hierarchy, e.g., object and
 * background (for pixel), image category (for img-class), etc.
 * @implements iftBovwDictEstimator()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - ret_real_cent: Used to be passed to iftBovwUnsupOPFDictEstimator
 * - clust_samp_perc: Used to be passed to iftBovwUnsupOPFDictEstimator
 * - knn_graph_n_neighb: Used to be passed to iftBovwUnsupOPFDictEstimator
 * - known_n_groups: Whether or not we know the number of groups (str)
 * - n_groups: Number of groups (int)
 * - local_feats_order_hrchy: Order hierarchy for all the local features (iftIntMatrix*)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018, updated Oct 9, 2019
 */
void* iftBovwSupOneLevelOrderHrchyOPFDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding);

/**
 * @brief Supervised two-level order-hierarchy OPF dictionary estimator for BoVW (e.g. img-class and position,
 * img-class and pixel, etc)
 * It computes an OPF dictionary estimator for each ordering combination at a certain two-level order-hierarchy, e.g.,
 * image category and position (for img-class and position), image category and object/background (for img-class and pixel), etc.
 * @implements iftBovwDictEstimator()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - ret_real_cent: Used to be passed to iftBovwUnsupOPFDictEstimator
 * - clust_samp_perc: Used to be passed to iftBovwUnsupOPFDictEstimator
 * - knn_graph_n_neighb: Used to be passed to iftBovwUnsupOPFDictEstimator
 * - known_n_groups: Whether or not we know the number of groups (str)
 * - n_groups: Number of groups (int)
 * - local_feats_order_hrchy: Order hierarchy for all the local features (iftIntMatrix*)
 *
 * @author Cesar Castelo
 * @date Feb 15, 2019, updated Oct 9, 2019
 */
void* iftBovwSupTwoLevelOrderHrchyOPFDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding);

/**
 * @brief Supervised one-level order-hierarchy Manual dictionary estimator for BoVW (e.g. pixel, img-class, position, etc)
 * It reads the centroids from a JSON file that contains the necessary information for each ordering at a certain
 * one-level order-hierarchy, e.g., object and background (for pixel), image category (for img-class), etc.
 * @implements iftBovwDictEstimator()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - sel_cent_json_filename: The name of the JSON file containing the necessary information
 *   If the order hierarchy is by position, the JSON file must contain these keys:
 *     {hierarchy_type, n_positions, position_1 (containing {x, y, dataset, selected_samples (int list)}, ..., position_n}
 *   If the order hierarchy is by img-class, the JSON file must contain these keys:
 *     {hierarchy_type, n_classes, class_1 (containing {dataset, selected_samples (int list)}, ..., class_n}
 * - local_feats_order_hrchy: Order hierarchy for all the local features (iftIntMatrix*)
 *
 * @author Cesar Castelo
 * @date Feb 15, 2019, updated Oct 9, 2019
 */
void* iftBovwSupOneLevelOrderHrchyManualDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding);

/**
 * @brief Supervised two-level order-hierarchy Manual dictionary estimator for BoVW (e.g. img-class and position,
 * img-class and pixel, etc)
 * It reads the centroids from a JSON file that contains the necessary information for each ordering at a certain
 * two-level order-hierarchy, e.g., image category and position (for img-class and position), image category and
 * object/background (for img-class and pixel), etc.
 * @implements iftBovwDictEstimator()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - sel_cent_json_filename: The name of the JSON file containing the necessary information
 *   If the order hierarchy is by img-class and position, the JSON file must contain these keys:
 *     {hierarchy_type, n_classes, n_positions, class_1_position_1 (containing {x, y, dataset, selected_samples (int list)}, ..., class_m_position_n}
 * - local_feats_order_hrchy: Order hierarchy for all the local features (iftIntMatrix*)
 *
 * @author Cesar Castelo
 * @date Feb 25, 2019, updated Oct 9, 2019
 */
void* iftBovwSupTwoLevelOrderHrchyManualDictEstimator(iftDataSet* localFeats, iftDict *funcParams, bool isFisherVectorsEncoding);

/**
 * @brief Hard clustering assignment function
 * @implements iftBovwCodFunc()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - sim_func: Similarity function to be used for the coding process (int)
 * - norm_type: Normalization to be used (int)
 * - pool_func: Pooling function to be used for dictionary size reduction (int)
 * - activ_func: Activation function to be used (int)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018, updated Oct 9, 2019
 */
iftFeatures* iftBovwHardAsgmtCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams);

/**
 * @brief Soft clustering assignment function
 * @implements iftBovwCodFunc()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - sim_func: Similarity function to be used for the coding process (int)
 * - norm_type: Normalization to be used (int)
 * - pool_func: Pooling function to be used for dictionary size reduction (int)
 * - activ_func: Activation function to be used (int)
 * - n_near_neighb: Number of neighbors to be used for the soft assignment (int)
 * - wgt_func: Weight function to be used for the soft assignment (str)
 *
 * @author Cesar Castelo
 * @date Jan 04, 2018, updated Oct 9, 2019
 */
iftFeatures* iftBovwSoftAsgmtCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams);

/**
 * @brief Hard clustering assignment function (with batch processing).
 * The dataset with the local feat vectors must be composed as follows:
 *   < n1 local_feat_vects for the n1 patches of the image1, n2 local_feat_vects for the n2 patches of the image2, etc >
 * The feature vector that is returned is composed as follows:
 *   < final_feat_vect for the image1, feat_vect for the image2, etc >
 * being that each final_feat_vect has the size of the visual dictionary
 * @implements iftBovwCodFunc()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - sim_func: Similarity function to be used for the coding process (int)
 * - norm_type: Normalization to be used (int)
 * - pool_func: Pooling function to be used for dictionary size reduction (int)
 * - activ_func: Activation function to be used (int)
 * - n_imgs_in_batch: Number of images in the batch
 *
 * @author Cesar Castelo
 * @date Feb 06, 2019, updated Oct 9, 2019
 */
iftFeatures* iftBovwHardAsgmtBatchCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams);

/**
 * @brief Soft clustering assignment function (with batch processing).
 * The dataset with the local feat vectors must be composed as follows:
 *   < n1 local_feat_vects for the n1 patches of the image1, n2 local_feat_vects for the n2 patches of the image2, etc >
 * The feature vector that is returned is composed as follows:
 *   < final_feat_vect for the image1, feat_vect for the image2, etc >
 * being that each final_feat_vect has the size of the visual dictionary
 * @implements iftBovwCodFunc()
 * @warning The parameters that should be sent to this function in funcParams are:
 * - sim_func: Similarity function to be used for the coding process (int)
 * - norm_type: Normalization to be used (int)
 * - pool_func: Pooling function to be used for dictionary size reduction (int)
 * - activ_func: Activation function to be used (int)
 * - n_near_neighb: Number of neighbors to be used for the soft assignment (int)
 * - wgt_func: Weight function to be used for the soft assignment (str)
 * - n_imgs_in_batch: Number of images in the batch
 *
 * @author Cesar Castelo
 * @date Feb 06, 2019, updated Oct 9, 2019
 */
iftFeatures* iftBovwSoftAsgmtBatchCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams);

/**
 * @brief Second-order statistics coding function (Fisher Vectors).
 * @implements iftFVCodFunc()
 * @warning This function requires the dictionary to be represented as a clustering model
 * @warning The parameters that should be sent to this function in funcParams are:
 * - dict_estim_id: Dictionary estimator method ID
 * - n_clust_models: Number of clustering models (only if using a supervised dictionary estimation method)
 * - n_clust_models_1: Number of clustering models in hierarchy 1 (only if using a supervised dictionary estimation method of two hierarchies)
 * - n_clust_models_2: Number of clustering models in hierarchy 2 (only if using a supervised dictionary estimation method of two hierarchies)
 *
 * @author Cesar Castelo
 * @date Oct 1, 2019
 */
iftFeatures* iftBovw2ndOrderStatsFisherVectorsCodFunc(void* dict, iftDataSet* localFeats, iftDict *funcParams);

/* ---------------------------------------------------------------------------------------------------------------------*/
/*                                                Complementary functions                                               */
/* ---------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Verifies if a given coding function produces a fisher vector encoding.
 * @param codFuncId Coding function ID
 * @return Whether it produces a fisher vector enconding
 * @author Cesar Castelo
 * @date Oct 10, 2019
 */
bool iftBovwIsFisherVectorsEncoding(iftBovwCodFuncId codFuncId);

/**
 * @brief Gets the number of visual words in a BoVW structure
 * @param bovw Input BoVW structure
 * @return Number of visual words
 * @author Cesar Castelo
 * @date Oct 10, 2019
 */
int iftBovwGetNumbVisualWords(iftBagOfVisualWords *bovw);

/**
 * @brief Gets the number of features per visual word in a BoVW structure
 * @param bovw Input BoVW structure
 * @return Number of features per visual word
 * @author Cesar Castelo
 * @date Oct 11, 2019
 */
int iftBovwGetNumbFeatsPerWord(iftBagOfVisualWords *bovw);

/**
 * @brief Computes the number of global features that will be produced
 * @param bovw Input BoVW structure
 * @return Number of global features
 * @author Cesar Castelo
 * @date Oct 10, 2019
 */
int iftBovwComputeNumbFeatsForCodFunc(iftBagOfVisualWords *bovw);

/**
 * @brief Create the fileset reference data for a given set of regions of interest.
 * It returns the fileset and also saves the image patches in the given patchDirPath
 * @param img Image to extract the patches to be saved
 * @param roiArray Regions of interest inside the image
 * @param imgFilePtr Pointer to the image file
 * @param patchDirPath Directory where the image patches will be saved
 * @author Cesar Castelo
 * @date Abr 23, 2018
 */
iftFileSet *iftBoVWCreateFileSetRefDataFromImgROIs(iftImage *img, iftRoiArray *roiArray, iftFile *imgFilePtr, char *patchDirPath);

/**
 * @brief Compute the number of groups per ordering combination for a specific hierarchy
 * @param localFeatsOrderHrchy Integer matrix that contains the entire order hierarchy (one or more levels)
 * @param nSamples Number of samples in the dataset
 * @param nGroups Total number of groups
 * @author Cesar Castelo
 * @date Feb 02, 2018
 */
iftIntArray* iftBovwComputeNumbGroupsPerOrderHrchy(iftIntMatrix *localFeatsOrderHrchy, int nSamples, int nGroups);

/**
 * @brief Compute a matrix with the ordering corresponding to a specific hierarchy, for one given image.
 * This method calls the method iftBovwComputeOrderHrchyArray for all the hierarchies present
 * @param localFeats Feature vectors representing the image patches
 * @param orderHrchies Order hierarchies
 * @param dictLearnFileset Fileset for dictionary learning
 * @param dictLearnMasksFileset Fileset with dictionary learning masks
 * @param imgId Image ID
 * @author Cesar Castelo
 * @date Feb 01, 2018
 */
iftIntMatrix *iftBovwComputeOrderHrchiesMatrixForImg(iftDataSet *localFeats, iftIntArray *orderHrchies, iftFileSet *dictLearnFileset,
                                                        iftFileSet *dictLearnMasksFileset, int imgId);

/**
 * @brief Compute the array with the ordering corresponding to a specific hierarchy
 * @param feats Feature vectors representing the image patches
 * @param hrchyId Hierarchy Id
 * @param nParams Number of variable parameters
 * @param ... Variable parameters
 * @author Cesar Castelo
 * @date Feb 01, 2018
 */
int *iftBovwComputeOrderHrchyArray(iftDataSet *feats, int hrchyId, int nParams, ...);

/**
 * @brief Sample a dataset based on the given order hierarchy
 * @param Z Dataset to be sampled
 * @param localFeatsOrderHrchy Matrix that contains the order hierarchy for the dataset
 * @param chosenLabels The specific combination of labels that will be used for the sampling
 * @date Feb 02, 2018
 * @author Cesar Castelo
 */
iftDataSet *iftBovwSampleDataSetByOrderHrchy(iftDataSet *Z, iftIntMatrix *localFeatsOrderHrchy, int *chosenLabels);

/**
 * @brief Prints messages regarding the chosen methods for each part of BoVW
 * @param intPointDetec Chosen interest point detection method (IFT_NIL to not print)
 * @param localFeatExtr Chosen local features extraction method (IFT_NIL to not print)
 * @param dictEstim Chosen dictionary estimation method (IFT_NIL to not print)
 * @param codFunc Chosen coding method (IFT_NIL to not print)
 * @date Feb 21, 2018
 * @author Cesar Castelo
 */
void iftBovwPrintChosenMethods(iftBovwIntPointDetectorId intPointDetec, iftBovwLocalFeatExtractorId localFeatExtr, iftBovwDictEstimatorId dictEstim, iftBovwCodFuncId codFunc);
char *iftBovwIntPointDetName(iftBovwIntPointDetectorId method, bool fullName);
char *iftBovwLocalFeatExtrName(iftBovwLocalFeatExtractorId method, bool fullName);
char *iftBovwDictEstimName(iftBovwDictEstimatorId method, bool fullName);
char *iftBovwCodFuncName(iftBovwCodFuncId method, bool fullName);

/**
 * @brief Saves an image with the interest point detected from it
 * @param img Original image
 * @param roiArray Set of interest points to be added to the image
 * @param imgFilePtr Pointer to the image file
 * @date July 16, 2018
 * @author Cesar Castelo
 */
void iftBovwSaveImgWithIntPoints(iftImage *img, iftRoiArray *roiArray, iftFile *imgFilePtr);

/**
 * @brief Computes the batch size (number of images) to be used to execute the coding functions that use batch processing.
 * It returns the batch size when we are usign either GPU or CPU
 * @param nVisWords Number of visual words
 * @param localFeatsDim Dimensionality of the local feat vectors/visual words
 * @param nPatches Number of patches per image
 * @param maxMemUsePerc Maximum GPU/CPU memory use percentage
 * @date Feb 07, 2019
 * @author Cesar Castelo
 */
int iftBovwComputeBatchSize(int nVisWords, int localFeatsDim, int nPatches, float maxMemUsePerc);

/**
 * @brief Converts a BovwWgtFunc name (str) into its corresponding int value
 * @author Cesar Castelo
 * @date Mar 19, 2019
 */
int iftBovwWgtFuncStrToWgtFunc(char *wgtFuncStr);

/**
 * @brief Creates a file containing the interest points (roiArray)
 * @author Cesar Castelo
 * @date Nov 29, 2019
 */
void iftBoVWCreateIntPointsFile(iftRoiArray *roiArray, iftDict *funcParams);

/**
 * @brief Creates a bank of random kernels according to the given image set
 * @param bovw Bag of Visual Words structure
 * @param dictLearnFileset Fileset containing the images for dictionary learning
 * @date Nov 26, 2019
 * @author Cesar Castelo
 */
void iftBovwCreateRandomKernels(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset);

/**
 * @brief Creates several banks of random kernels (one for each layer) according to the given image set and
 * network architecture (stored in the BoVW struct)
 * @param bovw Bag of Visual Words structure
 * @param dictLearnFileset Fileset containing the images for dictionary learning
 * @date Nov 26, 2019
 * @author Cesar Castelo
 */
void iftBovwCreateRandomKernelsMultiLayer(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset);

/**
 * @brief Creates a joint mask for interest point detection from a set of image masks
 * @param bovw Bag of Visual Words structure
 * @param dictLearnMasksFileset Fileset containing the image masks for dictionary learning
 * @date Nov 26, 2019
 * @author Cesar Castelo
 */
void iftBovwCreateJointSamplingMask(iftBagOfVisualWords *bovw, iftFileSet **dictLearnMasksFileset);

/**
 * @brief Creates class-specific masks for interest point detection using ISF
 * @param bovw Bag of Visual Words structure
 * @param dictLearnFileset Fileset containing the images for dictionary learning
 * @param dictLearnMasksFileset Fileset containing the image masks for dictionary learning
 * @date Nov 26, 2019
 * @author Cesar Castelo
 */
void iftBovwCreateSamplingMasksPerClassISF(iftBagOfVisualWords *bovw, iftFileSet *dictLearnFileset, iftFileSet *dictLearnMasksFileset);

/**
 * @brief Creates the directory to save the patches extracted from the interest points
 * @param bovw Bag of Visual Words structure
 * @date Nov 26, 2019
 * @author Cesar Castelo
 */
char *iftBovwCreateDirectoryForPatches(iftBagOfVisualWords *bovw);

/**
 * @brief Performs local batch centralization to a set of images given a set of parameters in the BoVW structure.
 * It saves the centralized images in the specified output folder.
 * @param bovw Bag of Visual Words structure containing all the params needed to compute the features
 * @param imgFileset Fileset with the images to extract the features from
 * 
 * @date Apr 30, 2019
 * @author Cesar Castelo
 */
void iftBovwApplyLocalBatchCentralization(iftBagOfVisualWords *bovw, iftFileSet *imgFileset);

/**
 * @brief Creates a dataset containing the local feature vectors extracted from a mimage that represents the features
 * extracted with some convolutional method (filtered images).
 * @param mimg Multiband image
 * @param img Original image
 * @param imgFilePtr Pointer to the file that contains the original image
 * @param roiArray Array of ROIs that were extracted from the image
 * @param resizeImage Whether or not to resize the resulting filtered image (str).
 *        If the image is not resized, the interest points previously detected are no longer used (all the positions from the reduced image are used instead).
 * @param vectConstr Feature vector construction method (str): 'pixel_vals_in_patch','mean_val_in_patch'.
 * 
 * @date Apr 23, 2019
 * @author Cesar Castelo
 */
iftDataSet *iftBovwBuildLocalFeatureVectorsFromConvolutionalMImage(iftMImage *mimg, iftImage *img, iftFile *imgFilePtr, iftRoiArray *roiArray, char *resizeImage, char *vectConstr);

/**
 * @brief Creates a dataset containing the BIC feature vectors extracted from a mimage that represents the features
 * extracted with some convolutional method (filtered images).
 * @param mimg Multiband image
 * @param img Original image
 * @param roiArray Array of ROIs that were extracted from the image
 * @param resizeImage Whether or not to resize the resulting filtered image (str).
 *        If the image is not resized, the interest points previously detected are no longer used (all the positions from the reduced image are used instead).
 * @param nBinsPerBand: Number of bins to extract BIC features (int).
 * 
 * @date Jun 18, 2019
 * @author Cesar Castelo
 */
iftDataSet *iftBovwExtractBICFeatureVectorsFromConvolutionalMImage(iftMImage *mimg, iftImage *img, iftRoiArray *roiArray, char *resizeImage, int nBinsPerBand);

/**
 * @brief Computes multilayer convolutional features from a set of images according to a given network architecture
 * indicated in the BoVW structure. It uses batch processing to speed-up the process.
 * It saves the resulting filtered images in the specified output folder.
 * @param bovw Bag of Visual Words structure containing all the params needed to compute the features
 * @param imgFileset Fileset with the images to extract the features from
 * 
 * @date Apr 30, 2019
 * @author Cesar Castelo
 */
void iftBovwComputeMultiLayerConvolutionalFeaturesBatch(iftBagOfVisualWords *bovw, iftFileSet *imgFileset);

/**
 * @brief Computes the cluster model params necessary for Fisher Vectors encoding
 * @param Z Dataset containing the groups
 * @return means Output matrix containing the means for each group (passed by reference)
 * @return covariances Output array of matrices containing the covariance matrices for each group (passed by reference)
 * @return covariancesInv Output array of matrices containing the covariance inverse matrices for each group (passed by reference)
 * @return covariancesDet Output array of matrices containing the covariance matrices determinant for each group (passed by reference)
 * @return weights Output float array containing the weights for each group (passed by reference)
 * 
 * @date Oct 10, 2019
 * @author Cesar Castelo
 */
void iftBovwComputeClusterModelParamsFV(iftDataSet *Z, iftMatrix **means, iftMatrix ***covariances, iftMatrix ***covariancesInv,
                                        iftFloatArray **covariancesDet, iftFloatArray **weights);

/**
 * @brief Create the feature vector for one sample using the 2nd order statistics given a set of clustering params and moments
 * @param localFeats Dataset with the local feats
 * @param means Matrix containing the means per group
 * @param covariances Array of matrices containing the covariance matrices per group
 * @param covariancesInv Output array of matrices containing the covariance inverse matrices for each group (passed by reference)
 * @param covariancesDet Output array of matrices containing the covariance matrices determinant for each group (passed by reference)
 * @param weights Float array containing the weights per group
 * @return The feature vector
 * 
 * @date Oct 10, 2019
 * @author Cesar Castelo
 */
iftFeatures *iftBovwComputeFeatVect2ndOrderStatsFV(iftDataSet *localFeats, iftMatrix *means, iftMatrix **covariances,
                                                    iftMatrix **covariancesInv, iftFloatArray *covariancesDet,
                                                    iftFloatArray *weights);


/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif //IFT_IFTBAGOFVISUALWORDS_H
