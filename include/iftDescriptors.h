#ifndef IFT_DESCRIPTORS_H_
#define IFT_DESCRIPTORS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/stat.h>

#include "iftCommon.h"
#include "iftImage.h"
#include "iftMImage.h"
#include "iftAdjacency.h"
#include "iftRadiometric.h"

/**
 * @brief Feature vector
 * @details [long description]
 * @author Adan Echemendia
 * @date Jun 10, 2016
 */
typedef struct ift_features {

    /** List of feature values */
    float *val;
    /** Size of the feature vector */
    int    n;
} iftFeatures;

typedef struct ift_hog {
    iftVoxel cellSize;
    iftVoxel blockSize;
    iftVoxel blockStride;
    int nbins;
} iftHog;

/**
 * @brief Creates a feature vector
 *
 * @param n Size of the feature vector
 * @return A feature vector
 */
iftFeatures *iftCreateFeatures(int n);

/**
 * @brief Destroys a feature vector
 *
 * @param feat Pointer to pointer to a feature vector
 */
void         iftDestroyFeatures(iftFeatures **feat);

/**
 * @brief Reads a set of features from a file and builds a feature vector
 * @details [long description]
 *
 * @param filename File path
 * @return A feature vector
 */
iftFeatures *iftReadFeatures(const char *filename);

/**
 * @brief Writes a feature vector to a file
 *
 * @param feats Feature vector
 * @param filename File path
 */
void         iftWriteFeatures(const iftFeatures *feats, const char *filename);

/**
 * @brief Joins an array of feature vectors into one feature vector
 * @param array Array of feature vectors
 * @param n Size of the array
 * @return A joint feature vector
 * @author Cesar Castelo
 * @date Oct 10, 2019
 */
iftFeatures *iftJoinFeatures(iftFeatures **array, int n);

/**
 * @brief [brief description]
 * @details [long description]
 *
 * @param img [description]
 * @param A [description]
 *
 * @return [description]
 */
iftImage *iftLocalBinaryPattern(const iftImage *img, const iftAdjRel *A);

/**
 * @brief [brief description]
 * @details [long description]
 *
 * @param img [description]
 * @param A [description]
 *
 * @return [description]
 */
iftFeatures *iftExtractLBP(iftImage *img, iftAdjRel *A);


/**
 * @brief Codes a 3D-LBP from a 3D image.
 * 
 * Codes a 3D-LBP from a 3D image by using the algorithm LBP-TOP proposed in [1].
 * For each orientation (XY, ZX, YX), it computes a 2D LBP for each corresponding slice into
 * a 3D image.
 * We use 8-neighborhood to compute the 2D LBP values.
 * The function returns an array with 3 images that will be used to extract the feature vectors (@see iftExtract3DLBPTOPFeats)
 * [0] = 3D LBP (XY)
 * [1] = 3D LBP (ZX)
 * [2] = 3D LBP (YZ)
 * 
 * [1] Zhao, Guoying, and Matti Pietikainen. "Dynamic texture recognition using local binary patterns
 * with an application to facial expressions." IEEE transactions on pattern analysis and
 * machine intelligence 29.6 (2007): 915-928.
 * 
 * @param img Input 3D Image.
 * @return Returned array of 3D-LBP images.
 */
iftImageArray *ift3DLBPTOP(const iftImage *img);


/**
 * @brief Extracts the LBP features from an input 3D Image.
 * 
 * This function is based on the algorithm LBP-TOP proposed in [1].
 * Given the 3 LBP images, one for each orientation, coded by the function @see ift3DLBPTOP,
 * the resulting feature vector is the concatenation of the histogram of each LBP image.
 * 
 * If a mask is passed, the histograms are only computed inside it.
 * 
 * @param img Input 3D image.
 * @param mask (Optional) Mask that indicates the voxels to be considered during feature extraction.
 * @param n_bins Number of bins of each one of the 3 LBP histograms. Ex: 256
 * @param normalize_histograms Normalize histogram features.
 * @return Resulting LBP feature vector.
 * 
 * @author Samuel Martins
 * @date Aug 27, 2018
 */
iftFeatures *iftExtract3DLBPTOPFeats(const iftImage *img, const iftImage *mask,
                                     int n_bins, bool normalize_histograms);

/**
 * @brief Extracts the LBP features from an input 3D Image for each object of a label image.
 * 
 * This function is based on the algorithm LBP-TOP proposed in [1].
 * It extracts a LBP-TOP for each object of a label image.
 * The return is a matrix (max_label + 1, 3 * n_bins), where max_label the highest label value of label_img.
 * Therefore, row [i] corresponds to the feature vector for object with label i.
 * 
 * @param img Input 3D image.
 * @param label_img Label image with objects of interest.
 * @param n_bins Number of bins of each one of the 3 LBP histograms. Ex: 256
 * @param normalize_histograms Normalize histogram features.
 * @return iftMatrix* (n_objs + 1, 3 * n_bins) Matrix with the LBP feats for each object.
 *         Row [i] contains the LBP feats for object with label i.
 * 
 * @author Samuel Martins
 * @date Aug 27, 2018
 */
//! swig(newobject)
iftMatrix *iftExtract3DLBPTOPFeatsForLabels(const iftImage *img, const iftImage *label_img,
                                            int n_bins, bool normalize_histograms);

/**
 * @brief Codes a 3D-LBP from a 3D image using the method VLBP proposed in [1].
 * 
 * [1] Zhao, Guoying, and Matti Pietikainen. "Dynamic texture recognition using local binary patterns
 * with an application to facial expressions." IEEE transactions on pattern analysis and
 * machine intelligence 29.6 (2007): 915-928.
 * 
 * @param img Input 3D image.
 * @return Resulting LBP image.
 * 
 * @author Samuel Martins
 * @date Aug 28, 2018
 */
iftImage *iftVLBP(const iftImage *img);


/**
 * @brief Extracts the LBP features from an input 3D Image using the method VLBP proposed in [1].
 * 
 * [1] Zhao, Guoying, and Matti Pietikainen. "Dynamic texture recognition using local binary patterns
 * with an application to facial expressions." IEEE transactions on pattern analysis and
 * machine intelligence 29.6 (2007): 915-928.
 * 
 * @param img Input 3D image.
 * @param mask Mask that indicates the voxels to be considered during feature extraction.
 * @param n_bins Number of bins of the LBP histogram. If <= 0, it uses the default number 2ˆ14.
 * @param normalize_hist Normalize the LBP histogram.
 * @return Resulting LBP feature vector.
 * 
 * @author Samuel Martins
 * @date Aug 28, 2018
 */
iftFeatures *iftExtractVLBPFeats(const iftImage *img, const iftImage *mask,
                                 int n_bins, bool normalize_hist);

/**
 * @brief Extracts the LBP features from an input 3D Image, for each object of a label image,
 * using the method VLBP proposed in [1].
 * 
 * [1] Zhao, Guoying, and Matti Pietikainen. "Dynamic texture recognition using local binary patterns
 * with an application to facial expressions." IEEE transactions on pattern analysis and
 * machine intelligence 29.6 (2007): 915-928.
 * 
 * @param img Input 3D image.
 * @param label_img Label image with objects of interest.
 * @param n_bins Number of bins of the LBP histogram. If <= 0, use the default number 2ˆ14.
 * @param normalize_hist Normalize the LBP histogram.
 * @return iftMatrix* (n_objs + 1, n_bins) Matrix with the LBP feats for each object.
 *         Row [i] contains the LBP feats for object with label i.
 * 
 * @author Samuel Martins
 * @date Dec 3, 2019
 */
//! swig(newobject)
iftMatrix *iftExtractVLBPFeatsForLabels(const iftImage *img, const iftImage *label_img,
                                        int n_bins, bool normalize_hist);

/**
 * @brief Builds a feature vector with the normalized brightness values of an image
 * @details [long description]
 *
 * @param img Target image
 * @return Feature vector with the normalized brightness values of an image
 */
iftFeatures *iftExtractBrightness(iftImage *img);

/**
 * @brief Builds a feature vector with the normalized brightness values of the gradient magnitude image of image <br>img</br>
 * @details [long description]
 *
 * @param img Target image
 * @param A Adjacency relation
 *
 * @return Feature vector with the normalized brightness values of the gradient magnitude image of image <br>img</br>
 */
iftFeatures *iftExtractGradient(iftImage *img, iftAdjRel *A);


/**
 * @brief extract BIC features for a mimage. It extracts the BIC features for each band and then joins the histograms of all bands
 * @author Cesar Castelo
 * @date Jun 13, 2019
 * @details
 *
 * @param mimg input mimage
 * @param bins_per_band number of bins per band
 *
 * @return [description]
 */
iftFeatures *iftMExtractBIC(iftMImage *mimg, int bins_per_band);

/**
 * @brief extract BIC features performing the quantization process by clustering in 1D (e.g. a grayscale image)
 * @author Cesar Castelo
 * @date Nov 28, 2017
 * @details
 *
 * @param img input image in the YCbCr color space
 * @param bins_per_band number of bins per band
 *
 * @return [description]
 */
 iftFeatures *iftExtractBICClustering1D(iftImage *img, int bins_per_band);


 /**
 * @brief Extract the BIC descriptor for an image. If a binary mask is passed (!= NULL),
 * BIC is computed for the object of it.
 * 
 * Given a number of bins per channel n_bins_per_channel, the input image
 * is firstly quantized to n_bins_per_channel.
 * Then, a histogram for the borders and the interior of the quantized image
 * are computed.
 * The number of bins n_bins of these histograms corresponds to
 * n_bins = n_bins_per_channel for gray images, and
 * n_bins =  (n_bins_per_channel ˆ3) for color images.
 * 
 * Both histograms are normalized to [0, 255].
 * Finally, a log2 is applied to each histogram bins so that the resulting
 * log2 histograms are combined to form the resulting feature vector
 * with 2 * n_bins.
 * 
 * If a binary mask is passed, BIC is only computed for the object of it.
 * If such binary mask is actually a multi-label image, this function will
 * binary it then considering all objects as a single one.
 * 
 * @param img Image.
 * @param bin_mask Binary mask. If != NULL, BIC is computed for the object of it.
 * @param n_bins_per_channel Number of bins per channel. 
 * @return BIC features for the image.
 */
 iftFeatures *iftExtractBIC(const iftImage *img, const iftImage *bin_mask,
                            int n_bins_per_channel);

 /**
 * @brief Computes a BIC feature vector for each object of a Label Image.
 * 
 * This function returns a Matrix (n_objs + 1) x (2 * bins), where each
 * row [i] contains the BIC feats for the object label i.
 * Label 0 consists to the Background.
 * 
 * The input image is firstly quantized with <n_bins_per_channel> per channel/band.
 * Thus, a quantized gray image will have <n_bins_per_channel> possible values
 * whereas a quantized color image will have (n_bins_per_channel ^ 3) possible values.
 * This number of possible values corresponds to the final number of bins <bins>
 * of each histogram.
 * Therefore, the final feature vector has (2 * bins).
 * 
 * PS: Resulting histogram values are normalized to [0, norm_val], such that
 * norm_val is the normalization value (maximum integer value) of the input image.
 * E.g.: 255 for 8-bit image, 4095 for 12-bit image, etc.
 * 
 * @param img (Gray/Color) Image used for BIC.
 * @param label_img Label Image with the target objects.
 * @param n_bins_per_channel Number of bins per image band/channel.
 * @return Matrix with the feature vector for each object. The matrix shape is
 *         (n_objs + 1) x (2 * bins), such that each row [i] consists of the
 *         feature vector for the object label i.
 * 
 * @author Samuel Martins
 * @date Sep 30, 2019
 */
 //! swig(newobject)
 iftMatrix *iftExtractBICForLabels(const iftImage *img, const iftImage *label_img,
                                   int n_bins_per_channel);

 /**
 * @brief Computes the Manhattan distance between two feature vectors
 * @details [long description]
 *
 * @param feat1 Feature vector one
 * @param feat2 Feature vector two
 *
 * @return Manhattan distance value between the two feature vectors
 */
 float iftFeatDistL1(iftFeatures *feat1, iftFeatures *feat2);

 /**
 * @brief Builds a feature vector with all values of a multi-band image
 * @details [long description]
 *
 * @param img Target image
 * @return Feature vector with all values of a multi-band image
 */
 iftFeatures *iftMImageToFeatures(iftMImage *img);

 /**
 * @brief Writes a feature vector to a file. The first saved value is the true label of the sample.
 * @details [long description]
 *
 * @param features Feature vector of the sample
 * @param truelabel True label of the sample
 * @param fp File path
 */
 void iftWriteFeaturesInFile(iftFeatures *features, int truelabel, FILE *fp);

 /**
 * @brief [brief description]
 * @details [long description]
 *
 * @param img [description]
 * @param u1 [description]
 * @param u2 [description]
 * @return [description]
 */
 iftFeatures *iftIntensityProfile(iftImage *img, iftPoint u1, iftPoint u2);

 /**
 * @brief [brief description]
 * @details [long description]
 *
 * @param img [description]
 * @param bins_per_band [description]
 * @param normalization_value [description]
 * @return [description]
 */
 iftImage *iftQuantize(const iftImage *img, int bins_per_band, int normalization_value);

 /**
 * @brief quantize image by clustering in 1D (e.g. for a grayscale image)
 * @details [long description]
 * @author Cesar Castelo
 * @date Nov 28, 2017
 *
 * @param img input image
 * @param bins_per_band number of bins per band
 * @param normalization_value normalization value
 * @param knn_create_method method to be used to create the knn graph for OPF clustering
 * @return [description]
 */
 iftImage *iftQuantizeByClustering1D(iftImage *img, int bins_per_band);

 /**
 * @brief Creates and defines the HOG feature extractor parameters.
 * @param cellSize The cell size considered in the histogram computation.
 * @param blockSize The block size defines the number of cells to be considered in the normalization.
 * @param blockStride The block stride defines the overlap between blocks in the normalization.
 * @author Peixinho
 * @date April, 2016
 */
 iftHog *iftCreateHog2D(int cellSize, int blockSize, int blockStride, int nbins);

 /**
 * @brief Destroys a HOG
 *
 * @param Pointer to pointer to a HOG
 */
 void iftDestroyHog(iftHog **hog);

 /**
 * @brief Extracts the Hog descriptor of an image.
 * @author Peixinho
 * @date April, 2016
 */
 iftFeatures *iftExtractHOG2D(const iftHog *hog, iftImage *img);

#ifdef __cplusplus
}
#endif

#endif
