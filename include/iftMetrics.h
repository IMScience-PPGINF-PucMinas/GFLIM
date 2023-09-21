#ifndef IFT_METRICS_H
#define IFT_METRICS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftDataSet.h"


/**
 * @brief Type of the Metric Functions: Distance or Similarity Functions
 * @author Samuel Martins
 * @date Mar 23, 2016
 * @ingroup Metrics
 */
typedef enum {
    IFT_SIMILARITY_FUNCTION,
    IFT_DISTANCE_FUNCTION
} iftMetricType;


  
//True positives, false positives, false negatives, true negatives.
typedef struct ift_error_classification {
	int tp;
	int fp;
	int fn;
	int tn;
} iftErrorClassification;

float iftPrecisionGivenErrors(iftErrorClassification* error);
float iftRecallGivenErrors(iftErrorClassification* error);

float iftAccuracyGivenErrors(iftErrorClassification* error);
float iftFScoreGivenErrors(iftErrorClassification* error);
float iftProbabilityRatioGivenErrors(iftErrorClassification* error);

/**
 * Computes the median value for a float array.
 * @warning In case of an even size array, this method computes the average of both median values
 * @author Peixinho
 * @date Jul, 2017
 * @param array
 * @param length
 * @return median of array
 */
float iftMedianFloatArray(float* array, int length);

/**
 * Computes the weighted median value for a float array.
 * @warning In case of an even size array, this method computes the average of both median values
 * @author Peixinho
 * @date Jul, 2017
 * @param array
 * @param length
 * @param weight weight values (must be greater than zero)
 * @return median of array
 */
float iftWeightedMedianFloatArray (float *array, int length, float* weight);

/**
 * Computes the weighted average value for a float array.
 * @warning In case of an even size array, this method computes the average of both median values
 * @author Peixinho
 * @date Jul, 2017
 * @param array
 * @param length
 * @param weight weight values (must be greater than zero)
 * @return median of array
 */
float iftWeightedMeanFloatArray (float *array, int length, float* weight);

/**
 * Computes the median value for an int array.
 * @warning In case of an even size array, this method computes the rounded average of both median values
 * @author Peixinho
 * @date Jul, 2017
 * @param array
 * @param length
 * @return median of array
 */
int iftMedianIntArray(int* array, int length);

int iftMinIndexFloatArray(float *array, int length);
int iftMaxIndexFloatArray(float *array, int length);
float iftMinFloatArray(float* array, int length);
float iftMaxFloatArray(float* array, int length);
double iftSumDoubleArray (double *array, int length);
float iftSumFloatArray (float *array, int length);
int iftSumIntArray (int *array, int length);
float iftMeanFloatArray (float *array, int length);
void iftMultScalarFloatArray(float* array, int length, float scalar);
float iftMeanFloatArrayDiffZero(float *array, int length);
float iftVarianceFloatArray (float *array, int length);
float iftStddevFloatArray (float *array, int length);
bool iftCompareFloatArrays (float *array1, float *array2, int length);
bool iftCompareIntArrays (int *array1, int *array2, int length);

/**
 * Computes a given metric in a dataset
 * @param Z Dataset
 * @param metric The metric to be computed (full name)
 * @return score
 * @author Cesar Castelo
 * @date Aug, 2018
 */
float iftComputeMetricInDataSetByName(iftDataSet *Z, char *metric);
/**
 * Computes the cohen kappa score on dataset.
 * @param Z
 * @return score
 * @author Peixinho
 * @date Aug, 2017
 */
//! swig()
float iftCohenKappaScore(iftDataSet* Z);

float iftFscore(iftDataSet *Z, int PositiveClass);
float iftFscoreMulticlass(iftDataSet *Z);
float iftTruePositives(iftDataSet *Z);
float iftSkewedTruePositives(iftDataSet *Z);
float iftNormAccuracy(iftDataSet* Z);
float iftClassifAccuracy(iftDataSet *Z);
float iftOPFAccuracy(iftDataSet *Z, bool printAccByClass);
iftFloatArray *iftTruePositivesByClass(iftDataSet *Z);

float iftMeasureERR(float *negatives, int total_neg, float *positives, int total_pos, float threshold);
void  iftMeasureFarFrr(iftDataSet *Z, float *FAR, float *FRR);
iftErrorClassification iftMeasureErrorClassification(iftDataSet *Ztest);
float iftWeightedAccuracyBinaryClasses(iftErrorClassification cm);
float iftMeasureMCC(iftErrorClassification cm);

float iftCosineDistance(float *f1, float *f2, int n);
float iftCosineDistance2(float *f1, float *f2, int n);

float iftCorrCoef(float *x, float *y, int n);
float iftCorrCoefDist(float *x, float *y, int n);
float iftMICI(float *x, float *y, int n);


/**
 * @brief It Finds the best Sample Score from a Score Matrix, returning the indice of the sample in the matrix.
 *
 * Each row i=p and column j=p from the Score Matrix corresponds to the sample.\n
 * The Best Score is from the Sample whose sum of all scores is MAXIMUM (for SIMILARITY functions) or MINIMUM (for DISTANCE functions).\n
 *
 * @param score_matrix Score Matrix that will used to figure out the best Sample Score.
 * @param metric_type Type of the Score Function that generated the considered Score Matrix.
 * @param out_best_score Reference of the variable which will store the score of the Most Similar Sample. If NULL, nothing is returned to it.
 * @return The index of the Most Similar Sample.
 *
 * @author Samuel Martins
 * @date Nov 07, 2017
 * @ingroup Metrics
 */
int iftBestSampleScoreSum(const iftMatrix *score_matrix, iftMetricType metric_type, double *out_best_score);


/**
 * @brief Computes the Kullback-Leibler divergence between two distributions.
 * 
 * @note The distributions are normalized before computation.
 * @note The distributions must have the same size.
 * 
 * KL divergence is defined when the distributions have ONLY NONZERO entries.
 * To circumvent that, we use the following function as used by scipy.special.rel_entr [1]:
 * 
 * rel_entr(Pi, Qi) = Pi * log (Pi / Qi), if Pi > 0,  Qi > 0
 *                  = 0                 , if Pi == 0, Qi >= 0
 *                  = inf               , otherwise
 * 
 * [1] https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rel_entr.html#scipy.special.rel_entr
 * 
 * 
 * @param P Defines the first (discrete) distribution.
 * @param Q (Discrete) distribution against which the KL diverse is computed.
 * @return  The calculated KL divergence.
 * 
 * @author Samuel Martins
 * @date Aug 21, 2018
 */
float iftKullbackLeiblerDivergence(const iftFloatArray *P, const iftFloatArray *Q);

double iftMImageMahalanobis(const iftMImage *mimg, const iftMatrix *M, int p, int q);

double iftMImageMahalanobisDouble(const iftMImage *mimg, const double *M, int p, int q);

#ifdef __cplusplus
}
#endif

#endif
