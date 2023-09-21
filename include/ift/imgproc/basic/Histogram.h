//
// Created by samuel on 18/12/18.
//

#ifndef IFT_HISTOGRAM_H
#define IFT_HISTOGRAM_H


#include "ift/core/dtypes/BasicDataTypes.h"
#include "iftImage.h"
#include "iftMImage.h"


/**
 * @brief Struct for a histogram of float values.
 * @author Samuel Martins
 * @date Dec 18, 2018
 */
//! swig(extend = iftHistogramExt.i, destroyer = iftDestroyHist)
typedef struct ift_hist {
    /** Number of bins */
    long nbins;
    /** Histogram array */
    float *val;
} iftHist;



/**
 * @brief Creates a Histogram of <nbins> bins.
 * @author Samuel Martins
 * @date Dec 18, 2018
 */
//! swig(newobject)
iftHist *iftCreateHist(long nbins);


/**
 * @brief Destroys a Histogram.
 * @author Samuel Martins
 * @date Dec 18, 2018
 */
void iftDestroyHist(iftHist **hist);


/**
 * @brief Prints a histogram on console.
 * @author Samuka
 * @date Jan 3, 2017
 */
//! swig()
void iftPrintHist(const iftHist *hist);

/**
 * @brief Computes a Histogram for a Gray Image.
 * 
 * The function expects the maximum value max_val of the image for histogram computation. That is, the image range: [0, max_val].
 * The bisize is simply computed by: binsize = (max_val + 1) / nbins
 * 
 * Common Usage:
 * (i)  max_val = iftMaximumValue(img) and nbins = max_val + 1
 * (ii) max_val = iftNormalizationValue(iftMaximumValue(img)) and nbins = max_val + 1 ===> image range [0, 2^bits - 1]
 * * 
 * @warning If the mask is not NULL, it must have the same domain and voxel size from the input image.
 *
 * @param img Input gray image.
 * @param mask Optional mask image. If it is not NULL, the non-zero mask elements mark the array elements counted in the histogram.
 * @param nbins Number of bins from the Output Histogram.
 * @param max_val Maximum value considered during histogram computation. The range [0, max_val] is considered to compute
 *                the bin size.
 * @param normalize If true, the histogram is normalized.
 * @return The computed histogram.
 *
 * @author Samuel Martins
 * @date Dec 19, 2018
 */
//! swig(newobject)
iftHist *iftCalcGrayImageHist(const iftImage *img, const iftImage *mask, long nbins, int max_val, bool normalize);

/**
 * @brief Computes a histogram for each object from an input gray image.
 *
 * The function expects the maximum value max_val of the image for histogram computation. That is, the image range: [0, max_val].
 * The bisize is simply computed by: binsize = (max_val + 1) / nbins.
 * 
 * @note The number of histograms is 1 + n_objs, where n_objs is the number of objects (maximum value) of the label image.
 * Thus, the histogram[x] corresponds to the histogram for the object with label x.
 *
 * @param img       The input grayscale image.
 * @param label_img Label image with the target objects.
 * @param nbins     The number of bins into which the histogram should be divided.
 * @param max_val Maximum value considered during histogram computation. The range [0, max_val] is considered to compute
 *                the bin size.
 * @param normalize If true, normalizes the histogram by dividing it by the number of pixels in the binary mask.
 * @param           n_hists_out Reference to return the number os computed histograms. If NULL, it is ignored.
 * @return          An array of histograms, one for each label.
 *
 * @author Samuel Martins
 * @date Aug 24, 2018
 */
//! swig(newobject)
iftHist **iftCalcGrayImageHistForLabels(const iftImage *img, const iftImage *label_img, int nbins, int max_val, bool normalize,
                                        int *n_hists_out);


/**
 * @brief Computes a Histogram for the RGB Values from a Color Image.
 *
 * @warning If the mask is not NULL, it must have the same domain and voxel size from the input image.
 *
 * @param img Input Color Image.
 * @param mask Optional mask image. If it is not NULL, the non-zero mask elements mark the array elements counted in the histogram.
 * @param nbins Number of bins from the Output Histogram.
 * @param normalize If true, the histogram is normalized.
 * @return The computed histogram from the RGB values.
 */
//! swig(newobject)
iftHist *iftCalcColorImageHist(const iftImage *img, const iftImage *mask, int nbins, bool normalize);


/**
 * @brief Computes the accumulate histogram from a histogram.
 * @param hist Input histogram.
 * @return Accumulate histogram.
 */
//! swig(newobject)
iftHist *iftCalcAccHist(const iftHist *hist);


/**
 * @brief Normalizes a histogram between 0 and 1
 * @author Ruppert
 * @date Aug 24, 2018
 */
//! swig(newobject)
iftHist *iftNormalizeHist(const iftHist *src);


/**
 * @brief Computes the Mode a histogram.
 * @param hist Histogram.
 * @param exclude_zero If true, the zero bin is not considered during the computation.
 * @return Histogram mode.
 */
//! swig(newobject)
int iftHistMode(const iftHist *hist, bool exclude_zero);


/**
 * @brief Computes the mean value of a histogram.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param hist The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The mean value
 *
 */
//! swig(newobject)
double iftHistMean(const iftHist *hist, bool exclude_zero);


/**
 * @brief Computes the median value of a histogram.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param hist The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The median value
 *
 */
//! swig(newobject)
int iftHistMedian(const iftHist *hist, bool exclude_zero);


/**
 * @brief Gets the Arg (Bin) Min of a histogram. For ties, the first arg min value found is selected.
 * @author Samuel Martins
 * @date Jun 2, 2016
 *
 * Ex: For a histogram with 256 values/bins of a Gray Image, the value/bin (in [0, 255]) with the lowest incidence is selected.
 *
 * @param hist The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The arg (bin) min value.
 *
 */
//! swig(newobject)
int iftHistArgMin(const iftHist *hist, bool exclude_zero);


/**
 * @brief Gets the Arg (Bin) Max of a histogram. For ties, the first arg max value found is selected.
 * @author Samuel Martins
 * @date Jun 2, 2016
 *
 * Ex: For a histogram with 256 values/bins of a Gray Image, the value/bin (in [0, 255]) with the highest incidence is selected.
 *
 * @param hist The input histogram
 * @param exclude_zero If true, the 0 bin will not be counted.
 * @return The arg (bin) max value.
 *
 */
//! swig(newobject)
int iftHistArgMax(const iftHist *hist, bool exclude_zero);


/**
 * @brief Return the first bin (argument) whose value is greater than <thres>. If no bin has value greater than <thes>,
 * it returns -1.
 *
 * @param hist Histogram for checking.
 * @param thres Value/tresholding used to find the first bin greater than it.
 * @param exclude_zero If true, the bin 0 is not considered.
 * @return The first bin whose value is greater than <thres>, or -1 if no bin is greater thatn <thres>.
 *
 * @author Samuka Martins
 * @date Nov 21, 2018
 */
//! swig(newobject)
int iftHistArgGreaterThan(const iftHist *hist, float thres, bool exclude_zero);


/**
 * @brief Adds the bins of two histograms.
 *
 * @author Thiago Vallin Spina
 * @date Feb 21, 2016
 *
 * @param hist1 The first histogram.
 * @param hist2 The second histogram.
 * @return The added histogram.
 */
//! swig(newobject)
iftHist *iftAddHists(const iftHist *hist1, const iftHist *hist2);


/**
 * @brief Adds the histogram src to dst.
 * @author Samuel Martins
 * @author Thiago Vallin Spina
 * @date APr 14, 2016
 *
 * @param src The source histogram.
 * @param dst The destination histogram.
 */
//! swig(newobject)
void iftAddHistsInPlace(const iftHist *src, iftHist *dst);

#endif //IFT_HISTOGRAM_H
