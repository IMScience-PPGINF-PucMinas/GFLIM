#include "ift/imgproc/basic/Histogram.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/LongArray.h"
#include "ift/core/io/Stream.h"
#include "iftCommon.h"
#include "iftMemory.h"


/*********************** PRIVATE FUNCTIONS ***********************/
iftHist *_iftCalcHist(const int *val, long n, const int *mask, long nbins, int max_val, bool normalize) {
    if (val == NULL)
        iftError("Input Integer Array is NULL", "_iftCalcHist");
    if (nbins <= 0)
        iftError("Invalid nbins: %d <= 0", "_iftCalcHist", nbins);
    if (max_val < 0)
        iftError("Invalid negative range max: %d", "_iftCalcHist", max_val);

    double binsize = iftMax((max_val + 1.0) / nbins,1);

    long n_elements = 0;
    iftHist *hist = iftCreateHist(nbins);
    
    if (mask) {
        for (int p = 0; p < n; p++) {
            if (mask[p]) {
                int bucket = (int) (val[p] / binsize);
                hist->val[bucket]++;
                n_elements++;
            }
        }
    }
    else {
        n_elements = n;
        for (int p = 0; p < n; p++) {
            int bucket = (int)(val[p] / binsize);
            hist->val[bucket]++;
        }
    }

    if (normalize) {
        n_elements = iftMax(n_elements, 1);

        #pragma omp parallel for
        for (int b = 0; b < hist->nbins; b++) {
            hist->val[b] /= n_elements;
        }
    }
    
    return hist;
}


/*********************** PUBLIC FUNCTIONS ***********************/
iftHist *iftCreateHist(long nbins) {
    iftHist *hist = (iftHist *) iftAlloc(1, sizeof(iftHist));
    
    hist->nbins = nbins;
    hist->val   = iftAllocFloatArray(nbins);
    
    return hist;
}


void iftDestroyHist(iftHist **hist) {
    if (hist) {
        iftHist *aux = *hist;
        
        if (aux != NULL) {
            iftFree(aux->val);
            iftFree(aux);
            *hist = NULL;
        }
    }
}


void iftPrintHist(const iftHist *hist) {
    printf("hist = [");
    for (int b = 0; b < (hist->nbins - 1); b++)
        printf("%.4lf, ", hist->val[b]);
    
    printf("%.4lf]\n", hist->val[hist->nbins - 1]);
}


iftHist *iftCalcGrayImageHist(const iftImage *img, const iftImage *mask, long nbins, int max_val, bool normalize) {
    if (img == NULL)
        iftError("Image is NULL", "iftCalcGrayImageHist");	
    if (iftIsColorImage(img))
        iftWarning("Image is colored. Only the Y channel will be used", "iftCalcGrayImageHist");
    
    if (mask) {
        iftVerifyImages(img, mask, "iftCalcGrayImageHist");	
        return _iftCalcHist(img->val, img->n, mask->val, nbins, max_val, normalize);
    }
    else {
        return _iftCalcHist(img->val, img->n, NULL, nbins, max_val, normalize);
    }
}

iftHist **iftCalcGrayImageHistForLabels(const iftImage *img, const iftImage *label_img, int nbins, int max_val, bool normalize,
                                        int *n_hists_out) {
    iftVerifyImageDomains(img, label_img, "iftCalcGrayImageHistForLabels");

    int n_labels = iftMaximumValue(label_img);
    int n_hists = n_labels + 1; // includes a "Zero-Object" to make the computation easier

    double binsize = (max_val + 1.0) / nbins;

    iftLongArray *obj_volumes = iftCreateLongArray(n_hists);
    iftHist **hists = iftAlloc(n_hists, sizeof(iftHist *));

    #pragma omp parallel for
    for (int label = 0; label < n_hists; label++)
        hists[label] = iftCreateHist(nbins);

    for (int p = 0; p < img->n; p++) {
        int label = label_img->val[p];
        long bin = (long) (img->val[p] / binsize);

        hists[label]->val[bin]++;
        obj_volumes->val[label]++;
    }

    if (normalize) {
        #pragma omp parallel for
        for (int label = 0; label <= n_labels; label++) {
            for (int bin = 0; bin < hists[label]->nbins; bin++) {
                hists[label]->val[bin] /= iftMax(obj_volumes->val[label], 1);
            }
        }
    }

    iftDestroyLongArray(&obj_volumes);

    if (n_hists_out) {
        *n_hists_out = n_hists;
    }

    return hists;
}


iftHist *iftCalcColorImageHist(const iftImage *img, const iftImage *mask, int nbins, bool normalize) {
    iftHist *hist;
    int p, b;
    iftColor RGB, YCbCr;
    long normalization_value = iftNormalizationValue(iftMaximumValue(img));
    int xsize = (pow((double)nbins,(1.0/3.0))+0.5);
    nbins = xsize*xsize*xsize;    
    float binsize = (normalization_value+1) / (float)xsize;
    int xysize = xsize*xsize;
    int nelems = 0;
    
    
    if (!iftIsColorImage(img)||(iftMinimumValue(img)<0))
        iftError("It is not a valid color image", "iftCalcColorImageHist");
    
    if (xsize < 1)
        iftError("Insufficient number of bins", "iftCalcColorImageHist");
    else
    if (xsize > normalization_value+1)
        iftError("Excessive number of bins", "iftCalcColorImageHist");
    
    hist = iftCreateHist(nbins);
    for (p=0; p < img->n; p++){
      if (mask->val[p]!=0){
	nelems++;
        YCbCr.val[0] = img->val[p];
        YCbCr.val[1] = img->Cb[p];
        YCbCr.val[2] = img->Cr[p];
        RGB          = iftYCbCrtoRGB(YCbCr,normalization_value);
        b = (int)( ((float)RGB.val[0]) / binsize )
            + xsize*(int)( ((float)RGB.val[1]) / binsize )
            + xysize*(int)( ((float)RGB.val[2]) / binsize );
        hist->val[b]++;
      }
    }
    
    if (normalize)
        for (b=0; b < hist->nbins; b++) {
            hist->val[b] /= nelems;
        }
    
    return(hist);
}


iftHist *iftCalcAccHist(const iftHist *hist) {
    iftHist *acchist = iftCreateHist(hist->nbins);
    
    acchist->val[0] = hist->val[0];
    
    for (int i = 1; i < hist->nbins; i++)
        acchist->val[i] = acchist->val[i - 1] + hist->val[i];
    
    return acchist;
}


iftHist *iftNormalizeHist(const iftHist *src){
    if (src == NULL)
        iftError("Histogram is NULL", "iftNormalizeHist");
    iftHist *result = iftCreateHist(src->nbins);
    
    float sum=0;
    for (int b = 0; b < src->nbins; b++)
        sum += src->val[b];
    for (int b = 0; b < src->nbins; b++)
        result->val[b]= src->val[b]/(float)sum;
    return(result);
}



int iftHistMode(const iftHist *hist, bool exclude_zero)
{
    int i, i0=0, mode;
    
    if (exclude_zero)
        i0 = 1;
    
    mode=i0;
    for (i=i0+1; i < hist->nbins; i++) {
        if (hist->val[i] > hist->val[mode])
            mode = i;
    }
    
    return(mode);
}


double iftHistMean(const iftHist *hist, bool exclude_zero) {
    double mean = 0.0, total_freq = 0.0;
    
    int i0 = 0;
    if (exclude_zero)
        i0 = 1;

    #pragma omp parallel for reduction(+:mean), reduction(+:total_freq)
    for (int i = i0; i < hist->nbins; i++) {
        mean       += hist->val[i]*i;
        total_freq += hist->val[i];
    }
    
    return (mean / iftMax(1, total_freq));
}


int iftHistMedian(const iftHist *hist, bool exclude_zero)
{
    int i=0, i0=0, median;
    double acc = 0.0, total_freq = 0.0;
    
    if (exclude_zero)
        i0 = 1;
    
    for (i=i0; i < hist->nbins; i++) {
        total_freq += hist->val[i];
    }
    
    median = 0;
    for (i=i0; i < hist->nbins; i++) {
        if(hist->val[i] > 0) {
            if(acc <= total_freq / 2.0) {
                median = i;
            }
            acc += hist->val[i];
        }
    }
    
    return median;
}


int iftHistArgMin(const iftHist *hist, bool exclude_zero) {
    if (hist == NULL)
        iftError("Histogram is NULL", "iftHistArgMin");
    
    double min  = IFT_INFINITY_DBL;
    int arg_min = -1;
    
    int start_bin = 0;
    if (exclude_zero)
        start_bin = 1;
    
    for (int bin = start_bin; bin < hist->nbins; bin++) {
        if (hist->val[bin] < min) {
            min     = hist->val[bin];
            arg_min = bin;
        }
    }
    
    return arg_min;
}


int iftHistArgMax(const iftHist *hist, bool exclude_zero) {
    if (hist == NULL)
        iftError("Histogram is NULL", "iftHistArgMax");
    
    double max  = IFT_INFINITY_DBL_NEG;
    int arg_max = -1;
    
    int start_bin = 0;
    if (exclude_zero)
        start_bin = 1;
    
    for (int bin = start_bin; bin < hist->nbins; bin++) {
        if (hist->val[bin] > max) {
            max     = hist->val[bin];
            arg_max = bin;
        }
    }
    
    return arg_max;
}


int iftHistArgGreaterThan(const iftHist *hist, float thres, bool exclude_zero) {
    int start_bin = (exclude_zero) ? 1 : 0;
    
    for (int b = start_bin; b < hist->nbins; b++)
        if (hist->val[b] > thres)
            return b;
    
    return -1;
}


iftHist *iftAddHists(const iftHist *hist1, const iftHist *hist2) {
    int i;
    iftHist *result = NULL;
    
    if(hist1->nbins != hist2->nbins)
        iftError("Histograms must have the same size!", "iftAddHists");
    
    result = iftCreateHist(hist1->nbins);
    
    for(i = 0; i < result->nbins; i++)
        result->val[i] = hist1->val[i] + hist2->val[i];
    
    return result;
}


void iftAddHistsInPlace(const iftHist *src, iftHist *dst) {
    if (src->nbins != dst->nbins)
        iftError("Histograms must have the same size! (%d, %d)", "iftAddHistsInPlace", src->nbins, dst->nbins);

    #pragma omp parallel for
    for (int i = 0; i < dst->nbins; i++)
        dst->val[i] += src->val[i];
}

