#include "ift/core/dtypes/SaliencyArray.h"
#include "ift/core/io/NumPy.h"
#include "iftCommon.h"


iftSaliencyArray *iftCreateSaliencyArray(long n) {
    iftSaliencyArray *darr = (iftSaliencyArray *) iftAlloc(1, sizeof(iftSaliencyArray));

    darr->n   = n;
    darr->val = iftAllocFloatArray(n);

    return darr;
}

void iftDestroySaliencyArray(iftSaliencyArray **darr) {
    if (darr != NULL && *darr != NULL) {
        iftSaliencyArray *darr_aux = *darr;

        if (darr_aux->val != NULL)
            iftFree(darr_aux->val);
        iftFree(darr_aux);
        *darr = NULL;
    }
}

iftSaliencyArray *iftCopySaliencyArray(iftSaliencyArray *to_be_copied){
    iftSaliencyArray *copy = iftCreateSaliencyArray(to_be_copied->n);
    for(int i = 0; i < copy->n; i++)
        copy->val[i] = to_be_copied->val[i];

    return copy;
}

iftSaliencyArray *iftSaliencyPriorFromImage(iftImage *label_img, iftImage *saliency_image){
    iftSaliencyArray *region_prior;
    iftIntArray *region_sizes;

    region_sizes = iftCountLabelSpels(label_img);
    int n_superpixels = (int)region_sizes->n -1;
    region_prior = iftCreateSaliencyArray(n_superpixels);

    float max_value = 0, min_value = IFT_INFINITY_FLT;

    /*------- Compute the mean saliency of each region -------*/
    for(int p = 0; p < label_img->n; p++)
        region_prior->val[label_img->val[p] - 1] += (float)saliency_image->val[p];
    for(int s = 0; s < n_superpixels; s++)
        region_prior->val[s] /= (float)region_sizes->val[s + 1];


    /*--------- Normalize it between 0-1 ----------------*/
    for(int s = 0; s < n_superpixels; s++){
        if(region_prior->val[s] > max_value)
            max_value = region_prior->val[s];
        if(region_prior->val[s] < min_value)
            min_value = region_prior->val[s];
    }

    if(max_value == min_value)
        for(int s = 0; s < n_superpixels; s++)
            region_prior->val[s] = 1;
    else
        for(int s = 0; s < n_superpixels; s++)
            region_prior->val[s] = (region_prior->val[s] - min_value) / (max_value - min_value);
    iftDestroyIntArray(&region_sizes);

    region_prior->superpixel_image = label_img;

    return region_prior;
}

iftImage *iftSaliencyImageFromArray(iftSaliencyArray *priors, iftImage *label_img){
    iftImage *saliency_image = iftCreateImageFromImage(label_img);
    for(int p = 0; p < saliency_image->n; p++)
        saliency_image->val[p] = (int)(255 * priors->val[label_img->val[p]-1]);

    return saliency_image;
}