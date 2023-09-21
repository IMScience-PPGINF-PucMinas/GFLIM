//
// Created by Samuel Martins on 10/12/18.
//
#include "iftImage.h"

#ifndef IFT_SALIENCYARRAY_H
#define IFT_SALIENCYARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ift_sal_array{
    float *val; //The saliency values for each superpixel
    long n; //The number of superpixels
    iftImage *superpixel_image; //The superpixel map used to compute the saliency values
}iftSaliencyArray;



/**
 * @brief Creates an iftSaliencyArray.
 * @author Leonardo de Melo
 * @date Feb 10, 2022
 * @ingroup Memory
 */
iftSaliencyArray *iftCreateSaliencyArray(long n);


/**
 * @brief Destroys an iftSaliencyArray.
 * @author Leonardo de Melo
 * @date Feb 10, 2022
 * @ingroup Memory
 */
void iftDestroySaliencyArray(iftSaliencyArray **darr);

/**
 * @brief Create a saliency array based on the mean saliency of a each label
 * @date November 05, 2019
 * @author Leo Melo
 * @note This works for both saliency and prior maps
 * @note The arrays starts at 0 and the label image starts at 1. Deduct the label_img val in 1 when finding its corresponding saliency
 *
 * @param  label_img Superpixel labeled image.
 * @param  saliency_image Pixel-wise saliency map.
 * @return      Mean saliency array for each label.
 *
 * @warning The saliency values returned are in the range [0-1].
 */
iftSaliencyArray *iftSaliencyPriorFromImage(iftImage *label_img, iftImage *saliency_image);

/**
 * @brief Create a saliency map based on a prior array and a label image
 * @date November 05, 2019
 * @author Leo Melo
 * @note This works for both saliency and prior maps
 * @note The arrays starts at 0 and the label image starts at 1. Deduct the label_img val in 1 when finding its corresponding saliency
 *
 * @param  priors Prior array where each element relates to the saliency of a label.
 * @param  label_img Superpixel labeled image.
 * @return      Mean saliency array for each label.
 *
 * @warning The saliency values returned are in the range [0-1].
 */
iftImage *iftSaliencyImageFromArray(iftSaliencyArray *priors, iftImage *label_img);

iftSaliencyArray *iftCopySaliencyArray(iftSaliencyArray *to_be_copied);


#ifdef __cplusplus
}
#endif

#endif //IFT_INTARRAY_H
