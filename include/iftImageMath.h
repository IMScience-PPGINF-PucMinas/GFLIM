/**
 * @file iftImageMath.h
 * @brief Image Math operations.
 *
 * @note Programs:
 * * 
 */

#ifndef IFT_IMAGEMATH_H
#define IFT_IMAGEMATH_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftFImage.h"


/**
 * @brief Template (function pointer) for a function that computes the subtraction (only the Y channel) between two images: img1 - img2.
 *
 * @param img1 First image.
 * @param img2 Second image.
 * @return Resulting image from the difference img1 - img2.
 *
 * @author Samuka Martins
 * @date Dec 12, 2018
 */
typedef iftImage* (*iftImageSubFunc)(const iftImage *img1, const iftImage *img2);


/** @brief Adds two Images (only the Y channel).*/
//! swig(newobject, stable)
iftImage *iftAdd(const iftImage *img1, const iftImage *img2);

/**
 * @brief Adds src and dst image (only the Y channel), saving the values in place in dst.
 * @param src Source Image.
 * @param dst Destination Image
 *
 * @author Samuka
 * @date Jan 16, 2017
 */
void iftAddInPlace(const iftImage *src, const iftImage *dst);


/**
 * @brief Add a scalar into all voxels of an image (only the Y channel).
 * @author Samuka
 * @date Nov 16, 2017
 */
iftImage *iftAddScalar(const iftImage *src, int scalar);


/**
 * @brief Add (in place) a scalar into all voxels of an image (only the Y channel).
 * @author Samuka
 * @date Nov 16, 2017
 */
void iftAddScalarInPlace(iftImage *dst, int scalar);




/**
 * @brief Computes the voxel-wise difference between two images (img1 - img2).
 * 
 * @note It is only performed on Y channel.
 * @note Note that the resulting image can have negative values.
 *
 * @param img1 First image.
 * @param img2 Second image.
 * @return Resulting image from the difference img1 - img2.
 */
//! swig(newobject, stable)
iftImage *iftSub(const iftImage *img1, const iftImage *img2);


/**
 * @brief Computes the voxel-wise difference between two images max(0, img1 - img2)
 * followed by ReLU (negative subtracted values become 0).
 * 
 * @note It is only performed on Y channel.
 *
 * @param img1 First image.
 * @param img2 Second image.
 * @return Resulting image from the difference max(0, img1 - img2).
 * 
 * @author Samuel Martins
 * @date Nov 26, 2019
 */
//! swig(newobject, stable)
iftImage *iftSubReLU(const iftImage *img1, const iftImage *img2);


/**
 * @brief Computes the absolute subtraction (only the Y channel) between two images.
 *
 * @note If the images are colored, the channels Cb and Cr from image <img1> are copied into output image.
 *
 * @author Samuel Martins
 * @date Aug 27, 2018
 */
//! swig(newobject, stable)
iftImage *iftAbsSub(const iftImage *img1, const iftImage *img2);


/**
 * @brief Computes the Intersection between two Label Images
 * (intersection between objects is true if they have the same label).
 * @param  label_img1 Label Image 1.
 * @param  label_img2 label_img 2.
 * @return            Resulting Intersection Image.
 *
 * @author Samuka Martins
 * @date Aug 10, 2017
 */
//! swig(newobject, stable)
iftImage *iftIntersec(const iftImage *label_img1, const iftImage *label_img2);

/**
 * @brief Applies a logical AND between two Images (only the Y channel).
 *
 * For each pixel from the images, this functions gets the minimum value (Y channel from YCbCr, regardless
 * if the image is gray or colored) from two images.
 */
//! swig(newobject, stable)
iftImage *iftAnd(const iftImage *img1, const iftImage *img2);


/**
 * @brief Applies a logical OR between two Images (only the Y channel).
 *
 * For each pixel from the images, this functions gets the maximum value (Y channel from YCbCr, regardless
 * if the image is gray or colored) from two images.
 */
//! swig(newobject, stable)
iftImage *iftOr(const iftImage *img1, const iftImage *img2);


/** @brief Multiplies two Images (only the Y channel).*/
//! swig(newobject, stable)
iftImage *iftMult(const iftImage *img1, const iftImage *img2);

/**
 * @brief Multiplies an integer image (only the Y Channel) by a float image. The images must have the same domain.
 * @param  img  Input Integer Image
 * @param  fimg Input Float Image.
 * @return      Resulting multiplied Image.
 *
 * @author Samuka Martins
 * @date Dec 26, 2017
 */
iftFImage *iftMultByFImage(const iftImage *img, const iftFImage *fimg);

/**
 * @brief Multiplies an image by an integer scalar.
 * @author Samuka Martins
 * @date Nov 30, 2017
 */
//! swig(newobject, stable)
iftImage *iftMultByIntScalar(const iftImage *img, int scalar);


/**
 * @brief Multiplies an image by a float scalar.
 * @author Samuka Martins
 * @date Nov 30, 2017
 */
iftFImage *iftMultByScalar(const iftImage *img, float scalar);




/** @brief Gets an Image with the absolute values (only the Y channel) from the Input Image.*/
//! swig(newobject, stable)
iftImage *iftAbs(const iftImage *img);


/**
 * @brief Computes the square subtraction (only the Y channel) between two images.
 * 
 * @note If the images are colored, the channels Cb and Cr from image <img1> are copied into output image.
 * 
 * @author Samuel Martins
 * @date Aug 27, 2018
 */
iftImage *iftSquareSub(const iftImage *img1, const iftImage *img2);



/**
 * @brief Apply a xor on img1 and img2 (if img1->val[p] != img2->val[p], return true)
 * @return      Xor Binary Image.
 * 
 * @author Samuel Martins
 * @date Oct 2, 2017
 */
//! swig(newobject, stable)
iftImage *iftXor(const iftImage *img1, const iftImage *img2);


/** @brief Gets the Complement (only the Y channel) from an Image.*/
//! swig(newobject, stable))
iftImage *iftComplement(const iftImage *img);


/**
 * @brief Applies the activation function ReLU (Rectified Linear Units f(x) = max(0, x)) in the
 * image values (channel Y) from <img>.
 * 
 * @note If the image is colored, the channels Cb and Cr are copied.
 * 
 * @author Samuel Martins
 * @date Sep 17, 2018
 */
//! swig(newobject)
iftImage *iftReLUImage(const iftImage *img);


/**
 * @brief Copies an image only where its mask is, i.e., only for the voxels (!= 0) from the mask.
 *
 * @note The Cb and Cr are also copied for Colored Images.
 * 
 * @param  img  Input Image.
 * @param  mask Mask that determines which voxels are copied.
 * @return      Copied Image.
 */
//! swig(newobject, stable)
iftImage *iftMask(const iftImage *img, const iftImage *mask);


/**
 * @brief Copies an image only where its object mask is, i.e., only for the voxels == label.
 *
 * @note The Cb and Cr are also copied for Colored Images.
 * 
 * @param  img   Input Image.
 * @param  mask  Mask with multiple objects, where the voxels with label <label> will be copied.
 * @param  label Label from the object to be masked.
 * @return       Copied Image.
 *
 * @author Samuel Martins
 * @date Jan 24, 2019
 */
iftImage *iftObjectMask(const iftImage *img, const iftImage *mask, int label);


/**
 * @brief Binarize a labeled image.
 * @author Samuka Martins
 * @date Aug 10, 2017
 */
//! swig(newobject, stable)
iftImage *iftBinarize(const iftImage *label_img);


/** @brief Adds a value for all voxel only in the Y channel */
//! swig(newobject, stable)
iftImage *iftAddValue(const iftImage *img, int val);


/** Transforms a Binary Image <b>bin</b> into a Binary Image where each label value has value 255 */
//! swig(newobject, stable)
iftImage *iftBinaryFrom1To255(const iftImage *bin);


/** Transforms a Binary Image <b>bin</b> with object values != 1 into a Binary Image with values 1 */
iftImage *iftBinaryFrom255To1(const iftImage *bin);


/** @brief Adds a Value in an Input Image only in the voxels of a Mask (!= 0) */
iftImage *iftAddValueInRegion(const iftImage *img, int value, const iftImage *mask);


/** @brief Applies a square root into an Image */
iftFImage *iftSQRT(iftImage const *img1);


/** @brief Applies a square root into a Float Image */
iftFImage *iftFSQRT(const iftFImage *img);


/**
 * @brief Linearly combines two images given img0[x]*alpha+(1-alpha)*img1[x]. The images must have the same
 * size and must be normalized between the same minimum and maximum values.
 *
 * @author Thiago Vallin Spina
 * @date Jan 23, 2016
 *
 * @param img0 First input image.
 * @param img1 Second input image.
 *
 * @return Linearly combined image.
 */
iftImage *iftImageLinearCombination(const iftImage *img0, const iftImage *img1, double alpha);


/**
 * @brief Computes the mean intensity in a given region of the image, as
 * provided by a mask. It uses the entire image when the mask is not
 * provided.
 * @author Alexandre Falcao   
 * @date   Nov 10th, 2016
 * @param  img  Input image.  
 * @param  mask Input region (or NULL).  
 * @return Mean intensity.
 */
//! swig(stable)
float iftMeanValue(const iftImage *img, iftImage *mask);


/**
 * @brief Computes the most frequent intensity (mode) in a given
 * region of the image, as provided by a mask. It uses the entire
 * image when the mask is not provided.  
 * @author Alexandre Falcao
 * @date   Nov 10th, 2016
 * @param  img Input image.  
 * @param  mask Input region (or NULL).
 * @return The most frequent intensity.
 */
int iftModeValue(const iftImage *img, iftImage *mask);


/**
 * @brief Computes the Mean Value of each object/region from a image.
 * 
 * @note The Cb and Cr channels are not considered.
 * 
 * @param img Gray Image.
 * @param label_img Label image with the objects/regions to be considered.
 * @return Array of mean values with range [0, n_objects], where n_objects is the label with the highest value.
 *         Thus, the mean of the object i is at index [i].
 * 
 * @author Samuel Martins
 * @date Sep 4, 2018
 */
iftFloatArray *iftMeanValueInRegions(const iftImage *img, const iftImage *label_img);


/**
 * @brief Computes the Mean Value and Standard Deviation of each object/region from a image.
 * 
 * @note The Cb and Cr channels are not considered.
 * 
 * @param img Gray Image.
 * @param label_img Label image with the objects/regions to be considered.
 * @param norm_val Normalization value (2^image_bits - 1) to normalize the resulting metrics to [0, 1].
 *                 Use 1.0 to ignore normalization. 
 * @return Matrix of shape (nrows, ncols) = (n_objects + 1, 2), where n_objects is the label with the highest value.
 *         The row [i] indicates the object i, such [i][0] is the mean value and [i][1] is the stdev of the object i.
 *           
 * 
 * @author Samuel Martins
 * @date Sep 4, 2018
 */
iftMatrix *iftMeanStdevValuesInRegions(const iftImage *img, const iftImage *label_img, float norm_val);


/**
 * @brief Rounds a float image to a integer one (floor values for decimals < 0.5, ceil values otherwise).
 * @param fimg Float Image to be rounded.
 * @return Rounded integer image
 *
 * @author Samuel Martins
 * @date Dec 5, 2018
 */
//! swig(newobject)
iftImage *iftRoundFImage(const iftFImage *fimg);


#ifdef __cplusplus
}
#endif

#endif
