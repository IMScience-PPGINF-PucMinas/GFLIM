#ifndef IFT_MIMAGE_H_
#define IFT_MIMAGE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftFImage.h"
#include "ift/core/dtypes/Color.h"
#include "ift/metriclearn/MetricLearnCore.h"
#include "iftAdjacency.h"
#include "iftMatrix.h"



#define iftMGetXCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define iftMGetYCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define iftMGetZCoord(s,p) ((p) / (((s)->xsize)*((s)->ysize)))
#define iftMGetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])
#define iftMDiagonalSize(s) (iftRound(sqrtf(s->xsize*s->xsize + s->ysize*s->ysize + s->zsize*s->zsize)))

/* Multiband image: iftMImage *img; img->val[p][b]  */

//! swig(extend = iftMImageExt.i, destroyer = iftDestroyMImage)
typedef struct ift_mimage {
    iftMatrix      *data; /* matrix of pixel features */
    float          **val; /* fast access for each pixel feature (data row) ex. mimg->val[pixel][band]*/
    int            xsize,ysize,zsize; /* image dimensions of each band */
    float          dx,dy,dz;  /* voxel size */
    int            *tby, *tbz; /* LUT to speed up index to/from coordinate conversions */
    unsigned long  n,m; /* number of voxels and number of bands */
} iftMImage;

//! swig(extend = iftMImageArrayExt.i, destroyer = iftDestroyMImageArray)
typedef struct ift_mimage_array {
  iftMImage **val;
  int n;
} iftMImageArray;

int         iftMXSize(iftMImage *img);
int         iftMYSize(iftMImage *img);
int         iftMZSize(iftMImage *img);

//! swig(newobject, stable)
iftVoxel    iftMGetVoxelCoord(const iftMImage *img, int p);

//! swig(newobject, stable)
iftMImage  *iftCreateMImage(int xsize,int ysize,int zsize, int nbands);

void        iftDestroyMImage(iftMImage **img);

iftMImage  *iftCopyMImage(iftMImage *img);

iftMImageArray *iftCreateMImageArray(long n);

void iftDestroyMImageArray(iftMImageArray **arr);

iftMImageArray  *iftCopyMImageArray(iftMImageArray *arr);

//! swig(newobject, stable)
char iftMValidVoxel(const iftMImage *img, iftVoxel v);

void        iftMCopyVoxelSize(const iftMImage *img1, iftMImage *img2);
/**
 * @brief Copies the voxel size from an iftImage to an iftMImage
 *
 * @author Thiago Vallin Spina
 *
 * @param img1 An iftImage.
 * @param img2 An iftMImage.
 */
void        iftMCopyVoxelSizeFromImage(const iftImage *img1, iftMImage *img2);
/**
 * @brief Copies the voxel size from an iftMImage to an iftImage
 *
 * @author Thiago Vallin Spina
 *
 * @param img1 An iftMImage.
 * @param img2 An iftImage.
 */
void        iftMCopyVoxelSizeToImage(const iftMImage *img1, iftImage *img2);

/**
 * @brief Converts an image into a multiband image. 
 *
 * @author Alexandre Falcão
 * 
 * @param img: an image.
 * @param color_space: the desired color space. 
 */
//! swig(newobject, stable)
iftMImage  *iftImageToMImage(const iftImage *img, char color_space); /* See options for color_space in iftColor.h */


//! swig(newobject, stable)
  iftImage   *iftMImageToImage(const iftMImage *img, int Imax, int band);


//! swig()
static inline bool iftIs3DMImage(const iftMImage *img) {
    return (img->zsize > 1);
}


//! swig(newobject, stable)
iftMImage   *iftReadMImage(const char *filename);

/**
 * @brief Read a numpy array as a MImage.
 * 
 * @param npy_path Pathname from the numpy array.
 * @return Loaded MImage.
 * 
 * @author Samuka Martins
 * @date Sep 25, 2019 
 */
//! swig(newobject, stable)
iftMImage *iftReadMImageFromNumPy(const char *npy_path, ...);


//! swig(stable)
void  	     iftWriteMImage(iftMImage *mimg, const char *filename);

/**
 * @brief Writes a MImage as a numpy array.
 * 
 * @param mimg MImage. 
 * @param npy_path Numpy array path.
 * 
 * @warning Voxel sizes are not written.
 * 
 * @author Samuka Martins
 * @date Sep 25, 2019
 */
//! swig(stable)
void iftWriteMImageAsNumPy(const iftMImage *mimg, const char *npy_path, ...);

//! swig()
void        iftWriteMImageBands(iftMImage *mimg, char *base_filename);

iftMImage  *iftMAddFrame(iftMImage *img, int bx, int by, int bz, float value);
iftMImage  *iftMRemFrame(iftMImage *fimg, int bx, int by, int bz);
void        iftSetMImage(iftMImage *img, float value);
iftImage   *iftEuclMImageBasins(iftMImage *img, iftAdjRel *A);

//! swig(newobject, stable)
iftImage   *iftMImageBasins(const iftMImage *img, iftAdjRel *A);
//! swig(newobject)
iftImage   *iftMImageGradient(const iftMImage *img, iftAdjRel *A, int Imax);
iftImage   *iftBorderProbImage(iftMImage *img);
iftImage   *iftRegionProbImage(iftMImage *img);
iftImage   *iftUniformProbImage(iftMImage *img);
void        iftMultMImageByScalar(iftMImage *Z, float scalar);

iftMImage  *iftGradientVector(iftImage *img, iftImage *mask, iftAdjRel *A);

  /* Voxel sampling methods that return a binary image with 0/1 value,
     by using the border information to avoid outliers. */

  iftImage *iftSelectNonBorderVoxels(iftMImage *img, iftImage *mask, int nsamples);

/**
 * @brief Does a grid sampling in the multi-image
 * @param img A multi-image
 * @param mask The mask that defines the region in img to be sampled
 * @param nsamples Desired number of samples
 * @return A mask image with the selected samples
 *
 */
//! swig(newobject, stable)
iftImage *iftGridSampling(iftMImage *img, iftImage *mask, int nsamples);

//! swig(newobject, stable)
iftImage *iftAltMixedSampling(iftMImage *img, iftImage *mask, int nsamples);

//! swig(newobject, stable)
iftImage *iftSelectNSamplesFromMask(iftMImage *img, iftImage *mask1, int nsamples);

  // If band is negative then the function will search for the maximum value among all bands
  float iftMMaximumValue(const iftMImage *img, int band);
  // If band is negative then the function will search for the minimum value among all bands
  float iftMMinimumValue(iftMImage *img, int band);

  /**
  * @brief Computes an iftImageTiles record given a multidimensional image and a number of tiles for each image dimension.
  *
  * The dimensions of the image tiles are equally sized to approximately divide the image into regular patches.
  * In practice, this function uses iftComputeBoundingBoxImageTiles with the entire image as a bounding box.
  *
  * @author Thiago V. Spina
  * @date Feb 26, 2015
  *
  * @param mimg    A multidimensional image that must be tiled..
  * @param ntiles_x The number of tiles in the X axis.
  * @param ntiles_y The number of tiles in the Y axis.
  * @param ntiles_z The number of tiles in the Z axis.
  * @return The allocated record.
  *
  * @warning The number of tiles may actually be less than ntiles_x*ntiles_y*ntiles_z since we prioritize the fact that
  * the patches must have approximately the same size. REFER TO tiles->ntiles_x/y/z when accessing the actual number of
  * computed tiles.
  *
  * @sa iftBoundingBoxImageTilesByEquallyDividingAxes
  */
  iftImageTiles *iftMImageTilesByEquallyDividingAxes(const iftMImage *mimg, int ntiles_x, int ntiles_y, int ntiles_z);


/**
 * @brief Extends a MImage with only one band with its Local Binary Pattern (LPB)
 */
iftMImage *iftExtendMImageByLBP(iftMImage *img, iftAdjRel *A, char normalize);

  /**
  * @brief Extends a multi-dimensional image with the color information of adjacent pixels and their own spatial coordinates.
  * @author Adan Echemendia
  * @date Nov 20, 2016
  * @param img The input multi-band image.
  * @param A The adjacency relation.
  * @param normalize_voxel_coord A flag indicating whether or not the voxel coordinates must be normalized before
   * adding them
  * @return A extended multidimensional image.
  *
  */
iftMImage *iftExtendMImageByAdjacencyAndVoxelCoord(const iftMImage *mimg, const iftAdjRel *A, char normalize_voxel_coord);

/**
 * @brief Extends a multi-dimensional image with the color information of adjacent pixels
 *
 * @author Adan Echemendia
 * @date Nov 20, 2016
 *
 * @param img The input multi-band image.
 * @param A The adjacency relation
 * @return An extended multidimensional image.
 *
 */
iftMImage *iftExtendMImageByAdjacency(iftMImage *mimg, iftAdjRel *A);

/**
 * @brief Concatenates an image at the end of an multi-band image
 *
 * @author Azael M Sousa
 * @date Sept 26, 2020
 *
 * @param mimg The input multi-band image.
 * @param img The image to be concatenated.
 * @return An extended multidimensional image.
 *
 */
iftMImage *iftExtendMImageByImage(iftMImage *mimg, iftImage *img);

/**
 * @brief Concatenates a multi-band image at the end of another multi-band image
 *
 * @author Ilan Silva
 * @date May 10, 2022
 *
 * @param mimg1 The input multi-band image.
 * @param mimg2 The multi-band image to be concatenated.
 * @return An extended multidimensional image.
 *
 */
iftMImage *iftExtendMImageByMImage(iftMImage *mimg1, iftMImage *mimg2);

/**
 * @brief Creates a label image where each pixel has the value of the corresponding tile in a <b>n_tiles</b> tile division
 * @author Adan Echemendia
 */
iftImage *iftMImageTilesToLabelImage(iftMImage *mimg, int n_tiles);

/**
 * @brief Computes the squared L2-norm distance between 2 multi-band pixels p and q.
 * @author Jordão Bragantini
 * @date Dez 10, 2018
 */
float iftMImageSqDist(const iftMImage *mimg, int p, int q);

/**
 * @brief Computes the distance between 2 multi-band pixels p and q.
 * @author Samuka Martins
 * @date Oct 26, 2017
 */
float iftMImageDist(const iftMImage *mimg, int p, int q);

/**
 * @brief Computes the depth of a Multiband image.
 * @param  mimg Multiband Image.
 * @return      Image depth
 *
 * @author Cesar Castelo
 * @date Dec 20, 2017
 */
int iftMImageDepth(iftMImage *mimg);

  /**
  * @brief Extends multidimensional image by adding the voxel coordinates. The coordinates can be normalized (normalize_voxel_coord = true).
  *
  * @param mimg      An input multidimensional image to be extended with the voxel coordinates.
  * @param normalize A flag that indicates whether or not the voxel coordinates must be normalized for addition.
  * @return          The extended multidimensional image.
  *
  * @author Alexandre Falcao
  * @date Nov 10, 2016
  *
  * @lastupdate Aug 8, 2018 by Alexandre Falcao
  *
  */
  iftMImage *iftExtendMImageByVoxelCoord(const iftMImage *mimg, const iftImage *label_img, const int n_adj_voxels, bool normalize);


/**
 * @brief Extract image features with the option of adding color of the adjacent voxels and their coordinates.
 * 
 * @param img         Input Image.
 * @param label_img   (Optional) Label Image that indicates the voxels (values != 0) that will be
 *                    considered in the feature extraction. If NULL, all voxels are considered.
 * @param A           (Optional) Adjacency relation that defines which neighbors from the voxel will be considered for feature extraction.
 *                    If NULL, only the own voxel is considered.
 * @param voxel_coord True indicates the voxel coordinates should be added as features and false indicates otherwise.
 * @return            It returns a multiband image with the extracte image features.
 * 
 * @date July 15th, 2018
 * @author Alexandre Falcao
 * 
 * @lastupdate Jul 17, 2018 by Samuel Martins
 */
//! swig(newobject, stable)
iftMImage *iftExtractImageFeatures(const iftImage *img, const iftImage *label_img, const iftAdjRel *A, bool voxel_coord);

iftMImage *iftTransformMImage(const iftMImage *input, const double *L, int L_nrow, int L_ncol);

//! swig(newobject)
iftMImage *iftKernelTransformMImage(const iftMImage *input, const iftDoubleMatrix *train_data,
                                    const iftDoubleMatrix *omega, iftKernelFunction *K,
                                    double param1, double param2, double scale);


//! swig(newobject)
iftMImage *iftColorQuantization(const iftMImage *input, float lower_bound, float upper_bound, int bins_per_band);

/**
 * @brief Converts an mimage between color spaces
 * @param img Input mimage
 * @param colSpaceIn Input color space
 * @param colSpaceOut Output color space
 *
 * @author Cesar Castelo
 * @date Mar 20, 2019
 * @ingroup MImage
 */
iftMImage *iftMConvertColorSpace(iftMImage* mimg, iftColorSpace colSpaceIn, iftColorSpace colSpaceOut);

//! swig(newobject)
iftMImage *iftMImageCentralize(const iftMImage *mimg);

//! swig(newobject)
iftMImage *iftMGaussianFilter(const iftMImage *mimg, float radius, float stddev);

//! swig(newobject)
iftMImage *iftBinarySetsNCA(const iftMImage *mimg, const iftSet *obj, const iftSet *bkg,
                            int d_out, int n_iters, double learn_rate);

/**
 * @brief Rescale features to unitary scale, notice that location of features is changed
 * @param mimg  Input mimage
 *
 * @author Jordao Bragantini
 * @date Sep, 13, 2019
 * @ingroup MImage
 */
//! swig()
void iftRescaleMImage(iftMImage *mimg, bool single_max);


/**
 * @brief Stack a set of images to build a multi-band image.
 * Each image is a band.
 * 
 * @param n Number of images to be stacked.
 * @param ... Images to be stacked.
 * @return iftMImage* Multi-band image.
 * 
 * @author Samuel Martins
 * @date Nov 20, 2019
 */
//! swig(newobject)
iftMImage *iftStackGrayImages(int n, ...);

//! swig(newobject, stable)
int iftMGetVoxelIndex_pyift(iftMImage *s, iftVoxel v); 

#ifdef __cplusplus
}
#endif

#endif
