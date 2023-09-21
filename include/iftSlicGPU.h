/**
 * @file iftSlic.h
 * @brief Definitions and functions about SLIC supervoxel technique.
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 *
 * @note Programs:
 * * @ref iftSuperpixelsBySlic.c (demo/Slic/iftSuperpixelsBySlic.c) = Generates the superpixels for a given Image using SLIC.
 */

#ifndef IFT_SLIC_H
#define IFT_SLIC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftImage.h"
#include "iftIGraph.h"
#include "iftMImage.h"
#include "iftSegmentation.h"


/**
 * @brief Initializes the Slic parameters.
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 * 
 * @param  grad          Gradient of the Image to be oversegmented that is used to move the initial
 *                       the centroid of each cluster.
 * @param  mimg          Multi-Band Image that contains the converted (3 bands) LAB Image
 *                       (if the YCbCr Image is colored) or the original Grayscale Image (1 band) 
 * @param  input_n_cells The input (approximate) number of superpixels in the oversegmented output image
 * @param  comp          Compactness: Relative Importance of color similarity and spatial proximity
 * @param  img_max_range Maximum range of the image to be oversegmented: e.g 255 (for 8 bits images)
 * @return               The Slic parameters.
 */
iftSlic *iftInitSlicGPU(const iftImage *grad, const iftMImage *mimg, int input_n_cells, double comp,
                     int img_max_range);


/**
 * @brief Destroys a Slic Parameter structure.
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 */
void iftDestroySlicGPU(iftSlic **slic);


/**
 * @brief Computes the (Slic) Distance between a voxel u and a cell's center.
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 *
 * @note If the Input Image is colored, 3 features (LAB) are considered.
 * @note Otherwise, only the first feature (gray) is considered.
 * @note Likewise, if the Input Image is 3D, the z-axis is used in the distance
 * 
 * @param  center              Center (double) coordinates.
 * @param  u                   Voxel u (double) coordinates.
 * @param  cfeats_center       Color Feature Vector of the cell's Center
 * @param  cfeats_u            Color Feature Vector of the Voxel u
 * @param  step                Side of the initial regular cell = (step) distance between two initial adjacency cell's center
 * @param  compactness         Compactness: Relative Importance of color similarity and spatial proximity
 * @param  spatial_dist_factor Pre-computed Distance factor used to multiply the spatial distance.
 * @param  is_3D_img           Tells if the input image is 3D or not.
 * @param  is_color_img        Tells if the input image is colored or not.
 * @return                     The Slic distance between the voxel u and the cell's center <center>
 */
double iftSlicDistanceGPU(iftPoint center, iftPoint u, iftFColor cfeats_center, iftFColor cfeats_u,
                       int step, double compactness, double spatial_dist_factor,
                       bool is_3D_img, bool is_color_img);


/**
 * @brief Update the centers (Kmeans-like) from all cells by averaging the features vectors of
 * all their voxels.
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 * 
 * @param slic        Slic Parameter structure that contains the cell's feature vectors and where
 *                    the new center are saved.
 * @param superpixels Superpixel Image.
 * @param mimg        Multi-Band Image that contains the converted (3 bands) LAB Image
 *                    (if the YCbCr Image is colored) or the original Grayscale Image (1 band)
 */
// void iftUpdateCellCentersGPU(iftSlic *slic, const iftImage *superpixels, const iftMImage *mimg);


void iftEnforceLabelConnectivityGPU(iftImage *superpixels, iftSlic *slic);



/**
 * @brief Generate Superpixels for a Image using Slic.
 * @authors Samuel Martins, John Vargas
 * @date May 24, 2016
 * 
 * @param  img               Image to be oversegmented.
 * @param  input_grad        Pre-computed Image Gradient from <b>img</b>. If NULL, the gradient image will be computed.
 * @param  input_n_cells     The input (approximate) number of superpixels in the oversegmented output image
 * @param  comp              Compactness: Relative Importance of color similarity and spatial proximity
 * @param  img_max_range     Maximum range of the image to be oversegmented: e.g 255 (for 8 bits images)
 * @param  out_blocks        Returns by Reference (if != NULL) the final number of clusters/cells.
 * @param  out_n_cells       Returns by Reference (if != NULL) the final number of clusters/cells.
 * @return                   The Superpixel Image.
 */
iftImage *iftGenerateSuperpixelsBySlicGPU(const iftImage *img, const iftImage *input_grad, int input_n_cells,
                                       double comp, int img_max_range, int *out_n_cells);


/**
 * @brief Divides the Input Image into blocks, generates superpixels for each block, and then re-enumerates
 * all block's superpixels, getting the final Superpixel Image.
 * @authors Samuel Martins, John Vargas
 * @date Jun 16, 2016
 * 
 * @param  img               Image to be oversegmented.
 * @param  input_grad        Pre-computed Image Gradient from <b>img</b>. If NULL, the gradient image will be computed.
 * @param  input_n_blocks    Expected Number of Blocks to divide the Image.
 * @param  input_n_cells     The input (approximate) number of superpixels in the oversegmented output image
 * @param  comp              Compactness: Relative Importance of color similarity and spatial proximity
 * @param  img_max_range     Maximum range of the image to be oversegmented: e.g 255 (for 8 bits images)
 * @param  out_blocks        Returns by Reference (if != NULL) the coordinates of the generated blocks.
 * @param  out_n_cells       Returns by Reference (if != NULL) the final number of clusters/cells.
 * @return                   The Superpixel Image.
 */
iftImage *iftGenerateSuperpixelsBySlicInBlocksGPU(const iftImage *img, const iftImage *input_grad, int n_blocks,
                                               int input_n_cells, double comp, int img_max_range,
                                               iftImageTiles **out_blocks, int *out_n_cells);


#ifdef __cplusplus
}
#endif

#endif // IFT_SLIC_H