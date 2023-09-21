/**
 * @file iftIGraph.h
 * @brief Image Graph - data structure and functions.
 * 
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 * 
 */

#ifndef IFT_IGRAPH_H_
#define IFT_IGRAPH_H_

#ifdef __cplusplus
extern "C" {
#endif

  /* TODO LIST: iftRegionIGraph (graph of superpixels from label
     image), iftSurfaceIGraph (graph of 2D/3D object surface), and all
     IFTs */ 

#include "ift/core/dtypes/DHeap.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "iftAdjacency.h"
#include "iftCommon.h"
#include "iftCompTree.h"
#include "iftImage.h"
#include "iftGraphics.h"
#include "iftMImage.h"
#include "iftSegmentation.h"

#define   COMPLETE 0 /* graph where all nodes are adjacent to each other */  
#define   EXPLICIT 1 /* graph with adjacency list of the nodes */
#define   IMPLICIT 2 /* graph with translation-invariant adjacency relation */

/* Criteria for iftWaterGray*/
#define   HEIGHT 0
#define   AREA   1
#define   VOLUME 2

/**
 * @brief A struct that represents a node in the ImageGraph.
 * 
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
typedef struct ift_inode {
  int     voxel;  /* voxel (or supervoxel representative) */
  float   weight; /* the weight of the node */
  iftSet *adj;    /* index list to adjacent nodes (explicit graphs only) */
} iftINode;


/**
 * @brief A struct that represents the ImageGraph containing a list of nodes and its features (three types of graphs can be defined 'COMPLETE', 'IMPLICIT' and 'EXPLICIT')
 * 
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(destroyer = iftDestroyIGraph)
typedef struct ift_igraph {
  iftINode   *node;    /* list of graph nodes */
  int         nnodes;  /* number of graph nodes */
  int         nfeats;  /* number of image features */
  iftImage   *index;   /* node index */
  float     **feat;    /* image features */
  int        *label, *marker, *root, *pred; /* forest annotation */
  double      *pvalue;  /* forest annotation */
  iftAdjRel  *A;       /* adjacency relation (implicit graphs only) */
  char        type;    /* COMPLETE, IMPLICIT, EXPLICIT */
} iftIGraph;


/**
 * @brief Create a graph with nodes indicated by mask values greater than 0.
 * 
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  iftIGraph *iftCompleteIGraph(iftMImage *img, iftImage *mask); 

/**
 * @brief Create a graph with nodes indicated by mask values greater than 0. 
 * The arcs are defined by the adjacency relation A.
 * 
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */

//! swig(newobject, stable)
  iftIGraph *iftImplicitIGraph(iftMImage *img, const iftImage *mask, iftAdjRel *A);


/**
 * @brief Create a graph with nodes indicated by mask values greater than 0. 
 * The arcs are defined by the adjacency relation A and stored with the nodes.
 * 
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(newobject, stable)
  iftIGraph *iftExplicitIGraph(const iftMImage *img, const iftImage *mask, const iftImage *label, iftAdjRel *A);


/**
 * @brief Create a graph with nodes defined by mask and arcs defined by the k 
 * closest neighbors within adjacency A.
 * 
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(newobject, stable)
  iftIGraph *iftSpatialKnnIGraph(iftMImage *img, iftImage *mask, iftAdjRel *A, int K);

/**
 * @brief Create a graph with nodes defined by mask and arcs defined by the threshold df.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(newobject, stable)
iftIGraph *iftSpatialIGraph(iftMImage *img, iftImage *mask, iftAdjRel *A, float df);

/**
 * @brief Create a graph with nodes indicated by mask and arcs defined by the k closest neighbors.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(newobject, stable)
  iftIGraph *iftKnnIGraph(iftMImage *img, iftImage *mask, int K);

/**
 * @brief Destroy image graph from memory.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */ 
  void       iftDestroyIGraph(iftIGraph **igraph);

/**
 * @brief Estimate domes of the pdf.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void       iftIGraphDomes(iftIGraph *igraph);

/**
 * @brief Complement of the domes.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void       iftIGraphBasins(iftIGraph *igraph);  

/**
 * @brief Get path value image.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(newobject, stable)
  iftImage *iftIGraphPathValue(iftIGraph *igraph);


/**
 * @brief Get image with node weights.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(newobject, stable)
  iftFImage *iftIGraphWeight(iftIGraph *igraph);

/**
 * @brief Get label image.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(newobject, stable)
  iftImage *iftIGraphLabel(iftIGraph *igraph);

/**
 * @brief Get root image.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(newobject, stable)
  iftImage *iftIGraphRoot(iftIGraph *igraph);
  
  /*
    old function : void iftIGraphClusters(iftIGraph *igraph, char label_propagation); 
    Compute clusters among the graph nodes and, for explicit graphs,
    propagate labels to the remaining voxels when
    label_propagation=1 */


/**
 * @brief Compute superpixels by using the ISF algorithm.
 * using mean color of the superpixel as feature vector for the root
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(stable)
  int iftIGraphISF_Mean(iftIGraph *igraph, iftImage *seeds, double alpha, double beta0, int niters);

/**
 * @brief Compute superpixels by using the ISF algorithm.
 * using color of the root pixel in the superpixel as feature vector for the root
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
//! swig(stable)
  int iftIGraphISF_Root(iftIGraph *igraph, iftImage *seeds, double alpha, double beta0, int niters);

/**
 * @brief Compute superpixels by using the OISF algorithm.
 * 
 * @details This function computes the superpixels using the Object-based iterative spanning forest segmentation
 *          method. The function computes taking into account a grayscale probability map <b>AS THE LAST FEATURE BAND</b>.
 *          Thus, by selecting the gamma value, is possible to force the superpixels to better resemble the map.
 * @author Felipe Belem
 * @date Dec 1, 2018
 * @ingroup ImageGraph
 *
 * @param[in] igraph - Image Graph
 * @param[in] seeds - Image with seeds depicted with values different than zero
 * @param[in] alpha - Compactness factor (of superpixels)
 * @param[in] beta - Superpixel's boundary adherence
 * @param[in] gamma - Probability map confidence factor
 * @param[in] iters - Number of iterations for the segmentation
 * @return Segmented image represented by an IGraph
 * 
 * @warning This method was not tested with explicit graphs!
 * @warning This method is not applicable (yet) for complete graphs!
 */
iftIGraph *iftIGraphOISF(iftIGraph *igraph, iftImage *seeds, double alpha, double beta, double gamma, int iters);

/**
 * @brief Generate Superpixels for a Image using ISF as IGraph.
 * @author Samuel Martins
 * @date Jun 17, 2016
 * 
 * @param  img            Image to be oversegmented.
 * @param  input_n_cells  The input (approximate) number of superpixels in the oversegmented output image.
 * @param  alpha          Alpha factor.
 * @param  beta           Beta factor.
 * @param  n_iters        Maximum Number of Iterations.
 * @param  smooth_n_iters Number of Iterations for Smoothness. If 0, nothing is done.
 * @return                The Superpixel Image Graph (IGraph).
 */
iftIGraph *iftGenerateIGraphSuperpixelsByISF(const iftImage *img, int input_n_cells, float alpha,
                                             float beta, int n_iters, int smooth_n_iters);


/**
 * @brief Generate Superpixels for a Image using ISF.
 * @author Samuel Martins
 * @date Jun 17, 2016
 * 
 * @param  img            Image to be oversegmented.
 * @param  input_n_cells  The input (approximate) number of superpixels in the oversegmented output image.
 * @param  alpha          Alpha factor.
 * @param  beta           Beta factor.
 * @param  n_iters        Maximum Number of Iterations.
 * @param  smooth_n_iters Number of Iterations for Smoothness. If 0, nothing is done.
 * @param  out_igraph     Returns by Reference (if != NULL) the image graph used with the supervoxels.
 * @return                The Superpixel Image.
 */
iftImage *iftGenerateSuperpixelsByISF(const iftImage *img, int input_n_cells, float alpha, float beta,
                                     int n_iters, int smooth_n_iters, iftIGraph **out_igraph);


/**
 * @brief Recompute seeds if their distances to the old seeds are
 * greater than computed thresholds in the feature and coordinate
 * spaces.
 * @author Alexandre Falcao
 * @date June 6th, 2016
 * @ingroup ImageGraph
 *
 */
void iftIGraphISFRecomputeSeeds(iftIGraph *igraph, int *seed, int nseeds, iftSet **trees_for_removal, iftSet **new_seeds, float *new_seeds_flag);

/**
 * @brief Recomputes ISF seeds using spatial information.
 * @author John Vargas-Munoz
 * @ingroup ImageGraph
 *
 * Exposed as public function by Thiago Vallin Spina.
 */
void iftIGraphISFRecomputeSeedsUsingSpatialInformation(iftIGraph *igraph, int *seed, int nseeds, iftSet **trees_for_removal,
                                                       iftSet **new_seeds, float *new_seeds_flag, float **seed_features);

/**
 * @brief Normalize feature vector.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftNormIGraphFeatures(iftIGraph *igraph);

/**
 * @brief Create a max tree from a graph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  iftCompTree *iftIGraphCreateMaxTree(iftIGraph *igraph);

/**
 * @brief Create a min tree from a graph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  iftCompTree *iftIGraphCreateMinTree(iftIGraph *igraph);

/**
 * @brief Return an igraph that is identical to the input igtraph, except for the 
 * features, which are simply allocated (using the nFeatures parameter).
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  iftIGraph * iftIGraphResetFeatureSpace(iftIGraph *igraph, int nFeatures);

/**
 * @brief Return the maximum value of the feature indexed by "feature" in the 
 * input igraph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  float iftIGraphMaximumFeatureValue (iftIGraph * igraph, int feature);

/**
 * @brief Return the maximum value of the weight indexed by in
 * the input igraph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  float iftIGraphMaximumWeight (iftIGraph * igraph);

/**
 * @brief Return the minimum value of the weight indexed by in
 * the input igraph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  float iftIGraphMinimumWeight (iftIGraph * igraph);

/**
 * @brief Propagate forest attributes by using the maximum arc weight between
 * adjacent nodes as clustering criterion and the path values of the
 * nodes as priority function to resolve ties.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphClusterVoxelsByMaxArcWeight(iftIGraph *igraph, uchar pvalue_order);


/**
 * @brief Computes the optimum-path forest by superior reconstruction with
 * root, predecessor, and label propagation. The marker of this
 * watershed transform from grayscale marker is created by adding
 * D+1 to the weight.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphWaterGrayByDepth(iftIGraph *igraph, int D);

/**
 * @brief Computes the optimum-path forest by superior reconstruction with
 * root, predecessor, and label propagation. The marker of this
 * watershed transform from grayscale marker is created by subtracting
 * H+1 from the weight.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphDualWaterGrayByHeight(iftIGraph *igraph, int H);

/**
 * @brief Computes the optimum-path forest by superior reconstruction with
 * root, predecessor, and label propagation. The marker of this
 * watershed transform from grayscale marker is created by using 
 * one of the next criterions: HEIGHT, AREA or VOLUME.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphDualWaterGray(iftIGraph *igraph, int criterion, int thres);

/**
 * @brief Computes the optimum-path forest by superior reconstruction with
 * root, predecessor, and label propagation. The marker of this
 * watershed transform from grayscale marker is created by using 
 * one of the next criterions : HEIGHT, AREA or VOLUME.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphWaterGray(iftIGraph *igraph, int criterion, int thres);


/**
 * @brief Computes the minima of the node weights.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  iftImage *iftIGraphWeightMinima(iftIGraph *igraph);

/**
 * @brief Computes the maxima of the node weights.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  iftImage *iftIGraphWeightMaxima(iftIGraph *igraph);

/**
 * @brief Copies the feature indicated by feature_index to 
 * the weight of the graph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */  
  void iftIGraphCopyFeatureToWeight(iftIGraph *igraph, int feature_index);

/**
 * @brief Compute the maximum arc weight in the graph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  float iftIGraphMaximumFeatureDist(iftIGraph *igraph);


/**
 * @brief Return the root of the forest by applying path compression on the
 * root map.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  int iftIGraphRootVoxel(iftIGraph *igraph, int q);

/**
 * @brief Enumerate the roots of the forest from 1 to N and propagate its
 * labels to the remaining voxels.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  int iftIGraphEnumerateRootsAndPropagateLabels(iftIGraph *igraph); 

/**
 * @brief Copy an image graph
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  iftIGraph *iftCopyIGraph(iftIGraph *igraph);


/**
 * @brief Compute a minimum spanning tree of an image graph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  iftIGraph *iftIGraphMST(iftIGraph *igraph);

/**
 * @brief Initializes data for performing differential Watershed.
 * @author Thiago Vallin Spina
 * @date Oct 28, 2016
 * @ingroup ImageGraph
 *
 * @param igraph Input iftIGraph.
 * @param pvalue Input array that will store path values.
 * @return The allocated priority queue.
 */
iftDHeap *iftIGraphInitDiffWatershed(iftIGraph *igraph, double *pvalue);

/**
 * @brief Resets data for performing differential Watershed.
 * @author Thiago Vallin Spina
 * @date Oct 28, 2016
 * @ingroup ImageGraph
 *
 * @param igraph Input iftIGraph.
 * @param pvalue Priority queue used during watershed.
 */
void iftIGraphResetWatershed(iftIGraph *igraph, iftDHeap *Q);

/**
 * @brief Performs the differential Watershed after initilization. @sa iftIGraphInitDiffWatershed
 * @author Thiago Vallin Spina
 * @date Oct 28, 2016
 * @ingroup ImageGraph
 *
 * @param igraph Input iftIGraph.
 * @param seeds The (new) seeds to be used for delineation.
 * @param trees_for_removal The set of tree roots to be removed.
 * @param Q The priority queue to be used.
 */
void iftIGraphDiffWatershed(iftIGraph *igraph, iftLabeledSet *seeds, iftSet *trees_for_removal, iftDHeap *Q);

/**
 * @brief Performs the differential Watershed after initilization. @sa iftIGraphInitDiffWatershed
 * @author Thiago Vallin Spina
 * @date Oct 28, 2016
 * @ingroup ImageGraph
 *
 * @param igraph Input iftIGraph.
 * @param seeds The (new) seeds to be used for delineation.
 * @param trees_for_removal The set of tree roots to be removed.
 * @param Q The priority queue to be used.
 */

void iftIGraphDiffOrientedWatershed(iftIGraph *igraph, iftImage *orien, iftLabeledSet *seeds, iftSet *trees_for_removal,
                                    iftDHeap *Q, char orientation);
/**
 * @brief Set weight of the nodes in the graph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphSetWeight(iftIGraph *igraph, iftImage *weight);

/**
 * @brief Set weight of the nodes in the graph for FImage input.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphSetFWeight(iftIGraph *igraph, iftFImage *weight);

/**
 * @brief Set weight of the nodes to perform regions smoothing.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphSetWeightForRegionSmoothing(iftIGraph *igraph, const iftImage *img);


/**
 * @brief Smooth regions in label map and fix label and root inconsistencies.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphSmoothRegions(iftIGraph *igraph, int num_smooth_iterations);

/**
 * @brief Fix label and root inconsistencies.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  void iftIGraphFixLabelRootMap(iftIGraph *igraph, iftSet **T);

/**
 * @brief Tree removal function to implement the differential version of ISF.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
  iftSet    *iftIGraphTreeRemoval(iftIGraph *igraph, iftSet **trees_for_removal, double *pvalue, double INITIAL_PATH_VALUE);


/**
 * @brief Re-enumerates (in place) a Superpixel Image, divided in blocks, where each block b has its own superpixels
 * enumerated from 1 to nb. Re-enumerates to 1..n.
 * @author Samuel Martins
 * @date Jun 17, 2016
 * 
 * @param super_img    Superpixel Image to be re-enumerated.
 * @param blocks       Coordinates from the blocks.
 * @param  out_n_cells Returns by Reference (if != NULL) the final number of clusters/cells.
 */
void iftRenumerateBlockSuperpixels(iftImage *super_img, const iftImageTiles *blocks, int *out_n_cells);


/**
 * @brief Writes the image <b>img</b> with the border from the superpixels <b>super_img</b>.
 * @author Samuel Martins
 * @date Jun 17, 2016
 * 
 * @note If <b>img</b> is 3D, it saves a zip file with all XY slices (axial) with the supervoxel's borders drawn,
 * and a 3D image, also inside in the zip, with the supervoxel's borders.
 * 
 * @param img       Image where the superpixel's borders will be drawn.
 * @param super_img Superpixel Image.
 * @param out_path  Output Pathname: *.[ppm,png] for 2D images, *.zip for 3D.
 */
void iftWriteSuperpixelBorders(const iftImage *img, const iftImage *super_img, const char *out_path);

//! swig(newobject)
iftImage *iftExtract_ISF_MIX_ROOT_Superpixels(iftImage *img, int nsuperpixels, float alpha, float beta, int niters, int smooth_niters);

//! swig(newobject)
iftImage *iftExtract_ISF_MIX_MEAN_Superpixels(iftImage *img, int nsuperpixels, float alpha, float beta, int niters, int smooth_niters);


/**
 * @brief Returns iftSet <b>S</b> containing superpixel center voxels from output igraph.
 * @author Gabriel Seabra
 * @date Dec 24, 2021
 * @ingroup ImageGraph
 * 
 * @param igraph       IGraph after running superpixel algorithm.
 */

//! swig(newobject)
iftSet *iftSuperpixelCenterSetFromIGraph(iftIGraph *igraph);

#ifdef __cplusplus
}
#endif

#endif

