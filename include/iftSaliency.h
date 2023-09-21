
#ifndef IFT_Saliency_H_
#define IFT_Saliency_H_


#ifdef __cplusplus
extern "C" {
#endif

#include "iftDataSet.h"
#include "iftImage.h"
#include "iftIGraph.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/SaliencyArray.h"

#define IFT_NO_QUERY 0
#define IFT_TOP_BORDER_QUERY 1
#define IFT_RIGHT_BORDER_QUERY 2
#define IFT_BOTTOM_BORDER_QUERY 3
#define IFT_LEFT_BORDER_QUERY 4
#define IFT_ALL_BORDERS_QUERY 5
#define IFT_SALIENCY_MAP_QUERY 6
#define IFT_SALIENCY_MAP_FOREGROUND_QUERY 6
#define IFT_SALIENCY_MAP_BACKGROUND_QUERY 7
#define IFT_SEED_QUERY 8
#define IFT_SEED_FOREGROUND_QUERY 8
#define IFT_SEED_BACKGROUND_QUERY 9
#define IFT_SALIENCY_MAP_JOINT_QUERY 10


#define ITSELF_BY_SALIENCY 1
#define ITSELF_BY_SCRIBBLES 2

#define SIGMA_SQUARED 0.4

#define square(x) ((x)*(x))

typedef struct ift_itself_priors{
    float center_variance;
    float color_variance;
    float white_variance;
    float red_yellow_variance;
    float focus_variance;
    float ellipse_variance;
    float foreground_distance_variance;
    float foreground_feature_variance;
}iftITSELFPriors;

/**
* @brief Array of float values.
* @author Thiago Vallin Spina
* @date Mar 4, 2016
* @ingroup Memory
*/
//! swig(destroyer = iftDestroyITSELFParameters)
typedef struct ift_itself_parameters{
    int itself_iterations;
    int number_superpixels;
    int obj_seeds_number;
    float oisf_gamma;
    int oisf_iterations;
    float oisf_beta;
    float oisf_alpha;
    float oisf_dist_penalty;
    float superpixel_increase_ratio;
    float query_importance;
    float normalization_value;
    float integration_lambda;
    int integration_iteration;
    int enhancement_type;
    iftITSELFPriors prior_params;
}iftITSELFParameters;

iftITSELFParameters *iftInitializeITSELFParametersByDefault();

/**
 * @brief Initializes the parameters for SESS with optimized values for U²Net.
 * @author Leonardo de Melo;
 * @date Feb 15, 2022
 */
//! swig(newobject)
iftITSELFParameters *iftInitializeITSELFParametersByDefaultU2Net();
//! swig(newobject)
iftITSELFParameters *iftInitializeITSELFParametersByDefaultScribbles();

//! swig(newobject)
void iftITSELFParametersSetParam(iftITSELFParameters *params, char *parameter_name, float value);

/**
 * @brief Creates an iftITSELFParameters structure.
 * @author Leonardo de Melo
 * @date Feb 15, 2022
 * @ingroup Memory
 */
iftITSELFParameters *iftCreateITSELFParameters();

/**
 * @brief Destroys an iftITSELFParameters structure.
 * @author Leonardo de Melo
 * @date Feb 15, 2022
 * @ingroup Memory
 */
void iftDestroyITSELFParameters(iftITSELFParameters **darr);
void iftCopyITSELFParameters(iftITSELFParameters *copied_to, iftITSELFParameters *to_copy);
void iftWriteITSELFParamsToFile(iftITSELFParameters *params, char *out_file);

typedef struct ift_saliency_graph {
    iftIntArray *region_sizes;
    iftFloatArray *region_sizes_prob;
    iftFloatArray **adjacency; // The query of each superpixel can be the same (e.g. border query) or different (e.g. adjacency)
    iftFloatArray **superpixel_mean_color;
    iftFloatArray **superpixel_mean_feat;
    int feature_number;
    int query_threshold;
    iftSaliencyArray *saliency; //The computed saliency for each superpixel
    iftIntArray *query; // a binary array stating whether a region (index) is part of the selected query
    int query_size; // a binary array stating whether a region (index) is part of the selected query
    iftSaliencyArray *combined_priors;
    iftSaliencyArray **prior_list;
    int npriors;
    float query_weight;
    int n; // Number of superpixels
    int nfeats; // Number of features on each superpixel
    int improved;
} iftSaliencyGraph;

/**
 * @brief Creates a saliency map by using a saliency map to identify query regions based on foreground
 * @date December 09, 2019
 * @author Leo Melo
 * @note The "saliency_graph_ref" argument can be NULL and the graph will be computed on-the-fly.
 * @note Referencing a pre-computed graph will contribute to reduce memory alloc/free operations.
 * @note To generate a color quantized image dataset, computeQuantizedDatasetKMeans can be used.
 * @note Multiple prior maps can be combined using iftCuboidSaliencyIntegration
 * @note Any pre-computed saliency map can be used for this method, therefore, it can be used as a saliency enhancer.
 * @note The arrays starts at 0 and the label image starts at 1. Deduct the label_img val in 1 when finding its corresponding saliency
 *
 * @param  saliency_graph_ref A pre-computed saliency graph.
 * @param  image_dataset An iftDataset with quantized colors.
 * @param  label_img A superpixel labeled image.
 * @param  adjacency_radius The radius used when defining the adjacency of regions
 * @param  query_importance The weight to balance the importance between query and adjacent regions
 * @param  query_image Pre-computed saliency image from which the foreground queries will derive.
 * @param  normalization_value Defines how fuzzy the map will be. A higher value will imply in a more binary-like saliency map.
 *
 * @return      The estimated saliency map.
 *
 */
iftSaliencyArray *iftGBSSingleForegroundMap(iftSaliencyGraph **saliency_graph_ref, iftMImage *lab, iftImage *label_img, float query_importance, iftImage *query_image, float normalization_value);


/**
 * @brief Creates a saliency map by using a saliency map to identify query regions based on background
 * @date December 09, 2019
 * @author Leo Melo
 * @note The "saliency_graph_ref" argument can be NULL and the graph will be computed on-the-fly.
 * @note Referencing a pre-computed graph will contribute to reduce memory alloc/free operations.
 * @note To generate a color quantized image dataset, computeQuantizedDatasetKMeans can be used.
 * @note Multiple prior maps can be combined using iftCuboidSaliencyIntegration
 * @note Any pre-computed saliency map can be used for this method, therefore, it can be used as a saliency enhancer.
 * @note The arrays starts at 0 and the label image starts at 1. Deduct the label_img val in 1 when finding its corresponding saliency
 *
 * @param  saliency_graph_ref A pre-computed saliency graph.
 * @param  image_dataset An iftDataset with quantized colors.
 * @param  label_img A superpixel labeled image.
 * @param  adjacency_radius The radius used when defining the adjacency of regions
 * @param  query_importance The weight to balance the importance between query and adjacent regions
 * @param  query_image Pre-computed saliency image from which the foreground queries will derive.
 * @param  normalization_value Defines how fuzzy the map will be. A higher value will imply in a more binary-like saliency map.
 *
 * @return      The estimated saliency map.
 *
 */
iftSaliencyArray *iftGBSSingleBackgroundMap(iftSaliencyGraph **saliency_graph_ref, iftMImage *lab, iftImage *label_img, float query_importance, iftImage *query_image, float normalization_value);

/**
 * @brief Creates or recomputes a saliency graph with all its initial attributes
 * @date December 02, 2020
 * @author Leo Melo
 * @note To generate a color quantized image dataset, computeQuantizedDatasetKMeans can be used
 * @note The arrays starts at 0 and the label image starts at 1. Deduct the label_img val in 1 when finding its corresponding saliency
 * @note The result is returned in the prev_saliency_graph parameter.
 *
 * @param  saliency_graph_ref The saliency graph to store the function's result. If you have one computed, pass it as reference to speed up the computation.
 * @param  image_dataset An iftDataset with quantized colors.
 * @param  label_img A superpixel labeled image.
 * @param  adjacency_radius The radius used when defining the adjacency of regions
 *
 * @warning This function do not outputs any saliency score. Check iftComputeSuperpixelGraphDissimilarity for a way to do so.
 */
void iftCreateOrUpdateSaliencyGraph(iftSaliencyGraph **prev_saliency_graph, iftMImage *lab, iftImage *label_img, float adjacency_radius, iftMImage *features);

/**
 * @brief Update which colors are adjacent to each other on an updated saliency graph
 * @date January 14, 2020
 * @author Leo Melo
 * @note This method is required if there is any alteration on the graph's superpixels or color matrix.
 *
 * @param saliency_graph The saliency graph to be updated
 *
 * @return      void.
 */
void iftSaliencyGraphUpdateAdjacentColors(iftSaliencyGraph **saliency_graph);

iftSaliencyArray *iftComputeSuperpixelGraphDissimilarityNew(iftSaliencyGraph *saliency_graph, float query_importance, iftSaliencyArray *old_saliency, float variance);

/**
 * @brief Auxiliary function called by iftComputeSuperpixelGraphDissimilarity when query_importance = 1 to reduce the number of computations.
 * @date December 09, 2019
 * @author Leo Melo
 * @note The saliency graph does not need a precomputed saliency score
 * @note The "query" argument can be set to NULL if saliency_graph->query contains a valid query array or when there are no query regions.
 * @note The arrays starts at 0 and the label image starts at 1. Deduct the label_img val in 1 when finding its corresponding saliency
 *
 * @param  saliency_graph The input saliency graph.
 * @param  query An array with 0s and 1s that define whether a superpixel is a query (1) or not (0)
 * @param  variance Defines how fuzzy the map will be. A higher value will imply in a more binary-like saliency map. It is the variance of the gaussian
 *
 * @return      The saliency score for each superpixel.
 *
 * @warning If you are calling this function multiple times in parallel, be sure to explicitly use the query argument to avoid problems
 */
iftSaliencyArray *iftComputeSuperpixelGraphDissimilarityQueryOnly(iftSaliencyGraph *saliency_graph, float query_importance, iftIntArray *query, float variance);

/**
 * @brief Sets the query of a saliency graph based on a selected query type.
 * @date December 09, 2019
 * @author Leo Melo
 * @note Check the top of this document to see all available query types (e.g. IFT_ALL_BORDERS_QUERY)
 * @note The result will be in the "query" attribute of the saliency graph. Regions will have value 1 if they were set as query and 0 otherwise
 * @note The arrays starts at 0 and the label image starts at 1. Deduct the label_img val in 1 when finding its corresponding saliency
 *
 * @param  label_img A superpixel labeled image.
 * @param  saliency_graph The input saliency graph.
 * @param  query_type The query type option (see note).
 * @param  border_radius The radius required for border queries.1
 * @param  query_image Query image required for: Saliency and Seed-based queries. It should be NULL otherwise.
 *
 * @warning Check whether or not your saliency type requires a query_image.
 */
void iftSetSaliencyQuery(iftImage *label_img, iftSaliencyGraph *saliency_graph, int query_type, int border_radius, iftImage *query_image);

/**
 * @brief Auxiliary function for iftDestroySaliencyGraph. Destroy the superpixel arrays of a saliency graph (avoid gcc warnings).
 * @date November 05, 2019
 * @author Leo Melo
 *
 * @param  saliency_graph The reference saliency graph to have its superpixel arrays destroyed.
 * @return      void.
 *
 */
void iftDestroySaliencyGraphSuperpixels(iftSaliencyGraph **saliency_graph);

/**
 * @brief Destroy all structs contained in a saliency graph.
 * @date November 05, 2019
 * @author Leo Melo
 *
 * @param  saliency_graph The reference saliency graph to be destroyed.
 * @return      void.
 *
 */
void iftDestroySaliencyGraph(iftSaliencyGraph **saliency_graph);

/**
 * @brief Compute the Euclidean Distance of two colors.
 * @date November 05, 2019
 * @author Leo Melo
 *
 * @param  saliency_graph The input saliency graph.
 * @return      The maximum distance between adjacent regions.
 */
double computeColorDistance(iftFColor color1, iftFColor color2);


double computeFeatureDistance(float *feat1, float *feat2, int feat_number);


/**
 * @brief Integrate saliency priors inside a saliency graph by looking at the foreground/background ratio given a cuboid adjacency (https://link.springer.com/content/pdf/10.1007%2Fs11263-017-1062-2.pdf)
 * @date December 18, 2019
 * @author Leo Melo
 * @note This works for both saliency and prior maps
 * @note The arrays starts at 0 and the label image starts at 1. Deduct the label_img val in 1 when finding its corresponding saliency
 * @note The cuboid adjacency is of size AxM - 1 pixels, where M is the amount of saliency maps and A is a given adjacency relation
 * @note This is modeled in the reference paper as a cellular automata and thus was implemented as such.
 * @note The result of the integration replaces the saliency_graph's combined_priors attribute.
 *
 * @param  ref_saliency_graph The reference saliency graph to integrate the priors from.
 * @param  A The flat adjacency relation used for the cuboid adjacency (a 4 adjacency is recommended by the author).
 * @param  label_img Superpixel labeled image.
 * @param  number_iterations The number of iterations used when integrating the maps.
 * @param  lambda The factor in which saliency will be increased or decreased (recommended value by the author = 0.04)
 *
 * @warning The saliency values returned are in the range [0-1].
 * @warning The saliency values of the input saliency maps should be in the range [0-1].
 */
void iftCuboidPriorIntegrationImplicit(iftSaliencyGraph **ref_saliency_graph, iftAdjRel *A, iftImage *label_img, iftITSELFParameters *params);

/**
 * @brief Integrate saliency maps by looking at the foreground/background ratio given a cuboid adjacency (https://link.springer.com/content/pdf/10.1007%2Fs11263-017-1062-2.pdf)
 * @date November 05, 2019
 * @author Leo Melo
 * @note This works for both saliency and prior maps
 * @note The arrays starts at 0 and the label image starts at 1. Deduct the label_img val in 1 when finding its corresponding saliency
 * @note The cuboid adjacency is of size AxM - 1 pixels, where M is the amount of saliency maps and A is a given adjacency relation
 * @note This is modeled in the reference paper as a cellular automata and thus was implemented as such.
 *
 * @param  maps An array of iftImage containing the saliency maps to be integrate.
 * @param  number_maps The number of maps contained in the maps array
 * @param  A The flat adjacency relation used for the cuboid adjacency (a 4 adjacency is recommended by the author)
 * @param  number_iterations The number of iterations used when integrating the maps.
 * @param  lambda The factor in which saliency will be increased or decreased (recommended value by the author = 0.04)
 * @return      Mean saliency array for each label.
 *
 * @warning The saliency values returned are in the range [0-1].
 * @warning The saliency values of the input saliency maps should be in the range [0-1].
 */
iftImage *iftCuboidSaliencyIntegration(iftImage **maps, int number_maps, iftAdjRel *A, int number_iterations, float lambda);

/**
 * @brief Write a superpixel based saliency map based on a float array and a label image
 * @date November 05, 2019
 * @author Leo Melo
 * @note This works for both saliency and prior maps
 *
 * @param  label_img Superpixel labeled image.
 * @param  saliency Saliency values of each superpixel [0-n].
 * @param  filename File path for output.
 * @param  heatmap 0-1 value to determine whether to print a grayscale (0) or a heatmap (1) image.
 *
 * @warning The saliency values should be in 0-1 range.
 */
void iftWriteSaliencyImage(iftImage *label_img, iftSaliencyArray *saliency, char *filename, bool heatmap);

iftImage *updateSuperpixels(iftImage *original_image, iftImage *saliency_img, iftITSELFParameters *params);

iftImage *computeAndUpdateSaliency(iftImage ***saliency_list_ref, int i, iftSaliencyGraph *saliency_graph, iftMImage *lab, iftMImage *mimg, iftImage *superpixel_img, iftImage *query_image, iftITSELFParameters *params);

iftImage *getFinalSaliency(iftImage **saliency_list, iftSaliencyGraph *saliency_graph, iftImage *superpixel_img, iftITSELFParameters *params);

iftImage *iftITSELF(iftImage *orig, iftImage *initial_saliency, iftMImage *features, iftITSELFParameters *params);

/**
 * @brief Improve saliency maps using SESS.
 * @author Leonardo de Melo
 * @date Feb, 2022
 * @param  original_image    Input color image
 * @param  initial_saliency Input saliency map
 * @return improved saliency map
 */

//! swig(newobject)
iftImage *iftSESS(iftImage *original_image, iftImage *initial_saliency, iftITSELFParameters *params, iftMImage *features);

// Additional Stuff that shouldn't be here
iftMImage *iftNormalizeByBand(iftMImage *feats);

iftMImage *iftExtendMImageByGrayObjSalMap(iftMImage *mimg, iftImage* objsm);

iftIGraph *iftInitOISFIGraph(iftImage *img, iftImage *mask, iftImage *objsm);

iftMImage *convertToLab(iftImage *image);

iftDataSet *computeQuantizedDatasetKMeans(iftMImage *mimg, int max_number_of_clusters, int *cluster_number,  int maxIterations, float minImprovement);

void iftWriteOverlay(iftImage* orig, iftImage *label, const char *filename);
void iftWriteLabAsRGB(iftMImage *lab, char *filename);



int iftConnCompCounting(const iftImage *label_img);

int iftIsEmptyMap(iftImage *saliency);


#ifdef __cplusplus
}
#endif


#endif //IFT_Saliency_H_