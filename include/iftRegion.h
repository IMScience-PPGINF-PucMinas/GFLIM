#ifndef IFTREGIONGRAPH_H_
#define IFTREGIONGRAPH_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/FHeap.h"
#include "ift/core/dtypes/LabeledSet.h"

#include "iftImage.h"
#include "iftCommon.h"
#include "iftImage.h"
#include "iftAdjacency.h"
#include "iftDataSet.h"
#include "iftMathMorph.h"

typedef struct{
		int node_type; /* default 0 */
	iftSet* adjacent;
} iftRegionGraphNode;

typedef struct{
	iftRegionGraphNode *node;
	iftDataSet *dataset;
	int         nnodes;

	int  *root;
	int  *pred;
	float      *pathval;
	iftFHeap   *heap;
} iftRegionGraph;

/**
 * @brief Creates a region graph
 */
iftRegionGraph* iftCreateRegionGraph(iftDataSet* dataset, int n_regions);

/**
 * @brief Destroys a region graph
 */
void            iftDestroyRegionGraph(iftRegionGraph** rg);

/**
 * @brief Creates a region graph from a label image
 * @param label_image
 * @param dataset
 * @param adjacency
 * @return
 */
iftRegionGraph* iftRegionGraphFromLabelImage(iftImage* label_image, iftDataSet* dataset, iftAdjRel* adjacency);

void            iftSuperpixelClassification(iftRegionGraph *graph, iftImage *label_image, iftLabeledSet *pixel_seeds);

void 			iftDiffSuperpixelClassification(iftRegionGraph *graph, iftImage *label_image, iftLabeledSet *new_seeds);

void			iftFastSuperpixelClassification(iftRegionGraph *graph, iftImage *label_image, iftLabeledSet *pixel_seeds);

void iftSuperpixelClassificationGeodesic(iftRegionGraph *graph, iftImage *label_image, iftLabeledSet *pixel_seeds, float beta);

/**
 * @brief Relabels a label image with n labels different from 0
 * @param labelled A labeled image
 * @param adj_rel Adjacency relation
 * @return An image relabeled from 1 to n leaving region 0 untouched
 */

//! swig(newobject, stable)
iftImage* iftRelabelRegions(iftImage* labelled, iftAdjRel* adj_rel);

//! swig(newobject)
iftLabeledSet* iftGeodesicCenters(const iftImage* label_image);

/**
 * @brief Gets the Voxel's Coordinates from the Geodesic Center from each Supervoxel.
 * @author Samuel Martins
 * @dataset Jun 28, 2016
 *
 * @warning The array of coordinates starts from label 1 to n, i.e, the position [0] corresponds to the
 * center coords from the supervoxel with label 1, and so on.
 *  
 * @param  super_img  	  Supervoxel Image with labels from 1 to n.
 * @param  out_n_clusters Returns by Reference (if != NULL) the number supervoxels (array size) from the image.
 * @return                The array of Voxel coordinates.
 */
iftVoxel *iftGeodesicCenterVoxelCoords(const iftImage *super_img, int *out_n_clusters);


iftLabeledSet* iftGeometricCenters(iftImage* label_image);

iftImage*       iftCreateRefinedLabelImage(iftImage* image, iftLabeledSet* seed, float spatial_radius, int volume_threshold, int steps, float vol_ratio);

//Given a label_image numbered from 1 to n and a dataset with n samples, sets the classes for the superpixels from the ground-truth image.
void iftSetSuperpixelClassesFromGroundTruth(iftImage* label_image, iftImage* groud_truth, iftDataSet* dataset);

//Creates an image where each superpixel has the grayscale value defined by its label in the dataset
iftImage* iftSuperpixelLabelImageFromDataset(iftImage *label_image, iftDataSet *dataset);

//Finds the regions that are adjacent to regions with different label in the /graph/'s dataset
iftBMap		*iftAdjacentToDifferentLabels(iftRegionGraph *graph);

//Binary tree node
typedef struct iftRegionHierarchyNode{
	struct iftRegionHierarchyNode *left;
	struct iftRegionHierarchyNode *right;

	int merging_time; //Indicates when the region was merged
	int region; //If this node is a leaf contains the region id. Otherwise contains NIL
    int id; // unique id starting from 1 used to identify leaves and non-leaves alike

	int xmax, ymax, zmax, xmin, ymin, zmin; //Bounding cube coordinates
} iftRegionHierarchyNode;

//Binary tree representing a hierarchical watershed segmentation
typedef struct{
	iftRegionHierarchyNode *root; //Root of the hierarchy (node with the highest merging_time)
	int nleaves; //Number of leaves in the tree (regions)

	iftSet** pixels_by_node; //List of pixels by (unmerged) region. Speeds up prunning.
	iftImage* label_image; //Label image representing the initial regions
} iftRegionHierarchy;

//Creates an internal region hierarchy node
iftRegionHierarchyNode* iftCreateRegionHierarchyNode(int merging_time, iftRegionHierarchyNode* left, iftRegionHierarchyNode *right, int node_id);

//Creates a region hierarchy leaf node
iftRegionHierarchyNode* iftCreateRegionHierarchyLeaf(iftImage* label_image, iftSet* pixels, int region, int node_id);

typedef float* (*iftMergingFun)(float *f1, float *f2, float *alpha, int n);

iftRegionHierarchy* iftCreateRegionHierarchy(iftDataSet *dataset, iftAdjRel* adjacency, iftMergingFun merge_func);

void iftDestroyRegionHierarchy(iftRegionHierarchy** rh);

//Creates a label_image with nregions given a region hierarchy
iftImage* iftFlattenRegionHierarchy(iftRegionHierarchy* rh, int nregions);

void iftRHDrawSubregions(iftRegionHierarchy *rh, iftRegionHierarchyNode *node, iftImage* label_image, int label);
int iftRHFlatten(iftRegionHierarchy *rh, iftRegionHierarchyNode* node, iftImage *label_image, int cutoff, int current_label);

//Eliminates regions with area below a given number of pixels
iftImage* iftEliminateRegionsByArea(iftDataSet* dataset, iftAdjRel *adj, int threshold);

//Returns a list containing every node in the hierarchy (nleaves + root->merging_time elements)
iftRegionHierarchyNode** iftGetRegionHierarchyNodes(iftRegionHierarchy *rh);

void iftSetRegionGraphNodeTypeForBorderSuperpixelsInTiles(iftRegionGraph *region_graph, iftImage *label, iftImage *tiles_img, iftAdjRel *A);

iftRegionHierarchy* iftMergeBorderRegions(iftDataSet *dataset, iftRegionGraph *region_graph, iftMergingFun merge_func, float threshold,iftRegionHierarchyNode ***output_hierarchy_nodes);

float* iftMergeLabColorMeanStdAndSizeFeatures(float *f1, float *f2, float *alpha, int n);

//! swig(newobject)
iftImage *iftJoinProbabilityRegions(iftMImage *prob, iftImage *label, int norm_val, bool decrement);

#ifdef __cplusplus
}
#endif

#endif // IFTREGIONGRAPH_H_
