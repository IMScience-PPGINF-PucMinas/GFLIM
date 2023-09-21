#ifndef IFT_MSF_H_
#define IFT_MSF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftDataSet.h"
#include "iftSort.h"
#include "iftRegion.h"

/* Minimum Spanning Tree */
typedef struct ift_mstnode {
  int   sample;    /* training sample in the original dataset Z, whose
		      weight becomes the maximum among the weights of
		      the arcs that include this node in the MST. */
  iftSet  *adj;    /* adjacent nodes */ 
  int maxarcadj;   /* adjacent node with the maximum arc weight */
  char color;
} iftMSTNode;

//! swig(destroyer = iftDestroyMST)
typedef struct ift_mst {
  iftMSTNode *node;     // node
  int         nnodes;   // number of nodes
  float       maxarcw;  // maximum arc weight in the tree
  float       minarcw;  // minimum arc weight in the tree
  iftDataSet *Z;        // Each graph node is a training sample in Z
} iftMST;

//! swig(newobject)
iftMST *iftCreateMST(iftDataSet *Z);

//Finds the minimum spanning tree for the complete graph defined implicitly by /samples/
iftMST *iftCreateMSTFromSet(iftDataSet *Z, iftSet* samples);

void    iftDestroyMST(iftMST **mstree);
void    iftNormalizeSampleWeightMST(iftMST *mstree);
int     iftSelectSupTrainSamplesByMST(iftDataSet *Z, float train_perc);
void    iftSortNodesByWeightMST(iftMST *mstree, int order);
iftSet *iftSelectSamplesForUserLabelingMST(iftMST *mstree, int n);

void iftSetNodeColors(iftMST *mst, iftSet *samples, char color);

//Finds the minimum spanning tree for the sugraph G' of /graph/ such that every edge in G'
//is incident on nodes with different labels on the dataset
iftMST *iftCreateSuperpixelActiveLearningMST(iftRegionGraph* graph);

#ifdef __cplusplus
}
#endif


#endif
