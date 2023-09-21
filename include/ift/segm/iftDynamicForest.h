#ifndef IFT_DYNAMIC_FOREST_H
#define IFT_DYNAMIC_FOREST_H

#include "iftFImage.h"
#include "iftMImage.h"
#include "ift/core/dtypes/FHeap.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "ift/core/dtypes/IntQueue.h"
#include "ift/core/dtypes/IntArray.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ift_dynamic_tree {
    struct ift_dynamic_tree *prev;
    struct ift_dynamic_tree *next;
    int nnodes;
    int nfeats;
    float *sum_feat;
} iftDynamicTree, iftDynamicNode;

iftDynamicNode *iftInsertDynamicTree(iftDynamicTree **list, int nfeats);
void iftRemoveDynamicTree(iftDynamicTree **list, iftDynamicNode **node);
void iftDestroyDynamicTrees(iftDynamicTree **list);

typedef struct ift_dyn_forest
{
    iftFImage *cost;
    iftFHeap  *heap;
    iftImage  *label;
    iftImage  *marker;
    iftImage  *root;
    iftImage  *pred;
    iftMImage *mimg;
    iftAdjRel *A;

    iftDynamicTree *dynamic_roots;
    iftDynamicTree **tree_map;
} iftDynamicForest;

iftDynamicForest *iftCreateDynamicForest(iftMImage *mimg, iftAdjRel *A);
void iftDestroyDynamicForest(iftDynamicForest **forest);
void iftResetDynamicForest(iftDynamicForest *forest);

float iftDynamicArcWeight(iftDynamicForest *forest, int p, int q);

typedef enum {
    IFT_ADD,
    IFT_REM
} iftUpdateMode;

void iftUpdateDynamicTree(iftDynamicForest *forest, int p, iftUpdateMode mode);

iftSet *iftDynamicTreeRemoval(iftDynamicForest *fst, iftSet *trees_for_removal);
void iftDynamicSubtreeRemoval(iftDynamicForest *forest, int t);

void iftAddMarkersToDynamicForest(iftDynamicForest *forest, iftLabeledSet *seeds);
void iftRemMarkersFromDynamicForest(iftDynamicForest *forest, iftSet *removal_markers);

void iftDiffDynamicIFT(iftDynamicForest *forest);

//! swig(newobject)
iftImage *iftDynamicTreesObject(iftMImage *features, iftLabeledSet *seeds, iftAdjRel *A);

//! swig(newobject)
iftImage *iftDynamicTreesRoot(iftMImage *features, iftLabeledSet *seeds, iftAdjRel *A);

#ifdef __cplusplus
}
#endif

#endif //IFT_DYNAMIC_FOREST_H
