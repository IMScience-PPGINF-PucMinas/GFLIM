#ifndef IFT_IMAGEFOREST_H_
#define IFT_IMAGEFOREST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BMap.h"
#include "ift/core/dtypes/GQueue.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "iftCommon.h"
#include "iftImage.h"
#include "iftDataSet.h"
#include "iftAdjacency.h"
#include "iftSeeds.h"


//! swig(destroyer = iftDestroyImageForest, extend = iftImageForestExt.i)
typedef struct ift_imageforest {
    iftImage  *pathval;
    iftImage  *label;
    iftImage  *root;
    iftImage  *pred;
    iftImage  *marker;
    iftBMap   *processed;
    iftGQueue *Q;
    iftImage  *img;
    iftAdjRel *A;
} iftImageForest;

//! swig(newobject)
iftImageForest  *iftCreateImageForest(iftImage *img, iftAdjRel *A);

//! swig(newobject)
void             iftResetImageForest(iftImageForest *fst);
void             iftDestroyImageForest(iftImageForest **fst);
void             iftCompRemovalWithoutFrontier(iftImageForest *fst, iftLabeledSet *seed);
iftSet          *iftCompRemoval(iftImageForest *fst, iftLabeledSet *seed);
iftSet 	  *iftTreeRemoval(iftImageForest *fst, iftSet *trees_for_removal);
iftSet 	  *iftMarkerRemoval(iftImageForest *fst, iftSet *removal_markers);
iftImage          *iftSwitchTreeLabels(iftImageForest *fst, iftLabeledSet *seed, iftImage *label);
void             iftRelabelBorderTreesAsBkg(iftImageForest *fst);
void 		   iftPropagateLabelMarkerAndRootToSubtree(iftImageForest *fst, iftAdjRel *A, iftImage *new_label, iftImage *new_marker, int r);
int 		   iftPathLength(iftImage *pred, int p);

//! swig(newobject)
iftDataSet *iftImageForestToDataSet(iftImageForest* fst, iftAdjRel *A);

#ifdef __cplusplus
}
#endif

#endif


