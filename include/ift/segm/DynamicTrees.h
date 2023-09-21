//
// Created by jordao on 05/06/18.
//

#ifndef IFT_IFTDYNAMICSET_H
#define IFT_IFTDYNAMICSET_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftImage.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "iftMImage.h"
#include "iftDataSet.h"
#include "iftMetrics.h"

//! swig(extend = iftDynamicSetExt.i, destroyer = iftDestroyDynamicSet)
typedef struct ift_dynamic_set
{
    double *mean;
    iftSet *begin;
    iftSet *end;
    int dim;
    int size;
} iftDynamicSet;

//! swig(extend = iftDTForestExt.i, destroyer = iftDestroyDTForest)
typedef struct ift_dt_forest
{
    iftImage *label;
    iftFImage *cost;
    iftImage *root;
    iftImage *pred;
    iftImage *order;
    iftImage *delay;
    iftDynamicSet **dyn_trees;
} iftDTForest;


static inline void iftInsertDynamicSet(iftDynamicSet *S, const iftMImage *mimg, int p)
{
    if (S->size) {
        iftSet *a = (iftSet*) iftAlloc(1, sizeof *a);
        if (!a)
            iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertDynamicSet");
        a->elem = p;
        S->end->next = a;
        S->end = a;
    } else {
        S->begin = (iftSet*) iftAlloc(1, sizeof *S->begin);
        S->begin->elem = p;
        S->end = S->begin;
    }

    S->size += 1;
    for (int i = 0; i < S->dim; i++) {
        S->mean[i] += (mimg->val[p][i] - S->mean[i]) / S->size;
    }
}


iftDynamicSet *iftCreateDynamicSet(int dim);


//! swig(newobject)
iftDTForest *iftCreateDTForest(const iftMImage *mimg, const iftAdjRel *A, const iftLabeledSet *seeds,
                               float delta, float gamma);


void iftDestroyDynamicSet(iftDynamicSet **S);


void iftDestroyDTForest(iftDTForest **forest);


static inline double iftDistDynamicSetMImage(const iftDynamicSet *S, const iftMImage *mimg, int p)
{
    double dist = 0;
    for (int i = 0; i < S->dim; i++) {
        dist += (S->mean[i] - mimg->val[p][i]) * (S->mean[i] - mimg->val[p][i]);
    }

    return dist;
}

//! swig(newobject)
iftImage *iftDynamicSetObjectPolicy(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, bool use_dist);

//! swig(newobject)
iftImage *iftDynamicSetRootPolicy(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist);

//! swig(newobject)
iftImage *iftDynamicSetRootPolicyInMask(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, iftImage *mask);

//! swig(newobject)
iftImage *iftDynamicSetMinRootPolicy(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist);

//! swig(newobject)
iftImage *iftDynamicSetWithCluster(iftMImage *mimg, iftImage *cluster, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist);

//! swig(newobject)
iftImage *iftDynamicSetRootEnhanced(iftMImage *mimg, iftImage *objmap, iftAdjRel *A, iftLabeledSet *seeds, int h, float alpha, bool use_dist);

//! swig(newobject)
iftImage *iftDynamicSetMinRootEnhanced(iftMImage *mimg, iftImage *objmap, iftAdjRel *A, iftLabeledSet *seeds, int h, float alpha, bool use_dist);

/**
 * @param mimg          multi-band original image
 * @param A             adjacency relation
 * @param seeds         labeled seeds nodes
 * @param delta         plato penalization height
 * @param gamma         neighbor distance scaling parameter, must be > 0.0
 * @param objmap        (optional) saliency image
 * @param alpha         (optional) saliency weight, is set to 1.0 if objmap is not given
 * @return label mapping
 */

//! swig(newobject, stable)
iftImage *iftDynTreeRoot(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int delta,
                         float gamma, iftImage *objmap, float alpha);

/**
 * @param mimg          multi-band original image
 * @param A             adjacency relation
 * @param seeds         labeled seeds nodes
 * @param delta         plato penalization height
 * @param gamma         neighbor distance scaling parameter, must be > 0.0
 * @param objmap        (optional) saliency image
 * @param alpha         (optional) saliency weight, is set to 1.0 if objmap is not given
 * @return label mapping
 */

//! swig(newobject, stable)
iftImage *iftDynTreeClosestRoot(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int delta,
                                float gamma, iftImage *objmap, float alpha);


iftImage *iftDynTreeLearned(iftMImage *mimg, iftMatrix *M, iftAdjRel *A, iftLabeledSet *seeds,
                            int delta, float gamma);


typedef struct ift_dynamic_set_CIARP
{
    double *mean;
    int size;
} iftDynamicSet_CIARP;

static inline void iftInsertDynamicSet_CIARP(iftDynamicSet_CIARP *S, iftMImage *mimg, int p)
{
    S->mean[0] += (mimg->val[p][0] - S->mean[0]) / (S->size + 1);
    S->mean[1] += (mimg->val[p][1] - S->mean[1]) / (S->size + 1);
    S->mean[2] += (mimg->val[p][2] - S->mean[2]) / (S->size + 1);

    ++S->size;
}

static inline double iftDistDynamicSetMImage_CIARP(iftDynamicSet_CIARP *S, iftMImage *mimg, int p)
{
    double dist = (S->mean[0] - mimg->val[p][0]) * (S->mean[0] - mimg->val[p][0]) +
                  (S->mean[1] - mimg->val[p][1]) * (S->mean[1] - mimg->val[p][1]) +
                  (S->mean[2] - mimg->val[p][2]) * (S->mean[2] - mimg->val[p][2]);

    return dist;
}

iftDynamicSet_CIARP *iftCreateDynamicSet_CIARP(void);
void iftDestroyDynamicSet_CIARP(iftDynamicSet_CIARP **S);
iftImage *iftDynamicSetObjectPolicy_CIARP(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, bool use_dist);
iftImage *iftDynamicSetRootPolicy_CIARP(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, bool use_dist);
iftImage *iftDynamicSetMinRootPolicy_CIARP(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, bool use_dist);


//! swig(newobject)
iftMImage *iftDTRootWeightsMap(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist);

#ifdef __cplusplus
}
#endif

#endif //IFT_IFTDYNAMICSET_H
