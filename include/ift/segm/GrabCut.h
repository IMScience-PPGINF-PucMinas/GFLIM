//
// Created by jbragantini on 07/12/18.
//

#ifndef IFT_GRABCUT_H
#define IFT_GRABCUT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftImage.h"
#include "iftMImage.h"
#include "iftDataSet.h"
#include "iftKmeans.h"
#include "iftMaxflow.h"

typedef enum ift_node_label {
    SURE_BKG = 0, SURE_OBJ = 1, MAYBE_BKG = 2, MAYBE_OBJ = 3
} iftNodeLabel;

typedef struct ift_gmm {
    double *coefs;
    double *mean;
    double *cov;

    double *invCov;
    double *covDet;

    double *sums;
    double *prods;

    int *grp_size;
    int total_size;

    int n_grp;
} iftGMM;

//! swig(newobject)
iftImage *iftUncertMapToLabel(const iftImage *map);

//! swig(newobject)
iftImage *iftLabelToUncertMap(const iftImage *label, const iftAdjRel *A);

//! swig(newobject)
iftImage *iftMaybeForeground(const iftImage *label);

//! swig(newobject)
iftImage *iftGrabCut(const iftMImage *mimg, const iftImage *regions, double beta, int n_iters);

//! swig(newobject)
iftFImage *iftGMMDataSetDist(iftDataSet *train, const iftMImage *mimg, int n_comps);


#ifdef __cplusplus
}
#endif

#endif //IFT_GRABCUT_H
