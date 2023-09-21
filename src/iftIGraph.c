#include "iftIGraph.h"

#include "ift/core/io/Dir.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"


/*------------------------------- PRIVATE ----------------------------------*/

float     *iftIGraphWeightMarkerByHeight(iftIGraph *igraph, int H);
float     *iftIGraphWeightMarkerByDepth(iftIGraph *igraph, int D);
void       iftIGraphInfRec(iftIGraph *igraph, float *marker);
void       iftIGraphSupRec(iftIGraph *igraph, float *marker);
iftIGraph *iftMImageToIGraph(const iftMImage *img, const iftImage *mask);
float     *iftIGraphAreaClose(iftIGraph *igraph, int area_thres);
float     *iftIGraphAreaOpen(iftIGraph *igraph, int area_thres);
float     *iftIGraphVolumeClose(iftIGraph *igraph, int volume_thres);
float     *iftIGraphVolumeOpen(iftIGraph *igraph, int volume_thres);
int       *iftIGraphSuperpixelCenters(iftIGraph *igraph, int *seed, int nseeds);

/**
 * @brief Remove tree in the image graph to run the differential IFT.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
iftSet *iftIGraphTreeRemoval(iftIGraph *igraph, iftSet **trees_for_removal, double *pvalue, double INITIAL_PATH_VALUE)
{
    int        i, p, q, r, s, t;
    iftVoxel   u, v;
    iftAdjRel *A = igraph->A;
    iftSet    *Frontier = NULL, *adj = NULL;
    iftBMap   *inFrontier = iftCreateBMap(igraph->nnodes);
    iftImage  *index = igraph->index;
    iftSet    *T1 = NULL, *T2 = NULL;

    switch(igraph->type) {

        case IMPLICIT:

            /* Remove all marked trees and find the frontier voxels
           afterwards. */

            while (*trees_for_removal != NULL){
                s = iftRemoveSet(trees_for_removal);
                p = igraph->node[s].voxel;
                r = igraph->root[p];

                if (pvalue[index->val[r]] != INITIAL_PATH_VALUE){ /* tree not marked yet */
                    igraph->pvalue[r] = pvalue[index->val[r]] = INITIAL_PATH_VALUE; /* mark removed root */
                    igraph->pred[r]   = IFT_NIL;
                    iftInsertSet(&T1, r);
                    while (T1 != NULL){
                        p = iftRemoveSet(&T1);
                        iftInsertSet(&T2, p); /* compute in T2 the union of removed trees */
                        u = iftGetVoxelCoord(index, p);
                        for (i = 1; i < A->n; i++){
                            v = iftGetAdjacentVoxel(A, u, i);
                            if (iftValidVoxel(index, v)){
                                q   = iftGetVoxelIndex(index, v);
                                t   = index->val[q];
                                if ((t != IFT_NIL) && (pvalue[t] != INITIAL_PATH_VALUE)){ /* q has not been removed */
                                    if (igraph->pred[q] == p){ /* q belongs to the tree under removal */
                                        iftInsertSet(&T1, q);
                                        pvalue[t]         = igraph->pvalue[q] = INITIAL_PATH_VALUE; /* mark removed node */
                                        igraph->pred[q]   = IFT_NIL;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /* Find the frontier voxels of non-removed trees */

            while (T2 != NULL){
                p = iftRemoveSet(&T2);
                u = iftGetVoxelCoord(index, p);
                for (i = 1; i < A->n; i++){
                    v = iftGetAdjacentVoxel(A, u, i);
                    if (iftValidVoxel(index, v)){
                        q   = iftGetVoxelIndex(index, v);
                        t   = index->val[q];
                        if ((t != IFT_NIL) && (pvalue[t] != INITIAL_PATH_VALUE)){ /* q is a frontier node */
                            if (iftBMapValue(inFrontier, t) == 0){ /* t has not been inserted in the frontier yet */
                                iftInsertSet(&Frontier, t);
                                iftBMapSet1(inFrontier, t);
                            }
                        }
                    }
                }
            }

            break;

        case EXPLICIT:

            /* Remove all marked trees and find the frontier voxels
           afterwards. */

            while (*trees_for_removal != NULL){
                s = iftRemoveSet(trees_for_removal);
                p = igraph->node[s].voxel;
                r = igraph->root[p];

                if (pvalue[index->val[r]] != INITIAL_PATH_VALUE){ /* tree not marked yet */
                    pvalue[index->val[r]] = igraph->pvalue[r] = INITIAL_PATH_VALUE; /* mark removed root */
                    igraph->pred[r]     = IFT_NIL;
                    iftInsertSet(&T1, r);
                    while (T1 != NULL){
                        p = iftRemoveSet(&T1);
                        s = index->val[p];
                        iftInsertSet(&T2, p); /* compute in T2 the union of removed trees */
                        for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
                            t   = adj->elem;
                            q   = igraph->node[t].voxel;
                            if (pvalue[t] != INITIAL_PATH_VALUE){ /* q has not been removed */
                                if (igraph->pred[q] == p){ /* q belongs to the tree under removal */
                                    iftInsertSet(&T1, q);
                                    pvalue[t]       = igraph->pvalue[q] = INITIAL_PATH_VALUE; /* mark removed node */
                                    igraph->pred[q] = IFT_NIL;
                                }
                            }
                        }
                    }
                }
            }


            /* Find the frontier voxels of non-removed trees */

            while (T2 != NULL){
                p = iftRemoveSet(&T2);
                s = igraph->index->val[p];
                for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
                    t   = adj->elem;
                    if (pvalue[t] != INITIAL_PATH_VALUE){ /* t is a frontier node */
                        if (iftBMapValue(inFrontier, t) == 0){  /* t has not been inserted in the frontier yet */
                            iftInsertSet(&Frontier, t);
                            iftBMapSet1(inFrontier, t);
                        }
                    }
                }
            }

            break;

        case COMPLETE:

            /* Remove all marked trees and find the frontier voxels
           afterwards. */

            while (*trees_for_removal != NULL){
                s = iftRemoveSet(trees_for_removal);
                p = igraph->node[s].voxel;
                r = igraph->root[p];

                if (pvalue[index->val[r]] != INITIAL_PATH_VALUE){ /* tree not marked yet */
                    igraph->pvalue[r] = pvalue[index->val[r]] = INITIAL_PATH_VALUE; /* mark removed root */
                    igraph->pred[r]   = IFT_NIL;
                    iftInsertSet(&T1, r);
                    while (T1 != NULL){
                        p = iftRemoveSet(&T1);
                        s = index->val[p];
                        iftInsertSet(&T2, p); /* compute in T2 the union of removed trees */

                        for (t=0; t < igraph->nnodes; t++) {
                            if (t != s ) {
                                q = igraph->node[t].voxel;
                                if (pvalue[t] != INITIAL_PATH_VALUE){ /* q has not been removed */
                                    if (igraph->pred[q] == p){ /* q belongs to the tree under removal */
                                        iftInsertSet(&T1, q);
                                        pvalue[t] = igraph->pvalue[q] = INITIAL_PATH_VALUE; /* mark removed node */
                                        igraph->pred[q]   = IFT_NIL;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /* Find the frontier voxels of non-removed trees */

            while (T2 != NULL){
                p = iftRemoveSet(&T2);
                s = index->val[p];
                for (t=0; t < igraph->nnodes; t++) {
                    if (t != s ) {
                        if (pvalue[t] != INITIAL_PATH_VALUE){ /* t is a frontier node */
                            if (iftBMapValue(inFrontier, t) == 0){ /* t has not been inserted in the frontier yet */
                                iftInsertSet(&Frontier, t);
                                iftBMapSet1(inFrontier, t);
                            }
                        }
                    }
                }
            }

        default:
            iftError("Invalid type of image graph", "iftIGraphTreeRemoval");
    }

    iftDestroyBMap(&inFrontier);

    return (Frontier);
}

/**
 * @brief Remove subtreef a given node s in the image graph during the
 * differential IFT with a non-MI function.
 * @author Alexandre Falcao
 * @date July 14th, 2017
 * @ingroup ImageGraph
 *
 */

void iftIGraphSubTreeRemoval(iftIGraph *igraph, int s, double *pvalue, double INITIAL_PATH_VALUE, iftDHeap *Q)
{
    int        i, p, q, t;
    iftVoxel   u, v;
    iftAdjRel *A = igraph->A;
    iftSet    *Frontier = NULL, *adj = NULL, *Subtree = NULL;
    iftImage  *index = igraph->index;

    switch(igraph->type) {

        case IMPLICIT:

            /* Reinitialize voxels (nodes) of the subtree of s to be
               reconquered and compute the frontier nodes (voxels) */

            iftInsertSet(&Subtree,s);

            while (Subtree != NULL){
                s = iftRemoveSet(&Subtree);
                p = igraph->node[s].voxel;
                u = iftGetVoxelCoord(index,p);

                if (Q->color[s]==IFT_GRAY)
                    iftRemoveDHeapElem(Q,s);


                igraph->pvalue[p] = pvalue[s] = INITIAL_PATH_VALUE;
                igraph->pred[p]   = IFT_NIL;

                for (i = 1; i < A->n; i++){
                    v = iftGetAdjacentVoxel(A, u, i);
                    if (iftValidVoxel(index, v)){
                        q   = iftGetVoxelIndex(index, v);
                        t   = index->val[q];
                        if (igraph->pred[q]==p)
                            iftInsertSet(&Subtree,t);
                        else{ /* consider t as a candidate to be a frontier node */
                            iftInsertSet(&Frontier,t);
                        }
                    }
                }
            }

            break;

        case EXPLICIT:


            /* Reinitialize voxels (nodes) of the subtree of s to be
           reconquered and compute the frontier nodes (voxels) */

            iftInsertSet(&Subtree,s);

            while (Subtree != NULL){
                s = iftRemoveSet(&Subtree);
                p = igraph->node[s].voxel;

                if (Q->color[s]==IFT_GRAY)
                    iftRemoveDHeapElem(Q,s);

                igraph->pvalue[p] = pvalue[s] = INITIAL_PATH_VALUE;
                igraph->pred[p]   = IFT_NIL;

                for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
                    t = adj->elem;
                    q = igraph->node[t].voxel;

                    if (igraph->pred[q]==p)
                        iftInsertSet(&Subtree,t);
                    else{ /* consider t as a candidate to be a frontier node */
                        iftInsertSet(&Frontier,t);
                    }
                }
            }

            break;

        case COMPLETE:

            /* Reinitialize voxels (nodes) of the subtree of s to be
           reconquered and compute the frontier nodes (voxels) */

            iftInsertSet(&Subtree,s);

            while (Subtree != NULL){
                s = iftRemoveSet(&Subtree);
                p = igraph->node[s].voxel;
                u = iftGetVoxelCoord(index,p);

                if (Q->color[s]==IFT_GRAY)
                    iftRemoveDHeapElem(Q,s);

                igraph->pvalue[p] = pvalue[s] = INITIAL_PATH_VALUE;
                igraph->pred[p]   = IFT_NIL;

                for (t=0; t < igraph->nnodes; t++) {
                    if (t != s ) {
                        q = igraph->node[t].voxel;

                        if (igraph->pred[q]==p)
                            iftInsertSet(&Subtree,t);
                        else{ /* consider t as a candidate to be a frontier node */
                            iftInsertSet(&Frontier,t);
                        }
                    }
                }
            }

            break;

        default:
            iftError("Invalid type of image graph", "iftIGraphSubtreeRemoval");
    }

    /* Identify the real frontier nodes and insert them in Queue to
       continue the DIFT */

    while (Frontier != NULL){
        s = iftRemoveSet(&Frontier);
        p = igraph->node[s].voxel;
        if (igraph->label[p] != 0){
            if (Q->color[s] == IFT_GRAY)
                iftGoUpDHeap(Q, Q->pos[s]);
            else
                iftInsertDHeap(Q,s);
        }
    }

}


/**
 * @brief Obtain the superpixel centers (the closest node to the center in the feature space)
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
int *iftIGraphSuperpixelCenters(iftIGraph *igraph, int *seed, int nseeds)
{
    int    i, j, p, q, s, *center;
    float  **feat,  *nelems, dist1, dist2;

    /* compute average feature vector for each superpixel */

    feat   = (float **)iftAlloc(nseeds,sizeof(float *));
    nelems = iftAllocFloatArray(nseeds);
    center = iftAllocIntArray(nseeds);
    for (i=0; i < nseeds; i++){
        feat[i]   = iftAllocFloatArray(igraph->nfeats);
        center[i] = seed[i];
    }
    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        int label = igraph->label[igraph->root[p]];
        if (label > 0) {
            i = label-1;
            nelems[i]++;
            for (j=0; j < igraph->nfeats; j++)
                feat[i][j] += igraph->feat[p][j];
        }
    }

    for (i=0; i < nseeds; i++) {
        for (j=0; j < igraph->nfeats; j++)
	  feat[i][j] /= nelems[i];
    }

    /* compute the closest node to each superpixel center */

    for (s=0; s < igraph->nnodes; s++) {
        p     = igraph->node[s].voxel;
        int label = igraph->label[igraph->root[p]];
        if (label > 0) {
            i     = label-1;
            q     = igraph->node[center[i]].voxel;
            dist1 = iftFeatDistance(feat[i],igraph->feat[q],igraph->nfeats);
            dist2 = iftFeatDistance(feat[i],igraph->feat[p],igraph->nfeats);
            if (dist2 < dist1)
                center[i]=s;
        }
    }

    for (i=0; i < nseeds; i++) {
        iftFree(feat[i]);
    }
    iftFree(feat);
    iftFree(nelems);

    return(center);
}

/**
 * @brief Obtain markers by volume closing.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
float *iftIGraphVolumeClose(iftIGraph *igraph, int volume_thres)
{
    iftCompTree *ctree=NULL;
    int i,s;
    float *marker;
    int *level=NULL;


    ctree = iftIGraphCreateMinTree(igraph);
    iftCumSize(ctree,ctree->root);
    level = iftAllocIntArray(ctree->numnodes);
    for (i=0; i < ctree->numnodes; i++)
        level[i]=ctree->node[i].level;
    for (i=0; i < ctree->numnodes; i++)
        if (ctree->node[i].numsons==0)
            level[i]=iftVolumeLevel(ctree,level,i,volume_thres,0);

    marker = iftAllocFloatArray(igraph->nnodes);

    for (s=0; s < igraph->nnodes; s++) {
        marker[s]=level[ctree->cmap->val[s]] + 1; // It guarantees that marker is above weight
    }

    iftDestroyCompTree(&ctree);
    iftFree(level);

    return(marker);
}

/**
 * @brief Obtain markers by volume opening.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
float *iftIGraphVolumeOpen(iftIGraph *igraph, int volume_thres)
{
    iftCompTree *ctree=NULL;
    int i,s;
    float *marker;
    int *level=NULL;

    ctree = iftIGraphCreateMaxTree(igraph);
    iftCumSize(ctree,ctree->root);
    level = iftAllocIntArray(ctree->numnodes);
    for (i=0; i < ctree->numnodes; i++){
        level[i]=ctree->node[i].level;
    }

    for (i=0; i < ctree->numnodes; i++)
        if (ctree->node[i].numsons==0)
            level[i]=iftVolumeLevel(ctree,level,i,volume_thres,0);

    marker = iftAllocFloatArray(igraph->nnodes);

    for (s=0; s < igraph->nnodes; s++) {
        marker[s]=level[ctree->cmap->val[s]]; // It guarantees that marker
        // is below weight
        igraph->node[s].weight = marker[s] + 1;
    }

    iftDestroyCompTree(&ctree);
    iftFree(level);

    return(marker);
}

/**
 * @brief Obtain markers by area closing.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
float *iftIGraphAreaClose(iftIGraph *igraph, int area_thres)
{
    iftCompTree *ctree=NULL;
    int i,s;
    float *marker;
    int *level=NULL;

    ctree = iftIGraphCreateMinTree(igraph);
    iftCumSize(ctree,ctree->root);
    level = iftAllocIntArray(ctree->numnodes);
    for (i=0; i < ctree->numnodes; i++)
        level[i]=ctree->node[i].level;

    for (i=0; i < ctree->numnodes; i++)
        if (ctree->node[i].numsons==0)
            level[i]=iftAreaLevel(ctree,level,i,area_thres);

    marker = iftAllocFloatArray(igraph->nnodes);

    for (s=0; s < igraph->nnodes; s++) {
        marker[s]=level[ctree->cmap->val[s]] + 1; // it guarantees that marker is above weight
    }
    iftDestroyCompTree(&ctree);
    iftFree(level);

    return(marker);
}

/**
 * @brief Obtain markers by area opening.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
float *iftIGraphAreaOpen(iftIGraph *igraph, int area_thres)
{
    iftCompTree *ctree=NULL;
    int i,s;
    float *marker;
    int *level=NULL;

    ctree = iftIGraphCreateMaxTree(igraph);
    iftCumSize(ctree,ctree->root);
    level = iftAllocIntArray(ctree->numnodes);
    for (i=0; i < ctree->numnodes; i++)
        level[i]=ctree->node[i].level;

    for (i=0; i < ctree->numnodes; i++)
        if (ctree->node[i].numsons==0)
            level[i]=iftAreaLevel(ctree,level,i,area_thres);

    marker = iftAllocFloatArray(igraph->nnodes);

    for (s=0; s < igraph->nnodes; s++) {
        marker[s]=level[ctree->cmap->val[s]];   // it guarantees that marker is below weight
        igraph->node[s].weight = marker[s] + 1;
    }
    iftDestroyCompTree(&ctree);
    iftFree(level);

    return(marker);
}


/**
 * @brief Create node weight marker at level (node weight - H - 1). Node
 * weight is modified to avoid negative numbers in the
 * marker. Afterwards, weight and path value must be subtracted of
 * (H+1). See iftDualWaterGrayByHeight.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
float *iftIGraphWeightMarkerByHeight(iftIGraph *igraph, int H)
{
    float *marker = iftAllocFloatArray(igraph->nnodes);
    int s;

    for (s=0; s < igraph->nnodes; s++) {
        marker[s]  = igraph->node[s].weight;
        igraph->node[s].weight = marker[s] + (H+1);
    }

    return(marker);
}


/**
 * @brief Create node weight marker at level (node weight + D).
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
float *iftIGraphWeightMarkerByDepth(iftIGraph *igraph, int D)
{
    float *marker = iftAllocFloatArray(igraph->nnodes);
    int s;

    for (s=0; s < igraph->nnodes; s++) {
        marker[s]  = igraph->node[s].weight + (D+1);
    }

    return(marker);
}


/**
 * @brief Compute inferior reconstruction between weight and marker with
 * label, root, predecessor propagation among nodes. Note that
 * marker[s] < igraph->node[s].weight is needed for inferior
 * reconstruction with label propagation.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
void iftIGraphInfRec(iftIGraph *igraph, float *marker)
{
    double     tmp, *pvalue;
    int        s, t, i, p, q, l=1;
    iftSet    *adj;
    iftVoxel   u, v;
    iftDHeap  *Q;

    /* check weight marker */

    for (s=0; s < igraph->nnodes; s++) {
        if (marker[s] >= igraph->node[s].weight)
            iftError("Weight marker must be less than weight", "iftIGraphInfRec");
    }

    /* path initialization */

    pvalue = iftAllocDoubleArray(igraph->nnodes);
    Q      = iftCreateDHeap(igraph->nnodes, pvalue);
    iftSetRemovalPolicyDHeap(Q, IFT_MAXVALUE);

    for (s=0; s < igraph->nnodes; s++) {
        p                  = igraph->node[s].voxel;
        pvalue[s]          = marker[s];
        igraph->pred[p]    = IFT_NIL;
        igraph->root[p]    = p;
        iftInsertDHeap(Q,s);
    }

    /* optimum-path forest computation */

    switch (igraph->type) {

        case IMPLICIT:

            while (!iftEmptyDHeap(Q)) {
                s = iftRemoveDHeap(Q);
                p = igraph->node[s].voxel;
                if (igraph->pred[p] == IFT_NIL){
                    pvalue[s] = igraph->node[s].weight;
                    igraph->label[p]    = l;
                    l++;
                }
                igraph->pvalue[p] = pvalue[s];
                u = iftGetVoxelCoord(igraph->index,p);
                for (i=1; i < igraph->A->n; i++) {
                    v = iftGetAdjacentVoxel(igraph->A,u,i);
                    if (iftValidVoxel(igraph->index,v)){
                        q = iftGetVoxelIndex(igraph->index,v);
                        t = igraph->index->val[q];
                        if ((t != IFT_NIL) && (Q->color[t] != IFT_BLACK)){
                            tmp = iftMin(pvalue[s], igraph->node[t].weight);
                            if (tmp > pvalue[t]) {
                                pvalue[t]          = tmp;
                                igraph->label[q]   = igraph->label[p];
                                igraph->root[q]    = igraph->root[p];
                                igraph->pred[q]    = p;
                                iftGoUpDHeap(Q, Q->pos[t]);
                            }
                        }
                    }
                }
            }

            break;

        case EXPLICIT:

            /* optimum-path forest computation */

            while (!iftEmptyDHeap(Q)) {
                s = iftRemoveDHeap(Q);
                p = igraph->node[s].voxel;
                if (igraph->pred[p] == IFT_NIL){
                    pvalue[s] = igraph->node[s].weight;
                    igraph->label[p]    = l;
                    l++;
                }
                igraph->pvalue[p]  = pvalue[s];
                for (adj = igraph->node[s].adj; adj != NULL; adj = adj->next) {
                    t = adj->elem;
                    q = igraph->node[t].voxel;
                    if (Q->color[t] != IFT_BLACK){
                        tmp = iftMin(pvalue[s], igraph->node[t].weight);
                        if (tmp > pvalue[t]) {
                            pvalue[t]          = tmp;
                            igraph->label[q]   = igraph->label[p];
                            igraph->root[q]    = igraph->root[p];
                            igraph->pred[q]    = p;
                            iftGoUpDHeap(Q, Q->pos[t]);
                        }
                    }
                }
            }

            break;

        case COMPLETE:
        default:
            iftError("Invalid graph type", "iftIGraphInfRec");

    }

    iftDestroyDHeap(&Q);
    iftFree(pvalue);
}

/**
 * @brief Compute superior reconstruction between weight and marker with
 * label, root, predecessor propagation among nodes. Note that
 * marker[s] > igraph->node[s].weight is needed for superior
 * reconstruction with label propagation.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
void iftIGraphSupRec(iftIGraph *igraph, float *marker)
{
    double     tmp, *pvalue;
    int        s, t, i, p, q, l=1;
    iftSet    *adj;
    iftVoxel   u, v;
    iftDHeap  *Q;

    /* check weight marker */

    for (s=0; s < igraph->nnodes; s++) {
        if (marker[s] <= igraph->node[s].weight)
            iftError("Weight marker must be greater than weight", "iftIGraphInfRec");
    }

    /* path initialization */

    pvalue = iftAllocDoubleArray(igraph->nnodes);
    Q      = iftCreateDHeap(igraph->nnodes, pvalue);

    for (s=0; s < igraph->nnodes; s++) {
        p                  = igraph->node[s].voxel;
        pvalue[s]          = marker[s];
        igraph->pred[p]    = IFT_NIL;
        igraph->root[p]    = p;
        iftInsertDHeap(Q,s);
    }

    /* optimum-path forest computation */

    switch (igraph->type) {

        case IMPLICIT:

            while (!iftEmptyDHeap(Q)) {
                s = iftRemoveDHeap(Q);
                p = igraph->node[s].voxel;
                if (igraph->pred[p] == IFT_NIL){
                    pvalue[s] = igraph->node[s].weight;
                    igraph->label[p]    = l;
                    l++;
                }
                igraph->pvalue[p] = pvalue[s];
                u = iftGetVoxelCoord(igraph->index,p);
                for (i=1; i < igraph->A->n; i++) {
                    v = iftGetAdjacentVoxel(igraph->A,u,i);
                    if (iftValidVoxel(igraph->index,v)){
                        q = iftGetVoxelIndex(igraph->index,v);
                        t = igraph->index->val[q];
                        if ((t != IFT_NIL) && (Q->color[t] != IFT_BLACK)){
                            tmp = iftMax(pvalue[s], igraph->node[t].weight);
                            if (tmp < pvalue[t]) {
                                pvalue[t]          = tmp;
                                igraph->label[q]   = igraph->label[p];
                                igraph->root[q]    = igraph->root[p];
                                igraph->pred[q]    = p;
                                iftGoUpDHeap(Q, Q->pos[t]);
                            }
                        }
                    }
                }
            }

            break;

        case EXPLICIT:

            /* optimum-path forest computation */

            while (!iftEmptyDHeap(Q)) {
                s = iftRemoveDHeap(Q);
                p = igraph->node[s].voxel;
                if (igraph->pred[p] == IFT_NIL){
                    pvalue[s] = igraph->node[s].weight;
                    igraph->label[p]    = l;
                    l++;
                }
                igraph->pvalue[p]  = pvalue[s];
                for (adj = igraph->node[s].adj; adj != NULL; adj = adj->next) {
                    t = adj->elem;
                    q = igraph->node[t].voxel;
                    if (Q->color[t] != IFT_BLACK){
                        tmp = iftMax(pvalue[s], igraph->node[t].weight);
                        if (tmp < pvalue[t]) {
                            pvalue[t]          = tmp;
                            igraph->label[q]   = igraph->label[p];
                            igraph->root[q]    = igraph->root[p];
                            igraph->pred[q]    = p;
                            iftGoUpDHeap(Q, Q->pos[t]);
                        }
                    }
                }
            }

            break;

        case COMPLETE:
        default:
            iftError("Invalid graph type", "iftIGraphSupRec");

    }

    iftDestroyDHeap(&Q);
    iftFree(pvalue);
}


/**
 * @brief Convert MImage into an image graph.
 * @author Alexandre Falcao
 * @date May 12, 2015
 * @ingroup ImageGraph
 *
 */
iftIGraph *iftMImageToIGraph(const iftMImage *img, const iftImage *mask)
{
    iftIGraph *igraph = (iftIGraph *)iftAlloc(1,sizeof(iftIGraph));
    int        p, i;

    if (mask != NULL)    
      igraph->nnodes  = iftNumberOfElements(mask);
    else
      igraph->nnodes  = img->n;
	
    igraph->node    = (iftINode *)iftAlloc(igraph->nnodes,sizeof(iftINode));
    igraph->index   = iftCreateImage(img->xsize, img->ysize, img->zsize);
    igraph->nfeats  = img->m;
    igraph->feat    = (float **)iftAlloc(img->n,sizeof(float *));

    iftCopyVoxelSize(img, igraph->index);
    for (p=0; p < img->n; p++) {
        igraph->feat[p] = iftAllocFloatArray(igraph->nfeats);
        for (i=0; i < img->m; i++)
            igraph->feat[p][i] = img->val[p][i];
    }
    igraph->label   = iftAllocIntArray(img->n);
    igraph->marker  = iftAllocIntArray(img->n);
    igraph->root    = iftAllocIntArray(img->n);
    igraph->pred    = iftAllocIntArray(img->n);
    igraph->pvalue  = iftAllocDoubleArray(img->n);

    if (mask != NULL){
      for (p=0, i=0; p < img->n; p++) {
        igraph->index->val[p]     = IFT_NIL;
        if (mask->val[p]>0){
	  igraph->node[i].adj     = NULL;
	  igraph->node[i].voxel   = p;
	  igraph->node[i].weight  = 0.0;
	  igraph->index->val[p]   = i;
	  i++;
        }
      }
    } else {
      for (p=0, i=0; p < img->n; p++) {
	igraph->node[i].adj     = NULL;
	igraph->node[i].voxel   = p;
	igraph->node[i].weight  = 0.0;
	igraph->index->val[p]   = i;
	i++;
      }
      
    }
    
    return(igraph);
}

/*------------------------------- PUBLIC ------------------------------------*/

iftIGraph *iftCompleteIGraph(iftMImage *img, iftImage *mask)
{
    iftIGraph *igraph = iftMImageToIGraph(img,mask);

    igraph->A       = NULL;
    igraph->type    = COMPLETE;

    return(igraph);
}

iftIGraph *iftImplicitIGraph(iftMImage *img, const iftImage *mask, iftAdjRel *A)
{
    iftIGraph *igraph = iftMImageToIGraph(img,mask);

    igraph->A       = iftCopyAdjacency(A);
    igraph->type    = IMPLICIT;

    return(igraph);
}

iftIGraph *iftExplicitIGraph(const iftMImage *img, const iftImage *mask, const iftImage *label, iftAdjRel *A)
{
    iftIGraph *igraph = iftMImageToIGraph(img, mask);
    int        i, s, t, p, q;
    iftVoxel   u, v;

    igraph->type = EXPLICIT;

    if ((mask == NULL)||(igraph->nnodes == img->n)) {

      /* label is ignored, voxels are nodes and arcs connect adjacent
	 voxels according to A */
    
      for (p=0; p < img->n; p++) {
	s = igraph->index->val[p];
        u = iftMGetVoxelCoord(img,p);
        for (i=1; i < A->n; i++) {
	  v = iftGetAdjacentVoxel(A,u,i);
	  if (iftMValidVoxel(img,v)) {
	    q = iftMGetVoxelIndex(img,v);
	    t = igraph->index->val[q];
	    iftInsertSet(&igraph->node[s].adj,t);
	  }
	}
      }
      
    } else { 

      if (label == NULL) {
	
	/* voxels in the mask are the nodes and it is expected that
	   among adjacent voxels, according to A, there are nodes to
	   form the arcs. One can use it to create a graph on the
	   surface of a 3D object, for instance. */
	
	for (p=0; p < img->n; p++) {
	  s = igraph->index->val[p];
	  if (s != IFT_NIL){
	    u = iftMGetVoxelCoord(img,p);
	    for (i=1; i < A->n; i++) {
	      v = iftGetAdjacentVoxel(A,u,i);
	      if (iftMValidVoxel(img,v)) {
		q = iftMGetVoxelIndex(img,v);
		t = igraph->index->val[q];
		if (t != IFT_NIL)
		  iftInsertSet(&igraph->node[s].adj,t);
	      }
	    }
	  }
	}
	
      } else {
	
	/* voxels in the mask are the nodes, but the arcs are formed
	   between nodes that fall in adjacent regions, according to
	   A, of distinct labels. This can be used to create a
	   superpixel graph, for instance. */

	int nsuperpixels = iftMaximumValue(label);
	int *root = iftAllocIntArray(nsuperpixels+1);
	for (int i=0; i <= nsuperpixels; i++)
	  root[i] = IFT_NIL;

	for (s=0; s < igraph->nnodes; s++){
	  p = igraph->node[s].voxel;
	  root[label->val[p]]=s;
	}
	
	for (p=0; p < img->n; p++) {
	  u = iftMGetVoxelCoord(img,p);
	  for (i=1; i < A->n; i++) {
	    v = iftGetAdjacentVoxel(A,u,i);
	    if (iftMValidVoxel(img,v)) {
	      q = iftMGetVoxelIndex(img,v);
	      if (label->val[p]!=label->val[q]){
		s = root[label->val[p]];
		t = root[label->val[q]];
		iftUnionSetElem(&igraph->node[s].adj,t);
	      }
	    }
	  }
	}
      }
    }
    
    return(igraph);
}

iftIGraph *iftKnnIGraph(iftMImage *img, iftImage *mask, int K)
{
    iftIGraph *igraph = iftMImageToIGraph(img,mask);
    int        p, s, t, q, k, nn[K+1], j;
    float      d[K+1];

    igraph->type  = EXPLICIT;
    igraph->A     = NULL;

    if (K > igraph->nnodes - 2)
        K = igraph->nnodes - 2;

    for (s=0; s < igraph->nnodes; s++) {
        p     = igraph->node[s].voxel;
        d[0]  = IFT_INFINITY_FLT; nn[0] = IFT_NIL;
        k     = 1;
        for (t=0; t < igraph->nnodes; t++) {
            if (s != t) {
                q     = igraph->node[t].voxel;
                d[k]  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
                j     = k;
                nn[k] = t;
                while ((j > 0)&&(d[j] < d[j-1])){
                    iftSwap(d[j], d[j-1]);
                    iftSwap(nn[j], nn[j-1]);
                    j--;
                }
                if (k < K) k++;
            }
        }
        for (j=k-1; j >= 0; j--)  {
            iftInsertSet(&igraph->node[s].adj,nn[j]);
        }
    }


    return(igraph);
}

iftIGraph *iftSpatialKnnIGraph(iftMImage *img, iftImage *mask, iftAdjRel *A, int K)
{
    iftIGraph *igraph = iftMImageToIGraph(img,mask);
    int        p, s, t, i, q, k, nn[K+1], j;
    iftVoxel   u, v;
    float      d[K+1];

    igraph->type  = EXPLICIT;
    igraph->A     = A;

    if (K > (A->n-2))
        K = A->n-2;

    for (s=0; s < igraph->nnodes; s++) {
        p     = igraph->node[s].voxel;
        u     = iftGetVoxelCoord(mask,p);
        d[0]  = IFT_INFINITY_FLT; nn[0] = IFT_NIL;
        k     = 1;
        for (i=1; i < A->n; i++) {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(mask,v)){
                q  = iftGetVoxelIndex(mask,v);
                if (igraph->index->val[q] != IFT_NIL){
                    t  = igraph->index->val[q];
                    d[k]  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
                    nn[k] = t;
                    j     = k;
                    while ((j > 0)&&(d[j] < d[j-1])){
                        iftSwap(d[j], d[j-1]);
                        iftSwap(nn[j], nn[j-1]);
                        j--;
                    }
                    if (k < K) k++;
                }
            }
        }
        for (j=k-1; j >= 0; j--)  {
            if (nn[j] != IFT_NIL)
                iftInsertSet(&igraph->node[s].adj,nn[j]);
        }
    }

    return(igraph);
}

iftIGraph *iftSpatialIGraph(iftMImage *img, iftImage *mask, iftAdjRel *A, float df)
{
    iftIGraph *igraph = iftMImageToIGraph(img,mask);
    int        p, s, t, i, q;
    iftVoxel   u, v;
    float      dist;

    igraph->type  = EXPLICIT;
    igraph->A     = A;


    for (s=0; s < igraph->nnodes; s++) {
        p     = igraph->node[s].voxel;
        u     = iftGetVoxelCoord(mask,p);

        for (i=1; i < A->n; i++) {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(mask,v)){
                q  = iftGetVoxelIndex(mask,v);
                if (igraph->index->val[q] != IFT_NIL){
                    t     = igraph->index->val[q];
                    dist  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
                    if (dist <= df)
                        iftInsertSet(&igraph->node[s].adj,t);
                }
            }
        }
    }

    return(igraph);
}



void iftDestroyIGraph(iftIGraph **igraph)
{
    if(igraph != NULL && *igraph != NULL) {
        iftIGraph *aux = *igraph;
        int p, i;

        for (i = 0; i < aux->nnodes; i++) {
            if (aux->node[i].adj != NULL)
                iftDestroySet(&aux->node[i].adj);
        }
        for (p = 0; p < aux->index->n; p++) {
            iftFree(aux->feat[p]);
        }
        iftFree(aux->feat);
        iftFree(aux->label);
        iftFree(aux->marker);
        iftFree(aux->root);
        iftFree(aux->pred);
        iftFree(aux->pvalue);

        if (aux->type == IMPLICIT)
            iftDestroyAdjRel(&aux->A);

        iftFree(aux->node);
        iftDestroyImage(&aux->index);
        iftFree(aux);
        (*igraph) = NULL;
    }
}



iftImage *iftIGraphPathValue(iftIGraph *igraph)
{
    int p;
    iftImage *pvalue = iftCreateImage(igraph->index->xsize,igraph->index->ysize,igraph->index->zsize);

    for (p=0; p < pvalue->n; p++) {
        pvalue->val[p] = (int) igraph->pvalue[p];
    }
    iftCopyVoxelSize(igraph->index, pvalue);

    return(pvalue);
}

iftFImage *iftIGraphWeight(iftIGraph *igraph)
{
    int p,s;
    iftFImage *weight = iftCreateFImage(igraph->index->xsize,igraph->index->ysize,igraph->index->zsize);

    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        weight->val[p] = igraph->node[s].weight;
    }
    iftCopyVoxelSize(igraph->index, weight);

    return(weight);
}

iftImage *iftIGraphLabel(iftIGraph *igraph)
{
    int p;
    iftImage *label = iftCreateImage(igraph->index->xsize,igraph->index->ysize,igraph->index->zsize);

    for (p=0; p < label->n; p++) {
        label->val[p] = igraph->label[p];
    }

    iftCopyVoxelSize(igraph->index, label);

    return(label);
}

iftImage *iftIGraphRoot(iftIGraph *igraph)
{
    int p;
    iftImage *root = iftCreateImage(igraph->index->xsize,igraph->index->ysize,igraph->index->zsize);

    for (p=0; p < root->n; p++) {
        root->val[p] = igraph->root[p];
    }

    iftCopyVoxelSize(igraph->index, root);

    return(root);
}

void iftIGraphFixLabelRootMap(iftIGraph *igraph, iftSet **T) {
    iftSet *T2=NULL, *adj=NULL;
    iftLabeledSet *T1=NULL;
    iftVoxel u,v;
    int p, q, r, label, s, t;

    switch (igraph->type) {

        case IMPLICIT:

            /* Select from T only nodes p that present at least one son, whose
               label is inconsistent with respect to the label of p, and
               insert it in T1 */

            while((*T) != NULL) {
                p = iftRemoveSet(T);
                u = iftGetVoxelCoord(igraph->index, p);

                for (int i = 1; i < igraph->A->n; ++i) {
                    v = iftGetAdjacentVoxel(igraph->A, u, i);
                    if (iftValidVoxel(igraph->index, v)){
                        q = iftGetVoxelIndex(igraph->index, v);
                        if ( (igraph->pred[q] == p) &&
                             (igraph->label[q]!=igraph->label[p]) ) {
                            iftInsertLabeledSet(&T1, p, igraph->label[p]);
                            break;
                        }
                    }
                }
            }

            /* The subtree of each node r in T1 must be painted by the label
               of r. However, if there are multiple nodes in T1 that belong to
               a same tree of the forest, the label of the node r, which is
               the closest to the root of the tree, precedes the label of the
               other nodes. */

            while(T1 != NULL) {

                r = iftRemoveLabeledSet(&T1, &label);

                if (igraph->label[r] == label) { /* The node r has not been painted
					  yet by any other node in T1
					  with higher precedence than
					  r. In this case, the label of r
					  propagates to its subtree. This
					  is just an optimization. */
                    iftInsertSet(&T2, r);
                    while (T2 != NULL) {
                        p = iftRemoveSet(&T2);
                        u = iftGetVoxelCoord(igraph->index, p);

                        for (int i = 1; i < igraph->A->n; ++i) {
                            v = iftGetAdjacentVoxel(igraph->A, u, i);
                            if (iftValidVoxel(igraph->index, v)) {
                                q = iftGetVoxelIndex(igraph->index, v);
                                if (igraph->pred[q] == p) { /* correct label and root maps */
                                    iftInsertSet(&T2, q);
                                    igraph->label[q] = igraph->label[p];
                                    igraph->root[q]  = igraph->root[p];
                                }
                            }
                        }
                    }
                }
            }
            break;

        case EXPLICIT:

            /* Select from T only nodes p that present at least one son, whose
               label is inconsistent with respect to the label of p, and
               insert it in T1 */

            while((*T) != NULL) {
                p = iftRemoveSet(T);
                s = igraph->index->val[p];

                for (adj = igraph->node[s].adj; adj != NULL; adj = adj->next){
                    t = adj->elem;
                    q = igraph->node[t].voxel;
                    if ( (igraph->pred[q] == p) &&
                         (igraph->label[q]!=igraph->label[p]) ) {
                        iftInsertLabeledSet(&T1, p, igraph->label[p]);
                        break;
                    }
                }
            }

            /* The subtree of each node r in T1 must be painted by the label
               of r. However, if there are multiple nodes in T1 that belong to
               a same tree of the forest, the label of the node r, which is
               the closest to the root of the tree, precedes the label of the
               other nodes. */

            while(T1 != NULL) {

                r = iftRemoveLabeledSet(&T1, &label);

                if (igraph->label[r] == label) { /* The node r has not been painted
					  yet by any other node in T1
					  with higher precedence than
					  r. In this case, the label of r
					  propagates to its subtree. This
					  is just an optimization. */
                    iftInsertSet(&T2, r);
                    while (T2 != NULL) {
                        p = iftRemoveSet(&T2);
                        s = igraph->index->val[p];
                        for (adj = igraph->node[s].adj; adj != NULL; adj = adj->next){
                            t = adj->elem;
                            q = igraph->node[t].voxel;

                            if (igraph->pred[q] == p) { /* correct label and root maps */
                                iftInsertSet(&T2, q);
                                igraph->label[q] = igraph->label[p];
                                igraph->root[q]  = igraph->root[p];
                            }
                        }
                    }
                }
            }

            break;

        case COMPLETE:

            /* Select from T only nodes p that present at least one son, whose
               label is inconsistent with respect to the label of p, and
               insert it in T1 */

            while((*T) != NULL) {
                p = iftRemoveSet(T);
                s = igraph->index->val[p];

                for (t = 0; t < igraph->nnodes; t++){
                    if (s != t) {
                        q = igraph->node[t].voxel;
                        if ( (igraph->pred[q] == p) &&
                             (igraph->label[q]!=igraph->label[p]) ) {
                            iftInsertLabeledSet(&T1, p, igraph->label[p]);
                            break;
                        }
                    }
                }
            }

            /* The subtree of each node r in T1 must be painted by the label
               of r. However, if there are multiple nodes in T1 that belong to
               a same tree of the forest, the label of the node r, which is
               the closest to the root of the tree, precedes the label of the
               other nodes. */

            while(T1 != NULL) {

                r = iftRemoveLabeledSet(&T1, &label);

                if (igraph->label[r] == label) { /* The node r has not been painted
					  yet by any other node in T1
					  with higher precedence than
					  r. In this case, the label of r
					  propagates to its subtree. This
					  is just an optimization. */
                    iftInsertSet(&T2, r);
                    while (T2 != NULL) {
                        p = iftRemoveSet(&T2);
                        s = igraph->index->val[p];
                        for (t = 0; t < igraph->nnodes; t++){
                            if (s != t) {
                                q = igraph->node[t].voxel;
                                if (igraph->pred[q] == p) { /* correct label and root maps */
                                    iftInsertSet(&T2, q);
                                    igraph->label[q] = igraph->label[p];
                                    igraph->root[q]  = igraph->root[p];
                                }
                            }
                        }
                    }
                }
            }

            break;

        default:
            iftError("Invalid igraph type", "iftIGraphFixLabelRootMap");
    }

    iftDestroyLabeledSet(&T1);
    iftDestroySet(&T2);
}

int iftIGraphISF_Root(iftIGraph *igraph, iftImage *seeds, double alpha, double beta, int niters)
{
    double      tmp;
    int        r, s, t, i, p, q, it, *seed, nseeds;
    float      new_seeds_flag;
    iftVoxel   u, v;
    iftDHeap  *Q;
    double    *pvalue = iftAllocDoubleArray(igraph->nnodes);
    iftSet    *adj=NULL, *S=NULL, *new_seeds=NULL, *frontier_nodes=NULL, *trees_for_removal=NULL;
    //  iftSet    *T =NULL; /* This is used for the old differential version */
    //timer *t1, *t2;


    new_seeds_flag = 1.0;

    /* Initial seed set and trivial path initialization to infinity
       cost */

    Q      = iftCreateDHeap(igraph->nnodes, pvalue);

    nseeds = 0;
    for (s=0; s < igraph->nnodes; s++) {
        p               = igraph->node[s].voxel;
        pvalue[s]       = igraph->pvalue[p] = IFT_INFINITY_DBL;
        igraph->pred[p] = IFT_NIL;
        if (seeds->val[p]!=0){
            iftInsertSet(&new_seeds,s);
            nseeds++;
        }
    }

    seed    = iftAllocIntArray(nseeds);
    S       = new_seeds; i = 0;
    while (S != NULL) {
        seed[i] = S->elem;
        p       = igraph->node[seed[i]].voxel;
        igraph->label[p] = i+1;
        i++; S = S->next;
    }


    /* differential optimum-path forest computation */

    //  for (I=0; (I < niters)&&(new_seeds_flag >= 0.1); I++) {  // If 90% of the seed do not change much (satisfy both the criteria in LAB and XY) we stop
    for (it=0; (it < niters); it++) {


        // printf("iteration %d\n",it+1);
        //t1 = iftTic();

        if (trees_for_removal != NULL)
            frontier_nodes = iftIGraphTreeRemoval(igraph, &trees_for_removal, pvalue, IFT_INFINITY_DBL);

        while (new_seeds != NULL) { /* insert seeds in the priority queue Q with cost zero */
            s = iftRemoveSet(&new_seeds);
            p = igraph->node[s].voxel;
	    if (igraph->label[p] > 0){ /* we must avoid removed seeds
					  from previous iteration */
	      pvalue[s] = igraph->pvalue[p] = 0;
	      igraph->root[p] = p;
	      igraph->pred[p] = IFT_NIL;
	      iftInsertDHeap(Q,s);
	    }
        }

        while (frontier_nodes != NULL) { /* insert frontier nodes in Q to represent the previous seeds */
            s = iftRemoveSet(&frontier_nodes);
            if (Q->color[s] == IFT_WHITE){
                iftInsertDHeap(Q,s);
            }
        }


        switch(igraph->type) {

            case IMPLICIT:

                while (!iftEmptyDHeap(Q)) {
                    s = iftRemoveDHeap(Q);
                    p = igraph->node[s].voxel;
                    r = igraph->root[p];
                    igraph->pvalue[p] = pvalue[s];
                    u = iftGetVoxelCoord(igraph->index,p);
                    for (i=1; i < igraph->A->n; i++) {
                        v = iftGetAdjacentVoxel(igraph->A,u,i);
                        if (iftValidVoxel(igraph->index,v)){
                            q   = iftGetVoxelIndex(igraph->index,v);
                            t   = igraph->index->val[q];
                            if ((t != IFT_NIL) && (Q->color[t] != IFT_BLACK)){
                                tmp = pvalue[s] + pow(alpha*(double)iftFeatDistance(igraph->feat[r],igraph->feat[q],igraph->nfeats),beta) + (double)iftVoxelDistance(u,v);

                                if (tmp < pvalue[t]){
                                    pvalue[t]            = tmp;

                                    /* /\* This is used for the old differential version *\/ */
                                    /* if (igraph->label[p] != igraph->label[q]){ /\* voxels */
                                    /* 						that */
                                    /* 						have */
                                    /* 						changed */
                                    /* 						labels *\/ */
                                    /*   iftInsertSet(&T, q); */
                                    /* } */

                                    igraph->root[q]      = igraph->root[p];
                                    igraph->label[q]     = igraph->label[p];
                                    igraph->pred[q]      = p;

                                    if (Q->color[t] == IFT_GRAY){
                                        iftGoUpDHeap(Q, Q->pos[t]);
                                    } else {
                                        iftInsertDHeap(Q,t);
                                    }
                                }	else {
                                    if (igraph->pred[q] == p){
                                        if (tmp > pvalue[t]) {
                                            iftIGraphSubTreeRemoval(igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                                        } else { /* tmp == pvalue[t] */
                                            if ((igraph->label[q] != igraph->label[p])&&(igraph->label[q]!=0))
                                                iftIGraphSubTreeRemoval(igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                                        }
                                    }
                                }

                            }
                        }
                    }
                }



                break;

            case EXPLICIT:

                while (!iftEmptyDHeap(Q)) {
                    s = iftRemoveDHeap(Q);
                    p = igraph->node[s].voxel;
                    r = igraph->root[p];
                    u = iftGetVoxelCoord(igraph->index,p);
                    igraph->pvalue[p] = pvalue[s];

                    for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
                        t   = adj->elem;
                        if (Q->color[t] != IFT_BLACK){
                            q   = igraph->node[t].voxel;
                            v   = iftGetVoxelCoord(igraph->index,q);
                            tmp = pvalue[s] + pow(alpha*(double)iftFeatDistance(igraph->feat[r],igraph->feat[q],igraph->nfeats),beta) + (double)iftVoxelDistance(u,v);

                            if (tmp < pvalue[t]){
                                pvalue[t]            = tmp;

                                /* /\* This is used for the old differential version *\/ */
                                /* if (igraph->label[p] != igraph->label[q]){ /\* voxels */
                                /* 						that */
                                /* 						have */
                                /* 						changed */
                                /* 						labels *\/ */
                                /*   iftInsertSet(&T, q); */
                                /* } */

                                igraph->root[q]      = igraph->root[p];
                                igraph->label[q]     = igraph->label[p];
                                igraph->pred[q]      = p;

                                if (Q->color[t] == IFT_GRAY){
                                    iftGoUpDHeap(Q, Q->pos[t]);
                                } else {
                                    iftInsertDHeap(Q,t);
                                }
                            } else {
                                if (igraph->pred[q] == p){
                                    if (tmp > pvalue[t]) {
                                        iftIGraphSubTreeRemoval(igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                                    } else { /* tmp == pvalue[t] */
                                        if ((igraph->label[q] != igraph->label[p]))
                                            iftIGraphSubTreeRemoval(igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                                    }
                                }
                            }
                        }
                    }
                }

                break;

            case COMPLETE:
                iftError("Not implemented for complete graphs", "iftIGraphISF_Root");
        }


        iftResetDHeap(Q);

        /* This is for the old differential version */
        /* if (it>0) */
        /*     iftIGraphFixLabelRootMap(igraph, &T); */
        /* iftDestroySet(&T); */


        /* End of comment */

        /* Recompute new seeds */

        iftIGraphISFRecomputeSeeds(igraph, seed, nseeds, &trees_for_removal, &new_seeds, &new_seeds_flag);

        //t2 = iftToc();
        // iftPrintCompTime(t1,t2,"Computational time for iteration %d",it+1);

        /* Uncomment this block and comment the one above for
           non-differential IFT */

        /* iftResetDHeap(Q); */
        /* iftDestroySet(&trees_for_removal); */
        /* iftDestroySet(&new_seeds); */

        /* for (s=0; s < igraph->nnodes; s++) { */
        /*   p = igraph->node[s].voxel; */
        /*   pvalue[s]       = IFT_INFINITY_DBL; */
        /* } */
        /* for (int i = 0; i < nseeds; ++i) { */
        /*   s = seed[i]; */
        /*   iftInsertSet(&new_seeds,s); */
        /*   p = igraph->node[s].voxel; */
        /*   igraph->label[p] = i+1; */
        /* } */


        /* End of comment */

    }


    iftDestroySet(&adj);
    iftDestroySet(&S);
    iftDestroySet(&new_seeds);
    iftDestroySet(&frontier_nodes);
    iftDestroySet(&trees_for_removal);
    iftDestroyDHeap(&Q);
    iftFree(pvalue);
    iftFree(seed);

    return it;
}

iftIGraph *iftIGraphOISF
(iftIGraph *igraph, iftImage *seeds, double alpha, double beta, double gamma, int iters)
{
  // Variables
  int nseeds, seed_count;
  double tmp;
  int *seed;
  float *dummy;
  double *pvalue;
  iftDHeap *Q;
  iftSet *S, *new_seeds, *frontier_nodes, *trees_rm;
  iftIGraph *out_igraph;

  // Init
  S = NULL; new_seeds = NULL; frontier_nodes = NULL;
  trees_rm = NULL;
  out_igraph = iftCopyIGraph(igraph);

  // Assign
  nseeds = 0;
  pvalue = iftAllocDoubleArray(out_igraph->nnodes);
  Q = iftCreateDHeap(out_igraph->nnodes, pvalue);
  dummy = iftAllocFloatArray(nseeds);

  // Trivial paths
  for(int s = 0; s < out_igraph->nnodes; s++) 
  {
    // Variables
    int p;

    // Assign   
    p = out_igraph->node[s].voxel;
    pvalue[s] = IFT_INFINITY_DBL;
    out_igraph->pvalue[p] = IFT_INFINITY_DBL;
    out_igraph->pred[p] = IFT_NIL;
    
    // If it happens to be a seed
    if(seeds->val[p] != 0) {
        // Save it for later
        iftInsertSet(&new_seeds,s);
        nseeds++;
    }
  }

  // Assign II
  seed_count = 0;
  S = new_seeds; 
  seed = iftAllocIntArray(nseeds);

  // Process all seeds
  while (S != NULL) {
    // Variables
    int p;

    // Assign
    seed[seed_count] = S->elem;
    p = out_igraph->node[seed[seed_count]].voxel;
    out_igraph->label[p] = seed_count + 1;

    seed_count++; 
    S = S->next;
  }

  // Object-based Iterative Spanning Forest
  for(int it = 0; it < iters; it++) {

    // Remove trees which are marked for update (DIFT)
    if (trees_rm != NULL) { 
      frontier_nodes = iftIGraphTreeRemoval(out_igraph, &trees_rm, pvalue, IFT_INFINITY_DBL);
    }

    // Set trivial paths for seeds
    while(new_seeds != NULL)  {
      // Variables
      int s, p;

      // Assign
      s = iftRemoveSet(&new_seeds);
      p = out_igraph->node[s].voxel;  

      // If its reachable
      if(out_igraph->label[p] > 0) { 
        // Set trivial paths
        pvalue[s] = 0;
        out_igraph->pvalue[p] = 0;
        out_igraph->root[p] = p;
        out_igraph->pred[p] = IFT_NIL;

        iftInsertDHeap(Q,s);
      }
    }

    // Insert the nodes from the removed trees (DIFT)
    while(frontier_nodes != NULL) {
      // Variables
      int s;

      // Assign
      s = iftRemoveSet(&frontier_nodes);
      
      if (Q->color[s] == IFT_WHITE) iftInsertDHeap(Q,s);
    } 

    // What kind of iftIGraph is it?
    switch(out_igraph->type) {

      case IMPLICIT:
        // Image Foresting Transform
        while (!iftEmptyDHeap(Q)) {
          // Variables
          int s, p, r;
          iftVoxel u;

          // Assign
          s = iftRemoveDHeap(Q);
          p = out_igraph->node[s].voxel;
          r = out_igraph->root[p];
          out_igraph->pvalue[p] = pvalue[s];
          u = iftGetVoxelCoord(out_igraph->index,p);

          // For every adjacent node
          for(int i = 1; i < out_igraph->A->n; i++) {
            // Variables
            iftVoxel v;

            // Assign
            v = iftGetAdjacentVoxel(out_igraph->A,u,i);

            // If it is in the domain
            if(iftValidVoxel(out_igraph->index,v)) {
              // Variables
              int q, t;

              // Assign
              q   = iftGetVoxelIndex(out_igraph->index,v);
              t   = out_igraph->index->val[q];

              // If it is reachable and it is in the queue
              if( t != IFT_NIL && Q->color[t] != IFT_BLACK ) {
                // Variables
                double color_dist, spat_dist, obj_dist;

                // Compute distances
                color_dist = (double)iftFeatDistance(out_igraph->feat[r], out_igraph->feat[q], out_igraph->nfeats-1);
                spat_dist = (double)iftVoxelDistance(u,v);
                obj_dist = (double)(fabs((out_igraph->feat[r][out_igraph->nfeats-1] - out_igraph->feat[q][out_igraph->nfeats-1])));

                // Sum all distances
                tmp = pow( alpha*color_dist*pow(gamma, obj_dist) +gamma*obj_dist, beta);
                tmp += spat_dist;
                tmp += pvalue[s];

                // If it is a path with lower cost
                if(tmp < pvalue[t]) {              
                  // Assign
                  pvalue[t] = tmp;
                  out_igraph->root[q] = out_igraph->root[p];
                  out_igraph->label[q] = out_igraph->label[p];
                  out_igraph->pred[q] = p;

                  // Update if its in the queue
                  if (Q->color[t] == IFT_GRAY) iftGoUpDHeap(Q, Q->pos[t]);
                  else iftInsertDHeap(Q,t); // Insert it otherwise
                } 
                else {
                  // If the predecessor of 'q' is 'p' (DIFT)
                  if (out_igraph->pred[q] == p)
                  {
                    if (tmp > pvalue[t]) iftIGraphSubTreeRemoval(out_igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                    else if( out_igraph->label[q] != out_igraph->label[p] && out_igraph->label[q] != 0 ) {
                      iftIGraphSubTreeRemoval(out_igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                    }
                  }
                }
              } 
            }
          }
        }
        break;

      case EXPLICIT: 
        iftWarning("CAUTION! Needs testing for explicit graphs!", "iftIGraphOISF" );
        // Image Foresting Transform
        while (!iftEmptyDHeap(Q)) {
          // Variables
          int s, p, r;
          iftVoxel u;
          iftSet *adj;

          // Init
          adj = NULL;

          // Assign
          s = iftRemoveDHeap(Q);
          p = out_igraph->node[s].voxel;
          r = out_igraph->root[p];
          out_igraph->pvalue[p] = pvalue[s];
          u = iftGetVoxelCoord(out_igraph->index,p);

          // For every adjacent node
          for(adj = out_igraph->node[s].adj; adj != NULL; adj = adj->next) {
            // Variables
            int t;

            // Assign
            t = adj->elem;

            // If it is in the queue
            if(Q->color[t] != IFT_BLACK) {
              // Variables
              int q;
              double color_dist, spat_dist, obj_dist;
              iftVoxel v;

              // Assign
              q   = out_igraph->node[t].voxel;
              v = iftGetVoxelCoord(out_igraph->index,q);

              // Compute distances
              color_dist = (double)iftFeatDistance(out_igraph->feat[r], out_igraph->feat[q], out_igraph->nfeats-1);
              spat_dist = (double)iftVoxelDistance(u,v);
              obj_dist = (double)(fabs((out_igraph->feat[r][out_igraph->nfeats-1] - out_igraph->feat[q][out_igraph->nfeats-1])));

              // Sum all distances
              tmp = pow( alpha*color_dist*pow(gamma, obj_dist) +gamma*obj_dist, beta);
              tmp += spat_dist;
              tmp += pvalue[s];

              // If it is a path with lower cost
              if(tmp < pvalue[t]) {              
                // Assign
                pvalue[t] = tmp;
                out_igraph->root[q] = out_igraph->root[p];
                out_igraph->label[q] = out_igraph->label[p];
                out_igraph->pred[q] = p;

                // Update if its in the queue
                if (Q->color[t] == IFT_GRAY) iftGoUpDHeap(Q, Q->pos[t]);
                else iftInsertDHeap(Q,t); // Insert it otherwise
              } 
              else {
                // If the predecessor of 'q' is 'p' (DIFT)
                if (out_igraph->pred[q] == p)
                {
                  if (tmp > pvalue[t]) iftIGraphSubTreeRemoval(out_igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                  else if( out_igraph->label[q] != out_igraph->label[p] && out_igraph->label[q] != 0 ) {
                    iftIGraphSubTreeRemoval(out_igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                  }
                }
              }
            }
          }
        }

        break;
      case COMPLETE:
          iftError("Not implemented for complete graphs", "iftIGraphOISF");
        break;
      default:
          iftError("Unknown type of graph", "iftIGraphOISF");
        break;
    }

     

    // Reset properties
    iftResetDHeap(Q);
    
    if( iters > 1 ) iftIGraphISFRecomputeSeeds(out_igraph, seed, nseeds, &trees_rm, &new_seeds, dummy);
  }

  // Free
  iftDestroySet(&S);
  iftDestroySet(&new_seeds);
  iftDestroySet(&frontier_nodes);
  iftDestroySet(&trees_rm);
  iftDestroyDHeap(&Q);
  iftFree(pvalue);
  iftFree(seed);
  iftFree(dummy);

  return (out_igraph);
}


int *iftIGraphSuperpixelCentersUsingSpatialInformation(iftIGraph *igraph, int *seed, int nseeds, float **seed_features)
{
    int    i, j, p, q, s, *center, ndim;
    float  **feat,  *nelems, dist1, dist2, *featp, *featq;
    iftVoxel u, v;
    ndim = 3; // X,Y and Z

    /* compute average feature vector for each superpixel */

    featp = iftAllocFloatArray(ndim);
    featq = iftAllocFloatArray(ndim);
    
    feat   = (float **)iftAlloc(nseeds,sizeof(float *));
    nelems = iftAllocFloatArray(nseeds);
    center = iftAllocIntArray(nseeds);
    for (i=0; i < nseeds; i++){
        feat[i] = iftAllocFloatArray(ndim);
        center[i] = seed[i];
        /* update seed features NEW */
        for (j = 0; j < igraph->nfeats; ++j)
            seed_features[i][j] = 0;
    }

    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        i = igraph->label[igraph->root[p]]-1;

        nelems[i]++;

        u = iftGetVoxelCoord(igraph->index,p);

        feat[i][0] += u.x;
        feat[i][1] += u.y;
        feat[i][2] += u.z;
        /* update seed features NEW */
        for (j = 0; j < igraph->nfeats; ++j)
            seed_features[i][j] += igraph->feat[p][j];
    }

    for (i=0; i < nseeds; i++) {
        for (j=0; j < ndim; j++)
            feat[i][j] /= nelems[i];
        /* update seed features NEW */
        for (j = 0; j < igraph->nfeats; ++j)
            seed_features[i][j] /= nelems[i];
    }

    /* compute the closest node to each superpixel center */

    for (s=0; s < igraph->nnodes; s++) {
        p     = igraph->node[s].voxel;
        i     = igraph->label[igraph->root[p]]-1;
        q     = igraph->node[center[i]].voxel;

        u = iftGetVoxelCoord(igraph->index,p);
        v = iftGetVoxelCoord(igraph->index,q);

        featp[0] = u.x;
        featp[1] = u.y;
        featp[2] = u.z;

        featq[0] = v.x;
        featq[1] = v.y;
        featq[2] = v.z;

        dist1 = iftFeatDistance(feat[i],featq,ndim);
        dist2 = iftFeatDistance(feat[i],featp,ndim);
        if (dist2 < dist1)
            center[i]=s;
    }

    for (i=0; i < nseeds; i++) {
        iftFree(feat[i]);
    }
    iftFree(feat);
    iftFree(nelems);
    iftFree(featp);
    iftFree(featq);

    return(center);
}

int iftIGraphISF_Mean(iftIGraph *igraph, iftImage *seeds, double alpha, double beta, int niters)
{
    double      tmp;
    int        s, t, i, p, q, it, *seed, nseeds, j, index_seed;
    float      new_seeds_flag;
    iftVoxel   u, v;
    iftDHeap  *Q;
    double    *pvalue = iftAllocDoubleArray(igraph->nnodes);
    iftSet    *adj=NULL, *S=NULL, *new_seeds=NULL, *frontier_nodes=NULL, *trees_for_removal=NULL;
    /* iftSet    *T =NULL; /\* This is part of the old version of the DIFT *\/ */
    float     **seed_features; // NEW

    //    timer *t1, *t2;


    new_seeds_flag = 1.0;

    /* Initial seed set and trivial path initialization to infinity
       cost */

    Q      = iftCreateDHeap(igraph->nnodes, pvalue);

    nseeds = 0;
    for (s=0; s < igraph->nnodes; s++) {
        p               = igraph->node[s].voxel;
        pvalue[s]       = igraph->pvalue[p] = IFT_INFINITY_DBL;
        igraph->pred[p] = IFT_NIL;
        if (seeds->val[p]!=0){
            iftInsertSet(&new_seeds,s);
            nseeds++;
        }
    }

    /* Alloc seed features NEW */
    seed_features = (float**) iftAlloc(nseeds ,sizeof(float*));
    for (i = 0; i < nseeds; ++i)
        seed_features[i] = iftAllocFloatArray(igraph->nfeats);

    seed  = iftAllocIntArray(nseeds);
    S     = new_seeds; i = 0;
    while (S != NULL) {
        seed[i] = S->elem;
        p       = igraph->node[seed[i]].voxel;
        igraph->label[p] = i+1;

        // set initial seed feature NEW
        for (j = 0; j < igraph->nfeats; ++j)
            seed_features[i][j] = igraph->feat[p][j];

        i++; S = S->next;
    }


    /* differential optimum-path forest computation */
    for (it=0; (it < niters); it++) {

      //        printf("iteration %d\n",it+1);
      //        t1 = iftTic();

        if (trees_for_removal != NULL)
            frontier_nodes = iftIGraphTreeRemoval(igraph, &trees_for_removal, pvalue, IFT_INFINITY_DBL);

        while (new_seeds != NULL) { /* insert seeds in the priority queue Q with cost zero */
            s = iftRemoveSet(&new_seeds);
            p = igraph->node[s].voxel;
	    if (igraph->label[p] != 0){ /* we must avoid removed seeds
					  from previous iteration */
	      pvalue[s] = igraph->pvalue[p] = 0;
	      igraph->root[p] = p;
	      igraph->pred[p] = IFT_NIL;
	      iftInsertDHeap(Q,s);
	    }
	}
	
        while (frontier_nodes != NULL) { /* insert frontier nodes in Q to represent the previous seeds */
            s = iftRemoveSet(&frontier_nodes);
            if (Q->color[s] == IFT_WHITE)
                iftInsertDHeap(Q,s);
        }

        switch(igraph->type) {

            case IMPLICIT:


                while (!iftEmptyDHeap(Q)) {
                    s = iftRemoveDHeap(Q);
                    p = igraph->node[s].voxel;
                    //r = igraph->root[p];
                    index_seed = igraph->label[p] - 1;
                    //iftVoxel w = iftGetVoxelCoord(igraph->index,r);
                    igraph->pvalue[p] = pvalue[s];
                    u = iftGetVoxelCoord(igraph->index,p);
                    for (i=1; i < igraph->A->n; i++) {
                        v = iftGetAdjacentVoxel(igraph->A,u,i);
                        if (iftValidVoxel(igraph->index,v)){
                            q   = iftGetVoxelIndex(igraph->index,v);
                            t   = igraph->index->val[q];
                            if ((t != IFT_NIL) && (Q->color[t] != IFT_BLACK)){

			      tmp = pvalue[s] + pow(alpha*(double)iftFeatDistance(seed_features[index_seed],igraph->feat[q],igraph->nfeats),beta) + (double)iftVoxelDistance(u,v);

                                if (tmp < pvalue[t]){
                                    pvalue[t]            = tmp;

                                    /* /\* This is used for the old differential version *\/ */
                                    /* if (igraph->label[p] != igraph->label[q]){ /\* voxels */
                                    /* 						that */
                                    /* 						have */
                                    /* 						changed */
                                    /* 						labels *\/ */
                                    /*   iftInsertSet(&T, q); */
                                    /* } */

                                    igraph->root[q]      = igraph->root[p];
                                    igraph->label[q]     = igraph->label[p];
                                    igraph->pred[q]      = p;

                                    if (Q->color[t] == IFT_GRAY){
                                        iftGoUpDHeap(Q, Q->pos[t]);
                                    } else {
                                        iftInsertDHeap(Q,t);
                                    }
                                }	else {
                                    if (igraph->pred[q] == p){
                                        if (tmp > pvalue[t]) {
                                            iftIGraphSubTreeRemoval(igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                                            if (Q->color[s] != IFT_GRAY)
                                              iftInsertDHeap(Q, s);
                                        } else { /* tmp == pvalue[t] */
                                            if ((igraph->label[q] != igraph->label[p])&&(igraph->label[q]!=0)){
                                                iftIGraphSubTreeRemoval(igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                                                if (Q->color[s] != IFT_GRAY)
                                                  iftInsertDHeap(Q, s);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                break;

            case EXPLICIT:

                while (!iftEmptyDHeap(Q)) {
                    s = iftRemoveDHeap(Q);
                    p = igraph->node[s].voxel;
                    //r = igraph->root[p];
                    index_seed = igraph->label[p] - 1;
                    u = iftGetVoxelCoord(igraph->index,p);
                    igraph->pvalue[p] = pvalue[s];

                    for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
                        t   = adj->elem;
                        if (Q->color[t] != IFT_BLACK){
                            q   = igraph->node[t].voxel;
                            v   = iftGetVoxelCoord(igraph->index,q);
                            tmp = pvalue[s] + pow(alpha*(double)iftFeatDistance(seed_features[index_seed],igraph->feat[q],igraph->nfeats),beta) + (double)iftVoxelDistance(u,v);


                            if (tmp < pvalue[t]){
                                pvalue[t]            = tmp;

                                /* /\* This is used for the old differential version *\/ */
                                /* if (igraph->label[p] != igraph->label[q]){ /\* voxels */
                                /* 						that */
                                /* 						have */
                                /* 						changed */
                                /* 						labels *\/ */
                                /*   iftInsertSet(&T, q); */
                                /* } */

                                igraph->root[q]      = igraph->root[p];
                                igraph->label[q]     = igraph->label[p];
                                igraph->pred[q]      = p;

                                if (Q->color[t] == IFT_GRAY){
                                    iftGoUpDHeap(Q, Q->pos[t]);
                                } else {
                                    iftInsertDHeap(Q,t);
                                }
                            }	else {
                                if (igraph->pred[q] == p){
                                    if (tmp > pvalue[t]) {
                                        iftIGraphSubTreeRemoval(igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                                        if (Q->color[s] != IFT_GRAY)
                                          iftInsertDHeap(Q, s);
                                    } else { /* tmp == pvalue[t] */
                                        if (igraph->label[q] != igraph->label[p]){
                                            iftIGraphSubTreeRemoval(igraph,t,pvalue,IFT_INFINITY_DBL,Q);
                                            if (Q->color[s] != IFT_GRAY)
                                              iftInsertDHeap(Q, s);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                break;

            case COMPLETE:
                iftError("Not implemented for complete graphs", "iftIGraphISF_Mean");
        }


        iftResetDHeap(Q);

        /* /\* This is part of the old version of the DIFT *\/ */

        /* if (it>0) */
        /*     iftIGraphFixLabelRootMap(igraph, &T); */
        /* iftDestroySet(&T); */


        /* End of comment */

        /* Recompute new seeds */

        iftIGraphISFRecomputeSeedsUsingSpatialInformation(igraph, seed, nseeds, &trees_for_removal, &new_seeds, &new_seeds_flag, seed_features);


	//        t2 = iftToc();
	//        iftPrintCompTime(t1,t2,"Computational time for iteration %d",it+1);

        /* Uncomment this block and comment the one above for
           non-differential IFT */
        /*
        iftResetDHeap(Q);
        iftDestroySet(&trees_for_removal);
        iftDestroySet(&new_seeds);

        for (s=0; s < igraph->nnodes; s++) {
           p = igraph->node[s].voxel;
           pvalue[s]       = IFT_INFINITY_DBL;
        }
        for (int i = 0; i < nseeds; ++i) {
           s = seed[i];
           iftInsertSet(&new_seeds,s);
           p = igraph->node[s].voxel;
           igraph->label[p] = i+1;
        }
        */

        /* End of comment */

    }

    /* Free seed_features NEW */
    for (i = 0; i < nseeds; ++i)
        iftFree(seed_features[i]);
    iftFree(seed_features);

    iftDestroySet(&adj);
    iftDestroySet(&S);
    iftDestroySet(&new_seeds);
    iftDestroySet(&frontier_nodes);
    iftDestroySet(&trees_for_removal);
    iftDestroyDHeap(&Q);
    iftFree(pvalue);
    iftFree(seed);

    return it;
}


iftIGraph *iftGenerateIGraphSuperpixelsByISF(const iftImage *img, int input_n_cells, float alpha,
                                             float beta, int n_iters, int smooth_n_iters) {
    if (img == NULL)
        iftError("Image is NULL", "iftGenerateIGraphSuperpixelsByISF");

    iftAdjRel *A = NULL;
    if (iftIs3DImage(img))
        A = iftSpheric(1.0);
    else A = iftCircular(1.0);

    iftMImage *mimg = NULL;
    if (iftIsColorImage(img))
        mimg = iftImageToMImage(img, LABNorm_CSPACE);
    else mimg = iftImageToMImage(img, GRAY_CSPACE);


    iftImage *mask    = iftSelectImageDomain(mimg->xsize, mimg->ysize, mimg->zsize);
    iftIGraph *igraph = iftImplicitIGraph(mimg, mask, A); // minima of a basins manifold in that domain

    iftImage *seeds = iftAltMixedSampling(mimg, mask, input_n_cells);
    iftDestroyMImage(&mimg);
    iftDestroyImage(&mask);

    int final_n_iters = iftIGraphISF_Mean(igraph, seeds, alpha, beta, n_iters);
    printf("--> Final Number of Iterations (ISF):  %d\n", final_n_iters);

    // Smooth regions in the label map of igraph
    if (smooth_n_iters > 0) {
        iftIGraphSetWeightForRegionSmoothing(igraph, img);
        iftIGraphSmoothRegions(igraph, smooth_n_iters);
    }

    iftDestroyAdjRel(&A);
    iftDestroyImage(&seeds);

    return igraph;
}


iftImage *iftGenerateSuperpixelsByISF(const iftImage *img, int input_n_cells, float alpha, float beta,
                                      int n_iters, int smooth_n_iters, iftIGraph **out_igraph) {
    if (img == NULL)
        iftError("Image is NULL", "iftGenerateSuperpixelsByISF");

    iftIGraph *igraph = iftGenerateIGraphSuperpixelsByISF(img, input_n_cells, alpha, beta, n_iters, smooth_n_iters);

    // gets the superpixel image from the igraph
    iftImage *super_img = iftIGraphLabel(igraph);

    if (out_igraph == NULL)
        iftDestroyIGraph(&igraph);
    else *out_igraph = igraph;

    return super_img;
}


void iftIGraphISFRecomputeSeeds(iftIGraph *igraph, int *seed, int nseeds, iftSet **trees_for_removal, iftSet **new_seeds, float *new_seeds_flag)
{
    int      i, *center, p, q, s;
    iftVoxel u, v;
    float    distColor, distVoxel, distColorThres, distVoxelThres;

    /* Compute superpixel centers (i.e., the closest node to the center in the feature space) */

    center = iftIGraphSuperpixelCenters(igraph, seed, nseeds);


    /* Estimate distColorThres */

    distColorThres = 0.0;  distVoxelThres = 0.0;
    for (s=0; s < igraph->nnodes; s++) {
        p     = igraph->node[s].voxel;
        i     = igraph->label[igraph->root[p]]-1;
        q     = igraph->node[seed[i]].voxel;
        u = iftGetVoxelCoord(igraph->index,p);
        v = iftGetVoxelCoord(igraph->index,q);
        distColor = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
        distColorThres += distColor;
        distVoxel = iftVoxelDistance(u,v);
        distVoxelThres += distVoxel;
    }
    distColorThres /= igraph->nnodes;
    distVoxelThres /= igraph->nnodes;
    distColorThres = sqrtf(distColorThres);
    distVoxelThres = sqrtf(distVoxelThres);
    //distColorThres = 0.70*distColorThres; // converges faster
    //distVoxelThres = 0.70*distVoxelThres;

    /* Verify if the centers can be new seeds */

    *new_seeds_flag = 0;
    for (i=0; i < nseeds; i++) {
        p = igraph->node[seed[i]].voxel;
        q = igraph->node[center[i]].voxel;
        u = iftGetVoxelCoord(igraph->index,p);
        v = iftGetVoxelCoord(igraph->index,q);
        distColor = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
        distVoxel = iftVoxelDistance(u,v);

        if ((distColor > distColorThres)||(distVoxel > distVoxelThres)){
            seed[i] = center[i];
            iftInsertSet(new_seeds,center[i]);
            iftInsertSet(trees_for_removal,seed[i]);
            *new_seeds_flag+= 1;
        }
    }
    *new_seeds_flag /= (float)nseeds;

    iftFree(center);

}

void iftIGraphISFRecomputeSeedsUsingSpatialInformation(iftIGraph *igraph, int *seed, int nseeds, iftSet **trees_for_removal,
                                                       iftSet **new_seeds, float *new_seeds_flag, float **seed_features) {
    int      i, *center, p, q, s;
    iftVoxel u, v;
    float    distColor, distVoxel, distColorThres, distVoxelThres;
    //int    *area;
    //int     min_area;

    /* Compute superpixel centers (i.e., the closest node to the center in the feature space) */

    center = iftIGraphSuperpixelCentersUsingSpatialInformation(igraph, seed, nseeds, seed_features);

    //min_area  =  (int)((igraph->nnodes / (float)nseeds) / 5.0);
    //area = iftAllocIntArray(nseeds);

    /* Estimate distColorThres */

    distColorThres = 0.0;  distVoxelThres = 0.0;
    for (s=0; s < igraph->nnodes; s++) {
        p     = igraph->node[s].voxel;
        i     = igraph->label[igraph->root[p]]-1;
        q     = igraph->node[seed[i]].voxel;
        u = iftGetVoxelCoord(igraph->index,p);
        v = iftGetVoxelCoord(igraph->index,q);
        distColor = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
        distColorThres += distColor;
        distVoxel = iftVoxelDistance(u,v);
        distVoxelThres += distVoxel;
        //area[i]++;
    }
    distColorThres /= igraph->nnodes;
    distVoxelThres /= igraph->nnodes;
    distColorThres = sqrtf(distColorThres);
    distVoxelThres = sqrtf(distVoxelThres);

    //printf("distColorThres %f distVoxelThres %f\n", distColorThres, distVoxelThres);

    /* Verify if the centers can be new seeds */
    //distColorThres = 5.0;
    //distVoxelThres = 2.0;

    *new_seeds_flag = 0;
    for (i=0; i < nseeds; i++) {
        p = igraph->node[seed[i]].voxel;
        q = igraph->node[center[i]].voxel;
        u = iftGetVoxelCoord(igraph->index,p);
        v = iftGetVoxelCoord(igraph->index,q);
        distColor = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
        distVoxel = iftVoxelDistance(u,v);
        //if ((distColor > distColorThres) || (distVoxel > distVoxelThres) || area[i] < min_area) {
        if ((distColor > distColorThres) || (distVoxel > distVoxelThres)) {
            seed[i] = center[i];
            iftInsertSet(new_seeds,center[i]);
            iftInsertSet(trees_for_removal,seed[i]);
            *new_seeds_flag+= 1;
        }
    }
    *new_seeds_flag /= (float)nseeds;

    //iftFree(area);
    iftFree(center);
}

void iftNormIGraphFeatures(iftIGraph *igraph)
{
    int       p,s,i;
    double    module;

    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        module = 0.0;
        for (i=0; i < igraph->nfeats; i++){
            module += igraph->feat[p][i]*igraph->feat[p][i];
        }
        module = sqrt(module);
        if (module > IFT_EPSILON)
            for (i=0; i < igraph->nfeats; i++){
                igraph->feat[p][i] /= module;
            }
    }
}


/* Propagate forest attributes by using the maximum arc weight between
   adjacent nodes as clustering criterion and the path values of the
   nodes as priority function to resolve ties */

void iftIGraphClusterVoxelsByMaxArcWeight(iftIGraph *igraph, uchar pvalue_order)
{
    iftDHeap  *Q      = iftCreateDHeap(igraph->nnodes, igraph->pvalue);
    int       *node   = iftAllocIntArray(igraph->nnodes);
    float      dist;
    int        s,t,p,q;
    float     *max_arc_weight = iftAllocFloatArray(igraph->nnodes);
    iftSet    *adj;

    if (igraph->type != EXPLICIT)
        iftError("Graph must be explicit", "iftIGraphClusterVoxelsByMaxArcWeight");

    if (pvalue_order == IFT_DECREASING)
        iftSetRemovalPolicyDHeap(Q, IFT_MAXVALUE);

    /* Compute maximum arc-weight for each node and insert nodes in the
       priority queue for sorting */

    for (s=0; s < igraph->nnodes; s++) {
        max_arc_weight[s]=0.0;
        iftInsertDHeap(Q,s);
        p = igraph->node[s].voxel;
        for (adj = igraph->node[s].adj; adj != NULL; adj = adj->next) {
            t = adj->elem;
            q = igraph->node[t].voxel;
            dist  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
            if (dist > max_arc_weight[s])
                max_arc_weight[s]=dist;
        }
    }

    /* Sort nodes by path value */

    t = 0;
    while (!iftEmptyDHeap(Q)) {
        s = iftRemoveDHeap(Q);
        node[t]=s; t++;
    }


    /* Propagate forest attributes */

    for (q=0; q < (igraph->index->n); q++) {
        if (igraph->index->val[q] == IFT_NIL) { /* not in the graph */
            for (t=0; t < igraph->nnodes; t++) {
                s      = node[t];
                p     = igraph->node[s].voxel;
                dist  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
                if (dist <= max_arc_weight[s]){
                    igraph->label[q]   = igraph->label[p];
                    igraph->root[q]    = igraph->root[p];
                    igraph->pred[q]    = p;
                    igraph->pvalue[q]  = igraph->pvalue[p];
                    break;
                }
            }
        }
    }

    iftDestroyDHeap(&Q);
    iftFree(max_arc_weight);
    iftFree(node);

}


iftCompTree *iftIGraphCreateMaxTree(iftIGraph *igraph) {
    iftCompTree *ctree = (iftCompTree *) iftAlloc(1, sizeof(iftCompTree));
    iftImage *dad, *cmap, *tmp;
    int s, rs, t, rt, r, i, p, q;
    int Imax = (int) iftIGraphMaximumWeight(igraph);
    iftGQueue *Q;
    iftVoxel u, v;
    int *nsons = NULL;
    int *size = NULL;
    iftSet *adj;
    int *weight;
    iftAdjRel *A = igraph->A;

    dad = iftCreateImage(igraph->nnodes, 1, 1);
    size = iftAllocIntArray(igraph->nnodes);
    weight = iftAllocIntArray(igraph->nnodes);
    ctree->cmap = iftCreateImage(igraph->nnodes, 1, 1);
    cmap = ctree->cmap;
    ctree->root = IFT_NIL; /* Tree is empty */

    //Copy node weights as integers.
    for (s = 0; s < igraph->nnodes; s++) {
        weight[s] = (int) igraph->node[s].weight;
    }

    Q = iftCreateGQueue(Imax + 1, igraph->nnodes, weight);
    iftSetRemovalPolicy(Q, IFT_MAXVALUE);
    iftSetTieBreak(Q, LIFOBREAK);

    for (s = 0; s < cmap->n; s++) {
        dad->val[s] = IFT_NIL;
        cmap->val[s] = s;
        size[s] = 1;
        iftInsertGQueue(&Q, s);
    }

    switch (igraph->type) {

        case IMPLICIT:

            while (!iftEmptyGQueue(Q)) {
                s = iftRemoveGQueue(Q);
                p = igraph->node[s].voxel;
                rs = iftRepresentative(cmap, s);
                u = iftGetVoxelCoord(igraph->index, p);

                for (i = 1; i < A->n; i++) {
                    v = iftGetAdjacentVoxel(A, u, i);
                    if (iftValidVoxel(igraph->index, v)) {
                        q = iftGetVoxelIndex(igraph->index, v);
                        if (igraph->index->val[q] != IFT_NIL) {
                            t = igraph->index->val[q];
                            if (weight[s] == weight[t])   /* propagate on component */
                            {
                                if (Q->L.elem[t].color == IFT_GRAY) {
                                    iftRemoveGQueueElem(Q, t);
                                    cmap->val[t] = rs;
                                    size[rs] = size[rs] + 1;
                                    iftInsertGQueue(&Q, t);
                                }
                            }
                            else {
                                if (weight[s] < weight[t]) /* find current dad of rt */
                                {
                                    rt = iftRepresentative(cmap, t);
                                    r = iftAncestor(dad, cmap, rt);
                                    if (r == IFT_NIL)   /* rs is dad of the rt */
                                    {
                                        dad->val[rt] = rs;
                                    }
                                    else {
                                        if (weight[r] == weight[rs])  /* merge components */
                                        {
                                            if (r != rs) {
                                                if (size[rs] <= size[r]) {
                                                    cmap->val[rs] = r;
                                                    size[r] = size[r] + size[rs];
                                                    rs = r;
                                                }
                                                else {
                                                    cmap->val[r] = rs;
                                                    size[rs] = size[rs] + size[r];
                                                }
                                            }
                                        }
                                        else     /* weight[s] > weight[t] */
                                        {
                                            dad->val[r] = rs; /* rs is dad of r */
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            iftFree(size);
            iftDestroyGQueue(&Q);

            break;

        case EXPLICIT:

            while (!iftEmptyGQueue(Q)) {
                s = iftRemoveGQueue(Q);
                rs = iftRepresentative(cmap, s);
                for (adj = igraph->node[s].adj; adj != NULL; adj = adj->next) {
                    t = adj->elem;
                    if (weight[s] == weight[t])   /* propagate on component */
                    {
                        if (Q->L.elem[t].color == IFT_GRAY) {
                            iftRemoveGQueueElem(Q, t);
                            cmap->val[t] = rs;
                            size[rs] = size[rs] + 1;
                            iftInsertGQueue(&Q, t);
                        }
                    }
                    else {
                        if (weight[s] < weight[t]) /* find current dad of rt */
                        {
                            rt = iftRepresentative(cmap, t);
                            r = iftAncestor(dad, cmap, rt);
                            if (r == IFT_NIL)   /* rs is dad of the rt */
                            {
                                dad->val[rt] = rs;
                            }
                            else {
                                if (weight[r] == weight[rs])  /* merge components */
                                {
                                    if (r != rs) {
                                        if (size[rs] <= size[r]) {
                                            cmap->val[rs] = r;
                                            size[r] = size[r] + size[rs];
                                            rs = r;
                                        }
                                        else {
                                            cmap->val[r] = rs;
                                            size[rs] = size[rs] + size[r];
                                        }
                                    }
                                }
                                else     /* weight[s] > weight[t] */
                                {
                                    dad->val[r] = rs; /* rs is dad of r */
                                }
                            }
                        }
                    }
                }
            }

            iftFree(size);
            iftDestroyGQueue(&Q);
            break;

        case COMPLETE:
        default:
            iftError("Invalid graph type for tree creation", "iftIGraphCreateMaxTree");
    }

    /* Compress cmap map and count number of nodes */

    ctree->numnodes = 0;
    for (s = 0; s < cmap->n; s++) {
        if (dad->val[cmap->val[s]] != IFT_NIL)
            r = cmap->val[s];
        cmap->val[s] = iftRepresentative(cmap, s);

        if (cmap->val[s] == s)
            ctree->numnodes++;
    }

    /* Create and initialize nodes of the MaxTree. */

    ctree->node = (iftCompTreeNode *) iftAlloc(ctree->numnodes, sizeof(iftCompTreeNode));
    tmp = iftCreateImage(igraph->nnodes, 1, 1);
    for (s = 0; s < cmap->n; s++) {
        tmp->val[s] = IFT_NIL;
    }

    i = 0;
    for (s = 0; s < cmap->n; s++) {
        if (cmap->val[s] == s) {
            ctree->node[i].level = weight[s];
            ctree->node[i].comp = s;
            tmp->val[s] = i;
            ctree->node[i].dad = IFT_NIL;
            ctree->node[i].son = NULL;
            ctree->node[i].numsons = 0;
            ctree->node[i].size = 0;
            i++;
        }
    }

    iftFree(weight);

    /* Make the component map to point back to the maxtree. */

    for (s = 0; s < tmp->n; s++) {
        if (tmp->val[s] == IFT_NIL)
            tmp->val[s] = tmp->val[cmap->val[s]];
    }

    for (s = 0; s < cmap->n; s++) {
        cmap->val[s] = tmp->val[s];
    }
    iftDestroyImage(&tmp);

    /* Copy dad information to the maxtree and find its root */

    for (i = 0; i < ctree->numnodes; i++) {
        if (dad->val[ctree->node[i].comp] != IFT_NIL)
            ctree->node[i].dad = cmap->val[dad->val[ctree->node[i].comp]];
        else {
            ctree->node[i].dad = IFT_NIL;
            ctree->root = i;
        }
    }
    iftDestroyImage(&dad);

    /* Copy son information to the maxtree */

    nsons = iftAllocIntArray(ctree->numnodes);
    for (i = 0; i < ctree->numnodes; i++) {
        s = ctree->node[i].dad;
        if (s != IFT_NIL) {
            nsons[s]++;
        }
    }
    for (i = 0; i < ctree->numnodes; i++) {
        if (nsons[i] != 0) {
            ctree->node[i].son = iftAllocIntArray(nsons[i]);
        }
    }
    iftFree(nsons);

    for (i = 0; i < ctree->numnodes; i++) {
        s = ctree->node[i].dad;
        if (s != IFT_NIL) {
            ctree->node[s].son[ctree->node[s].numsons] = i;
            ctree->node[s].numsons++;
        }
    }

    /* Compute size of each node */

    for (s = 0; s < cmap->n; s++)
        ctree->node[cmap->val[s]].size++;


    return (ctree);
}


iftCompTree *iftIGraphCreateMinTree(iftIGraph *igraph) {
    iftCompTree *ctree = (iftCompTree *) iftAlloc(1, sizeof(iftCompTree));
    iftImage *dad, *cmap, *tmp;
    int s, rs, t, rt, r, i, p, q;
    int Imax = (int) iftIGraphMaximumWeight(igraph);
    iftGQueue *Q;
    iftVoxel u, v;
    int *nsons = NULL;
    int *size = NULL;
    iftSet *adj;
    int *weight;
    iftAdjRel *A = igraph->A;

    dad = iftCreateImage(igraph->nnodes, 1, 1);
    size = iftAllocIntArray(igraph->nnodes);
    weight = iftAllocIntArray(igraph->nnodes);
    ctree->cmap = iftCreateImage(igraph->nnodes, 1, 1);
    cmap = ctree->cmap;
    ctree->root = IFT_NIL; /* Tree is empty */

    //Copy node weights as integers.
    for (s = 0; s < igraph->nnodes; s++) {
        weight[s] = (int) igraph->node[s].weight;
    }

    Q = iftCreateGQueue(Imax + 1, igraph->nnodes, weight);
    iftSetTieBreak(Q, LIFOBREAK);

    for (s = 0; s < cmap->n; s++) {
        dad->val[s] = IFT_NIL;
        cmap->val[s] = s;
        size[s] = 1;
        iftInsertGQueue(&Q, s);
    }

    switch (igraph->type) {

        case IMPLICIT:

            while (!iftEmptyGQueue(Q)) {
                s = iftRemoveGQueue(Q);
                p = igraph->node[s].voxel;
                rs = iftRepresentative(cmap, s);
                u = iftGetVoxelCoord(igraph->index, p);

                for (i = 1; i < A->n; i++) {
                    v = iftGetAdjacentVoxel(A, u, i);
                    if (iftValidVoxel(igraph->index, v)) {
                        q = iftGetVoxelIndex(igraph->index, v);
                        if (igraph->index->val[q] != IFT_NIL) {
                            t = igraph->index->val[q];
                            if (weight[s] == weight[t])   /* propagate on component */
                            {
                                if (Q->L.elem[t].color == IFT_GRAY) {
                                    iftRemoveGQueueElem(Q, t);
                                    cmap->val[t] = rs;
                                    size[rs] = size[rs] + 1;
                                    iftInsertGQueue(&Q, t);
                                }
                            }
                            else {
                                if (weight[s] > weight[t]) /* find current dad of rt */
                                {
                                    rt = iftRepresentative(cmap, t);
                                    r = iftAncestor(dad, cmap, rt);
                                    if (r == IFT_NIL)   /* rs is dad of the rt */
                                    {
                                        dad->val[rt] = rs;
                                    }
                                    else {
                                        if (weight[r] == weight[rs])  /* merge components */
                                        {
                                            if (r != rs) {
                                                if (size[rs] <= size[r]) {
                                                    cmap->val[rs] = r;
                                                    size[r] = size[r] + size[rs];
                                                    rs = r;
                                                }
                                                else {
                                                    cmap->val[r] = rs;
                                                    size[rs] = size[rs] + size[r];
                                                }
                                            }
                                        }
                                        else     /* weight[s] < weight[t] */
                                        {
                                            dad->val[r] = rs; /* rs is dad of r */
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            iftFree(size);
            iftDestroyGQueue(&Q);

            break;

        case EXPLICIT:

            while (!iftEmptyGQueue(Q)) {
                s = iftRemoveGQueue(Q);
                rs = iftRepresentative(cmap, s);
                for (adj = igraph->node[s].adj; adj != NULL; adj = adj->next) {
                    t = adj->elem;
                    if (weight[s] == weight[t])   /* propagate on component */
                    {
                        if (Q->L.elem[t].color == IFT_GRAY) {
                            iftRemoveGQueueElem(Q, t);
                            cmap->val[t] = rs;
                            size[rs] = size[rs] + 1;
                            iftInsertGQueue(&Q, t);
                        }
                    }
                    else {
                        if (weight[s] > weight[t]) /* find current dad of rt */
                        {
                            rt = iftRepresentative(cmap, t);
                            r = iftAncestor(dad, cmap, rt);
                            if (r == IFT_NIL)   /* rs is dad of the rt */
                            {
                                dad->val[rt] = rs;
                            }
                            else {
                                if (weight[r] == weight[rs])  /* merge components */
                                {
                                    if (r != rs) {
                                        if (size[rs] <= size[r]) {
                                            cmap->val[rs] = r;
                                            size[r] = size[r] + size[rs];
                                            rs = r;
                                        }
                                        else {
                                            cmap->val[r] = rs;
                                            size[rs] = size[rs] + size[r];
                                        }
                                    }
                                }
                                else     /* weight[s] < weight[t] */
                                {
                                    dad->val[r] = rs; /* rs is dad of r */
                                }
                            }
                        }
                    }
                }
            }

            iftFree(size);
            iftDestroyGQueue(&Q);
            break;

        case COMPLETE:
        default:
            iftError("Invalid graph type for tree creation", "iftIGraphCreateMaxTree");
    }

    /* Compress cmap map and count number of nodes */

    ctree->numnodes = 0;
    for (s = 0; s < cmap->n; s++) {
        if (dad->val[cmap->val[s]] != IFT_NIL)
            r = cmap->val[s];
        cmap->val[s] = iftRepresentative(cmap, s);

        if (cmap->val[s] == s)
            ctree->numnodes++;
    }

    /* Create and initialize nodes of the MinTree. */

    ctree->node = (iftCompTreeNode *) iftAlloc(ctree->numnodes, sizeof(iftCompTreeNode));
    tmp = iftCreateImage(igraph->nnodes, 1, 1);
    for (s = 0; s < cmap->n; s++) {
        tmp->val[s] = IFT_NIL;
    }

    i = 0;
    for (s = 0; s < cmap->n; s++) {
        if (cmap->val[s] == s) {
            ctree->node[i].level = weight[s];
            ctree->node[i].comp = s;
            tmp->val[s] = i;
            ctree->node[i].dad = IFT_NIL;
            ctree->node[i].son = NULL;
            ctree->node[i].numsons = 0;
            ctree->node[i].size = 0;
            i++;
        }
    }

    iftFree(weight);

    /* Make the component map to point back to the mintree. */

    for (s = 0; s < tmp->n; s++) {
        if (tmp->val[s] == IFT_NIL)
            tmp->val[s] = tmp->val[cmap->val[s]];
    }

    for (s = 0; s < cmap->n; s++) {
        cmap->val[s] = tmp->val[s];
    }
    iftDestroyImage(&tmp);

    /* Copy dad information to the mintree and find its root */

    for (i = 0; i < ctree->numnodes; i++) {
        if (dad->val[ctree->node[i].comp] != IFT_NIL)
            ctree->node[i].dad = cmap->val[dad->val[ctree->node[i].comp]];
        else {
            ctree->node[i].dad = IFT_NIL;
            ctree->root = i;
        }
    }
    iftDestroyImage(&dad);

    /* Copy son information to the mintree */

    nsons = iftAllocIntArray(ctree->numnodes);
    for (i = 0; i < ctree->numnodes; i++) {
        s = ctree->node[i].dad;
        if (s != IFT_NIL) {
            nsons[s]++;
        }
    }
    for (i = 0; i < ctree->numnodes; i++) {
        if (nsons[i] != 0) {
            ctree->node[i].son = iftAllocIntArray(nsons[i]);
        }
    }
    iftFree(nsons);

    for (i = 0; i < ctree->numnodes; i++) {
        s = ctree->node[i].dad;
        if (s != IFT_NIL) {
            ctree->node[s].son[ctree->node[s].numsons] = i;
            ctree->node[s].numsons++;
        }
    }

    /* Compute size of each node */

    for (s = 0; s < cmap->n; s++)
        ctree->node[cmap->val[s]].size++;


    return (ctree);
}


iftIGraph * iftIGraphResetFeatureSpace(iftIGraph *igraph, int nFeatures) {
    iftIGraph *new_igraph = (iftIGraph *)iftAlloc(1,sizeof(iftIGraph));
    int        p;

    new_igraph->nnodes  = igraph->nnodes;
    new_igraph->node    = (iftINode *)iftAlloc(igraph->nnodes,sizeof(iftINode));
    new_igraph->index   = iftCopyImage(igraph->index);
    new_igraph->nfeats  = nFeatures;

    //Check: Really igraph->nnodes?
    new_igraph->feat    = (float **)iftAlloc(igraph->nnodes,sizeof(float *));
    for (p=0; p < igraph->nnodes; p++) {
        new_igraph->feat[p] = iftAllocFloatArray(nFeatures);
    }
    //Check: Really igraph->nnodes?

    //Number of elements: nnodes?
    iftCopyIntArray(new_igraph->label, igraph->label, igraph->nnodes);
    iftCopyIntArray(new_igraph->marker, igraph->marker, igraph->nnodes);
    iftCopyIntArray(new_igraph->root, igraph->root, igraph->nnodes);
    iftCopyIntArray(new_igraph->pred, igraph->pred, igraph->nnodes);
    iftCopyDoubleArray(new_igraph->pvalue, igraph->pvalue, igraph->nnodes);

    return new_igraph;
}

float iftIGraphMaximumFeatureValue (iftIGraph * igraph, int feature) {
    float max;
    int s, p;

    max = IFT_INFINITY_FLT_NEG;
    for (s = 0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        if (igraph->feat[p][feature] > max) {
            max = igraph->feat[p][feature];
        }
    }

    return max;
}

float iftIGraphMaximumWeight (iftIGraph * igraph) {
    float max;
    int s;

    max = IFT_INFINITY_FLT_NEG;
    for (s = 0; s < igraph->nnodes; s++) {
        if (igraph->node[s].weight > max) {
            max = igraph->node[s].weight;
        }
    }

    return(max);
}

float iftIGraphMinimumWeight (iftIGraph * igraph) {
    float min;
    int s;

    min = IFT_INFINITY_FLT;
    for (s = 0; s < igraph->nnodes; s++) {
        if (igraph->node[s].weight < min) {
            min = igraph->node[s].weight;
        }
    }

    return(min);
}

void iftIGraphDualWaterGray(iftIGraph *igraph, int criterion, int thres)
{
    int s,p;
    float * marker = NULL;
    int value_to_be_restored = 0; //Used to restore the weight and path value maps to
    //the value they had before the marker creation, in
    //the case of the HEIGHT criterion.

    switch (criterion) {
        case HEIGHT:
            marker = iftIGraphWeightMarkerByHeight(igraph,thres);
            value_to_be_restored = thres+1;

            break;

        case AREA:
            marker = iftIGraphAreaOpen(igraph, thres);
            value_to_be_restored = 1;

            break;

        case VOLUME:
            marker = iftIGraphVolumeOpen(igraph, thres);
            value_to_be_restored = 1;

            break;

        default:
            iftError("Invalid criterion", "iftIGraphDualWaterGray");
            break;
    }

    iftIGraphInfRec(igraph,marker);


    for (s=0; s < igraph->nnodes; s++) {
        igraph->node[s].weight -= value_to_be_restored;
        p = igraph->node[s].voxel;
        igraph->pvalue[p]      -= value_to_be_restored;
    }
    iftFree(marker);

}

void iftIGraphWaterGray(iftIGraph *igraph, int criterion, int thres)
{
    float *marker = NULL;

    switch (criterion) {
        case HEIGHT: /*This case actually refers to the image depth, which can be
                 interpreted as a negative height.*/

            marker = iftIGraphWeightMarkerByDepth(igraph,thres);
            break;

        case AREA:
            marker = iftIGraphAreaClose(igraph, thres);

            break;

        case VOLUME:
            marker = iftIGraphVolumeClose(igraph, thres);

            break;

        default:
            iftError("Invalid criterion", "iftIGraphWaterGray");
            break;
    }

    iftIGraphSupRec(igraph,marker);
    iftFree(marker);
}

iftImage *iftIGraphWeightMinima(iftIGraph *igraph)
{
    iftImage *mask=iftCreateImage(igraph->index->xsize,igraph->index->ysize,igraph->index->zsize);
    int s,p;

    //iftIGraphWaterGrayByDepth(igraph,1);
    iftIGraphWaterGray(igraph, HEIGHT, 1);
    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        if (igraph->root[p]==p)
            mask->val[p]=255;
    }

    return(mask);
}

iftImage *iftIGraphWeightMaxima(iftIGraph *igraph)
{
    iftImage *mask=iftCreateImage(igraph->index->xsize,igraph->index->ysize,igraph->index->zsize);
    int s,p;

    iftIGraphDualWaterGray(igraph, HEIGHT, 1);
    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        if (igraph->root[p]==p)
            mask->val[p]=255;
    }

    return(mask);
}

void iftIGraphCopyFeatureToWeight(iftIGraph *igraph, int feature_index)
{
    int s, p;

    if ((feature_index < 0)||(feature_index > igraph->nfeats-1))
        iftError("Invalid feature", "iftIGraphCopyFeatureToWeight");

    for (s = 0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        igraph->node[s].weight = igraph->feat[p][feature_index];
    }

}

float iftIGraphMaximumFeatureDist(iftIGraph *igraph)
{
    float      maxdist= IFT_INFINITY_FLT_NEG,dist;
    int        i, p, q, s, t;
    iftAdjRel *A = igraph->A;
    iftVoxel   u, v;
    iftImage  *index = igraph->index;
    iftSet    *adj;

    switch (igraph->type) {

        case IMPLICIT:

            for (s = 0; s < igraph->nnodes; s++) {
                p = igraph->node[s].voxel;
                u = iftGetVoxelCoord(index,p);
                for (i=1; i < A->n; i++) {
                    v = iftGetAdjacentVoxel(A,u,i);
                    if (iftValidVoxel(index,v)){
                        q = iftGetVoxelIndex(index,v);
                        t = index->val[q];
                        if (t != IFT_NIL){
                            dist  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
                            if (dist > maxdist)
                                maxdist = dist;
                        }
                    }
                }
            }

            break;

        case EXPLICIT:

            for (s = 0; s < igraph->nnodes; s++) {
                p = igraph->node[s].voxel;
                for (adj=igraph->node[s].adj; adj != NULL; adj = adj->next) {
                    t = adj->elem;
                    q = igraph->node[t].voxel;
                    dist  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
                    if (dist > maxdist)
                        maxdist = dist;
                }
            }

            break;

        case COMPLETE:

            for (s = 0; s < igraph->nnodes; s++) {
                p = igraph->node[s].voxel;
                for (t = s+1; t < igraph->nnodes; t++) {
                    q = igraph->node[t].voxel;
                    dist  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
                    if (dist > maxdist)
                        maxdist = dist;
                }
            }

            break;

        default:
            iftError("Invadid graph type", "iftIGraphMaximumFeatureDist");
    }

    return(maxdist);
}

int iftIGraphRootVoxel(iftIGraph *igraph, int q)
{
    if (igraph->root[q]==q)
        return(q);
    else
        return(igraph->root[q]=iftIGraphRootVoxel(igraph, igraph->root[q]));
}

int iftIGraphEnumerateRootsAndPropagateLabels(iftIGraph *igraph)
{
    int nroots=0, i, s, p;

    for (s=0, i=1, nroots=0; s < igraph->nnodes; s++)
    {
        p = igraph->node[s].voxel;
        if (igraph->root[p]==p) {
            igraph->label[p]=i;
            i++; nroots++;
        }
    }

    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        if (igraph->root[p]!=p) {
            igraph->label[p] = igraph->label[igraph->root[p]];
        }
    }

    return(nroots);
}

iftIGraph *iftIGraphMST(iftIGraph *igraph)
{
    int        s, t, i, p, q;
    iftVoxel   u, v;
    float      dist;
    iftSet    *adj;
    iftIGraph *mst = iftCopyIGraph(igraph);
    iftDHeap  *Q   = iftCreateDHeap(mst->nnodes, mst->pvalue);

    if (igraph->type==IMPLICIT)
        iftDestroyAdjRel(&mst->A);
    if (igraph->type==EXPLICIT)
        for (s=0; s < mst->nnodes; s++) {
            iftDestroySet(&mst->node[s].adj);
        }
    mst->type = EXPLICIT;

    for (s=0; s < mst->nnodes; s++) {
        p=mst->node[s].voxel;
        mst->pvalue[p] = IFT_INFINITY_DBL;
        mst->label[p]  = mst->marker[p] = 0;
        mst->root[p]   = p;
        mst->pred[p]   = IFT_NIL;
    }
    p=mst->node[0].voxel;
    mst->pvalue[p]=0;
    iftInsertDHeap(Q,0);

    switch (igraph->type) {

        case IMPLICIT:

            while(!iftEmptyDHeap(Q)){

                s  = iftRemoveDHeap(Q);
                p  = mst->node[s].voxel;
                u  = iftGetVoxelCoord(mst->index,p);

                if (mst->pred[p] != IFT_NIL){
                    q = mst->pred[p];
                    t = mst->index->val[q];
                    iftInsertSet(&mst->node[t].adj,s);
                    iftInsertSet(&mst->node[s].adj,t);
                }

                for (i=1; i < igraph->A->n; i++) {
                    v = iftGetAdjacentVoxel(igraph->A,u,i);
                    if (iftValidVoxel(mst->index,v)){
                        q = iftGetVoxelIndex(mst->index,v);
                        t = mst->index->val[q];
                        if ((t != IFT_NIL) && (Q->color[t] != IFT_BLACK)){
                            dist  =
                                    iftFeatDistance(mst->feat[p],mst->feat[q],mst->nfeats);
                            if (dist < mst->pvalue[q]){
                                mst->pvalue[q] = dist;
                                mst->pred[q]   = p;
                                if(Q->color[t] == IFT_WHITE)
                                    iftInsertDHeap(Q, t);
                                else
                                    iftGoUpDHeap(Q, Q->pos[t]);
                            }
                        }
                    }
                }
            }

            break;

        case EXPLICIT:

            while(!iftEmptyDHeap(Q)){

                s  = iftRemoveDHeap(Q);
                p  = mst->node[s].voxel;

                if (mst->pred[p] != IFT_NIL){
                    q = mst->pred[p];
                    t = mst->index->val[q];
                    iftInsertSet(&mst->node[t].adj,s);
                    iftInsertSet(&mst->node[s].adj,t);
                }

                for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
                    t = adj->elem;
                    q = mst->node[t].voxel;
                    if (Q->color[t] != IFT_BLACK){
                        dist  =
                                iftFeatDistance(mst->feat[p],mst->feat[q],mst->nfeats);
                        if (dist < mst->pvalue[q]){
                            mst->pvalue[q] = dist;
                            mst->pred[q]   = p;
                            if(Q->color[t] == IFT_WHITE)
                                iftInsertDHeap(Q, t);
                            else
                                iftGoUpDHeap(Q, Q->pos[t]);
                        }
                    }
                }
            }

            break;

        case COMPLETE:

            while(!iftEmptyDHeap(Q)){

                s  = iftRemoveDHeap(Q);
                p  = mst->node[s].voxel;

                if (mst->pred[p] != IFT_NIL){
                    q = mst->pred[p];
                    t = mst->index->val[q];
                    iftInsertSet(&mst->node[t].adj,s);
                    iftInsertSet(&mst->node[s].adj,t);
                }

                for (t=0; t < igraph->nnodes; t++) {
                    if (Q->color[t] != IFT_BLACK){
                        q = igraph->node[t].voxel;
                        dist  =
                                iftFeatDistance(mst->feat[p],mst->feat[q],mst->nfeats);
                        if (dist < mst->pvalue[q]){
                            mst->pvalue[q] = dist;
                            mst->pred[q]   = p;
                            if(Q->color[t] == IFT_WHITE)
                                iftInsertDHeap(Q, t);
                            else
                                iftGoUpDHeap(Q, Q->pos[t]);
                        }
                    }
                }
            }

            break;

    }

    iftDestroyDHeap(&Q);

    return(mst);
}

iftDHeap *iftIGraphInitDiffWatershed(iftIGraph *igraph, double *pvalue) {
    iftDHeap *Q      = iftCreateDHeap(igraph->nnodes, pvalue);

    iftIGraphResetWatershed(igraph, Q);

    return Q;
}

void iftIGraphResetWatershed(iftIGraph *igraph, iftDHeap *Q) {
    double *pvalue = Q->value;

    for(int s = 0; s < igraph->nnodes; s++) {
        int p           = igraph->node[s].voxel;
        pvalue[s]       = igraph->pvalue[p] = IFT_INFINITY_DBL;
        igraph->pred[p] = IFT_NIL;
        igraph->root[p] = p;
        igraph->label[p] = 0;
    }

    iftResetDHeap(Q);
}

void iftIGraphDiffWatershed(iftIGraph *igraph, iftLabeledSet *seeds, iftSet *trees_for_removal, iftDHeap *Q)
{
    double      tmp;
    int        s, t, i, p, q;
    iftVoxel   u, v;
    double    *pvalue = Q->value;
    iftSet    *adj=NULL, *frontier_nodes=NULL;
    iftLabeledSet *S = NULL;

    if (trees_for_removal != NULL)
        frontier_nodes = iftIGraphTreeRemoval(igraph, &trees_for_removal, pvalue, IFT_INFINITY_DBL);

    S = seeds;
    while (S != NULL) { /* insert seeds in the priority queue Q with cost zero */
        s = S->elem;
        p = igraph->node[s].voxel;
        pvalue[s] = igraph->pvalue[p] = 0;
        igraph->root[p] = p;
        igraph->pred[p] = IFT_NIL;
        igraph->label[p] = S->label;

        iftInsertDHeap(Q,s);

        S = S->next;
    }

    while (frontier_nodes != NULL) { /* insert frontier nodes in Q to represent the previous seeds */
        s = iftRemoveSet(&frontier_nodes);

        if (Q->color[s] == IFT_WHITE){
            iftInsertDHeap(Q,s);
        }
    }

    switch(igraph->type) {

        case IMPLICIT:


            while (!iftEmptyDHeap(Q)) {
                s = iftRemoveDHeap(Q);
                p = igraph->node[s].voxel;

                igraph->pvalue[p] = pvalue[s];
                u = iftGetVoxelCoord(igraph->index,p);
                for (i=1; i < igraph->A->n; i++) {
                    v = iftGetAdjacentVoxel(igraph->A,u,i);
                    if (iftValidVoxel(igraph->index,v)){
                        q   = iftGetVoxelIndex(igraph->index,v);

                        t   = igraph->index->val[q];

//                        if ((t != IFT_NIL) && (Q->color[t] != IFT_BLACK)) {
                        if (t != IFT_NIL) {

                            tmp = iftMax(pvalue[s], igraph->node[t].weight);

                            if (tmp < pvalue[t] || igraph->pred[q] == p)
                            {
                                pvalue[t]            = tmp;

                                igraph->root[q]      = igraph->root[p];
                                igraph->label[q]     = igraph->label[p];
                                igraph->pred[q]      = p;

                                if (Q->color[t] == IFT_GRAY){
                                    iftGoUpDHeap(Q, Q->pos[t]);
                                } else {
                                    iftInsertDHeap(Q,t);
                                }
                            }
                        }
                    }
                }
            }

            break;

        case EXPLICIT:

            while (!iftEmptyDHeap(Q)) {
                s = iftRemoveDHeap(Q);
                p = igraph->node[s].voxel;

                igraph->pvalue[p] = pvalue[s];

                for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
                    t   = adj->elem;
                    q   = igraph->node[t].voxel;

//                    if(Q->color[t] != IFT_BLACK) {
                    tmp = iftMax(pvalue[s], igraph->node[t].weight);

                    if (tmp < pvalue[t] || igraph->pred[q] == p) {
                        pvalue[t] = tmp;
                        igraph->root[q] = igraph->root[p];
                        igraph->label[q] = igraph->label[p];
                        igraph->pred[q] = p;
                        if (Q->color[t] == IFT_GRAY) {
                            iftGoUpDHeap(Q, Q->pos[t]);
                        } else {
                            iftInsertDHeap(Q, t);
                        }
                    }
//                    }
                }
            }
            break;

        case COMPLETE:
            iftError("Not implemented for complete graphs", "iftIGraphDiffWatershed");
    }

    iftResetDHeap(Q);

    iftDestroySet(&adj);
    iftDestroySet(&frontier_nodes);
}

void iftIGraphDiffOrientedWatershed(iftIGraph *igraph, iftImage *orien, iftLabeledSet *seeds, iftSet *trees_for_removal,
                                    iftDHeap *Q, char orientation)
{
    double      tmp;
    int        s, t, i, p, q;
    iftVoxel   u, v;
    double    *pvalue = Q->value;
    iftSet    *adj=NULL, *frontier_nodes=NULL, *T=NULL;
    iftLabeledSet *S = NULL;
    bool is_first_iteration = true;

    if (trees_for_removal != NULL)
        frontier_nodes = iftIGraphTreeRemoval(igraph, &trees_for_removal, pvalue, IFT_INFINITY_DBL);

    /* Checking if we are in the first iteration, in order to fix the label map when we are not */
    for(s = 0; s < igraph->nnodes && is_first_iteration; s++) {
        p           = igraph->node[s].voxel;
        is_first_iteration = igraph->pred[p] == IFT_NIL;
    }

    fprintf(stderr,"Is first oriented DIFT iteration: %d\n", is_first_iteration);

    S = seeds;
    while (S != NULL) { /* insert seeds in the priority queue Q with cost zero */
        s = S->elem;
        p = igraph->node[s].voxel;
        pvalue[s] = igraph->pvalue[p] = 0;
        igraph->root[p] = p;
        igraph->pred[p] = IFT_NIL;
        igraph->label[p] = S->label;
        iftInsertDHeap(Q,s);

        S = S->next;
    }

    while (frontier_nodes != NULL) { /* insert frontier nodes in Q to represent the previous seeds */
        s = iftRemoveSet(&frontier_nodes);

        if (Q->color[s] == IFT_WHITE){
            iftInsertDHeap(Q,s);
        }
    }


    switch(igraph->type) {

        case IMPLICIT:


            while (!iftEmptyDHeap(Q)) {
                s = iftRemoveDHeap(Q);
                p = igraph->node[s].voxel;

                igraph->pvalue[p] = pvalue[s];
                u = iftGetVoxelCoord(igraph->index,p);
                for (i=1; i < igraph->A->n; i++) {
                    v = iftGetAdjacentVoxel(igraph->A,u,i);
                    if (iftValidVoxel(igraph->index,v)){
                        q   = iftGetVoxelIndex(igraph->index,v);

                        t   = igraph->index->val[q];

                        if ((t != IFT_NIL) && (Q->color[t] != IFT_BLACK)) {

                            if (orientation == 1){ /* object brighter than background */
                                if (igraph->label[p]!=0){
                                    if (orien->val[p] > orien->val[q])
                                        tmp = iftMax(pvalue[s], igraph->node[t].weight* 1.5);
                                    else
                                        tmp = iftMax(pvalue[s], igraph->node[t].weight* 0.5);
                                } else {
                                    if (orien->val[p] < orien->val[q])
                                        tmp = iftMax(pvalue[s], igraph->node[t].weight* 1.5);
                                    else
                                        tmp = iftMax(pvalue[s], igraph->node[t].weight* 0.5);
                                }
                            } else { /* object darker than background */
                                if (igraph->label[p]!=0){
                                    if (orien->val[p] < orien->val[q])
                                        tmp = iftMax(pvalue[s], igraph->node[t].weight* 1.5);
                                    else
                                        tmp = iftMax(pvalue[s], igraph->node[t].weight* 0.5);
                                } else {
                                    if (orien->val[p] > orien->val[q])
                                        tmp = iftMax(pvalue[s], igraph->node[t].weight* 1.5);
                                    else
                                        tmp = iftMax(pvalue[s], igraph->node[t].weight* 0.5);
                                }
                            }

                            if (tmp < pvalue[t])
                            {
                                pvalue[t]            = tmp;

                                if (igraph->label[p] != igraph->label[q]){ /* voxels that have changed labels */
                                    iftInsertSet(&T, q);
                                }

                                igraph->root[q]      = igraph->root[p];
                                igraph->label[q]     = igraph->label[p];
                                igraph->pred[q]      = p;

                                if (Q->color[t] == IFT_GRAY){
                                    iftGoUpDHeap(Q, Q->pos[t]);
                                } else {
                                    iftInsertDHeap(Q,t);
                                }
                            }
                        }
                    }
                }
            }

            break;

        case EXPLICIT:

            while (!iftEmptyDHeap(Q)) {
                s = iftRemoveDHeap(Q);
                p = igraph->node[s].voxel;

                igraph->pvalue[p] = pvalue[s];

                for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
                    t = adj->elem;
                    if (Q->color[t] != IFT_BLACK) {
                        q = igraph->node[t].voxel;

                        if (orientation == 1) { /* object brighter than background */
                            if (igraph->label[p] != 0) {
                                if (orien->val[p] > orien->val[q])
                                    tmp = iftMax(pvalue[s], igraph->node[t].weight* 1.5);
                                else
                                    tmp = iftMax(pvalue[s], igraph->node[t].weight* 0.5);
                            } else {
                                if (orien->val[p] < orien->val[q])
                                    tmp = iftMax(pvalue[s], igraph->node[t].weight* 1.5);
                                else
                                    tmp = iftMax(pvalue[s], igraph->node[t].weight* 0.5);
                            }
                        } else { /* object darker than background */
                            if (igraph->label[p] != 0) {
                                if (orien->val[p] < orien->val[q])
                                    tmp = iftMax(pvalue[s], igraph->node[t].weight* 1.5);
                                else
                                    tmp = iftMax(pvalue[s], igraph->node[t].weight* 0.5);
                            } else {
                                if (orien->val[p] > orien->val[q])
                                    tmp = iftMax(pvalue[s], igraph->node[t].weight* 1.5);
                                else
                                    tmp = iftMax(pvalue[s], igraph->node[t].weight* 0.5);
                            }
                        }

                        if (tmp < pvalue[t]) {
                            pvalue[t] = tmp;

                            if (igraph->label[p] != igraph->label[q]){ /* voxels that have changed labels */
                                iftInsertSet(&T, q);
                            }

                            igraph->root[q] = igraph->root[p];
                            igraph->label[q] = igraph->label[p];
                            igraph->pred[q] = p;
                            if (Q->color[t] == IFT_GRAY) {
                                iftGoUpDHeap(Q, Q->pos[t]);
                            } else {
                                iftInsertDHeap(Q, t);
                            }
                        }
                    }
                }
            }
            break;

        case COMPLETE:
            iftError("Not implemented for complete graphs", "iftIGraphDiffOrientedWatershed");
    }

    iftResetDHeap(Q);

    /* Fixing the label map after the first iteration */
    if(!is_first_iteration)
        iftIGraphFixLabelRootMap(igraph, &T);

    iftDestroySet(&adj);
    iftDestroySet(&T);
    iftDestroySet(&frontier_nodes);
}



iftIGraph *iftCopyIGraph(iftIGraph *igraph)
{
    iftIGraph *igraph_copy = (iftIGraph *)iftAlloc(1,sizeof(iftIGraph));
    int        s, p, i;

    igraph_copy->nnodes  = igraph->nnodes;
    igraph_copy->nfeats  = igraph->nfeats;
    igraph_copy->type    = igraph->type;
    igraph_copy->node    = (iftINode *)iftAlloc(igraph->nnodes,sizeof(iftINode));
    igraph_copy->index   = iftCopyImage(igraph->index);
    if (igraph->A != NULL)
        igraph_copy->A       = iftCopyAdjacency(igraph->A);

    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        igraph_copy->node[s].voxel  = p;
        igraph_copy->node[s].weight = igraph->node[s].weight;
        if (igraph->node[s].adj != NULL)
            igraph_copy->node[s].adj = iftSetCopy(igraph->node[s].adj);
    }
    igraph_copy->feat    = (float **)iftAlloc(igraph->index->n,sizeof(float *));
    for (p=0; p < igraph->index->n; p++) {
        igraph_copy->feat[p] = iftAllocFloatArray(igraph->nfeats);
        for (i=0; i < igraph->nfeats; i++)
            igraph_copy->feat[p][i] = igraph->feat[p][i];
    }

    igraph_copy->label   = iftAllocIntArray(igraph->index->n);
    iftCopyIntArray(igraph_copy->label,igraph->label,igraph->index->n);
    igraph_copy->marker  = iftAllocIntArray(igraph->index->n);
    iftCopyIntArray(igraph_copy->marker,igraph->marker,igraph->index->n);
    igraph_copy->root    = iftAllocIntArray(igraph->index->n);
    iftCopyIntArray(igraph_copy->root,igraph->root,igraph->index->n);
    igraph_copy->pred    = iftAllocIntArray(igraph->index->n);
    iftCopyIntArray(igraph_copy->pred, igraph->pred,igraph->index->n);
    igraph_copy->pvalue  = iftAllocDoubleArray(igraph->index->n);
    iftCopyDoubleArray(igraph_copy->pvalue,igraph->pvalue,igraph->index->n);

    return(igraph_copy);
}

void iftIGraphDomes(iftIGraph *igraph)
{
    int      s, t, i, p, q;
    iftVoxel u, v;
    float    dist, maxval=0.0, minval= IFT_INFINITY_FLT;
    iftSet  *adj1, *adj2;
    char     is_neighbor;

    switch (igraph->type) {

        case IMPLICIT:
            for (s=0; s < igraph->nnodes; s++) {
                p = igraph->node[s].voxel;
                u = iftGetVoxelCoord(igraph->index,p);
                igraph->node[s].weight=0.0;
                for (i=1; i < igraph->A->n; i++) {
                    v = iftGetAdjacentVoxel(igraph->A,u,i);
                    if (iftValidVoxel(igraph->index,v)){
                        q = iftGetVoxelIndex(igraph->index,v);
                        if (igraph->index->val[q] != IFT_NIL){
                            dist  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
                            igraph->node[s].weight += dist;
                        }
                    }
                }
                if (igraph->node[s].weight > maxval)
                    maxval = igraph->node[s].weight;
                if (igraph->node[s].weight < minval)
                    minval = igraph->node[s].weight;
            }
            break;

        case EXPLICIT:

            for (s=0; s < igraph->nnodes; s++) {
                p = igraph->node[s].voxel;
                igraph->node[s].weight=0.0;
                for (adj1 = igraph->node[s].adj; adj1 != NULL; adj1 = adj1->next){
                    q     = igraph->node[adj1->elem].voxel;
                    dist  = iftFeatDistance(igraph->feat[p],igraph->feat[q],igraph->nfeats);
                    igraph->node[s].weight += dist;
                }
                if (igraph->node[s].weight > maxval)
                    maxval = igraph->node[s].weight;
                if (igraph->node[s].weight < minval)
                    minval = igraph->node[s].weight;
            }

            break;

        case COMPLETE:
        default:
            iftError("Invalid graph type for pdf estimation", "iftIGraphDomes");
    }

    if (maxval > minval) {
        for (s=0; s < igraph->nnodes; s++) {
            igraph->node[s].weight = (IFT_MAXWEIGHT - 1.0) *
                                     ((maxval - igraph->node[s].weight) / (maxval-minval)) + 1.0;
        }
    }


    /* make the graph symmetric on plateaus */

    if (igraph->type == EXPLICIT) {

        for (s=0; s < igraph->nnodes; s++) {
            for (adj1 = igraph->node[s].adj; adj1 != NULL; adj1 = adj1->next){
                t = adj1->elem;
                if (fabs(igraph->node[s].weight-igraph->node[t].weight) < IFT_EPSILON){
                    is_neighbor = 0;
                    for (adj2 = igraph->node[t].adj; (adj2 != NULL)&&(is_neighbor==0); adj2 = adj2->next){
                        if (s == adj2->elem){
                            is_neighbor = 1;
                        }
                    }
                    if (is_neighbor==0) {
                        iftInsertSet(&igraph->node[t].adj,s);
                    }
                }
            }
        }
    }

}

void iftIGraphBasins(iftIGraph *igraph)
{
    int      s;

    iftIGraphDomes(igraph);
    for (s=0; s < igraph->nnodes; s++) {
        igraph->node[s].weight = IFT_MAXWEIGHT - igraph->node[s].weight;
    }

}

void iftIGraphSetWeight(iftIGraph *igraph, iftImage *weight)
{

    for (int s=0; s < igraph->nnodes; s++) {
        int p = igraph->node[s].voxel;
        igraph->node[s].weight = weight->val[p];
    }

}


void iftIGraphSetFWeight(iftIGraph *igraph, iftFImage *weight)
{
    for (int s=0; s < igraph->nnodes; s++) {
        int p = igraph->node[s].voxel;
        igraph->node[s].weight = weight->val[p];
    }
}

void iftIGraphSetWeightForRegionSmoothing(iftIGraph *igraph, const iftImage *img)
{
    iftAdjRel *A;

    if (iftIs3DImage(img))
        A = iftSpheric(sqrtf(3.0));
    else
        A = iftCircular(sqrtf(2.0));

    iftImage  *grad      = iftImageGradientMagnitude(img,A);
    iftDestroyAdjRel(&A);

    iftFImage *weight    = iftSmoothWeightImage(grad,0.5);

    iftIGraphSetFWeight(igraph, weight);

    iftDestroyImage(&grad);
    iftDestroyFImage(&weight);
}


void iftIGraphSmoothRegions(iftIGraph *igraph, int num_smooth_iterations)
{
    iftImage  *prev_label,  *next_label;
    iftFImage *prev_weight, *next_weight, *norm_factor, *weight;
    float     *sum, max_membership;
    int        l, i, p, q, r, max_label, iter;
    iftVoxel   u, v;
    iftAdjRel *A = igraph->A;
    iftSet    *prev_frontier = NULL, *next_frontier = NULL, *S = NULL;
    iftBMap   *inFrontier;

    if (igraph->type != IMPLICIT)
        iftError("For implicit graphs only", "iftIGraphSmoothRegions");

    /* Initialization */

    prev_label  = iftIGraphLabel(igraph);
    next_label  = iftCopyImage(prev_label);
    weight      = iftIGraphWeight(igraph);
    norm_factor = iftWeightNormFactor(weight,A);
    inFrontier  = iftCreateBMap(prev_label->n);

    int prev_label_max_val = iftMaximumValue(prev_label);
    sum         = iftAllocFloatArray(prev_label_max_val + 1);
    prev_weight = iftCreateFImage(prev_label->xsize, prev_label->ysize, prev_label->zsize);
    next_weight = iftCreateFImage(next_label->xsize, next_label->ysize, next_label->zsize);
    prev_frontier = iftObjectBorderSet(prev_label, A);

    for (p = 0; p < prev_label->n; p++){
        prev_weight->val[p] = next_weight->val[p] = 1.0;
    }

    S = prev_frontier;
    while (S != NULL) {
        p = S->elem;
        iftBMapSet1(inFrontier,p);
        S = S->next;
    }

    /* Smooth frontier and reset its path values */

    for (iter = 0; iter < num_smooth_iterations; iter++)
    {
        while (prev_frontier != NULL)
        {
            p = iftRemoveSet(&prev_frontier);
            iftInsertSet(&next_frontier, p);
            u   = iftGetVoxelCoord(prev_label, p);

            for (l = 0; l <= prev_label_max_val; l++)
            {
                sum[l] = 0.0;
            }

            for (i = 1; i < A->n; i++)
            {
                v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(prev_label, v))
                {
                    q = iftGetVoxelIndex(prev_label, v);
                    sum[prev_label->val[q]] += prev_weight->val[q] * weight->val[q];
                    if (iftBMapValue(inFrontier, q) == 0) /* expand frontier */
                    {
                        if (igraph->pred[q] != IFT_NIL)
                        {
                            iftInsertSet(&next_frontier, q);
                            iftBMapSet1(inFrontier, q);
                        }
                    }
                }
            }

            for (l = 0; l <= prev_label_max_val; l++)
                sum[l]  = sum[l] / norm_factor->val[p];

            max_membership = IFT_INFINITY_FLT_NEG; max_label = IFT_NIL;
            for (l = 0; l <= prev_label_max_val; l++)
            {
                if (sum[l] > max_membership)
                {
                    max_membership = sum[l];
                    max_label      = l;
                }
            }
            next_label->val[p]  = max_label;
            next_weight->val[p] = sum[max_label];
        }

        prev_frontier = next_frontier;
        next_frontier = NULL;

        for (r = 0; r < prev_label->n; r++)
        {
            prev_weight->val[r] = next_weight->val[r];
            prev_label->val[r]  = next_label->val[r];
        }
    }

    iftFree(sum);
    iftDestroyFImage(&prev_weight);
    iftDestroyFImage(&next_weight);
    iftDestroyImage(&next_label);
    iftDestroyFImage(&weight);
    iftDestroyFImage(&norm_factor);
    iftDestroyBMap(&inFrontier);
    iftDestroySet(&prev_frontier);

    /* It fixes the label map, by eliminating the smallest regions and
       relabel them with the adjaceny labels */

    prev_label_max_val = iftMaximumValue(prev_label);
    next_label = iftSelectKLargestRegionsAndPropagateTheirLabels(prev_label, A, prev_label_max_val);
    for (p=0; p < next_label->n; p++)
        igraph->label[p]=next_label->val[p];

    iftDestroyImage(&next_label);
    iftDestroyImage(&prev_label);

}


void iftRenumerateBlockSuperpixels(iftImage *super_img, const iftImageTiles *blocks, int *out_n_cells) {
    if (super_img == NULL)
        iftError("Superpixel Image is NULL", "iftRenumerateBlockSuperpixels");
    if (blocks == NULL)
        iftError("Blocks is NULL", "iftRenumerateBlockSuperpixels");
    if (blocks->ntiles <= 0)
        iftError("Number of Tiles/Blocks %d <= 0", "iftRenumerateBlockSuperpixels", blocks->ntiles);

    //
    int n_cells = 0;


    for (int b = 0; b < blocks->ntiles; b++) {
        iftBoundingBox bb  = blocks->tile_coords[b];
        int n_cells_block = iftMaximumValueInRegion(super_img, bb);

        // re-enumerates the block's superpixels by adding the total number of labels until the last visited block
        iftVoxel v;
        for (v.z = bb.begin.z; v.z <= bb.end.z; v.z++) {
            for (v.y = bb.begin.y; v.y <= bb.end.y; v.y++) {
                for (v.x = bb.begin.x; v.x <= bb.end.x; v.x++) {
                    int p = iftGetVoxelIndex(super_img, v);
                    super_img->val[p] += n_cells;
                }
            }
        }
        n_cells += n_cells_block;
    }

    if (out_n_cells != NULL)
        *out_n_cells = n_cells;
}


void iftWriteSuperpixelBorders(const iftImage *img, const iftImage *super_img, const char *out_path) {
    if (img == NULL)
        iftError("Image is NULL", "iftWriteSuperpixelBorders");
    if (super_img == NULL)
        iftError("Superpixel Image is NULL", "iftWriteSuperpixelBorders");
    if (!iftIsDomainEqual(img, super_img))
        iftError("Image and Superpixel Image have different domain\n" \
                 "Img: %dx%dx%d\nSuperpixel: %dx%dx%d\n", "iftWriteSuperpixelBorders",
                 img->xsize, img->ysize, img->zsize, super_img->xsize, super_img->ysize, super_img->zsize);
    if (out_path == NULL)
        iftError("Output Pathname is NULL", "iftWriteSuperpixelBorders");

    int max_range      = iftNormalizationValue(iftMaximumValue(img));
    float color_factor = 255.0 / max_range;
    iftColor RGB       = {.val = {255, 0, 0}};
    iftColor YCbCr     = iftRGBtoYCbCr(RGB, 255);

    iftImage *copy_img = iftCopyImage(img);
    iftAdjRel *A       = iftSpheric(0.0);

    if (iftIs3DImage(copy_img)) {
        if (!iftEndsWith(out_path, ".zip"))
            iftWarning("Output Pathname should be *.zip, because the image is 3D", "iftWriteSuperpixelBorders");
        // if (iftIsColorImage(copy_img))
        //     iftError("Not supported 3D Color Images yet", "iftWriteSuperpixelBorders");

        char *tmp_dir = iftMakeTempDir("slices", NULL, NULL);

        for (int z = 0; z < img->zsize; z++) {
            char filename[32];
            sprintf(filename, "slice_%06d.png", z);
            char *out_slice_path = iftJoinPathnames(2, tmp_dir, filename);

            iftImage *slice_img       = iftGetXYSlice(img, z);
            iftSetCbCr(slice_img, (iftMaxImageRange(iftImageDepth(img))+1)/2);

            // converts the gray image to 8-bits
            for (int p = 0; p < slice_img->n; p++) {
                slice_img->val[p] *= color_factor; // it's in Gray 0..max_val -  the correct is to convert this value to YCbCr
            }

            iftImage *slice_super_img = iftGetXYSlice(super_img, z);
            iftImage *slice_bimg      = iftBorderImage(slice_super_img,1);

            iftDrawBorders(slice_img, slice_bimg, A, YCbCr, A);
            iftWriteImageByExt(slice_img, out_slice_path);

            iftFree(out_slice_path);
            iftDestroyImage(&slice_img);
            iftDestroyImage(&slice_super_img);
            iftDestroyImage(&slice_bimg);
        }

        // Draws the Supervoxel's border on 3D image
        char *border_path = iftJoinPathnames(2, tmp_dir, "img_with_borders_3D.scn");
        iftImage *bimg    = iftBorderImage(super_img,1);

        iftColor RGB_white   = {.val = {max_range, max_range, max_range}};
        iftColor YCbCr_white = iftRGBtoYCbCr(RGB_white, max_range); // this is the white color in YCbCr
        iftDrawBorders(copy_img, bimg, A, YCbCr_white, A);
        iftWriteImageByExt(copy_img, border_path);
        iftFree(border_path);

        // Writes the Supervoxel's borders on 3D image
        border_path = iftJoinPathnames(2, tmp_dir, "borders_3D.scn");
        // sets the label borders to 1 (binary)
        for (int p = 0; p < bimg->n; p++)
            bimg->val[p] = (bimg->val[p] != 0);

        iftWriteImageByExt(bimg, border_path);
        iftFree(border_path);


        iftZipDirContent(tmp_dir, out_path);

        iftRemoveDir(tmp_dir);
        iftDestroyImage(&bimg);
        iftFree(tmp_dir);
    }
    else {
        if ((!iftEndsWith(out_path, ".ppm")) && (!iftEndsWith(out_path, ".png")))
            iftWarning("Output Pathname should be *.[ppm,png], because the image is 2D", "iftWriteSuperpixelBorders");

        iftImage *bimg = iftBorderImage(super_img,1);
        iftDrawBorders(copy_img, bimg, A, YCbCr, A);
        iftWriteImageByExt(copy_img, out_path);

        iftDestroyImage(&bimg);
    }


    iftDestroyAdjRel(&A);
    iftDestroyImage(&copy_img);
}

iftImage *iftExtract_ISF_MIX_ROOT_Superpixels(iftImage *img, int nsuperpixels, float alpha, float beta, int niters, int smooth_niters)
{
    iftImage  *mask1, *seeds, *label;
    iftMImage *mimg;
    iftAdjRel *A;
    iftIGraph *igraph;
    timer     *t1=NULL,*t2=NULL;

    /* Compute ISF superpixels */
    if (iftIs3DImage(img)){
        A      = iftSpheric(1.0);
    } else {
        A      = iftCircular(1.0);
    }

    if (iftIsColorImage(img)){
        mimg   = iftImageToMImage(img,LABNorm_CSPACE);
    } else {
        mimg   = iftImageToMImage(img,GRAY_CSPACE);
    }

    mask1  = iftSelectImageDomain(mimg->xsize,mimg->ysize,mimg->zsize);

    /* minima of a basins manifold in that domain */
    igraph = iftImplicitIGraph(mimg,mask1,A);

    t1 = iftTic();
    /* seed sampling for ISF */
    seeds   = iftAltMixedSampling(mimg,mask1,nsuperpixels);

    iftDestroyImage(&mask1);
    iftDestroyMImage(&mimg);

    iftIGraphISF_Root(igraph,seeds,alpha,beta,niters);

    /* Smooth regions in the label map of igraph */
    if (smooth_niters > 0){
        iftIGraphSetWeightForRegionSmoothing(igraph, img);
        iftIGraphSmoothRegions(igraph, smooth_niters);
    }
    label   = iftIGraphLabel(igraph);
    t2 = iftToc();
    printf("ISF proc time im ms: %f\n", iftCompTime(t1,t2));

    iftDestroyImage(&seeds);
    iftDestroyIGraph(&igraph);
    iftDestroyAdjRel(&A);

    return label;
}


iftImage *iftExtract_ISF_MIX_MEAN_Superpixels(iftImage *img, int nsuperpixels, float alpha, float beta, int niters, int smooth_niters) {
    iftImage  *mask1, *seeds, *label;
    iftMImage *mimg;
    iftAdjRel *A;
    iftIGraph *igraph;
    timer     *t1=NULL,*t2=NULL;

    /* Compute ISF superpixels */
    if (iftIs3DImage(img)){
        A      = iftSpheric(1.0);
    } else {
        A      = iftCircular(1.0);
    }

    if (iftIsColorImage(img)){
        mimg   = iftImageToMImage(img,LABNorm_CSPACE);
    } else {
        mimg   = iftImageToMImage(img,GRAY_CSPACE);
    }

    mask1  = iftSelectImageDomain(mimg->xsize,mimg->ysize,mimg->zsize);

    /* minima of a basins manifold in that domain */
    igraph = iftImplicitIGraph(mimg,mask1,A);

    t1 = iftTic();
    /* seed sampling for ISF */
    seeds   = iftAltMixedSampling(mimg,mask1,nsuperpixels);

    iftDestroyImage(&mask1);
    iftDestroyMImage(&mimg);

    iftIGraphISF_Mean(igraph,seeds,alpha,beta,niters);

    /* Smooth regions in the label map of igraph */
    if (smooth_niters > 0){
        iftIGraphSetWeightForRegionSmoothing(igraph, img);
        iftIGraphSmoothRegions(igraph, smooth_niters);
    }
    label   = iftIGraphLabel(igraph);
    t2 = iftToc();
    printf("ISF proc time im ms: %f\n", iftCompTime(t1,t2));

    iftDestroyImage(&seeds);
    iftDestroyIGraph(&igraph);
    iftDestroyAdjRel(&A);


    return label;
}

iftSet *iftSuperpixelCenterSetFromIGraph(iftIGraph *igraph) {
    iftSet* S = NULL;
    for (int t = 0; t < igraph->nnodes; t++)
    {
        int p = igraph->node[t].voxel;
        int r = igraph->root[p];
        if (p == r)
        {
            iftInsertSet(&S, p);
        }
    }

    return S;
}

