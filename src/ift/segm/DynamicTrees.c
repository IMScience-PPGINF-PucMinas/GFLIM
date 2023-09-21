#include "ift/segm/DynamicTrees.h"

#include "ift/core/dtypes/FHeap.h"
#include "ift/core/dtypes/GQueue.h"
#include "ift/core/io/Stream.h"
#include "ift/metriclearn/MetricLearnCore.h"
#include "ift/metriclearn/LargeMargin.h"

iftDynamicSet *iftCreateDynamicSet(int dim)
{
    iftDynamicSet *S = iftAlloc(sizeof (*S), 1);
    if (!S)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDynamicSet");
    S->size = 0;
    S->dim = dim;
    S->mean = iftAlloc(dim, sizeof *S->mean);
    if (!S->mean)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDynamicSet");
    S->begin = NULL;
    S->end = NULL;

    return S;
}


iftDTForest *iftCreateDTForest(const iftMImage *mimg, const iftAdjRel *A,
                               const iftLabeledSet *seeds, float delta, float gamma)
{
    if (gamma < 0.0f)
        iftError("Gamma must be greater than 0.0", "iftDynTreeClosestRoot");

    if (delta < 0)
        iftError("Delta must be greater than 0", "iftDynTreeClosestRoot");

    iftDTForest *forest = iftAlloc(1, sizeof *forest);
    forest->label = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    forest->cost = iftCreateFImage(mimg->xsize, mimg->ysize, mimg->zsize);
    forest->root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    forest->pred = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    forest->order = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    forest->delay = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    forest->dyn_trees = iftAlloc(mimg->n, sizeof (*forest->dyn_trees));

    iftImage *label = forest->label;
    iftFImage *cost = forest->cost;
    iftImage *root = forest->root;
    iftImage *pred = forest->pred;

    iftImage *order = forest->order;
    iftImage *delay = forest->delay;

    iftDynamicSet **S = forest->dyn_trees;

    iftFHeap *Q = iftCreateFHeap(mimg->n, cost->val);

    for (int p = 0; p < mimg->n; p++)
    {
        cost->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
        pred->val[p] = IFT_NIL;
        delay->val[p] = 0;
    }

    for (const iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        cost->val[p] = delta;
        root->val[p] = p;
        iftInsertFHeap(Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
    }

    int count = 0;
    while (!iftEmptyFHeap(Q))
    {
        int p = iftRemoveFHeap(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        order->val[p] = count++;

        if (pred->val[p] != IFT_NIL)
            delay->val[p] = order->val[p] - order->val[ pred->val[p] ];

        if (root->val[p] == p) {
            cost->val[p] = 0;
        } else { /* roots have already been inserted */
            iftInsertDynamicSet(S[root->val[p] ], mimg, p);
        }

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->color[q] != IFT_BLACK)
                {
                    double arc_weight = iftDistDynamicSetMImage(S[root->val[p] ], mimg, q) +
                                        gamma * iftMImageSqDist(mimg, p, q);

                    float tmp = iftMax((float) arc_weight, cost->val[p]);

                    if (tmp < cost->val[q])
                    {

                        label->val[q] = label->val[p];
                        cost->val[q] = tmp;
                        root->val[q] = root->val[p];
                        pred->val[q] = p;

                        if (Q->color[q] == IFT_GRAY)
                            iftGoUpFHeap(Q, Q->pos[q]);
                        else
                            iftInsertFHeap(Q, q);
                    }
                }
            }
        }
    }

    iftDestroyFHeap(&Q);

    return forest;
}


void iftDestroyDynamicSet(iftDynamicSet **S)
{
    iftDynamicSet *aux = *S;
    if (aux != NULL)
    {
        iftDestroySet(&aux->begin);
        iftFree(aux->mean);
        iftFree(aux);
    }
    aux = NULL;
}

void iftDestroyDTForest(iftDTForest **forest)
{
    iftDTForest *aux = *forest;
    if (aux) {
        for (int i = 0; i < aux->label->n; i++) {
            if (aux->dyn_trees[i])
                iftDestroyDynamicSet(&aux->dyn_trees[i]);
        }
        iftDestroyImage(&aux->label);
        iftDestroyFImage(&aux->cost);
        iftDestroyImage(&aux->root);
        iftDestroyImage(&aux->pred);
        iftDestroyImage(&aux->order);
        iftDestroyImage(&aux->delay);
    }
}

iftImage *iftDynamicSetObjectPolicy(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, bool use_dist)
{
    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    int n_labels = iftNumberOfLabels(seeds);
    iftDynamicSet **S = iftAlloc(n_labels, sizeof (*S));

    int max_val = iftMMaximumValue(mimg, -1);
    iftGQueue *Q = iftCreateGQueue(max_val * max_val * 3, mimg->n, pathval->val);

    for (int i = 0; i < n_labels; i++)
        S[i] = iftCreateDynamicSet(mimg->m);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = 0;
        iftInsertGQueue(&Q, p);
        iftInsertDynamicSet(S[label->val[p] ], mimg, p);
    }

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (pathval->val[p] != 0) /* roots have already been inserted*/
            iftInsertDynamicSet(S[label->val[p]], mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    double arc_weight = iftDistDynamicSetMImage(S[ label->val[p] ], mimg, q);

                    int tmp;
                    if (use_dist) {
                        tmp = iftMax(pathval->val[p], arc_weight + iftMImageSqDist(mimg, p, q));
                    } else {
                        tmp = iftMax(pathval->val[p], arc_weight);
                    }

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    for (int i = 0; i < n_labels; i++)
        iftDestroyDynamicSet(&S[i]);
    iftFree(S);
    iftDestroyGQueue(&Q);

    return label;
}


iftImage *iftDynamicSetRootPolicy(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist)
{
    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftDynamicSet **S = iftAlloc(mimg->n, sizeof (*S));

    int max_val = iftMMaximumValue(mimg, -1);
    iftGQueue *Q = iftCreateGQueue(max_val * max_val * 3, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = h;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
    }

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (root->val[p] == p) {
            pathval->val[p] = 0;
        } else /* roots have already been inserted */
            iftInsertDynamicSet(S[root->val[p] ], mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK) {
                    double arc_weight = iftDistDynamicSetMImage(S[root->val[p]], mimg, q);

                    int tmp;
                    if (use_dist) {
                        tmp = iftMax(pathval->val[p], arc_weight + iftMImageSqDist(mimg, p, q));
                    } else {
                        tmp = iftMax(pathval->val[p], arc_weight);
                    }

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        root->val[q] = root->val[p];
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    iftDestroyImage(&root);
    for (iftLabeledSet *M = seeds; M != NULL; M = M->next) {
        int p = M->elem;
        if (S[p] != NULL) {
            iftDestroyDynamicSet(&S[p]);
        }
    }
    iftFree(S);
    iftDestroyGQueue(&Q);

    return label;
}

iftImage *iftDynamicSetRootPolicyInMask(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, iftImage *mask) {
    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftDynamicSet **S = iftAlloc(mimg->n, sizeof (*S));

    int max_val = iftMMaximumValue(mimg, -1);
    iftGQueue *Q = iftCreateGQueue(max_val * max_val * 3, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
        if (mask->val[p] == 0) {
            pathval->val[p] = 0;
        }
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = 0;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
    }

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (root->val[p] == p) {
            pathval->val[p] = 0;
        } else /* roots have already been inserted */
            iftInsertDynamicSet(S[root->val[p] ], mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK) {
                    double arc_weight = iftDistDynamicSetMImage(S[root->val[p]], mimg, q);

                    int tmp = iftMax(pathval->val[p], arc_weight);

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        root->val[q] = root->val[p];
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    iftDestroyImage(&root);
    for (iftLabeledSet *M = seeds; M != NULL; M = M->next) {
        int p = M->elem;
        if (S[p] != NULL) {
            iftDestroyDynamicSet(&S[p]);
        }
    }
    iftFree(S);
    iftDestroyGQueue(&Q);

    return label;
}

iftImage *iftDynamicSetMinRootPolicy(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist)
{
    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftDynamicSet **S = iftAlloc(mimg->n, sizeof (*S));

    int max_val = iftMMaximumValue(mimg, -1);
    iftGQueue *Q = iftCreateGQueue(max_val * max_val * 3, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = h;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
    }

    iftLabeledSet *subSeeds = NULL;
    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (root->val[p] == p)
        {
            pathval->val[p] = 0;
        }

        if (p != root->val[p]) /* roots have already been inserted */
            iftInsertDynamicSet(S[root->val[p] ], mimg, p);
        else
            iftInsertLabeledSet(&subSeeds, p, label->val[p]);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    double arc_weight = IFT_INFINITY_INT;
                    int r_min = -1;

                    for (iftLabeledSet *M = subSeeds; M != NULL; M = M->next)
                    {
                        if (M->label == label->val[p]) {
                            double maybe_arc = iftDistDynamicSetMImage(S[M->elem], mimg, q);
                            if (maybe_arc < arc_weight) {
                                arc_weight = maybe_arc;
                                r_min = M->elem;
                            }
                        }
                    }

                    int tmp;
                    if (use_dist) {
                        tmp = iftMax(pathval->val[p], arc_weight + iftMImageSqDist(mimg, p, q));
                    } else {
                        tmp = iftMax(pathval->val[p], arc_weight);
                    }

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        root->val[q] = r_min;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    iftDestroyImage(&root);
    for (iftLabeledSet *M = seeds; M != NULL; M = M->next) {
        int p = M->elem;
        if (S[p] != NULL)
            iftDestroyDynamicSet(&S[p]);
    }
    iftFree(S);
    iftDestroyGQueue(&Q);
    iftDestroyLabeledSet(&subSeeds);

    return label;
}


iftImage *iftDynamicSetWithCluster(iftMImage *mimg, iftImage *cluster, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist) {

    iftImage *label = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftDynamicSet **S = iftAlloc(mimg->n, sizeof(*S));

    int max_val = iftMMaximumValue(mimg, -1);
    iftGQueue *Q = iftCreateGQueue(max_val * max_val * 3, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = h;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
    }

    iftLabeledSet *subSeeds = NULL;
    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (root->val[p] == p) {
            pathval->val[p] = 0;
        }

        if (p != root->val[p]) {
            if (cluster->val[p] == cluster->val[root->val[p] ])
                iftInsertDynamicSet(S[root->val[p] ], mimg, p);
        } else {
            iftInsertLabeledSet(&subSeeds, p, label->val[p]);
        }

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    double arc_weight = IFT_INFINITY_INT;
                    int r_min = -1;

                    for (iftLabeledSet *M = subSeeds; M != NULL; M = M->next)
                    {
                        int t = M->elem;
                        if (M->label == label->val[p]) {
                            double maybe_arc = iftDistDynamicSetMImage(S[t], mimg, q);
                            if (maybe_arc < arc_weight) {
                                arc_weight = maybe_arc;
                                r_min = t;
                            }
                        }
                    }

                    int tmp;
                    if (use_dist) {
                        tmp = iftMax(pathval->val[p], arc_weight + iftMImageSqDist(mimg, p, q));
                    } else {
                        tmp = iftMax(pathval->val[p], arc_weight);
                    }

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        root->val[q] = r_min;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    iftDestroyImage(&root);
    for (iftLabeledSet *M = seeds; M != NULL; M = M->next) {
        int p = M->elem;
        if (S[p] != NULL)
            iftDestroyDynamicSet(&S[p]);
    }
    iftFree(S);
    iftDestroyGQueue(&Q);
    iftDestroyLabeledSet(&subSeeds);

    return label;
}


iftImage *iftDynamicSetRootEnhanced(iftMImage *mimg, iftImage *objmap, iftAdjRel *A, iftLabeledSet *seeds, int h, float alpha, bool use_dist)
{
    if (alpha < 0.0f || alpha > 1.0f) {
        iftError("Alpha must be between (0, 1)", "iftDynamicSetRootEnhanced");
    }

    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftDynamicSet **S = iftAlloc(mimg->n, sizeof (*S));

    int max_val = iftMMaximumValue(mimg, -1);
    iftGQueue *Q = iftCreateGQueue(max_val * max_val * 3, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = h;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
    }

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (root->val[p] == p)
        {
            pathval->val[p] = 0;
        }

        if (p != root->val[p]) /* roots have already been inserted */
            iftInsertDynamicSet(S[root->val[p] ], mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    /* squared distance because, dynamic dist is also squared */
					double arc_weight = (1 - alpha) * iftDistDynamicSetMImage(S[root->val[p] ], mimg, q) +
                            (objmap->val[q] - objmap->val[p]) * (objmap->val[q] - objmap->val[p]) * alpha;

                    int tmp;
                    if (use_dist) {
                        tmp = iftMax(pathval->val[p], arc_weight + (1 - alpha) * iftMImageSqDist(mimg, p, q));
                    } else {
                        tmp = iftMax(pathval->val[p], arc_weight);
                    }

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        root->val[q] = root->val[p];
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    iftDestroyImage(&root);
    for (iftLabeledSet *M = seeds; M != NULL; M = M->next) {
        int p = M->elem;
        if (S[p] != NULL) {
            iftDestroyDynamicSet(&S[p]);
        }
    }
    iftFree(S);
    iftDestroyGQueue(&Q);

    return label;
}



iftImage *iftDynamicSetMinRootEnhanced(iftMImage *mimg, iftImage *objmap, iftAdjRel *A, iftLabeledSet *seeds, int h, float alpha, bool use_dist)
{
    if (alpha < 0.0f || alpha > 1.0f) {
        iftError("Alpha must be between (0, 1)", "iftDynamicSetMinRootEnhanced");
    }

    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftDynamicSet **S = iftAlloc(mimg->n, sizeof (*S));

    int max_val = iftMMaximumValue(mimg, -1);
    iftGQueue *Q = iftCreateGQueue(max_val * max_val * 3, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = h;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
    }

    iftLabeledSet *subSeeds = NULL;
    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (root->val[p] == p)
        {
            pathval->val[p] = 0;
        }

        if (p != root->val[p]) /* roots have already been inserted */
            iftInsertDynamicSet(S[root->val[p] ], mimg, p);
        else
            iftInsertLabeledSet(&subSeeds, p, label->val[p]);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    double arc_weight = IFT_INFINITY_INT;
                    int r_min = -1;

                    for (iftLabeledSet *M = subSeeds; M != NULL; M = M->next)
                    {
                        if (M->label == label->val[p]) {
                            double maybe_arc = iftDistDynamicSetMImage(S[M->elem], mimg, q);
                            if (maybe_arc < arc_weight) {
                                arc_weight = maybe_arc;
                                r_min = M->elem;
                            }
                        }
                    }

                    /* squared distance because, dynamic dist is also squared */
					arc_weight = arc_weight * (1 - alpha) +
                            (objmap->val[p] - objmap->val[q]) * (objmap->val[p] - objmap->val[q]) * alpha;

                    int tmp;
                    if (use_dist) {
                        tmp = iftMax(pathval->val[p], arc_weight + (1 - alpha) * iftMImageSqDist(mimg, p, q));
                    } else {
                        tmp = iftMax(pathval->val[p], arc_weight);
                    }

					if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        root->val[q] = r_min;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    iftDestroyImage(&root);
    for (iftLabeledSet *M = seeds; M != NULL; M = M->next) {
        int p = M->elem;
        if (S[p] != NULL)
            iftDestroyDynamicSet(&S[p]);
    }
    iftFree(S);
    iftDestroyGQueue(&Q);
    iftDestroyLabeledSet(&subSeeds);

    return label;
}


iftImage *iftDynTreeRoot(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int delta,
                         float gamma, iftImage *objmap, float alpha)
{
    if (objmap == NULL)
        alpha = 0.0f;

    if (alpha < 0.0f || alpha > 1.0f)
        iftError("Alpha must be between (0, 1)", "iftDynTreeClosestRoot");

    if (gamma < 0.0f)
        iftError("Gamma must be greater than 0.0", "iftDynTreeClosestRoot");

    if (delta < 0)
        iftError("Delta must be greater than 0", "iftDynTreeClosestRoot");

    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftFImage *pathval = iftCreateFImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftDynamicSet **S = iftAlloc(mimg->n, sizeof (*S));

    iftFHeap *Q = iftCreateFHeap(mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = delta;
        root->val[p] = p;
        iftInsertFHeap(Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
    }

    while (!iftEmptyFHeap(Q))
    {
        int p = iftRemoveFHeap(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (root->val[p] == p) {
            pathval->val[p] = 0;
        } else { /* roots have already been inserted */
            iftInsertDynamicSet(S[root->val[p] ], mimg, p);
        }

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->color[q] != IFT_BLACK)
                {
                    double arc_weight = iftDistDynamicSetMImage(S[root->val[p] ], mimg, q);

                    arc_weight += gamma * iftMImageSqDist(mimg, p, q);

                    /* squared distance because, dyntree dist is also squared */
                    if (objmap != NULL)
                        arc_weight = (1.0f - alpha) * arc_weight +
                                     alpha * (objmap->val[p] - objmap->val[q]) * (objmap->val[p] - objmap->val[q]);

                    double tmp = iftMax(arc_weight, pathval->val[p]);

                    if (tmp < pathval->val[q])
                    {
                        label->val[q] = label->val[p];
                        pathval->val[q] = (float) tmp;
                        root->val[q] = root->val[p];

                        if (Q->color[q] == IFT_GRAY)
                            iftGoUpFHeap(Q, Q->pos[q]);
                        else
                            iftInsertFHeap(Q, q);
                    }
                }
            }
        }
    }

    iftDestroyFImage(&pathval);
    iftDestroyImage(&root);
    for (iftLabeledSet *M = seeds; M != NULL; M = M->next) {
        int p = M->elem;
        if (S[p] != NULL) {
            iftDestroyDynamicSet(&S[p]);
        }
    }
    iftFree(S);
    iftDestroyFHeap(&Q);

    return label;
}


iftImage *iftDynTreeClosestRoot(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int delta,
                                float gamma, iftImage *objmap, float alpha)
{
    if (objmap == NULL)
        alpha = 0.0f;

    if (alpha < 0.0f || alpha > 1.0f)
        iftError("Alpha must be between (0, 1)", "iftDynTreeClosestRoot");

    if (gamma < 0.0f)
        iftError("Gamma must be greater than 0.0", "iftDynTreeClosestRoot");

    if (delta < 0)
        iftError("Delta must be greater than 0", "iftDynTreeClosestRoot");

    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftFImage *pathval = iftCreateFImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftDynamicSet **S = iftAlloc(mimg->n, sizeof (*S));

    iftFHeap *Q = iftCreateFHeap(mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = delta;
        root->val[p] = p;
        iftInsertFHeap(Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
    }

    /* Subset of seed nodes for speed up when delta > 0 because not all seed nodes are initialized
     * on this subset, some may not be used to find closest root until all seed nodes have left the queue */
    iftLabeledSet *subSeeds = NULL;

    while (!iftEmptyFHeap(Q))
    {
        int p = iftRemoveFHeap(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (root->val[p] == p) {
            pathval->val[p] = 0;
            iftInsertLabeledSet(&subSeeds, p, label->val[p]);
        } else {
            iftInsertDynamicSet(S[root->val[p] ], mimg, p);
        }

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->color[q] != IFT_BLACK)
                {
                    double arc_weight = IFT_INFINITY_INT;
                    int r_min = -1;

                    for (iftLabeledSet *M = subSeeds; M != NULL; M = M->next)
                    {
                        if (M->label == label->val[p]) {
                            double maybe_arc = iftDistDynamicSetMImage(S[M->elem], mimg, q);
                            if (maybe_arc < arc_weight) {
                                arc_weight = maybe_arc;
                                r_min = M->elem;
                            }
                        }
                    }


                    arc_weight += gamma * iftMImageSqDist(mimg, p, q);

                    /* squared distance because, dyntree dist is also squared */
                    if (objmap != NULL)
                        arc_weight = (1.0f - alpha) * arc_weight +
                                     alpha * (objmap->val[p] - objmap->val[q]) * (objmap->val[p] - objmap->val[q]);


                    double tmp = iftMax(arc_weight, pathval->val[p]);

                    if (tmp < pathval->val[q])
                    {

                        label->val[q] = label->val[p];
                        pathval->val[q] = (float) tmp;
                        root->val[q] = r_min;

                        if (Q->color[q] == IFT_GRAY)
                            iftGoUpFHeap(Q, Q->pos[q]);
                        else
                            iftInsertFHeap(Q, q);
                    }
                }
            }
        }
    }

    iftDestroyFImage(&pathval);
    iftDestroyImage(&root);
    for (iftLabeledSet *M = seeds; M != NULL; M = M->next) {
        int p = M->elem;
        if (S[p] != NULL)
            iftDestroyDynamicSet(&S[p]);
    }
    iftFree(S);
    iftDestroyFHeap(&Q);
    iftDestroyLabeledSet(&subSeeds);

    return label;
}


iftDynamicSet_CIARP *iftCreateDynamicSet_CIARP(void)
{
    iftDynamicSet_CIARP *S = iftAlloc(sizeof (*S), 1);

    S->size = 0;

    S->mean = iftAlloc(3, sizeof *S->mean);

    return S;
}


void iftDestroyDynamicSet_CIARP(iftDynamicSet_CIARP **S)
{
    iftDynamicSet_CIARP *aux = *S;
    if (aux != NULL)
    {
        iftFree(aux->mean);
        iftFree(aux);
    }
    aux = NULL;
}


iftImage *iftDynamicSetObjectPolicy_CIARP(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, bool use_dist)
{
    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    int n_labels = iftNumberOfLabels(seeds);
    iftDynamicSet_CIARP **S = iftAlloc(n_labels, sizeof (*S));

    iftGQueue *Q = iftCreateGQueue(255 * 2, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
    }

	for (int i = 0; i < n_labels; i++)
		S[i] = iftCreateDynamicSet_CIARP();

	for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = 0;
        iftInsertGQueue(&Q, p);
        iftInsertDynamicSet_CIARP(S[label->val[p] ], mimg, p);
    }

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (pathval->val[p] != 0) /* roots have already been inserted*/
            iftInsertDynamicSet_CIARP(S[label->val[p]], mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {

                    double arc_weight = iftDistDynamicSetMImage_CIARP(S[ label->val[p] ], mimg, q);

                    int tmp;
                    if (use_dist)
                        tmp = iftMax(pathval->val[p], arc_weight + iftMImageSqDist(mimg, p, q));
                    else
                        tmp = iftMax(pathval->val[p], arc_weight);

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    for (int i = 0; i < n_labels; i++)
    	iftDestroyDynamicSet_CIARP(&S[i]);
    iftFree(S);
    iftDestroyGQueue(&Q);

    return label;
}


iftImage *iftDynamicSetRootPolicy_CIARP(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, bool use_dist)
{
    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *output = iftCreateColorImage(mimg->xsize, mimg->ysize, mimg->zsize, 8);

    iftDynamicSet_CIARP **S = iftAlloc(mimg->n, sizeof (*S));

    iftGQueue *Q = iftCreateGQueue(255 * 2, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = 0;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
        S[p] = iftCreateDynamicSet_CIARP();
        iftInsertDynamicSet_CIARP(S[p], mimg, p);

    }


    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (p != root->val[p]) /* roots have already been inserted */
            iftInsertDynamicSet_CIARP(S[root->val[p] ], mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    double arc_weight = iftDistDynamicSetMImage_CIARP(S[root->val[p] ], mimg, q);

                    int tmp;
                    if (use_dist)
                        tmp = iftMax(pathval->val[p], arc_weight + iftMImageDist(mimg, p, q));
                    else
                        tmp =  iftMax(pathval->val[p], arc_weight);

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        root->val[q] = root->val[p];
                        iftInsertGQueue(&Q, q);

                        output->val[q] = output->val[p];
                        output->Cb[q] = output->Cb[p];
                        output->Cr[q] = output->Cr[p];
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    iftDestroyImage(&root);
    for (int p = 0; p < mimg->n; p++) {
        if (S[p] != NULL)
            iftDestroyDynamicSet_CIARP(&S[p]);
    }
    iftFree(S);
    iftDestroyGQueue(&Q);

    return label;
}


iftImage *iftDynamicSetMinRootPolicy_CIARP(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, bool use_dist)
{
    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftDynamicSet_CIARP **S = iftAlloc(mimg->n, sizeof (*S));

    iftGQueue *Q = iftCreateGQueue(255 * 2, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = 0;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
        S[p] = iftCreateDynamicSet_CIARP();
        iftInsertDynamicSet_CIARP(S[p], mimg, p);
    }

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (p != root->val[p]) /* roots have already been inserted */
            iftInsertDynamicSet_CIARP(S[root->val[p] ], mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    double arc_weight = IFT_INFINITY_INT;
                    int r_min = -1;

                    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
                    {
                        if (M->label == label->val[p]) {
                            double maybe_arc = iftDistDynamicSetMImage_CIARP(S[M->elem], mimg, q);
                            if (maybe_arc < arc_weight) {
                                arc_weight = maybe_arc;
                                r_min = M->elem;
                            }
                        }
                    }

                    int tmp;
                    if (use_dist)
                        tmp = iftMax(pathval->val[p], arc_weight + iftMImageDist(mimg, p, q));
                    else
                        tmp =  iftMax(pathval->val[p], arc_weight);

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        root->val[q] = r_min;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    iftDestroyImage(&root);
    for (int p = 0; p < mimg->n; p++) {
        if (S[p] != NULL)
            iftDestroyDynamicSet_CIARP(&S[p]);
    }
    iftFree(S);
    iftDestroyGQueue(&Q);

    return label;
}


iftMImage *iftDTRootWeightsMap(iftMImage *mimg, iftAdjRel *A, iftLabeledSet *seeds, int h, bool use_dist) {
    iftImage *label   = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *pathval = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *root = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftMImage *weights_map = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, 1);

    iftDynamicSet **S = iftAlloc(mimg->n, sizeof (*S));

    int max_val = iftMMaximumValue(mimg, -1);
    iftGQueue *Q = iftCreateGQueue(max_val * max_val * 3, mimg->n, pathval->val);

    printf("TO AQUI CACETE\n");

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
        root->val[p] = IFT_NIL;
    }

    for (iftLabeledSet *M = seeds; M != NULL; M = M->next)
    {
        int p = M->elem;
        label->val[p] = M->label;
        pathval->val[p] = h;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
        S[p] = iftCreateDynamicSet(mimg->m);
        iftInsertDynamicSet(S[p], mimg, p);
        weights_map->val[p][0] = 0;
    }

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (root->val[p] == p) {
            pathval->val[p] = 0;
        } else /* roots have already been inserted */
            iftInsertDynamicSet(S[root->val[p] ], mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK) {
                    double arc_weight = iftDistDynamicSetMImage(S[root->val[p]], mimg, q);

                    int tmp;
                    if (use_dist) {
                        tmp = iftMax(pathval->val[p], arc_weight + iftMImageSqDist(mimg, p, q));
                    } else {
                        tmp = iftMax(pathval->val[p], arc_weight);
                    }

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q] = tmp;
                        root->val[q] = root->val[p];
                        iftInsertGQueue(&Q, q);
                        weights_map->val[q][0] = tmp;
                    }
                }
            }
        }
    }

    iftDestroyImage(&pathval);
    iftDestroyImage(&root);
    for (iftLabeledSet *M = seeds; M != NULL; M = M->next) {
        int p = M->elem;
        if (S[p] != NULL) {
            iftDestroyDynamicSet(&S[p]);
        }
    }
    iftFree(S);
    iftDestroyGQueue(&Q);
    iftDestroyImage(&label);

    return weights_map;
}