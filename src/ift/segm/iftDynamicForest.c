#include "ift/segm/iftDynamicForest.h"

iftDynamicNode *iftInsertDynamicTree(iftDynamicTree **list, int nfeats)  {
    iftDynamicNode *node = (iftDynamicNode *) iftAlloc(1, sizeof(iftDynamicNode));

    node->nfeats     = nfeats;
    node->nnodes     = 0;
    node->sum_feat   = iftAllocFloatArray(nfeats);
    node->prev       = NULL;
    node->next       = *list;

    if (*list)
        (*list)->prev = node;
    *list = node;

    return node;
}

void iftRemoveDynamicTree(iftDynamicTree **list, iftDynamicNode **node) {
    iftDynamicTree *next = (*node)->next;
    iftDynamicTree *prev = (*node)->prev;

    if (*list == *node) {
        *list = next;
    }

    if (prev != NULL) {
        prev->next = next;
    }

    if (next != NULL) {
        next->prev = prev;
    }

    iftFree((*node)->sum_feat);
    iftFree(*node);
    *node = NULL;
}

void iftDestroyDynamicTrees(iftDynamicTree **list) {
    iftDynamicTree  *node = *list;

    while (node) {
        iftDynamicTree *next = node->next;
        iftFree(node->sum_feat);
        iftFree(node);
        node = next;
    }

    *list = NULL;
}

iftDynamicForest *iftCreateDynamicForest(iftMImage *mimg, iftAdjRel *A) {
    iftDynamicForest *forest = (iftDynamicForest*) iftAlloc(1, sizeof *forest);

    forest->cost           = iftCreateFImage(mimg->xsize, mimg->ysize, mimg->zsize);
    forest->heap           = iftCreateFHeap(mimg->n, forest->cost->val);
    forest->label          = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    forest->marker         = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    forest->root           = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    forest->pred           = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    forest->dynamic_roots  = NULL;
    forest->tree_map       = (iftDynamicTree **) iftAlloc(mimg->n, sizeof (iftDynamicTree*));

    forest->mimg           = mimg;
    forest->A              = A;

    for (ulong p = 0; p < mimg->n; p++)
    {
        forest->cost->val[p] = IFT_INFINITY_FLT;
        forest->pred->val[p] = IFT_NIL;
        forest->root->val[p] = p;
    }

    return forest;
}

void iftDestroyDynamicForest(iftDynamicForest **forest) {
    if ((*forest) != NULL) {
        iftDestroyDynamicTrees(&(*forest)->dynamic_roots);

        iftFree((*forest)->tree_map);

        iftDestroyFImage(&(*forest)->cost);
        iftDestroyFHeap(&(*forest)->heap);
        iftDestroyImage(&(*forest)->label);
        iftDestroyImage(&(*forest)->marker);
        iftDestroyImage(&(*forest)->root);
        iftDestroyImage(&(*forest)->pred);

        iftFree(*forest);
        *forest = NULL;
    }
}

void iftResetDynamicForest(iftDynamicForest *forest) {
    iftResetFHeap(forest->heap);
    #pragma omp parallel for
    for (ulong p = 0; p < forest->mimg->n; p++)
    {
        iftDestroyDynamicTrees(&forest->dynamic_roots);
        forest->cost->val[p]             = IFT_INFINITY_FLT;
        forest->marker->val[p]           = 0;
        forest->root->val[p]             = p;
        forest->label->val[p]            = 0;
        forest->pred->val[p]             = IFT_NIL;
        forest->tree_map[p]              = NULL;
    }
}

float iftDynamicArcWeight(iftDynamicForest *forest, int p, int q) {
    iftMImage *mimg = forest->mimg;

    iftDynamicTree *node = forest->tree_map[forest->root->val[p]];

    float distance = 0;
    for (int i = 0; i < node->nfeats; i++) {
        float mean = node->sum_feat[i] / node->nnodes;
        distance += (mean - mimg->val[q][i]) * (mean - mimg->val[q][i]);
    }
    return distance;
}

void iftUpdateDynamicTree(iftDynamicForest *forest, int p, iftUpdateMode mode) {
    int r = forest->root->val[p];
    iftDynamicTree *node = forest->tree_map[r];

    if (mode == IFT_REM) {
        if (forest->tree_map[r] && forest->tree_map[p] == forest->tree_map[r]) {
            for (int i = 0; i < node->nfeats; i++)
                node->sum_feat[i] -= forest->mimg->val[p][i];
            node->nnodes--;
        }
        forest->tree_map[p] = NULL;
    } else if (mode == IFT_ADD) {
        for (int i = 0; i < node->nfeats; i++)
            node->sum_feat[i] += forest->mimg->val[p][i];
        node->nnodes++;
        forest->tree_map[p] = node;
    }
}

iftSet *
iftDynamicTreeRemoval(iftDynamicForest *fst, iftSet *trees_for_removal) {
    int        i, p, q, r, n = fst->mimg->n;
    float V0;
    iftVoxel   u, v;
    iftAdjRel *A = fst->A;
    iftSet   *Frontier = NULL;
    iftBMap  *inFrontier = iftCreateBMap(n);
    iftImage  *pred = fst->pred, *root = fst->root;
    iftFImage *cost = fst->cost;
    iftMImage *mimg = fst->mimg;
    iftSet   *T1 = NULL, *T2 = NULL;

    if (fst->heap->removal_policy == IFT_MINVALUE)
        V0 = IFT_INFINITY_FLT;
    else // MAXVALUE
        V0 = IFT_INFINITY_FLT_NEG;

    /* Remove all marked trees and find the frontier voxels
       afterwards. */

    while (trees_for_removal != NULL) {
        p = trees_for_removal->elem;

        if (cost->val[root->val[p]] != V0) { //tree not marked yet
            r = root->val[p];
            cost->val[r] = V0;
            pred->val[r] = IFT_NIL;
            iftRemoveDynamicTree(&fst->dynamic_roots, &fst->tree_map[r]);
            fst->tree_map[r] = NULL;
            iftInsertSet(&T1, r);
            while (T1 != NULL) {
                p = iftRemoveSet(&T1);
                iftInsertSet(&T2, p);

                u = iftMGetVoxelCoord(mimg, p);
                for (i = 1; i < A->n; i++) {
                    v = iftGetAdjacentVoxel(A, u, i);
                    if (iftMValidVoxel(mimg, v)) {
                        q   = iftMGetVoxelIndex(mimg, v);
                        if (cost->val[q] != V0) {
                            if (pred->val[q] == p) {
                                iftInsertSet(&T1, q);
                                iftUpdateDynamicTree(fst, q, IFT_REM);
                                cost->val[q] = V0; // mark removed node
                                pred->val[q]    = IFT_NIL;
                            }
                        }
                    }
                }
            }
        }

        trees_for_removal = trees_for_removal->next;
    }

    /* Find the frontier voxels of non-removed trees */

    while (T2 != NULL)
    {
        p = iftRemoveSet(&T2);
        u = iftMGetVoxelCoord(mimg, p);
        for (i = 1; i < A->n; i++)
        {
            v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v))
            {
                q   = iftMGetVoxelIndex(mimg, v);
                if (cost->val[q] != V0)
                {
                    if (iftBMapValue(inFrontier, q) == 0)
                    {
                        iftInsertSet(&Frontier, q);
                        iftBMapSet1(inFrontier, q);
                    }
                }
            }
        }
    }
    iftDestroyBMap(&inFrontier);

    return (Frontier);
}

void iftDynamicSubtreeRemoval(iftDynamicForest *forest, int t) {
    iftMImage *mimg = forest->mimg;
    iftImage  *pred = forest->pred;
    iftImage  *label = forest->label;
    iftImage  *root = forest->root;
    iftFImage *cost = forest->cost;
    iftAdjRel *A    = forest->A;
    iftFHeap  *heap = forest->heap;

    iftSet      *K = NULL;
    iftIntQueue *T = iftCreateIntQueue(mimg->n);
    iftIntArray *arr = iftCreateIntArray(mimg->n);

    iftInsertIntQueue(T, t);

    while (!iftIsIntQueueEmpty(T)) {
        int p;
        iftRemoveIntQueue(T, &p);

        if (heap->color[p] == IFT_GRAY) {
            iftRemoveFHeapElem(heap, p);
        }

        iftUpdateDynamicTree(forest, p, IFT_REM);
        heap->color[p]     = IFT_WHITE;
        pred->val[p]       = IFT_NIL;
        label->val[p]      = IFT_NIL;
        root->val[p]       = p;
        cost->val[p]       = IFT_INFINITY_FLT;

        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v)) {
                int q = iftMGetVoxelIndex(mimg, v);

                if (pred->val[q] == p) {
                    iftInsertIntQueue(T, q);
                } else {
                    if (arr->val[q] == 0 && cost->val[q] != IFT_INFINITY_FLT && root->val[q] != q) {
                        iftInsertSet(&K, q);
                        arr->val[q] = 1;
                    }
                }
            }
        }
    }

    while (K != NULL) {
        int p = iftRemoveSet(&K);
        arr->val[p] = 0;
        if (cost->val[p] != IFT_INFINITY_FLT) {
            if (heap->color[p] == IFT_GRAY) {
                iftGoUpFHeap(heap, heap->pos[p]);
            } else {
                iftInsertFHeap(heap, p);
            }
        }
    }

    iftDestroySet(&K);
    iftDestroyIntQueue(&T);
    iftDestroyIntArray(&arr);
}

void iftRemMarkersFromDynamicForest(iftDynamicForest *forest, iftSet *removal_markers) {
    iftSet    *Frontier = NULL;
    iftFHeap *heap = forest->heap;

    if (removal_markers != NULL)
    {
        Frontier = iftDynamicTreeRemoval(forest, removal_markers);
        while (Frontier != NULL) {
            int p = iftRemoveSet(&Frontier);
            /* p is also a seed voxel, but the priority is it as a seed. */
            if (heap->color[p] != IFT_GRAY) {
                iftInsertFHeap(heap, p);
            }
        }
    }
}

void iftAddMarkersToDynamicForest(iftDynamicForest *forest, iftLabeledSet *seeds) {
    iftFImage *cost = forest->cost;
    iftImage  *label = forest->label;
    iftImage  *pred = forest->pred;
    iftImage  *root = forest->root, *marker = forest->marker;
    iftFHeap  *heap = forest->heap;

    iftLabeledSet *M = seeds;
    while (M != NULL)
    {
        int p = M->elem;

        if (heap->color[p] == IFT_GRAY) {
            /* p is also a frontier voxel, but the priority is it as a seed. */
            iftRemoveFHeapElem(heap, p);
        }

        cost->val[p]    = 0;
        label->val[p]   = M->label;
        root->val[p]    = p;
        pred->val[p]    = IFT_NIL;
        marker->val[p]  = M->marker;
        iftInsertFHeap(heap, p);
        forest->tree_map[p] = iftInsertDynamicTree(&forest->dynamic_roots,
                                                   forest->mimg->m);
        iftUpdateDynamicTree(forest, p, IFT_ADD);

        M = M->next;
    }

}

void iftDiffDynamicIFT(iftDynamicForest *forest) {
    iftAdjRel *A     = forest->A;
    iftFImage *cost  = forest->cost;
    iftImage  *label = forest->label;
    iftImage  *pred  = forest->pred;
    iftImage  *root  = forest->root;
    iftMImage *mimg  = forest->mimg;
    iftFHeap  *heap  = forest->heap;

    /* Image Foresting Transform */
    while (!iftEmptyFHeap(heap))
    {
        int p      = iftRemoveFHeap(heap);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (cost->val[p] != 0) {/* it is not a seed voxel */
            iftUpdateDynamicTree(forest, p, IFT_ADD);
        }

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (heap->color[q] != IFT_BLACK) {
                    float arc_weight = iftDynamicArcWeight(forest, p, q);

                    float tmp = iftMax(cost->val[p], arc_weight);

                    if (tmp < cost->val[q]) {
                        iftUpdateDynamicTree(forest, q, IFT_REM);
                        pred->val[q]     = p;
                        root->val[q]     = root->val[p];
                        label->val[q]    = label->val[p];
                        cost->val[q]     = tmp;

                        if (heap->color[q] == IFT_WHITE) {
                            iftInsertFHeap(heap, q);
                        } else {
                            //iftRemoveFHeapElem(heap, q);
                            //iftInsertFHeap(heap, q);
                            iftGoUpFHeap(heap, heap->pos[q]);
                        }

                    } else if (pred->val[q] == p && label->val[p] != label->val[q]) {
                        iftDynamicSubtreeRemoval(forest, q);
                    }
                }
            }
        }
    }
    iftResetFHeap(heap);

}

float ift_dynamic_dist_object(iftDynamicForest *forest, int p, int q) {
    iftMImage *mimg = forest->mimg;

    iftDynamicTree *node = forest->tree_map[forest->label->val[p]];

    float distance = 0;
    for (int i = 0; i < node->nfeats; i++) {
        float mean = node->sum_feat[i] / node->nnodes;
        distance += (mean - mimg->val[q][i]) * (mean - mimg->val[q][i]);
    }
    return distance;
}

void ift_update_dynamic_trees_object(iftDynamicForest *forest, int p) {
    int l = forest->label->val[p];
    iftDynamicTree *node = forest->tree_map[l];

    for (int i = 0; i < node->nfeats; i++)
        node->sum_feat[i] += forest->mimg->val[p][i];
    node->nnodes++;
}

iftImage *iftDynamicTreesObject(iftMImage *features, iftLabeledSet *seeds, iftAdjRel *A) {
    iftDynamicForest *forest = iftCreateDynamicForest(features, A);

    // Dynamic Trees object

    iftFImage *cost = forest->cost;
    iftMImage *mimg = forest->mimg;
    iftImage  *label = forest->label;
    iftImage  *pred = forest->pred;
    iftImage  *root = forest->root, *marker = forest->marker;
    iftFHeap  *heap = forest->heap;

    iftLabeledSet *M = seeds;
    while (M != NULL)
    {
        int p = M->elem;

        if (heap->color[p] == IFT_GRAY) {
            /* p is also a frontier voxel, but the priority is it as a seed. */
            iftRemoveFHeapElem(heap, p);
        }

        cost->val[p]    = 0;
        label->val[p]   = M->label;
        root->val[p]    = p;
        pred->val[p]    = IFT_NIL;
        marker->val[p]  = M->marker;
        iftInsertFHeap(heap, p);
        if (forest->tree_map[M->label] == NULL)
            forest->tree_map[M->label] = iftInsertDynamicTree(&forest->dynamic_roots,
                                                   forest->mimg->m);

        ift_update_dynamic_trees_object(forest, p);

        M = M->next;
    }

    /* Image Foresting Transform */
    while (!iftEmptyFHeap(heap))
    {
        int p      = iftRemoveFHeap(heap);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        if (cost->val[p] != 0) {/* it is not a seed voxel */
            ift_update_dynamic_trees_object(forest, p);
        }

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (heap->color[q] != IFT_BLACK) {
                    float arc_weight = ift_dynamic_dist_object(forest, p, q);

                    float tmp = iftMax(cost->val[p], arc_weight);

                    if (tmp < cost->val[q]) {
                        pred->val[q]     = p;
                        root->val[q]     = root->val[p];
                        label->val[q]    = label->val[p];
                        cost->val[q]     = tmp;

                        if (heap->color[q] == IFT_WHITE) {
                            iftInsertFHeap(heap, q);
                        } else {
                            //iftRemoveFHeapElem(heap, q);
                            //iftInsertFHeap(heap, q);
                            iftGoUpFHeap(heap, heap->pos[q]);
                        }
                    }
                }
            }
        }
    }
    iftResetFHeap(heap);


    // End Dynamic Trees object

    iftImage *_label = iftCopyImage(forest->label);
    iftDestroyDynamicForest(&forest);
    return _label;
}

iftImage *iftDynamicTreesRoot(iftMImage *features, iftLabeledSet *seeds, iftAdjRel *A) {
    iftDynamicForest *forest = iftCreateDynamicForest(features, A);
    iftAddMarkersToDynamicForest(forest, seeds);
    iftDiffDynamicIFT(forest);
    iftImage *label = iftCopyImage(forest->label);
    iftDestroyDynamicForest(&forest);
    return label;
}
