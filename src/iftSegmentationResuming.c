//
// Created by tvspina on 1/12/16.
//

#include "iftSegmentationResuming.h"

#include "ift/core/io/Stream.h"


iftLabeledSet *iftIGraphISF_Resuming_Root(iftIGraph *igraph, const iftImage *input_label, iftDHeap *Q, iftImage *seeds,
                                          double alpha, double beta, double gamma, int niters);
/**
 * @brief Copies all seeds from the new seeds image and inserts them in both current_seeds and segmentation_seeds.
 * They are inserted into segmentation_seeds in such a way that segmentation_seeds holds the original order in which
 * the seeds were added to current_seeds, to ensure that iftDiffWatershed and iftWatershed behave the same.
 *
 * @author Thiago Vallin Spina
 * @date Jan 26, 2016
 *
 * @param max_marker_id Last marker id that was added.
 * @param new_seeds_image The image with newly added seeds.
 * @param current_seeds Returns the newly added seeds.
 * @param segmentation_seeds Returns all seeds from the first iteration until now in orderly fashion.
 */
void iftGetNewlyAddedSeeds(int max_marker_id, iftImage *new_seeds_image, iftImage *new_seeds_image_mk_id,
                           iftLabeledSet **new_seeds, iftLabeledSet **all_seeds);


ulong iftSelectSeedsWithinDistance(iftVoxel center_marker_voxel, const iftImage *gt_image, const iftImage *dist,
                                   const iftImage *root, const iftImage *error_components, iftImage *seed_image,
                                   iftImage *seed_image_mk_id, double marker_length, double min_border_distance,
                                   int seed_lb, int p, int new_seed_id);


void iftCopyNewSeeds(iftImage *seeds_image, iftImage *seeds_image_mk_id, iftImage *seeds_image_copy,
                     iftImage *seeds_image_mk_id_copy);

iftLabeledSet *iftSortSeedsByErrorComponentArea(iftLabeledSet *seeds, iftImage *label, iftImage *gt, iftAdjRel *A) {
    iftImage *error_components = NULL;
    iftHist *hist = NULL;
    iftLabeledSet *aux = NULL, *Sorted = NULL;
    int number_of_components, i, nseeds;
    int *seed_component_area = NULL, *index = NULL, *marker = NULL, *elem = NULL, *handicap = NULL, *seed_label = NULL;

    if(seeds != NULL) {

        fprintf(stderr, "Sorting seeds by error component area\n");

        error_components = iftRelabelSegmentationErrorComponents(gt, label, A);
        number_of_components = iftMaximumValue(error_components);

        nseeds = iftLabeledSetSize(seeds);
        seed_component_area = iftAllocIntArray(nseeds);
        index = iftAllocIntArray(nseeds);
        marker = iftAllocIntArray(nseeds);
        elem = iftAllocIntArray(nseeds);
        seed_label = iftAllocIntArray(nseeds);
        handicap = iftAllocIntArray(nseeds);

        // Computing the area of the error components
        int nbins = number_of_components + 1;
        hist = iftCalcGrayImageHist(error_components, NULL, nbins, number_of_components, 0);

        for(i = 0, aux = seeds; aux != NULL; aux = aux->next, i++) {
            index[i] = i;
            seed_component_area[i] = iftRound(hist->val[error_components->val[aux->elem]]);
            elem[i] = aux->elem;
            seed_label[i] = aux->label;
            marker[i] = aux->marker;
            handicap[i] = aux->handicap;
        }


        // Using a *stable* sort algorithm to ensure that the order of the seeds is preserved. FIFO
        // policy must be in effect when using iftGQueue to this end
        iftBucketSort(seed_component_area, index, nseeds, IFT_DECREASING);

        // Inserting the elemnts in IFT_INCREASING order of area to ensure that the LIFO property of iftLabeledSet
        // makes them be sorted in IFT_DECREASING order.
        for(i = nseeds-1; i >= 0; i--) {
            iftInsertLabeledSetMarkerAndHandicap(&Sorted, elem[index[i]], seed_label[index[i]], marker[index[i]], handicap[index[i]]);
        }
    
        iftDestroyHist(&hist);
        iftDestroyImage(&error_components);

        free(seed_component_area);
        free(index);
        free(marker);
        free(elem);
        free(seed_label);
        free(handicap);
    }

    return Sorted;
}

//iftLabeledSet *iftLabelToForestPixelRobot(iftImage *gradient, iftImage *label, iftAdjRel *A, int seeds_per_iteration,
//                                          float border_dilation_radius_for_candidate_seeds,
//                                          int min_safe_distance_to_border, int min_marker_radius, int max_marker_radius,
//                                          iftLabeledSet *optional_input_seeds, double stopping_threshold,
//                                          int secondary_stopping_threshold,
//                                          uchar (*iftStoppingCriterion)(iftImage *, iftImage *, double),
//                                          uchar (*iftSecondaryStoppingCriterion)(int, int, int))
//{
//    int j = 0, seeds_added;
//    int total_added_seeds = 0;
//    int number_of_objects, number_of_labels;
//    iftLabeledSet  *available_seeds = NULL, *current_seeds = NULL;
////    iftLabeledSet *sorted = NULL;
//    iftBMap *seeds_bmap            = NULL;
//    iftImage *seeds_image         = NULL;
//    iftImage *seeds_image_copy     = NULL;
//    iftImage *current_segmentation = NULL;
//    iftImageForest *forest         = NULL;
//    iftImage *gt_image             = NULL;
//    iftLabeledSet *all_segmentation_seeds = NULL;
//
//    gt_image = label;
//
//    number_of_objects   = iftMaximumValue(gt_image);
//    number_of_labels    = number_of_objects + 1;
//    seeds_per_iteration = number_of_labels*seeds_per_iteration;
//
//    forest = iftCreateImageForest(gradient, A);
//    seeds_bmap  = iftCreateBMap(gradient->n);
//
//    // Computing the seed set from the ground truth image (original erroneously segmented image) that is close to the
//    // objects' boundaries and sorted by gradient value
//    available_seeds  = iftBorderMarkersForPixelSegmentation(gradient, gt_image,
//                                                            border_dilation_radius_for_candidate_seeds);
//
//    // Performing an initial segmentation if the optinal input seed set is passed as a parameter
//    if(optional_input_seeds != NULL) {
//        fprintf(stderr,"Using optional input seeds to perform an initial segmentation\n");
//        iftDiffWatershed(forest, optional_input_seeds, NULL);
//        current_segmentation = forest->label;
//        seeds_image = iftSeedImageFromLabeledSet(optional_input_seeds, gt_image);
//        all_segmentation_seeds = iftCopyOrderedLabeledSet(optional_input_seeds);
//
//    } else {
//        seeds_image = iftCreateImage(gradient->xsize, gradient->ysize, gradient->zsize);
//        iftSetImage(seeds_image, IFT_NIL);
//    }
//
//    do
//    {
//        printf("Rebuilding iteration: %03d\n", ++j);
//
//        // Selecting new seeds according to the available ones for each object
////        sorted = iftSortSeedsByErrorComponentArea(available_seeds, current_segmentation, gt_image, A);
////        iftDestroyLabeledSet(&available_seeds);
////        available_seeds = sorted;
//
////        iftImage *tmp = iftCreateImage(seeds_image->xsize, seeds_image->ysize, seeds_image->zsize);
////
////        for(iftLabeledSet *aux = available_seeds; aux != NULL; aux = aux->next) {
////            if(seeds_image->val[aux->elem] < 0 && (current_segmentation == NULL || current_segmentation->val[aux->elem] != gt_image->val[aux->elem]))
////                tmp->val[aux->elem] = aux->label + 1;
////        }
////
////
////        iftWriteImageByExt(tmp, "available_seeds_%03d.scn", j);
////        iftDestroyImage(&tmp);
//
//        seeds_image_copy = iftCopyImage(seeds_image);
//        seeds_added      = iftMarkersFromMisclassifiedSeeds(seeds_image, available_seeds, seeds_bmap,
//                                                            seeds_per_iteration, number_of_labels, gt_image,
//                                                            current_segmentation, min_safe_distance_to_border,
//                                                            max_marker_radius, min_marker_radius);
//
//        total_added_seeds += seeds_added;
//        fprintf(stderr,"Seeds added %d\n", seeds_added);
//
//        //This produces only the new seeds added this iteration
//        for (int p = 0; p < seeds_image_copy->n; p++)
//        {
//            if (seeds_image_copy->val[p] == seeds_image->val[p])
//                seeds_image_copy->val[p] = IFT_NIL;
//            else
//                seeds_image_copy->val[p] = seeds_image->val[p];
//        }
//
//
//        iftGetNewlyAddedSeeds(total_added_seeds, seeds_image_copy,  &current_seeds, &all_segmentation_seeds);
//
//        // Segmenting the image by adding the new seed set
//        iftDiffWatershed(forest, current_seeds, NULL);
//        current_segmentation = forest->label;
//
//
////        // Preventing negative values
////        tmp = iftLinearStretch(seeds_image, iftMinimumValue(seeds_image), iftMaximumValue(seeds_image), 0, number_of_labels);
////        iftWriteImageByExt(tmp, "seeds_%03d.scn", j);
////        iftDestroyImage(&tmp);
////        iftWriteImageByExt(current_segmentation, "label_%03d.scn", j);
//
//        iftDestroyImage(&seeds_image_copy);
//        iftDestroyLabeledSet(&current_seeds);
//
//    }while(seeds_added != 0 && !iftStoppingCriterion(current_segmentation, gt_image, stopping_threshold)
//           && !iftSecondaryStoppingCriterion(j, total_added_seeds, secondary_stopping_threshold));
//
//    iftDestroyBMap(&seeds_bmap);
//    iftDestroyImage(&seeds_image);
//
//    iftDestroyLabeledSet(&available_seeds);
//    iftDestroyImageForest(&forest);
//
//    return all_segmentation_seeds;
//}

void iftGetNewlyAddedSeeds(int max_marker_id, iftImage *new_seeds_image, iftImage *new_seeds_image_mk_id,
                           iftLabeledSet **new_seeds, iftLabeledSet **all_seeds) {
    iftLabeledSet *inverted_seeds = NULL;
    iftLabeledSet *aux = NULL;

    // We do this weird concatenation of seeds to ensure that the final seed set's order reflects the
    // insertion of seeds performed for DIFT. This minimizes some FIFO/LIFO related issues with iftWatershed
    (*new_seeds) = iftLabeledSetFromSeedImage(new_seeds_image, false);

    // Assigning marker ids related to the iteration when they were inserted to facilitate interactive removal
    for(aux = (*new_seeds); aux != NULL; aux = aux->next) {
        // each label is given a unique marker id related to the iteration when it was added
        aux->marker = max_marker_id + new_seeds_image_mk_id->val[aux->elem];
    }

    inverted_seeds = iftReverseLabeledSet((*all_seeds));
    iftDestroyLabeledSet(all_seeds);
    (*all_seeds) = iftCopyOrderedLabeledSet((*new_seeds));
    iftConcatLabeledSet(all_seeds, &inverted_seeds);
    iftDestroyLabeledSet(&inverted_seeds);
}

void iftIGraphDiffISF_Resuming_Root(iftIGraph *igraph, const iftImage *input_label, iftLabeledSet *new_seeds,
                                    iftSet **trees_for_removal, iftDHeap *Q, double alpha, double beta, double gamma,
                                    bool allow_label_fixing) {
    double      tmp;
    int        s, t, i, p, q, r;
    iftVoxel   u, v;
    double    *pvalue = Q->value;
    iftSet    *adj = NULL, *frontier_nodes = NULL;
    iftSet    *T = NULL; /* Uncomment for the differential version */
    iftLabeledSet *S = NULL;
    bool is_first_iteration = true;

    /* Checking if we are in the first iteration, in order to fix the label map when we are not */
    for(s = 0; s < igraph->nnodes && is_first_iteration; s++) {
        p           = igraph->node[s].voxel;
        is_first_iteration = igraph->pred[p] == IFT_NIL;
    }

    if (trees_for_removal != NULL && *trees_for_removal != NULL)
        frontier_nodes = iftIGraphTreeRemoval(igraph, trees_for_removal, pvalue, IFT_INFINITY_DBL);

    S = new_seeds;
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
                r = igraph->root[p];
                igraph->pvalue[p] = pvalue[s];
                u = iftGetVoxelCoord(igraph->index,p);
                for (i=1; i < igraph->A->n; i++) {
                    v = iftGetAdjacentVoxel(igraph->A,u,i);
                    if (iftValidVoxel(igraph->index,v)){
                        q   = iftGetVoxelIndex(igraph->index,v);
                        t   = igraph->index->val[q];
                        if ((t != IFT_NIL) && (Q->color[t] != IFT_BLACK)){
                            double bndr = 1.0-iftDiracDeltaFunction(input_label->val[r]-input_label->val[q]);

                            tmp = pvalue[s] + pow(alpha*(double)iftFeatDistance(igraph->feat[r],igraph->feat[q],igraph->nfeats)*pow(gamma, bndr)
                                                  + gamma*bndr, beta) + (double)iftVoxelDistance(u,v);

                            if (tmp < pvalue[t])
                            {
                                pvalue[t]            = tmp;
                                /* Uncomment for the differential version */
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
                r = igraph->root[p];
                u = iftGetVoxelCoord(igraph->index,p);
                igraph->pvalue[p] = pvalue[s];

                for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
                    t   = adj->elem;
                    if (Q->color[t] != IFT_BLACK){
                        q   = igraph->node[t].voxel;
                        v   = iftGetVoxelCoord(igraph->index,q);

                        double bndr = 1.0-iftDiracDeltaFunction(input_label->val[r]-input_label->val[q]);

                        tmp = pvalue[s] + pow(alpha*(double)iftFeatDistance(igraph->feat[r],igraph->feat[q],igraph->nfeats)*pow(gamma, bndr)
                                              + gamma*bndr, beta) + (double)iftVoxelDistance(u,v);

                        if (tmp < pvalue[t]) {

                            /* Uncomment for the differential version */
                            if (igraph->label[p] != igraph->label[q]){ /* voxels
	      						    that have
	      						    changed
	      						    labels */
                                iftInsertSet(&T, q);
                            }

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
            break;

        case COMPLETE:
            iftError("Not implemented for complete graphs", "iftIGraphISF_Root");
    }


    /* Uncomment this block and comment the one below for differential
       IFT */

    iftResetDHeap(Q);

    if (allow_label_fixing && !is_first_iteration)
        iftIGraphFixLabelRootMap(igraph, &T);
    iftDestroySet(&T);


    /* End of comment */
}

iftLabeledSet *iftIGraphISF_Resuming_Root(iftIGraph *igraph, const iftImage *input_label, iftDHeap *Q, iftImage *seeds,
                                          double alpha, double beta, double gamma, int niters) {
    int        s, i, p, I, *seed, nseeds;
//    int j, index_seed; /* Uncomment for the differential version */
    float      new_seeds_flag;
    double    *pvalue = Q->value;
    iftSet    *adj = NULL, *S = NULL, *new_seeds = NULL, *frontier_nodes = NULL, *trees_for_removal = NULL;
    iftLabeledSet *new_seeds_labeled = NULL;
    /* Uncomment for the differential version */

    timer *t1, *t2;


    new_seeds_flag = 1.0;

    /* Initial seed set and trivial path initialization to infinity
       cost */

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

        // The marker ID can be the same as the temporary supervoxel label
        iftInsertLabeledSetMarkerAndHandicap(&new_seeds_labeled, S->elem, igraph->label[p], igraph->label[p], 0);

        i++; S = S->next;
    }

    /* differential optimum-path forest computation */
    for (I=0; (I < niters); I++) {

        printf("iteration %d\n",I+1);
        t1 = iftTic();

        iftIGraphDiffISF_Resuming_Root(igraph, input_label, new_seeds_labeled, &trees_for_removal, Q, alpha, beta,
                                       gamma, true);

        /* Recompute new seeds */
        iftDestroySet(&new_seeds);

        iftIGraphISFRecomputeSeeds(igraph, seed, nseeds, &trees_for_removal, &new_seeds, &new_seeds_flag);

        /* Uncomment for differential version */
        iftDestroyLabeledSet(&new_seeds_labeled);

        for(S = new_seeds; S != NULL; S = S->next) {
            s = S->elem;
            p = igraph->node[s].voxel;

            // The marker ID can be the same as the temporary supervoxel label
            iftInsertLabeledSetMarkerAndHandicap(&new_seeds_labeled, S->elem, igraph->label[p], igraph->label[p], 0);
        }
        /* End of comment */

        t2 = iftToc();
        iftPrintCompTime(t1,t2,"Computational time for iteration %d",I+1);

        /* Uncomment this block and comment the one above for
           non-differential IFT */

//        iftResetDHeap(Q);
//        iftDestroySet(&trees_for_removal);
//        iftDestroySet(&new_seeds);
//
//        for (s=0; s < igraph->nnodes; s++) {
//            p = igraph->node[s].voxel;
//            pvalue[s]       = IFT_INFINITY_DBL;
//        }
//        for (int i = 0; i < nseeds; ++i) {
//            s = seed[i];
//            iftInsertSet(&new_seeds,s);
//            p = igraph->node[s].voxel;
//            igraph->label[p] = i+1;
//        }


        /* End of comment */

    }

    iftDestroyLabeledSet(&new_seeds_labeled);
    iftDestroySet(&adj);
    iftDestroySet(&S);
    iftDestroySet(&new_seeds);
    iftDestroySet(&frontier_nodes);
    iftDestroySet(&trees_for_removal);
    free(seed);

    iftLabeledSet *resuming_seeds = NULL;

    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;

        // Copying the seed set, when a seed is encountered
        if(igraph->pred[p] == IFT_NIL) {
            // The final seed label refers to the label in the presegmentation, while the marker ID refers to the
            // temporarily assigned superpixel label
            iftInsertLabeledSetMarkerAndHandicap(&resuming_seeds, p, input_label->val[p], igraph->label[p], 0);
        }

        // Copying the input label from the roots to all nodes
        igraph->label[p] = input_label->val[igraph->root[p]];
    }

    return resuming_seeds;
}

iftLabeledSet *iftLabelToForestISF_Root(const iftImage *img, const iftImage *input_label, double alpha, double beta,
                                        double gamma, int niters, int *nseeds, iftIGraph **igraph_out, iftDHeap **Q_out) {
    iftLabeledSet *resuming_seeds = NULL;
    iftImage *seeds, *mask1;
    iftMImage *mimg;
    iftAdjRel *A;
    iftIGraph *igraph = NULL;
    iftDHeap *Q = NULL;
    double *pvalue = NULL;
    iftBoundingBox bb;
    iftVoxel v;
    int nsuperpixels, bb_n;

    /* Compute ISF superpixels */
    if (iftIs3DImage(img)){
        A      = iftSpheric(1.0);
    } else {
        A      = iftCircular(1.0);
    }

    if (iftIsColorImage(img)){
        mimg   = iftImageToMImage(img, LABNorm_CSPACE);
    } else {
        mimg   = iftImageToMImage(img, GRAY_CSPACE);
    }

    mask1  = iftSelectImageDomain(mimg->xsize,mimg->ysize,mimg->zsize);

    /* minima of a basins manifold in that domain */
    igraph = iftImplicitIGraph(mimg,mask1,A);

    int xsize, ysize, zsize;
    bb = iftMinBoundingBox(input_label, NULL);

    xsize = bb.end.x - bb.begin.x + 1;
    ysize = bb.end.y - bb.begin.y + 1;
    zsize = bb.end.z - bb.begin.z + 1;

    bb_n = xsize*ysize*zsize;

//    bb.begin.x = iftMax(0, iftRound(bb.begin.x - xsize*0.1));
//    bb.begin.y = iftMax(0, iftRound(bb.begin.y - ysize*0.1));
//    bb.begin.z = iftMax(0, iftRound(bb.begin.z - zsize*0.1));
//
//    bb.end.x = iftMin(img->xsize-1, iftRound(bb.end.x + xsize*0.1));
//    bb.end.y = iftMin(img->ysize-1, iftRound(bb.end.y + ysize*0.1));
//    bb.end.z = iftMin(img->zsize-1, iftRound(bb.end.z + zsize*0.1));



    int obj_size = 0;
    for(v.z = 0; v.z < img->zsize; v.z++) {
        for(v.y = 0; v.y < img->ysize; v.y++) {
            for (v.x = 0; v.x < img->xsize; v.x++) {
                int p = iftGetVoxelIndex(mask1, v);
//                // Zeroing everything outside the bounding box
//                if(!(v.x >= bb.begin.x && v.x <= bb.end.x &&
//                        v.y >= bb.begin.y && v.y <= bb.end.y &&
//                        v.z >= bb.begin.z && v.z <= bb.end.z)) {
//                    mask1->val[p] = 0;
//                }
                if(input_label->val[p] > 0) {
                    obj_size++;
                }
            }
        }
    }

//    nsuperpixels = iftRound((obj_size) / (200 * gamma));

    // Adding extra supervoxel seeds outside the bounding box of the object. The original formulation
    // considered obj_size/(200*gamma) to determine a fixed number of supervoxel seeds only in the object's
    // bounding box (expanded by a small factor, see Barreto et. al. 2017, ISBI)
    nsuperpixels = iftRound((obj_size)/ (200*gamma))*(1.0 + (img->n - bb_n) / ((float)img->n));

    fprintf(stderr,"Number of selected superpixels: %d\n", nsuperpixels);

    /* seed sampling for ISF */
//    seeds   = iftGridSampling(mimg,mask1,nsuperpixels);
    seeds   = iftAltMixedSampling(mimg,mask1,nsuperpixels);

    *nseeds = iftNumberOfElements(seeds);

    iftDestroyImage(&mask1);
    iftDestroyMImage(&mimg);

    pvalue = iftAllocDoubleArray(igraph->nnodes);
    Q = iftCreateDHeap(igraph->nnodes, pvalue);

    resuming_seeds = iftIGraphISF_Resuming_Root(igraph, input_label, Q, seeds, alpha, beta, gamma, niters);

    if(igraph_out == NULL) {
        iftDestroyIGraph(&igraph);
    } else {
        *igraph_out = igraph;
    }

    if(Q_out == NULL) {
        iftFree(Q->value);
        iftDestroyDHeap(&Q);
    } else {
        *Q_out = Q;
    }

    iftDestroyImage(&seeds);
    iftDestroyAdjRel(&A);

    return resuming_seeds;
}


//
//iftLabeledSet * iftIGraphDiffISF_Resuming_Mean(iftIGraph *igraph, iftImage *input_label, iftImage *seeds, double alpha,
//                                               double beta,
//                                               double gamma, int niters) {
//    int        s, i, p, I, *seed, nseeds, j;
//    float      new_seeds_flag;
//    iftDHeap  *Q;
//    double    *pvalue = iftAllocDoubleArray(igraph->nnodes);
//    iftSet    *adj=NULL, *S=NULL, *new_seeds=NULL, *frontier_nodes=NULL, *trees_for_removal=NULL;
//
//    float     **seed_features; // NEW
//
//    timer *t1, *t2;
//
//
//    new_seeds_flag = 1.0;
//
//    /* Initial seed set and trivial path initialization to infinity
//       cost */
//
//    Q      = iftCreateDHeap(igraph->nnodes, pvalue);
//
//    nseeds = 0;
//    for (s=0; s < igraph->nnodes; s++) {
//        p               = igraph->node[s].voxel;
//        pvalue[s]       = igraph->pvalue[p] = IFT_INFINITY_DBL;
//        igraph->pred[p] = IFT_NIL;
//        if (seeds->val[p]!=0){
//            iftInsertSet(&new_seeds,s);
//            nseeds++;
//        }
//    }
//
//    /* Alloc seed features NEW */
//    seed_features = (float**) calloc(nseeds ,sizeof(float*));
//    for (i = 0; i < nseeds; ++i)
//        seed_features[i] = iftAllocFloatArray(igraph->nfeats);
//
//    seed  = iftAllocIntArray(nseeds);
//    S     = new_seeds; i = 0;
//    while (S != NULL) {
//        seed[i] = S->elem;
//        p       = igraph->node[seed[i]].voxel;
//        igraph->label[p] = i+1;
//
//        // set initial seed feature NEW
//        for (j = 0; j < igraph->nfeats; ++j)
//            seed_features[i][j] = igraph->feat[p][j];
//
//        i++; S = S->next;
//    }
//
//
//    /* differential optimum-path forest computation */
//    for (I=0; (I < niters); I++) {
//
//        printf("iteration %d\n",I+1);
//        t1 = iftTic();
//
//        iftIGraphDiffISF_IFT_Resuming_Mean(igraph, input_label, Q, pvalue, new_seeds, trees_for_removal, seed_features,
//                                           alpha, beta, gamma, I == 0);
//
//        /* Recompute new seeds */
//
//        iftIGraphISFRecomputeSeedsUsingSpatialInformation(igraph, seed, nseeds, &trees_for_removal, &new_seeds, &new_seeds_flag, seed_features);
//
//        t2 = iftToc();
//        iftPrintCompTime(t1,t2,"Computational time for iteration %d",I+1);
//
//        /* Uncomment this block and comment the one above for
//           non-differential IFT */
//        /*
//        iftResetDHeap(Q);
//        iftDestroySet(&trees_for_removal);
//        iftDestroySet(&new_seeds);
//
//        for (s=0; s < igraph->nnodes; s++) {
//           p = igraph->node[s].voxel;
//           pvalue[s]       = IFT_INFINITY_DBL;
//        }
//        for (int i = 0; i < nseeds; ++i) {
//           s = seed[i];
//           iftInsertSet(&new_seeds,s);
//           p = igraph->node[s].voxel;
//           igraph->label[p] = i+1;
//        }
//        */
//
//        /* End of comment */
//
//    }
//
//    /* Free seed_features NEW */
//    for (i = 0; i < nseeds; ++i)
//        free(seed_features[i]);
//    free(seed_features);
//
//    iftDestroySet(&adj);
//    iftDestroySet(&S);
//    iftDestroySet(&new_seeds);
//    iftDestroySet(&frontier_nodes);
//    iftDestroySet(&trees_for_removal);
//    iftDestroyDHeap(&Q);
//    free(pvalue);
//    free(seed);
//
//
//    iftLabeledSet *resuming_seeds = NULL;
//
//    for (s=0; s < igraph->nnodes; s++) {
//        p = igraph->node[s].voxel;
//        if(igraph->pred[p] == IFT_NIL) {
//            iftInsertLabeledSet(&resuming_seeds, )
//        }
//    }
//}
//
//void iftIGraphDiffISF_IFT_Resuming_Mean(iftIGraph *igraph, const iftImage *input_label, iftDHeap *Q, double *pvalue,
//                                        iftSet *new_seeds, iftSet *trees_for_removal, float **seed_features,
//                                        double alpha, double beta, double gamma, bool is_first_iteration) {
//    double      tmp;
//    int        r, s, t, i, p, q, index_seed;
//    iftVoxel    u,v;
//    iftSet    *adj=NULL, *frontier_nodes=NULL;
//    iftSet    *T = NULL;
//
//    if (trees_for_removal != NULL)
//        frontier_nodes = iftIGraphTreeRemoval(igraph, &trees_for_removal, pvalue, IFT_INFINITY_DBL);
//
//    while (new_seeds != NULL) { /* insert seeds in the priority queue Q with cost zero */
//        s = iftRemoveSet(&new_seeds);
//        p = igraph->node[s].voxel;
//        pvalue[s] = igraph->pvalue[p] = 0;
//        igraph->root[p] = p;
//        igraph->pred[p] = IFT_NIL;
//        iftInsertDHeap(Q,s);
//    }
//
//    while (frontier_nodes != NULL) { /* insert frontier nodes in Q to represent the previous seeds */
//        s = iftRemoveSet(&frontier_nodes);
//
//        if (Q->color[s] == IFT_WHITE)
//            iftInsertDHeap(Q,s);
//    }
//
//    switch(igraph->type) {
//
//        case IMPLICIT:
//
//
//            while (!iftEmptyDHeap(Q)) {
//                s = iftRemoveDHeap(Q);
//                p = igraph->node[s].voxel;
//                r = igraph->root[p];
//                index_seed = igraph->label[p] - 1;
//
//                igraph->pvalue[p] = pvalue[s];
//                u = iftGetVoxelCoord(igraph->index,p);
//                for (i=1; i < igraph->A->n; i++) {
//                    v = iftGetAdjacentVoxel(igraph->A,u,i);
//                    if (iftValidVoxel(igraph->index,v)){
//                        q   = iftGetVoxelIndex(igraph->index,v);
//                        t   = igraph->index->val[q];
//                        if ((t != IFT_NIL) && (Q->color[t] != IFT_BLACK)){
//                            long bndr = labs(input_label->val[r]-input_label->val[q]);
//
//                            tmp = pvalue[s] + pow(alpha*(double)iftFeatDistance(seed_features[index_seed],igraph->feat[q],igraph->nfeats)*pow(gamma, bndr)
//                                                  + gamma*bndr, beta) + (double)iftVoxelDistance(u,v);
//
//                            if (tmp < pvalue[t])
//                            {
//                                pvalue[t]            = tmp;
///* Uncomment for the differential version */
//                                if (igraph->label[p] != igraph->label[q]){ /* voxels that have changed labels */
//                                    iftInsertSet(&T, q);
//                                }
///* end Uncomment */
//
//                                igraph->root[q]      = igraph->root[p];
//                                igraph->label[q]     = igraph->label[p];
//                                igraph->pred[q]      = p;
//                                if (Q->color[t] == IFT_GRAY){
//                                    iftGoUpDHeap(Q, Q->pos[t]);
//                                } else {
//                                    iftInsertDHeap(Q,t);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//
//            break;
//
//        case EXPLICIT:
//
//            while (!iftEmptyDHeap(Q)) {
//                s = iftRemoveDHeap(Q);
//                p = igraph->node[s].voxel;
//                r = igraph->root[p];
//                index_seed = igraph->label[p] - 1;
//                u = iftGetVoxelCoord(igraph->index,p);
//                igraph->pvalue[p] = pvalue[s];
//
//                for (adj=igraph->node[s].adj; adj != NULL; adj=adj->next) {
//                    t   = adj->elem;
//                    if (Q->color[t] != IFT_BLACK){
//                        q   = igraph->node[t].voxel;
//                        v   = iftGetVoxelCoord(igraph->index,q);
//
//                        long bndr = labs(input_label->val[r]-input_label->val[q]);
//
//                        tmp = pvalue[s] + pow(alpha*(double)iftFeatDistance(seed_features[index_seed],igraph->feat[q],igraph->nfeats)*pow(gamma, bndr)
//                                              + gamma*bndr, beta) + (double)iftVoxelDistance(u,v);
//                        if (tmp < pvalue[t]) {
///* Uncomment for the differential version */
//                            if (igraph->label[p] != igraph->label[q]){ /* voxels that have changed labels */
//                                iftInsertSet(&T, q);
//                            }
///* end Uncomment */
//                            pvalue[t]            = tmp;
//                            igraph->root[q]      = igraph->root[p];
//                            igraph->label[q]     = igraph->label[p];
//                            igraph->pred[q]      = p;
//                            if (Q->color[t] == IFT_GRAY){
//                                iftGoUpDHeap(Q, Q->pos[t]);
//                            } else {
//                                iftInsertDHeap(Q,t);
//                            }
//                        }
//                    }
//                }
//            }
//            break;
//
//        case COMPLETE:
//            iftError("Not implemented for complete graphs","iftIGraphISF_Mean");
//    }
//
//
///* Uncomment this block and comment the one below for differential
//   IFT */
//
//    iftResetDHeap(Q);
//    if (!is_first_iteration)
//        iftIGraphFixLabelRootMap(igraph, &T);
//    iftDestroySet(&T);
//
//}

//iftLabeledSet *iftLabelToForestPixelRobotByErrorComponent(iftImage *gradient, iftImage *label, iftAdjRel *A, int seeds_per_iteration,
//                                          float border_dilation_radius_for_candidate_seeds,
//                                          int min_safe_distance_to_border, int min_marker_radius, int max_marker_radius,
//                                          double stopping_threshold,
//                                          uchar (*iftStoppingCriterion)(iftImage *, iftImage *, double))
//{
//    int p, j = 0, seeds_added;
//    int number_of_objects;
//    iftLabeledSet  *available_seeds = NULL, *current_seeds = NULL;
//    iftBMap *seeds_bmap            = NULL;
//    iftImage  *seeds_image         = NULL;
//    iftImage *seeds_image_copy     = NULL;
//    iftImage *current_segmentation = NULL;
//    iftImageForest *forest         = NULL;
//    iftImage *gt_image             = NULL;
//    iftLabeledSet *segmentation_seeds = NULL;
//    iftLabeledSet *inverted_seeds = NULL;
//
//    gt_image = label;
//
//    number_of_objects = iftMaximumValue(gt_image);
//    seeds_per_iteration = (number_of_objects+1)*seeds_per_iteration;
//
//    forest = iftCreateImageForest(gradient, A);
//    seeds_bmap  = iftCreateBMap(gradient->n);
//    seeds_image = iftCreateImage(gradient->xsize, gradient->ysize, gradient->zsize);
//    iftSetImage(seeds_image, IFT_NIL);
//
//    seeds_added = IFT_NIL;
//
//    // Computing the seed set from the ground truth image (original erroneously segmented image) that is close to the
//    // objects' boundaries and sorted by gradient value
//    available_seeds  = iftBorderMarkersForPixelSegmentation(gradient, gt_image,
//                                                            border_dilation_radius_for_candidate_seeds);
//
//    do
//    {
//        printf("Rebuilding iteration: %03d\n", ++j);
//
//        seeds_image_copy = iftCopyImage(seeds_image);
//        // Selecting new seeds according to the available ones for each object
//        seeds_added      = iftMarkersFromMisclassifiedSeeds(seeds_image, available_seeds, seeds_bmap, seeds_per_iteration,
//                                                            number_of_objects + 1, gt_image, current_segmentation,
//                                                            min_safe_distance_to_border, max_marker_radius,
//                                                            min_marker_radius);
//
//        //This produces only the new seeds added this iteration
//        for (p = 0; p < seeds_image_copy->n; p++)
//        {
//            if (seeds_image_copy->val[p] == seeds_image->val[p])
//                seeds_image_copy->val[p] = IFT_NIL;
//            else
//                seeds_image_copy->val[p] = seeds_image->val[p];
//        }
//
//
//        // We do this weird concatenation of seeds to ensure that the final seed set's order reflects the
//        // insertion of seeds performed for DIFT. This minimizes some FIFO/LIFO related issues with iftWatershed
//        current_seeds = iftLabeledSetFromSeedImage(seeds_image_copy, false);
//        inverted_seeds = iftReverseLabeledSet(segmentation_seeds);
//        iftDestroyLabeledSet(&segmentation_seeds);
//        segmentation_seeds = iftCopyOrderedLabeledSet(current_seeds);
//
//        iftConcatLabeledSet(&segmentation_seeds, &inverted_seeds);
//
//        // Segmenting the image by adding the new seed set
//        iftDiffWatershed(forest, current_seeds, NULL);
//        current_segmentation = forest->label;
//
////        iftImage *tmp = iftAddValue(seeds_image, 1);
////        iftWriteImage(forest->label, "label_%03d.scn", j);
////        iftWriteImage(tmp, "seeds_%03d.scn", j);
////        iftWriteImage(forest->img, "grad_%03d.scn", j);
////        iftDestroyImage(&tmp);
//
////        iftStopLabelReconstructionByASSD(current_segmentation, gt_image, stopping_threshold);
//
//        iftDestroyImage(&seeds_image_copy);
//        iftDestroyLabeledSet(&current_seeds);
//    }while(seeds_added != 0 && !iftStoppingCriterion(current_segmentation, gt_image, stopping_threshold));
//
//    iftDestroyBMap(&seeds_bmap);
//    iftDestroyImage(&seeds_image);
//
//    iftDestroyLabeledSet(&available_seeds);
//    iftDestroyImageForest(&forest);
//
//    return segmentation_seeds;
//}

uchar iftStopLabelReconstructionByDice(iftImage *label, iftImage *gt, double threshold) {
    int l, nobjects = 0;
    double min_dice = IFT_INFINITY_DBL;

    nobjects = iftMaximumValue(gt);

    iftDblArray *dices = iftDiceSimilarityMultiLabel(label, gt, nobjects);

    // Computing the minimum dices accuracy among all objects
    for(l = 1; l <= nobjects; l++) {
        min_dice = iftMin(min_dice, dices->val[l]);
        fprintf(stderr,"Dice for label %02d: %lf\n", l, dices->val[l]);
    }
    iftDestroyDblArray(&dices);

    fprintf(stderr,"Min dice: %lf\n", min_dice);

    // Stopping if the minimum Dice is above a given  threshold (i.e., we guarantee that all objects have a certain
    // accuracy)
    return min_dice >= threshold;
}

uchar iftStopLabelReconstructionByASSD(iftImage *label, iftImage *gt, double threshold) {
    int l, nobjects = 0;
    float max_assd = IFT_INFINITY_FLT_NEG;

    nobjects = iftMaximumValue(gt);

    iftDblArray *assd = iftASSDMultiLabel(label, gt, nobjects);

    // Computing the maximum ASSD error among all objects
    for(l = 1; l <= nobjects; l++) {
        max_assd = iftMax(max_assd, assd->val[l]);
        fprintf(stderr,"ASSD for label %02d: %lf\n", l, assd->val[l]);
    }
    iftDestroyDblArray(&assd);

    fprintf(stderr,"Max ASSD: %lf\n", max_assd);


    // Stopping if the maximum ASSD error is below a given threshold (i.e., we guarantee that all objects have a certain
    // accuracy)
    return max_assd <= threshold;
}


uchar iftStopLabelReconstructionByNumInterations(int cur_iteration, int cur_num_seeds, int threshold) {
    return cur_iteration >= threshold;
}

uchar iftStopLabelReconstructionByNumSeeds(int cur_iteration, int cur_num_seeds, int threshold) {
    return cur_num_seeds >= threshold;
}

/* EDT-based robot */

iftLabeledSet * iftLabelBorderSet(iftImage *label, iftAdjRel *A) {
    iftLabeledSet *border = NULL;
    iftVoxel u, v;
    int p, q;

    for(p = 0; p < label->n; p++) {
        u = iftGetVoxelCoord(label, p);

        for (int i = 1; i < A->n; i++) {
            v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(label, v)) {
                q = iftGetVoxelIndex(label, v);
                if (label->val[q] != label->val[p]) {
                    iftInsertLabeledSet(&border, p, label->val[p]);
                    break;
                }
            }
            else if (label->val[p] != 0) {
                iftInsertLabeledSet(&border, p, label->val[p]);
                break;
            }
        }
    }
    return border;
}



iftLabeledSet *iftBorderMarkersDistAndRootForPixelSegmentation(iftImage *grad_image, iftImage *gt_image,
                                                               double border_distance, iftImage **dist,
                                                               iftImage **root)
{
    int p, *grad = NULL;
    iftAdjRel *A1,*A2;
    iftSet    *B=NULL;
    iftGQueue *Q=NULL;
    iftLabeledSet *border = NULL, *S=NULL;

    /* Initialization */

    if (iftIs3DImage(grad_image)){
        A1 = iftSpheric(1.0);
        A2 = iftSpheric(sqrtf(3.0));
    }else
    {
        A1 = iftCircular(1.0);
        A2 = iftCircular(1.5);
    }

    /* Computing spels on the border of every label *including* the background */
    border = iftLabelBorderSet(gt_image, A1);
    /* Computing the squared Euclidean distance transform from the object set */
    B = iftLabeledSetElemsToSet(border);
    iftDestroyLabeledSet(&border);
    iftMultiLabelShellDistTransFromSet(B, gt_image, A2, IFT_BOTH, border_distance, dist, root);

    /* Sort seed in the decreasing order of gradient values (due to
       the LIFO nature of iftLabeledSet) */
    grad = iftAllocIntArray(grad_image->n);
    Q     = iftCreateGQueue(iftMaximumValue(grad_image)+1, grad_image->n, grad);
    iftSetRemovalPolicy(Q,MAXVALUE);

    iftDestroyAdjRel(&A2);

    if (iftIs3DImage(grad_image)){
        A2 = iftSpheric(1.0);
    }else
    {
        A2 = iftCircular(1.0);
    }

    while(B != NULL) {
        p = iftRemoveSet(&B);
        // Since the ground truth's border might be a little bit displaced, we compute the median gradient value around
        // the the as its weight, aiming to minimize it along the boundary at each iteration of the robot
        grad[p] = iftMedianValueInAdjacency(grad_image, p, A2);
        iftInsertGQueue(&Q, p);
    }

    while(!iftEmptyGQueue(Q)) {
        p=iftRemoveGQueue(Q);

        iftInsertLabeledSet(&S,p,gt_image->val[p]);
    }

    iftDestroyGQueue(&Q);
    free(grad);


    return(S);
}

ulong iftSelectSeedsBasedOnEDT(iftVoxel center_marker_voxel, iftImage *gt_image, iftImage *dist, iftImage *root,
                               iftImage *error_components, iftImage *seed_image, iftImage *seed_image_mk_id,
                               double marker_length, double min_border_distance, int new_seed_id,
                               bool select_seeds_on_axes_only) {
    int seed_lb, center;
    ulong nselected_seeds = 0;
    iftVoxel v_p, v_first, v_last;

    center = iftGetVoxelIndex(gt_image, center_marker_voxel);
    seed_lb = gt_image->val[center];


    if(select_seeds_on_axes_only) {
        v_first.x = v_last.x = center_marker_voxel.x;
        v_first.y = v_last.y = center_marker_voxel.y;
        v_first.z = v_last.z = center_marker_voxel.z;

        for(v_p.z = v_first.z; v_p.z <= v_last.z; v_p.z++) {
            for (v_p.y = 0; v_p.y <= dist->ysize - 1; v_p.y++) {
                for (v_p.x = 0; v_p.x <= dist->xsize - 1; v_p.x++) {

                    int p = iftGetVoxelIndex(dist, v_p);

                    nselected_seeds += iftSelectSeedsWithinDistance(center_marker_voxel, gt_image, dist, root,
                                                                    error_components, seed_image, seed_image_mk_id,
                                                                    marker_length, min_border_distance, seed_lb, p, new_seed_id);


                }
            }
        }

        for(v_p.y = v_first.y; v_p.y <= v_last.y; v_p.y++) {
            for (v_p.z = 0; v_p.z <= dist->zsize - 1; v_p.z++) {
                for (v_p.x = 0; v_p.x <= dist->xsize - 1; v_p.x++) {

                    int p = iftGetVoxelIndex(dist, v_p);

                    nselected_seeds += iftSelectSeedsWithinDistance(center_marker_voxel, gt_image, dist, root,
                                                                    error_components, seed_image, seed_image_mk_id,
                                                                    marker_length, min_border_distance, seed_lb, p, new_seed_id);


                }
            }
        }

        for(v_p.x = v_first.x; v_p.x <= v_last.x; v_p.x++) {
            for (v_p.z = 0; v_p.z <= dist->zsize - 1; v_p.z++) {
                for (v_p.y = 0; v_p.y <= dist->ysize - 1; v_p.y++) {

                    int p = iftGetVoxelIndex(dist, v_p);

                    nselected_seeds += iftSelectSeedsWithinDistance(center_marker_voxel, gt_image, dist, root,
                                                                    error_components, seed_image, seed_image_mk_id,
                                                                    marker_length, min_border_distance, seed_lb, p, new_seed_id);


                }
            }
        }
    } else {
        v_first.x = v_first.y = v_first.z = 0;
        v_last.x  = dist->xsize - 1;
        v_last.y  = dist->ysize - 1;
        v_last.z  = dist->zsize - 1;

        for(v_p.z = v_first.z; v_p.z <= v_last.z; v_p.z++) {
            for (v_p.y = v_first.y; v_p.y <= v_last.y; v_p.y++) {
                for (v_p.x = v_first.x; v_p.x <= v_last.x; v_p.x++) {

                    int p = iftGetVoxelIndex(dist, v_p);

                    nselected_seeds += iftSelectSeedsWithinDistance(center_marker_voxel, gt_image, dist, root,
                                                                    error_components, seed_image, seed_image_mk_id,
                                                                    marker_length, min_border_distance, seed_lb, p, new_seed_id);


                }
            }
        }
    }


    return nselected_seeds;
}

ulong iftSelectSeedsWithinDistance(iftVoxel center_marker_voxel, const iftImage *gt_image, const iftImage *dist,
                                   const iftImage *root, const iftImage *error_components, iftImage *seed_image,
                                   iftImage *seed_image_mk_id, double marker_length, double min_border_distance,
                                   int seed_lb, int p, int new_seed_id) {
    int r;
    iftVoxel v_r;
    ulong nselected_seeds = 0;
    double root_dist = 0.0;
    int center = iftGetVoxelIndex(dist, center_marker_voxel);

    /* If the p's EDT value is not infinity but is above the minimum allowed distance to the border then it is within
     * the maximum marker distance allowed and is a candidate seed voxel. We thus evalute it if it hasn't been
     * selected as seed before.
     */
    if (dist->val[p] != IFT_INFINITY_INT && dist->val[p] >= min_border_distance &&
        seed_image->val[p] == IFT_NIL) {
        /* We check if p's root has the same label of the given seed pixel that represents the current marker
         * being selected and then select p as seed if its root r is within the marker radius
         */
        r = root->val[p];

        if (gt_image->val[r] == seed_lb && error_components->val[center] == error_components->val[r]) {

            v_r = iftGetVoxelCoord(dist, r);

            root_dist = iftVoxelDistance(v_r, center_marker_voxel);

            // Marker length should be the end-to-end distance, hence we divide it by 2
            if (root_dist <= marker_length / 2.0) {
                seed_image->val[p] = gt_image->val[p];
                seed_image_mk_id->val[p] = new_seed_id;
                nselected_seeds++;
            }
        }
    }
    return nselected_seeds;
}

/**
 * @brief This function selects seed voxels on the border between foreground and background according to the euclidean
 * distance that each pixel has to it. It takes as input a voxel on the border of the object or the background, which
 * serves as the "center" of the marker. Then, it selects as seed voxels those whose root on the border is within a
 * given radius to the center voxel and belong to the same error component. Hence, the resulting marker follows the
 * contour instead of being a simple circle/sphere.
 *
 * @author Thiago Vallin Spina
 *
 * @param center_marker_voxel The voxel on the object's or background's border that has been selected as the marker's "center".
 * @param gt_image The ground truth image.
 * @param dist The EDT distance image. NOTE: This function selects ALL voxels with finite distance as seeds, hence, <dist>
 * should be created only until a given radius from the border to ensure that the markers will be small. Ideally, this
 * distance should also be <marker_radius>.
 * @param root The EDT's root image.
 * @param seed_image The seed image where the seeds will be placed on.
 * @param error_components The image with labeled error components.
 * @param marker_radius The radius used to select the roots around the center voxel.
 * @param min_border_distance A minimum distance to the border that must be respected in order to select seed voxels.
 *
 */
//ulong iftMarkRootsWithinRadius(iftVoxel center_marker_voxel, iftImage *gt_image, iftImage *dist, iftImage *root,
//                               iftImage *seed_image, iftImage *error_components, iftGQueue *Q, iftAdjRel *A,
//                               iftBMap *selected_roots, double marker_length) {
//    int i, seed_lb, p, r, center;
//    ulong nselected_seeds = 0;
//    double root_dist;
//    iftVoxel v;
//
//    iftResetGQueue(Q);
//
//    center = iftGetVoxelIndex(gt_image, center_marker_voxel);
//    seed_lb = gt_image->val[center];
//
//    Q->L.value[center] = 0;
//    iftInsertGQueue(&Q, center);
//
//    while(iftEmptyGQueue(Q)) {
//        p = iftRemoveGQueue(Q);
//        /* If the p's EDT value is not infinity but is above the minimum allowed distance to the border then it is within
//         * the maximum marker distance allowed and is a candidate seed voxel. We thus evalute it if it hasn't been
//         * selected as seed before.
//         */
//        if (dist->val[p] != INFINITY_INT && root->val[p] == p && iftBMapValue(selected_roots, ) seed_image->val[p] == IFT_NIL) {
//            /* We check if p's root has the same label of the given seed pixel that represents the current marker
//             * being selected and then select p as seed if its root r is within the marker radius
//             */
//            r = root->val[p];
//
//            if (gt_image->val[r] == seed_lb && error_components->val[center] == error_components->val[r]) {
//                v = iftGetVoxelCoord(dist, r);
//                root_dist = iftVoxelDistance(v, center_marker_voxel);
//
//                if (root_dist <= marker_length) {
//                    seed_image->val[p] = gt_image->val[p];
//                    nselected_seeds++;
//                }
//            }
//        }
//    }
//
//    return nselected_seeds;
//}

/**
 * @brief This function recursively removes seeds from a given set if they are within a given radius of a given voxel.
 *
 * @param Thiago Vallin Spina
 *
 * @param S Seed set.
 * @param center The center voxel.
 * @param label The center voxel's label
 * @param img The image in which the seed voxels have been selected.
 * @param radius The radius for removal.
 */
//void iftRemoveSeedsWithinRadius(iftLabeledSet **S, iftVoxel center, int label, iftImage *img, double radius) {
//    iftLabeledSet *aux = NULL;
//    if(S == NULL || *S == NULL) return; //(base case 1)
//
//    // Recursively removing the next elements that are within <radius> from <center>
//    iftRemoveSeedsWithinRadius(&(*S)->next, center, label, img, radius);
//
//    // Checking whether the head of the list must be removed as well (base case 2)
//    if((*S)->label == label && iftVoxelDistance(iftGetVoxelCoord(img, (*S)->elem), center) <= radius) {
//        aux = (*S)->next; // Saving the new head
//        free(*S); // Freeing the old one
//
//        *S = aux; // Updating the list's head
//    }
//}

void iftRemoveSeedsWithinRadius(iftLabeledSet **S, iftVoxel center, int label, iftImage *img, double radius) {
    iftLabeledSet *aux = NULL, *prev = NULL, *tmp = NULL;
//    if(S == NULL || *S == NULL) return; //(base case 1)

    // Recursively removing the next elements that are within <radius> from <center>
//    iftRemoveSeedsWithinRadius(&(*S)->next, center, label, img, radius);

    // Checking whether the head of the list must be removed as well (base case 2)
    prev = *S;
    for(aux = *S; aux != NULL; ) {
        if (aux->label == label && iftVoxelDistance(iftGetVoxelCoord(img, aux->elem), center) <= radius) {
            if(aux == *S) {
                aux = aux->next;
                free(*S);
                *S = aux;
                prev = *S;
            } else {
                tmp = aux;
                aux = aux->next;
                prev->next = aux;
                free(tmp);
            }
        }
        else
        {
            prev = aux;
            aux = aux->next;
        }
    }
}


iftSet **iftErrorComponentsPerSegmentationLabelSortedByArea(iftImage *gt, iftImage *error_components) {
    int p, lb, id, number_of_labels, number_of_components, *index = NULL, *area = NULL;
    int *selected_components = NULL;
    iftSet **component_ids_per_label = NULL, *aux = NULL;
    iftHist *hist = NULL;

    number_of_labels = iftMaximumValue(gt) + 1;
    number_of_components = iftMaximumValue(error_components);

    // Computing the area of the error components
    int nbins = number_of_components + 1;
    hist = iftCalcGrayImageHist(error_components, NULL, nbins, number_of_components, 0);

    // Allocating an iftSet array that will store the component ids per label
    component_ids_per_label = (iftSet**)calloc(number_of_labels, sizeof(iftSet*));
    selected_components = iftAllocIntArray(number_of_components);

    for(lb = 0; lb < number_of_labels; lb++)
        component_ids_per_label[lb] = NULL;

    //
    for(p = 0; p < gt->n; p++) {
        if(error_components->val[p] > 0 && !selected_components[error_components->val[p]-1]) {
            selected_components[error_components->val[p]-1] = 1;

            iftInsertSet(&component_ids_per_label[gt->val[p]], error_components->val[p]);
        }
    }

    // Sorting the components for each label by area
    for(lb = 0; lb < number_of_labels; lb++) {
        // Determining the number of components for the current label
        int n = iftSetSize(component_ids_per_label[lb]);

        if(n > 0) {
            area = iftAllocIntArray(n);
            index = iftAllocIntArray(n);

            aux = component_ids_per_label[lb];

            // Filling area and component id arrays
            for (id = 0; id < n; id++) {
                area[id] = iftRound(hist->val[aux->elem]);
                index[id] = aux->elem;

                aux = aux->next;
            }

            // Soring in increasing order since iftInsertSet will reverse the components and we want the components
            // to be sorted in decreasing area
            iftQuickSort(area, index, 0, n-1, IFT_INCREASING);

            iftDestroySet(&component_ids_per_label[lb]);

            // Inserting the components by decreasing area size
            for (id = 0; id < n; id++) {
                iftInsertSet(&component_ids_per_label[lb], index[id]);
            }
            free(area);
            free(index);
        }
    }

    free(selected_components);
    
    iftDestroyHist(&hist);

    return component_ids_per_label;
}

int iftMarkersFromMisclassifiedComponentsSortedByAreaAndSeedEDT(iftImage *seed_image, iftImage *seed_image_mk_id,
                                                                iftLabeledSet **all_seeds, iftImage *gt,
                                                                iftImage *label, iftImage *dist, iftImage *root,
                                                                int nseeds_per_object, double marker_length,
                                                                double min_border_distance, iftAdjRel *A,
                                                                bool select_seeds_on_axes_only) {
    int i, p, id, lb, number_of_labels = iftMaximumValue(gt) + 1, number_of_components;
    int total_seeds = 0, seed_added = 0;
    iftVoxel center_marker_voxel;
    iftLabeledSet *S = NULL, *Sorted_seeds = NULL, *Stmp = NULL;
    iftImage *error_components = NULL;
    iftLabeledSet **seeds_per_component = NULL;
    iftSet **components_per_label = NULL, *aux = NULL;

    // Computing error components
    error_components = iftRelabelSegmentationErrorComponents(gt, label, A);

    // For each label, we sort by area its corresponding error components
    components_per_label = iftErrorComponentsPerSegmentationLabelSortedByArea(gt, error_components);

    // We then determine, for each error component the seeds that are part of them
    number_of_components = iftMaximumValue(error_components);
    seeds_per_component = (iftLabeledSet**)calloc(number_of_components, sizeof(iftLabeledSet*));

    for(id = 0; id < number_of_components; id++)
        seeds_per_component[id] = NULL;



    for(S = *all_seeds; S != NULL; S = S->next) {
        int p = S->elem, misclassified;

        lb = gt->val[p];

        misclassified = ((label == NULL) || (lb != label->val[p]));

        // Separating all seeds per error component
        if (misclassified) {
            id = error_components->val[p] - 1;
            iftInsertLabeledSet(&seeds_per_component[id], p, lb); // This reverses the sorting order by gradient value
        }
    }

    // Resorting the seeds in increasing order of gradient
    for(id = 0; id < number_of_components; id++) {
        Stmp = iftReverseLabeledSet(seeds_per_component[id]);
        iftDestroyLabeledSet(&seeds_per_component[id]);
        seeds_per_component[id] = Stmp;
    }


    for(lb = 0; lb < number_of_labels; lb++) {

        Sorted_seeds = NULL;
        // For a given label, we interleave the seeds from each of its error components
        // to select them in order of gradient magnitude AND error component area
        do {
            seed_added = 0;
            // Going through the error components sorted by area
            for(aux = components_per_label[lb]; aux != NULL; aux = aux->next) {
                id = aux->elem - 1;
                // Adding the first seed from the error component
                if(seeds_per_component[id] != NULL) {
                    int foo;
                    iftInsertLabeledSet(&Sorted_seeds, seeds_per_component[id]->elem, seeds_per_component[id]->label);

                    // Removing the added seed for the current error component
                    iftRemoveLabeledSet(&seeds_per_component[id], &foo);
                    seed_added = 1;
                }
            }
        }while(seed_added);

        // As before, the sorted seeds are in reverse order, so we reverse them again
        Stmp = iftReverseLabeledSet(Sorted_seeds);
        iftDestroyLabeledSet(&Sorted_seeds);
        Sorted_seeds = Stmp;

        // Finally, we draw the seeds on the seed image in order of gradient magnitude and error component area
        for(i = 0; Sorted_seeds != NULL && i < nseeds_per_object; i++) {
            p = Sorted_seeds->elem;

            center_marker_voxel = iftGetVoxelCoord(seed_image, p);

            seed_added = iftSelectSeedsBasedOnEDT(center_marker_voxel, gt, dist, root, error_components, seed_image,
                                                  seed_image_mk_id, marker_length, min_border_distance, total_seeds + 1,
                                                  select_seeds_on_axes_only);

            if(seed_added > 0) {

                // We remove the seeds within the marker radius from the added seed center from both the sorted
                // seed set and the seed set containing all available seeds
                iftRemoveSeedsWithinRadius(all_seeds, center_marker_voxel, lb, gt, marker_length);
                iftRemoveSeedsWithinRadius(&Sorted_seeds, center_marker_voxel, lb, gt, marker_length);

                total_seeds++;
            }
        }

        iftDestroyLabeledSet(&Sorted_seeds);
    }

    for(lb = 0; lb < number_of_labels; lb++)
        iftDestroySet(&components_per_label[lb]);

    for(id = 0; id < number_of_components; id++)
        iftDestroyLabeledSet(&seeds_per_component[id]);

    free(seeds_per_component);
    free(components_per_label);

    iftDestroyImage(&error_components);

    return total_seeds;
}

//
//ulong iftMarkersFromMisclassifiedComponentsSortedByAreaAndCircularSeed(iftImage *seed_image, iftLabeledSet **all_seeds,
//                                                                       iftImage *gt, iftImage *label,
//                                                                       int nseeds_per_object, double max_marker_radius,
//                                                                       double min_marker_radius,
//                                                                       double min_border_distance, iftAdjRel *A) {
//    int i, p, id, lb, number_of_labels = iftMaximumValue(gt) + 1, number_of_components;
//    ulong total_seeds = 0, seed_added = 0;
//    iftVoxel center_marker_voxel;
//    iftLabeledSet *S = NULL, *Sorted_seeds = NULL, *Stmp = NULL;
//    iftImage *error_components = NULL;
//    iftLabeledSet **seeds_per_component = NULL;
//    iftSet **components_per_label = NULL, *aux = NULL;
//    iftBMap *used_seeds = NULL;
//    iftAdjRel *distance_border = NULL;
//
//    used_seeds = iftCreateBMap(seed_image->n);
//    // Computing error components
//    error_components = iftRelabelSegmentationErrorComponents(gt, label, A);
//
//    // For each label, we sort by area its corresponding error components
//    components_per_label = iftErrorComponentsPerSegmentationLabelSortedByArea(gt, error_components);
//
//    // We then determine, for each error component the seeds that are part of them
//    number_of_components = iftMaximumValue(error_components);
//    seeds_per_component = (iftLabeledSet**)calloc(number_of_components, sizeof(iftLabeledSet*));
//
//    if(iftIs3DImage(gt))
//        distance_border = iftSpheric((float)(min_border_distance + max_marker_radius));
//    else
//        distance_border = iftCircular((float)(min_border_distance + max_marker_radius));
//
//    for(id = 0; id < number_of_components; id++)
//        seeds_per_component[id] = NULL;
//
//    for(S = *all_seeds; S != NULL; S = S->next) {
//        int p = S->elem, misclassified;
//
//        lb = gt->val[p];
//
//        misclassified = ((label == NULL) || (lb != label->val[p]));
//
//        // Separating all seeds per error component
//        if (misclassified) {
//            id = error_components->val[p] - 1;
//            iftInsertLabeledSet(&seeds_per_component[id], p, lb); // This reverses the sorting order by gradient value
//        }
//    }
//
//    // Resorting the seeds in increasing order of gradient
//    for(id = 0; id < number_of_components; id++) {
//        Stmp = iftReverseLabeledSet(seeds_per_component[id]);
//        iftDestroyLabeledSet(&seeds_per_component[id]);
//        seeds_per_component[id] = Stmp;
//    }
//
//
//    for(lb = 0; lb < number_of_labels; lb++) {
//        seed_added = 1;
//
//        Sorted_seeds = NULL;
//        // For a given label, we interleave the seeds from each of its error components
//        // to select them in order of gradient magnitude AND error component area
//        do {
//            seed_added = 0;
//            // Going through the error components sorted by area
//            for(aux = components_per_label[lb]; aux != NULL; aux = aux->next) {
//                id = aux->elem - 1;
//                // Adding the first see from the error component
//                if(seeds_per_component[id] != NULL) {
//                    int foo;
//                    iftInsertLabeledSet(&Sorted_seeds, seeds_per_component[id]->elem, seeds_per_component[id]->label);
//
//                    // Removing the added seed for the current error component
//                    iftRemoveLabeledSet(&seeds_per_component[id], &foo);
//                    seed_added = 1;
//                }
//            }
//        }while(seed_added);
//
//        // As before, the sorted seeds are in reverse order, so we reverse them again
//        Stmp = iftReverseLabeledSet(Sorted_seeds);
//        iftDestroyLabeledSet(&Sorted_seeds);
//        Sorted_seeds = Stmp;
//
//        // Finally, we draw the seeds on the seed image in order of gradient magnitude and error component area
//        for(i = 0; Sorted_seeds != NULL && i < nseeds_per_object; i++) {
//            p = Sorted_seeds->elem;
//
//            center_marker_voxel = iftGetVoxelCoord(seed_image, p);
//
//            seed_added = iftSelectCircularRobotSeeds(seed_image, used_seeds, gt, min_border_distance, max_marker_radius,
//                                                     min_marker_radius, distance_border, p);
//
//            if(seed_added > 0) {
//                total_seeds++;
//            }
//        }
//
//        iftDestroyLabeledSet(&Sorted_seeds);
//    }
//
//    for(lb = 0; lb < number_of_labels; lb++)
//        iftDestroySet(&components_per_label[lb]);
//
//    for(id = 0; id < number_of_components; id++)
//        iftDestroyLabeledSet(&seeds_per_component[id]);
//
//    free(seeds_per_component);
//    free(components_per_label);
//
//    iftDestroyImage(&error_components);
//    iftDestroyBMap(&used_seeds);
//    iftDestroyAdjRel(&distance_border);
//
//    return total_seeds;
//}

ulong iftMarkersFromMisclassifiedComponentsSortedByAreaAndCircleSeeds(iftImage *seed_image, iftLabeledSet **all_seeds,
                                                                      iftImage *gt, iftImage *label, iftImage *dist,
                                                                      iftImage *root, int nseeds_per_object,
                                                                      double marker_length, double min_border_distance,
                                                                      iftAdjRel *A, bool select_seeds_on_axes_only) {
    int i, p, id, lb, number_of_labels = iftMaximumValue(gt) + 1, number_of_components;
    ulong total_seeds = 0, seed_added = 0;
    iftVoxel center_marker_voxel;
    iftLabeledSet *S = NULL, *Sorted_seeds = NULL, *Stmp = NULL;
    iftImage *error_components = NULL;
    iftLabeledSet **seeds_per_component = NULL;
    iftSet **components_per_label = NULL, *aux = NULL;

    // Computing error components
    error_components = iftRelabelSegmentationErrorComponents(gt, label, A);

    // For each label, we sort by area its corresponding error components
    components_per_label = iftErrorComponentsPerSegmentationLabelSortedByArea(gt, error_components);

    // We then determine, for each error component the seeds that are part of them
    number_of_components = iftMaximumValue(error_components);
    seeds_per_component = (iftLabeledSet**)calloc(number_of_components, sizeof(iftLabeledSet*));

    for(id = 0; id < number_of_components; id++)
        seeds_per_component[id] = NULL;



    for(S = *all_seeds; S != NULL; S = S->next) {
        int p = S->elem, misclassified;

        lb = gt->val[p];

        misclassified = ((label == NULL) || (lb != label->val[p]));

        // Separating all seeds per error component
        if (misclassified) {
            id = error_components->val[p] - 1;
            iftInsertLabeledSet(&seeds_per_component[id], p, lb); // This reverses the sorting order by gradient value
        }
    }

    // Resorting the seeds in increasing order of gradient
    for(id = 0; id < number_of_components; id++) {
        Stmp = iftReverseLabeledSet(seeds_per_component[id]);
        iftDestroyLabeledSet(&seeds_per_component[id]);
        seeds_per_component[id] = Stmp;
    }


    for(lb = 0; lb < number_of_labels; lb++) {
        seed_added = 1;

        Sorted_seeds = NULL;
        // For a given label, we interleave the seeds from each of its error components
        // to select them in order of gradient magnitude AND error component area
        do {
            seed_added = 0;
            // Going through the error components sorted by area
            for(aux = components_per_label[lb]; aux != NULL; aux = aux->next) {
                id = aux->elem - 1;
                // Adding the first see from the error component
                if(seeds_per_component[id] != NULL) {
                    int foo;
                    iftInsertLabeledSet(&Sorted_seeds, seeds_per_component[id]->elem, seeds_per_component[id]->label);

                    // Removing the added seed for the current error component
                    iftRemoveLabeledSet(&seeds_per_component[id], &foo);
                    seed_added = 1;
                }
            }
        }while(seed_added);

        // As before, the sorted seeds are in reverse order, so we reverse them again
        Stmp = iftReverseLabeledSet(Sorted_seeds);
        iftDestroyLabeledSet(&Sorted_seeds);
        Sorted_seeds = Stmp;

        // Finally, we draw the seeds on the seed image in order of gradient magnitude and error component area
        for(i = 0; Sorted_seeds != NULL && i < nseeds_per_object; i++) {
            p = Sorted_seeds->elem;

            center_marker_voxel = iftGetVoxelCoord(seed_image, p);

            seed_added = iftSelectSeedsBasedOnEDT(center_marker_voxel, gt, dist, root, error_components, seed_image,
                                                  NULL, marker_length, min_border_distance, 0,
                                                  select_seeds_on_axes_only);

            if(seed_added > 0) {

                // We remove the seeds within the marker radius from the added seed center from both the sorted
                // seed set and the seed set containing all available seeds
                iftRemoveSeedsWithinRadius(all_seeds, center_marker_voxel, lb, gt, marker_length);
                iftRemoveSeedsWithinRadius(&Sorted_seeds, center_marker_voxel, lb, gt, marker_length);

                total_seeds++;
            }
        }

        iftDestroyLabeledSet(&Sorted_seeds);
    }

    for(lb = 0; lb < number_of_labels; lb++)
        iftDestroySet(&components_per_label[lb]);

    for(id = 0; id < number_of_components; id++)
        iftDestroyLabeledSet(&seeds_per_component[id]);

    free(seeds_per_component);
    free(components_per_label);

    iftDestroyImage(&error_components);

    return total_seeds;
}


void iftDrawSeeds(iftImage *img, iftImage *seed_image, iftColorTable *cmap)
{
    int p;

    iftVerifyImageDomains(img, seed_image,"iftDrawSeeds");

    if (img->Cb==NULL)
        iftSetCbCr(img,128);

    for (p=0; p < img->n; p++) {
        if(seed_image->val[p] >= 0) {
            iftSetRGB(img, p, cmap->color[seed_image->val[p]].val[0],
                              cmap->color[seed_image->val[p]].val[1],
                              cmap->color[seed_image->val[p]].val[2], 255);
        }
    }
}

iftImage* iftSeedImageMkIdFromLabeledSet(iftLabeledSet* seeds, iftImage *image){
    iftImage* seed_image_mk_id = iftCreateImage(image->xsize,image->ysize,image->zsize);
    iftSetImage(seed_image_mk_id, -1);

    iftLabeledSet* i = seeds;

    while(i != NULL){
        image->val[i->elem] = i->marker;
        i = i->next;
    }
    return seed_image_mk_id;
}

iftLabeledSet *iftLabelToForestPixelRobotEDT(iftImage *gradient, iftImage *label, iftAdjRel *A, int nseeds_per_label_per_iteration,
                                             double min_safe_distance_to_border, double max_marker_width, double max_marker_length,
                                             iftLabeledSet *optional_input_seeds, double stopping_threshold,
                                             int secondary_stopping_threshold,
                                             uchar (*iftStoppingCriterion)(iftImage *, iftImage *, double),
                                             uchar (*iftSecondaryStoppingCriterion)(int, int, int), bool select_seeds_on_axes_only)
{
    int j = 0, seeds_added;
    int total_added_seeds = 0;
    iftLabeledSet  *available_seeds = NULL, *current_seeds = NULL;
    iftImage *seeds_image         = NULL;
    iftImage *seeds_image_mk_id    = NULL;
    iftImage *seeds_image_copy     = NULL;
    iftImage *seeds_image_mk_id_copy     = NULL;
    iftImage *current_segmentation = NULL;
    iftImageForest *forest         = NULL;
    iftImage *gt_image             = NULL;
    iftImage *dist = NULL, *root = NULL;
    iftLabeledSet *all_segmentation_seeds = NULL;
    timer *tic = NULL, *toc = NULL;

    iftMaximumValue(label);

    gt_image = label;

    forest = iftCreateImageForest(gradient, A);

    // Computing the seed set from the ground truth image (original erroneously segmented image) that is close to the
    // objects' boundaries and sorted by gradient value. Note that we add the minimum safe distance to the boundary
    // and the maximum marker width to ensure that the marker will be as wide as specified
    available_seeds  = iftBorderMarkersDistAndRootForPixelSegmentation(gradient, gt_image,
                                                                       min_safe_distance_to_border + max_marker_width,
                                                                       &dist, &root);

    // Performing an initial segmentation if the optinal input seed set is passed as a parameter
    if(optional_input_seeds != NULL) {
        iftDiffWatershed(forest, optional_input_seeds, NULL);
        current_segmentation = forest->label;
        seeds_image = iftSeedImageFromLabeledSet(optional_input_seeds, gt_image);
        seeds_image_mk_id = iftSeedImageMkIdFromLabeledSet(optional_input_seeds, gt_image);
        all_segmentation_seeds = iftCopyOrderedLabeledSet(optional_input_seeds);

    } else {
        seeds_image = iftCreateImage(gradient->xsize, gradient->ysize, gradient->zsize);
        seeds_image_mk_id = iftCreateImage(gradient->xsize, gradient->ysize, gradient->zsize);
        iftSetImage(seeds_image, IFT_NIL);
        iftSetImage(seeds_image_mk_id, IFT_NIL);
    }

    do
    {
        printf("Rebuilding iteration: %03d\n", ++j);

        tic = iftTic();
        seeds_image_copy = iftCopyImage(seeds_image);
        seeds_image_mk_id_copy = iftCopyImage(seeds_image_mk_id);
        // Selecting new seeds according to the available ones for each object in order of gradient and
        // error component area

        seeds_added = iftMarkersFromMisclassifiedComponentsSortedByAreaAndSeedEDT(seeds_image, seeds_image_mk_id, &available_seeds,
                                                                                  gt_image, current_segmentation, dist,
                                                                                  root, nseeds_per_label_per_iteration,
                                                                                  max_marker_length,
                                                                                  min_safe_distance_to_border, A,
                                                                                  select_seeds_on_axes_only);

        fprintf(stderr,"Seeds added %d\n", seeds_added);
        iftCopyNewSeeds(seeds_image, seeds_image_mk_id, seeds_image_copy, seeds_image_mk_id_copy);

        iftGetNewlyAddedSeeds(total_added_seeds, seeds_image_copy, seeds_image_mk_id_copy, &current_seeds, &all_segmentation_seeds);

        // Segmenting the image by adding the new seed set
        iftDiffWatershed(forest, current_seeds, NULL);

        current_segmentation = forest->label;

        toc = iftToc();

        iftPrintFormattedTime(iftCompTime(tic, toc));

        iftDestroyImage(&seeds_image_copy);
        iftDestroyImage(&seeds_image_mk_id_copy);
        iftDestroyLabeledSet(&current_seeds);

        total_added_seeds += seeds_added;
    }while(seeds_added != 0 && !iftStoppingCriterion(current_segmentation, gt_image, stopping_threshold)
           && !iftSecondaryStoppingCriterion(j, total_added_seeds, secondary_stopping_threshold));

    iftDestroyImage(&seeds_image);
    iftDestroyImage(&seeds_image_mk_id);
    iftDestroyImage(&dist);
    iftDestroyImage(&root);
    iftDestroyLabeledSet(&available_seeds);
    iftDestroyImageForest(&forest);


    return all_segmentation_seeds;
}

void iftCopyNewSeeds(iftImage *seeds_image, iftImage *seeds_image_mk_id, iftImage *seeds_image_copy,
                     iftImage *seeds_image_mk_id_copy) {//This produces only the new seeds added this iteration
    for (int p = 0; p < seeds_image_copy->n; p++)
    {
        if (seeds_image_copy->val[p] == seeds_image->val[p]) {
            seeds_image_copy->val[p] = IFT_NIL;
            seeds_image_mk_id_copy->val[p] = IFT_NIL;
        } else {
            seeds_image_copy->val[p] = seeds_image->val[p];
            seeds_image_mk_id_copy->val[p] = seeds_image_mk_id->val[p];
        }
    }
}

iftLabeledSet *iftLabelToForestPixelRobotEDTOrientedWatershed(iftImageForest *fst, iftImage *gt_image,
                                                              int nseeds_per_label_per_iteration,
                                                              double min_safe_distance_to_border,
                                                              double max_marker_width, double max_marker_length,
                                                              double gamma, iftLabeledSet *optional_input_seeds,
                                                              double stopping_threshold,
                                                              int secondary_stopping_threshold,
                                                              uchar (*iftStoppingCriterion)(iftImage *, iftImage *,
                                                                                            double),
                                                              uchar (*iftSecondaryStoppingCriterion)(int, int, int),
                                                              bool select_seeds_on_axes_only)
{
    int j = 0, seeds_added;
    int total_added_seeds = 0;
    iftLabeledSet  *available_seeds = NULL, *all_delineation_seeds = NULL, *current_seeds = NULL;
    iftImage  *seeds_image         = NULL;
    iftImage *seeds_image_mk_id    = NULL;
    iftImage *seeds_image_copy     = NULL;
    iftImage *seeds_image_mk_id_copy     = NULL;
    iftImage *current_segmentation = NULL;

    iftImage *dist = NULL, *root = NULL;
    iftLabeledSet *all_segmentation_seeds = NULL;
    timer *tic = NULL, *toc = NULL;

    iftMaximumValue(gt_image);

    // Computing the seed set from the ground truth image (original erroneously segmented image) that is close to the
    // objects' boundaries and sorted by fst->img value. Note that we add the minimum safe distance to the boundary
    // and the maximum marker width to ensure that the marker will be as wide as specified
    available_seeds  = iftBorderMarkersDistAndRootForPixelSegmentation(fst->img, gt_image,
                                                                       min_safe_distance_to_border + max_marker_width,
                                                                       &dist, &root);

    // Performing an initial segmentation if the optinal input seed set is passed as a parameter
    if(optional_input_seeds != NULL) {
        iftDiffOrientedWatershedResuming(fst, gt_image, optional_input_seeds, NULL, gamma);
        current_segmentation = fst->label;
        seeds_image = iftSeedImageFromLabeledSet(optional_input_seeds, gt_image);
        seeds_image_mk_id = iftSeedImageMkIdFromLabeledSet(optional_input_seeds, gt_image);
        all_segmentation_seeds = iftCopyOrderedLabeledSet(optional_input_seeds);

    } else {
        seeds_image = iftCreateImage(fst->img->xsize, fst->img->ysize, fst->img->zsize);
        seeds_image_mk_id = iftCreateImage(fst->img->xsize, fst->img->ysize, fst->img->zsize);
        iftSetImage(seeds_image, IFT_NIL);
        iftSetImage(seeds_image_mk_id, IFT_NIL);
    }

    iftColorTable *cmap = iftCreateColorTable(iftMaximumValue(gt_image)+1);

    do
    {
        printf("Rebuilding iteration: %03d\n", ++j);

        tic = iftTic();
        seeds_image_copy = iftCopyImage(seeds_image);
        seeds_image_mk_id_copy = iftCopyImage(seeds_image_mk_id);

        // Selecting new seeds according to the available ones for each object in order of fst->img and
        // error component area
        seeds_added = iftMarkersFromMisclassifiedComponentsSortedByAreaAndSeedEDT(seeds_image, seeds_image_mk_id, &available_seeds,
                                                                                  gt_image, current_segmentation, dist,
                                                                                  root, nseeds_per_label_per_iteration,
                                                                                  max_marker_length,
                                                                                  min_safe_distance_to_border, fst->A,
                                                                                  select_seeds_on_axes_only);

        fprintf(stderr,"Seeds added %d\n", seeds_added);
        iftCopyNewSeeds(seeds_image, seeds_image_mk_id, seeds_image_copy, seeds_image_mk_id_copy);

        iftGetNewlyAddedSeeds(total_added_seeds, seeds_image_copy, seeds_image_mk_id_copy, &current_seeds, &all_segmentation_seeds);

        iftConcatLabeledSet(&all_delineation_seeds, &current_seeds);

        iftDiffOrientedWatershedResuming(fst, gt_image, current_seeds, NULL, gamma);

        current_segmentation = fst->label;

        iftImage *tmp_img = iftLinearStretch(current_segmentation, iftMinimumValue(current_segmentation),
                                             iftMaximumValue(current_segmentation), 0, 255);
        iftWriteImageByExt(tmp_img, "seg_%04d.pgm", j);
        iftImage *seed_img = iftSeedImageFromLabeledSet(all_delineation_seeds, current_segmentation);
        iftDrawSeeds(tmp_img, seed_img, cmap);
        iftWriteImageByExt(tmp_img, "seeds_%04d.ppm", j);
        iftDestroyImage(&tmp_img);
        iftDestroyImage(&seed_img);

        toc = iftToc();

        iftPrintFormattedTime(iftCompTime(tic, toc));

        iftDestroyImage(&seeds_image_copy);
        iftDestroyImage(&seeds_image_mk_id_copy);
        iftDestroyLabeledSet(&current_seeds);
        total_added_seeds += seeds_added;

    }while(seeds_added != 0 && !iftStoppingCriterion(current_segmentation, gt_image, stopping_threshold)
           && !iftSecondaryStoppingCriterion(j, total_added_seeds, secondary_stopping_threshold));

    iftDestroyLabeledSet(&all_delineation_seeds);
    iftDestroyImage(&seeds_image);
    iftDestroyImage(&seeds_image_mk_id);
    iftDestroyImage(&dist);
    iftDestroyImage(&root);
    iftDestroyLabeledSet(&available_seeds);
    iftDestroyColorTable(&cmap);

    return all_segmentation_seeds;
}


iftLabeledSet *iftLabelToForestPixelRobotEDTIGraphOrientedWatershed(iftIGraph *igraph, iftImage *gt_image,
                                                                    int nseeds_per_label_per_iteration,
                                                                    double min_safe_distance_to_border,
                                                                    double max_marker_width, double max_marker_length,
                                                                    double gamma, iftLabeledSet *optional_input_seeds,
                                                                    double stopping_threshold,
                                                                    int secondary_stopping_threshold,
                                                                    uchar (*iftStoppingCriterion)(iftImage *,
                                                                                                  iftImage *, double),
                                                                    uchar (*iftSecondaryStoppingCriterion)(int, int,
                                                                                                           int),
                                                                    bool select_seeds_on_axes_only)
{
    int j = 0, seeds_added;
    int total_added_seeds = 0;
    iftLabeledSet  *available_seeds = NULL, *current_seeds = NULL, *all_delineation_seeds = NULL;
    iftImage  *seeds_image         = NULL;
    iftImage *seeds_image_mk_id    = NULL;
    iftImage *seeds_image_copy     = NULL;
    iftImage *seeds_image_mk_id_copy     = NULL;
    iftImage *current_segmentation = NULL;

    iftImage *dist = NULL, *root = NULL;
    iftImage *gradient           = NULL;
    iftFImage *weight            = NULL;
    iftLabeledSet *all_segmentation_seeds = NULL;
    iftDHeap *Q = NULL;
    double *pvalue = NULL;

    timer *tic = NULL, *toc = NULL;

    iftMaximumValue(gt_image);

    weight      = iftIGraphWeight(igraph);
    gradient    = iftFImageToImage(weight, IFT_MAXWEIGHT);
    pvalue      = iftAllocDoubleArray(igraph->nnodes);
    Q           = iftIGraphInitDiffWatershed(igraph, pvalue);

    // Computing the seed set from the ground truth image (original erroneously segmented image) that is close to the
    // objects' boundaries and sorted by gradient value. Note that we add the minimum safe distance to the boundary
    // and the maximum marker width to ensure that the marker will be as wide as specified
    available_seeds  = iftBorderMarkersDistAndRootForPixelSegmentation(gradient, gt_image,
                                                                       min_safe_distance_to_border + max_marker_width,
                                                                       &dist, &root);

    // Performing an initial segmentation if the optinal input seed set is passed as a parameter
    if(optional_input_seeds != NULL) {
        iftIGraphResumingDiffOrientedWatershed(igraph, gt_image, optional_input_seeds, NULL, Q, gamma);
        current_segmentation = iftIGraphLabel(igraph);

        seeds_image = iftSeedImageFromLabeledSet(optional_input_seeds, gt_image);
        all_segmentation_seeds = iftCopyOrderedLabeledSet(optional_input_seeds);

    } else {
        seeds_image = iftCreateImage(gradient->xsize, gradient->ysize, gradient->zsize);
        seeds_image_mk_id = iftCreateImage(gradient->xsize, gradient->ysize, gradient->zsize);
        iftSetImage(seeds_image, IFT_NIL);
    }

    iftColorTable *cmap = iftCreateColorTable(iftMaximumValue(gt_image)+1);

    iftImageForest *fst = iftCreateImageForest(gradient, igraph->A);
    iftDestroyGQueue(&fst->Q);
    fst->Q = iftCreateGQueue(iftRound(gamma*iftMaximumValue(gradient))+1, gradient->n, fst->pathval->val);

    do
    {
        printf("Rebuilding iteration: %03d\n", ++j);

        tic = iftTic();
        seeds_image_copy = iftCopyImage(seeds_image);
        seeds_image_mk_id_copy = iftCopyImage(seeds_image_mk_id);

        // Selecting new seeds according to the available ones for each object in order of gradient and
        // error component area
        seeds_added = iftMarkersFromMisclassifiedComponentsSortedByAreaAndSeedEDT(seeds_image, seeds_image_mk_id, &available_seeds,
                                                                                  gt_image, current_segmentation, dist,
                                                                                  root, nseeds_per_label_per_iteration,
                                                                                  max_marker_length,
                                                                                  min_safe_distance_to_border,
                                                                                  igraph->A,
                                                                                  select_seeds_on_axes_only);

        fprintf(stderr,"Seeds added %d\n", seeds_added);
        iftCopyNewSeeds(seeds_image, seeds_image_mk_id, seeds_image_copy, seeds_image_mk_id_copy);

        iftGetNewlyAddedSeeds(total_added_seeds, seeds_image_copy, seeds_image_mk_id_copy, &current_seeds, &all_segmentation_seeds);

        iftDestroyImage(&current_segmentation);

        iftConcatLabeledSet(&all_delineation_seeds, &current_seeds);

        /* Starting from scratch while the oriented DIFT version does not work */
        iftIGraphResetWatershed(igraph, Q);
        iftIGraphDiffWatershed(igraph, all_delineation_seeds, NULL, Q);

        iftIGraphResumingDiffOrientedWatershed(igraph, gt_image, all_delineation_seeds, NULL, Q, gamma);

//        iftIGraphResumingDiffOrientedWatershed(igraph, gt_image, current_seeds, NULL, Q, gamma);

        current_segmentation = iftIGraphLabel(igraph);

//          current_segmentation = iftOrientedWatershedResuming(gt_image, gradient, igraph->A, all_delineation_seeds, NULL, gamma);

        iftImage *tmp_img = iftLinearStretch(current_segmentation, iftMinimumValue(current_segmentation),
                                             iftMaximumValue(current_segmentation), 0, 255);
        iftWriteImageByExt(tmp_img, "seg_%04d.pgm", j);
        iftImage *seed_img = iftSeedImageFromLabeledSet(all_delineation_seeds, current_segmentation);
        iftDrawSeeds(tmp_img, seed_img, cmap);
        iftWriteImageByExt(tmp_img, "seeds_%04d.ppm", j);
        iftDestroyImage(&tmp_img);
        iftDestroyImage(&seed_img);

        toc = iftToc();

        iftPrintFormattedTime(iftCompTime(tic, toc));

        iftDestroyImage(&seeds_image_copy);
        iftDestroyImage(&seeds_image_mk_id_copy);
        iftDestroyLabeledSet(&current_seeds);

        total_added_seeds += seeds_added;

    }while(seeds_added != 0 && !iftStoppingCriterion(current_segmentation, gt_image, stopping_threshold)
           && !iftSecondaryStoppingCriterion(j, total_added_seeds, secondary_stopping_threshold));

    iftFree(pvalue);
    iftDestroyDHeap(&Q);
    iftDestroyImage(&gradient);
    iftDestroyFImage(&weight);
    iftDestroyImage(&seeds_image);
    iftDestroyImage(&seeds_image_mk_id);
    iftDestroyImage(&current_segmentation);
    iftDestroyImage(&dist);
    iftDestroyImage(&root);
    iftDestroyLabeledSet(&available_seeds);
    iftDestroyLabeledSet(&all_delineation_seeds);
    iftDestroyColorTable(&cmap);

    return all_segmentation_seeds;
}


//iftLabeledSet *iftLabelToForestGeodesicRobotByACC(iftImage *gradient, iftImage *label, iftAdjRel *A,
//                                                  int seeds_per_iteration, int min_distance_border, int max_marker_size,
//                                                  int min_marker_size, iftLabeledSet *optional_input_seeds,
//                                                  double stopping_threshold, int secondary_stopping_threshold,
//                                                  uchar (*iftStoppingCriterion)(iftImage *, iftImage *, double),
//                                                  uchar (*iftSecondaryStoppingCriterion)(int, int, int))
//{
//    int p, j = 0, seeds_added, number_of_seeds = 0;
//    iftLabeledSet  *available_seeds = NULL, *current_seeds = NULL;
//    iftBMap *seeds_bmap            = NULL;
//    iftImage *seeds_image         = NULL;
//    iftImage *seeds_image_copy     = NULL;
//    iftImage *current_segmentation = NULL;
//    iftImageForest *forest         = NULL;
//    iftImage *gt = NULL;
//    iftLabeledSet *all_segmentation_seeds = NULL;
//    int number_of_objects, number_of_labels;
//    int total_added_seeds = 0;
//
//    gt = label;
//
//    number_of_objects = iftMaximumValue(gt);
//    number_of_labels = number_of_objects + 1;
//    forest = iftCreateImageForest(gradient, A);
//    seeds_bmap  = iftCreateBMap(gradient->n);
//
//    // Performing an initial segmentation if the optinal input seed set is passed as a parameter
//    if(optional_input_seeds != NULL) {
//        iftDiffWatershed(forest, optional_input_seeds, NULL);
//        current_segmentation = forest->label;
//        seeds_image = iftSeedImageFromLabeledSet(optional_input_seeds, gt);
//        all_segmentation_seeds = iftCopyOrderedLabeledSet(optional_input_seeds);
//    } else {
//        seeds_image = iftCreateImage(gradient->xsize, gradient->ysize, gradient->zsize);
//        iftSetImage(seeds_image, IFT_NIL);
//    }
//
//    do
//    {
//        printf("Rebuilding iteration: %03d\n", ++j);
//        available_seeds  = iftGeodesicMarkersForSegmentation(gt, current_segmentation);
//        seeds_image_copy = iftCopyImage(seeds_image);
//        seeds_added      = iftMarkersFromMisclassifiedSeeds(seeds_image, available_seeds, seeds_bmap,
//                                                            seeds_per_iteration * (number_of_objects + 1),
//                                                            number_of_objects + 1, gt, current_segmentation,
//                                                            min_distance_border, max_marker_size, min_marker_size);
//
//        total_added_seeds += seeds_added;
//
//        fprintf(stderr,"Seeds added %d\n", seeds_added);
//
//        //This produces only the new seeds added this iteration
//        for (p = 0; p < seeds_image_copy->n; p++)
//        {
//            if (seeds_image_copy->val[p] == seeds_image->val[p])
//                seeds_image_copy->val[p] = IFT_NIL;
//            else
//                seeds_image_copy->val[p] = seeds_image->val[p];
//        }
//
//        iftGetNewlyAddedSeeds(j, number_of_labels, seeds_image_copy, NULL, &current_seeds, &all_segmentation_seeds);
//
//        iftDiffWatershed(forest, current_seeds, NULL);
//        current_segmentation = forest->label;
//
//        number_of_seeds += seeds_added;
//
////        iftWriteImageByExt(seeds_image, "seeds_%03d.scn", j);
////        iftWriteImageByExt(current_segmentation, "label_%03d.scn", j);
//
//        iftDestroyImage(&seeds_image_copy);
//        iftDestroyLabeledSet(&available_seeds);
//        iftDestroyLabeledSet(&current_seeds);
//    }while(seeds_added != 0 && !iftStoppingCriterion(gt, current_segmentation, stopping_threshold)
//           && !iftSecondaryStoppingCriterion(j, total_added_seeds, secondary_stopping_threshold));
//
//    iftDestroyBMap(&seeds_bmap);
//    iftDestroyImage(&seeds_image);
//
//    iftDestroyImageForest(&forest);
//
//    return all_segmentation_seeds;
//}

iftImage *iftRespSystemGradient(iftImage *img, iftAdjRel *A)
{
    iftImage *aux, *grad;

    /* Enhance dark regions of the respiratory system: good to assign
       arc weights for interactive lung and bronchi segmentation. */

    aux = iftCreateImage(img->xsize,img->ysize,img->zsize);
    for (int z=0; z < img->zsize; z++) {
        iftImage *slice = iftGetXYSlice(img,z);
        iftImage *cbas  = iftCloseBasins(slice, NULL, NULL);
        iftImage *res   = iftSub(cbas,slice);
        iftPutXYSlice(aux,res,z);
        iftDestroyImage(&res);
        iftDestroyImage(&cbas);
        iftDestroyImage(&slice);
    }

    grad = iftImageBasins(aux, A);

    iftDestroyImage(&aux);

    return grad;
}



void iftResumeImageSegmentation(iftImage *orig, iftImage *input_label, iftImage *basins, iftAdjRel *A, iftDict *json,
                                iftLabeledSet *optional_input_seeds, const char *out_dir, char *out_img_basename) {
    iftImage *reconstructed_label = NULL;
    iftImage *seed_image = NULL;
    iftLabeledSet *seeds = NULL;
    iftIGraph *igraph = NULL;
    iftImageForest *fst = NULL;

    char *out_seeds_img_file;
    char *out_reconstructed_label_file;
    char *out_basins_file;
    char *out_seeds_file_txt;
    char out_ext[10];
    char *tmp = NULL;

    int nseeds_per_label_per_iteration;
    double min_safe_distance_to_border;
    double max_marker_width;
    double max_marker_length;
    double stopping_threshold;
//    double min_marker_radius;
//    double max_marker_radius;
//    double border_dilation_radius;
    char *resuming_type = NULL;
    int measure_reconstruction_accuracy;
    int secondary_stopping_threshold;
    char *secondary_stopping_criterion = NULL;
    char *stopping_criterion = NULL;
    bool select_EDT_seeds_on_axes_only;
    double ori_watershed_gamma;

    uchar (*iftStoppingCriterion)(iftImage*, iftImage*, double) = NULL;
    uchar (*iftSecondaryStoppingCriterion)(int, int, int) = NULL;

    /* ISF-based resuming parameters */
    double alpha;
    double beta;
    double ISF_gamma;
//    int DinIFT_pred_length;
    int niters;
    int nsuperpixels_final;
    /* DynamicIFT */
    iftMImage *mimg = NULL;


    /* Reading configuration data */
    nseeds_per_label_per_iteration = iftGetLongValFromDict("nseeds_per_label_per_iteration", json);
    min_safe_distance_to_border = iftGetDblValFromDict("min_safe_distance_to_border", json);
    max_marker_width = iftGetDblValFromDict("max_marker_width", json);
    max_marker_length = iftGetDblValFromDict("max_marker_length", json);

    stopping_criterion = iftGetStrValFromDict("stopping_criterion", json);
    stopping_threshold = iftGetDblValFromDict("stopping_threshold", json);

    secondary_stopping_criterion = iftGetStrValFromDict("secondary_stopping_criterion", json);
    secondary_stopping_threshold = iftGetLongValFromDict("secondary_stopping_threshold", json);

//    min_marker_radius       = iftGetJDouble(json, "min_marker_radius");
//    max_marker_radius       = iftGetJDouble(json, "max_marker_radius");
//    border_dilation_radius  = iftGetJDouble(json, "border_dilation_radius");
    ori_watershed_gamma     = iftGetDblValFromDict("oriented_watershed_gamma", json);

    resuming_type           = iftGetStrValFromDict("resuming_type", json);


    measure_reconstruction_accuracy = iftGetLongValFromDict("measure_reconstruction_accuracy", json);
    select_EDT_seeds_on_axes_only   = iftGetBoolValFromDict("select_EDT_seeds_on_axes_only", json);

    alpha = iftGetDblValFromDict("ISF_alpha", json);
    beta = iftGetDblValFromDict("ISF_beta", json);
    ISF_gamma = iftGetDblValFromDict("ISF_gamma", json);
    niters = iftGetLongValFromDict("ISF_niters", json);

//    DinIFT_pred_length = iftGetLongValFromDict( "DynamicIFT_pred_length", json);


    /* Saving results */
    if(iftIs3DImage(orig))
        sprintf(out_ext, ".scn");
    else
        sprintf(out_ext, ".pgm");

    if(iftCompareStrings(stopping_criterion, "STOP_BY_DICE"))
        iftStoppingCriterion = iftStopLabelReconstructionByDice;
    else if(iftCompareStrings(stopping_criterion, "STOP_BY_ASSD"))
        iftStoppingCriterion = iftStopLabelReconstructionByASSD;
    else
        iftError(
                "Please choose between STOP_BY__NUM_ITERATIONS or STOP_BY_NUM_SEEDS for the primary stopping criterion",
                "main");


    if(iftCompareStrings(secondary_stopping_criterion, "STOP_BY_NUM_ITERATIONS"))
        iftSecondaryStoppingCriterion = iftStopLabelReconstructionByNumInterations;
    else if(iftCompareStrings(secondary_stopping_criterion, "STOP_BY_NUM_SEEDS"))
        iftSecondaryStoppingCriterion = iftStopLabelReconstructionByNumSeeds;
    else
        iftError(
                "Please choose between STOP_BY_NUM_ITERATIONS or STOP_BY_NUM_SEEDS for the secondary stopping criterion",
                "main");

    fprintf(stderr,"resuming_type %s\n", resuming_type);
    /* Resuming according to selected type */
    if(iftCompareStrings(resuming_type, "EDT_MARKER_ORIENTED")) {
        fst = iftCreateImageForest(basins, A);

        seeds = iftLabelToForestPixelRobotEDTOrientedWatershed(fst, input_label, nseeds_per_label_per_iteration,
                                                               min_safe_distance_to_border, max_marker_width,
                                                               max_marker_length, ori_watershed_gamma,
                                                               optional_input_seeds, stopping_threshold,
                                                               secondary_stopping_threshold, iftStoppingCriterion,
                                                               iftSecondaryStoppingCriterion,
                                                               select_EDT_seeds_on_axes_only);

    } else if(iftCompareStrings(resuming_type, "EDT_MARKER_ORIENTED_IGRAPH")) {
        iftImage  *mask1;

        if (iftIsColorImage(orig)){
            mimg   = iftImageToMImage(orig,LABNorm_CSPACE);
        } else {
            mimg   = iftImageToMImage(orig,GRAY_CSPACE);
        }

        mask1  = iftSelectImageDomain(mimg->xsize,mimg->ysize,mimg->zsize);

        /* minima of a basins manifold in that domain */
        igraph = iftImplicitIGraph(mimg,mask1,A);

        iftIGraphSetWeight(igraph, basins);

        fprintf(stderr,"Here\n");

        seeds = iftLabelToForestPixelRobotEDTIGraphOrientedWatershed(igraph, input_label,
                                                                     nseeds_per_label_per_iteration,
                                                                     min_safe_distance_to_border, max_marker_width,
                                                                     max_marker_length, ori_watershed_gamma,
                                                                     optional_input_seeds, stopping_threshold,
                                                                     secondary_stopping_threshold, iftStoppingCriterion,
                                                                     iftSecondaryStoppingCriterion,
                                                                     select_EDT_seeds_on_axes_only);

        iftDestroyImage(&mask1);
    } else if(iftCompareStrings(resuming_type, "EDT_MARKER")) {
        seeds = iftLabelToForestPixelRobotEDT(basins, input_label, A, nseeds_per_label_per_iteration,
                                              min_safe_distance_to_border, max_marker_width, max_marker_length,
                                              optional_input_seeds, stopping_threshold, secondary_stopping_threshold,
                                              iftStoppingCriterion, iftSecondaryStoppingCriterion, select_EDT_seeds_on_axes_only);

//    } else if(iftCompareStrings(resuming_type, "CIRCLE_PIXEL_MARKER")) {
//        seeds = iftLabelToForestPixelRobot(basins, input_label, A, nseeds_per_label_per_iteration,
//                                           border_dilation_radius, min_safe_distance_to_border, min_marker_radius,
//                                           max_marker_radius, optional_input_seeds, stopping_threshold, secondary_stopping_threshold,
//                                           iftStoppingCriterion, iftSecondaryStoppingCriterion);
//    } else if(iftCompareStrings(resuming_type, "GEODESIC_MARKER")) {
//        seeds = iftLabelToForestGeodesicRobotByACC(basins, input_label, A, nseeds_per_label_per_iteration,
//                                                   min_safe_distance_to_border, max_marker_radius, min_marker_radius,
//                                                   optional_input_seeds, stopping_threshold,
//                                                   secondary_stopping_threshold,
//                                                   iftStoppingCriterion, iftSecondaryStoppingCriterion);

    } else if(iftCompareStrings(resuming_type, "ISF")) {
        seeds = iftLabelToForestISF_Root(orig, input_label, alpha, beta, ISF_gamma, niters, &nsuperpixels_final,
                                         &igraph, NULL);
        reconstructed_label = iftIGraphLabel(igraph);

        iftDestroyIGraph(&igraph);
//    } else if(iftCompareStrings(resuming_type, "DYNAMIC_IFT")) {
//        if(iftIsColorImage(orig)){
//            mimg = iftImageToMImage(orig, LAB_CSPACE);
//        } else {
//            mimg = iftImageToMImage(orig, GRAY_CSPACE);
//        }
//
//        seeds = iftSeedsFromDynamicIFTRobot(mimg, basins, input_label, A, nseeds_per_label_per_iteration,
//                                            min_safe_distance_to_border, max_marker_width, max_marker_length,
//                                            optional_input_seeds, stopping_threshold,
//                                            secondary_stopping_threshold, alpha, DinIFT_pred_length, iftStoppingCriterion,
//                                            iftSecondaryStoppingCriterion, select_EDT_seeds_on_axes_only);

    } else {
        iftError(
                "Please choose among EDT_MARKER, EDT_MARKER_ORIENTED, EDT_MARKER_ORIENTED_IGRAPH, ISF and DYNAMIC_IFT for the resuming type",
                "main");
    }


    /* Performing one last segmentation to reconstruct the input label */
    seed_image = iftSeedImageFromLabeledSet(seeds, orig);
    int nlabels = iftMaximumValue(seed_image) + 1;

    for(int p = 0; p < seed_image->n; p++)
        seed_image->val[p] = (seed_image->val[p] < 0) ? 0 : iftRound((seed_image->val[p]+1)*255.0/nlabels);


    if(!iftCompareStrings(resuming_type, "ISF")) {
        if (iftCompareStrings(resuming_type, "EDT_MARKER_ORIENTED")) {
            reconstructed_label = iftCopyImage(fst->label);
        } else if(iftCompareStrings(resuming_type, "EDT_MARKER_ORIENTED_IGRAPH")) {
            reconstructed_label = iftIGraphLabel(igraph);
//        } else if(iftCompareStrings(resuming_type, "DYNAMIC_IFT")){
//            reconstructed_label = iftDynamicObjectDelineation(mimg,seeds, DinIFT_pred_length, NULL);
//            iftWriteSeeds2D("tmp/teste.txt", seeds, reconstructed_label);
        } else {
            reconstructed_label = iftWatershed(basins, A, seeds, NULL);
        }
    }

    /* Measuring accuracy one last time */
    if (measure_reconstruction_accuracy) {
        fprintf(stderr, "\nFinal reconstruction accuracy\n");

        iftStopLabelReconstructionByDice(reconstructed_label, input_label, stopping_threshold);
//        iftStopLabelReconstructionByASSD(reconstructed_label, input_label, stopping_threshold);
    }

    iftImage *tmp_img = iftLinearStretch(reconstructed_label, iftMinimumValue(reconstructed_label), iftMaximumValue(reconstructed_label),
                                         0, 255);
    iftDestroyImage(&reconstructed_label);
    reconstructed_label = tmp_img;

    tmp = iftConcatStrings(4, out_img_basename, "_", resuming_type, out_ext);
    out_reconstructed_label_file = iftJoinPathnames(2, out_dir, tmp);
    free(tmp);

    tmp = iftConcatStrings(5, out_img_basename, "_seeds", "_", resuming_type, out_ext);
    out_seeds_img_file = iftJoinPathnames(2, out_dir, tmp);
    free(tmp);

    tmp = iftConcatStrings(4, out_img_basename, "_", resuming_type, ".txt");
    out_seeds_file_txt = iftJoinPathnames(2, out_dir, tmp);
    free(tmp);

    tmp = iftConcatStrings(5, out_img_basename, "_basins", "_", resuming_type, out_ext);
    out_basins_file = iftJoinPathnames(2, out_dir, tmp);
    free(tmp);

    iftWriteSeeds(seeds, reconstructed_label, out_seeds_file_txt);
    iftWriteImageByExt(seed_image, out_seeds_img_file);
    iftWriteImageByExt(reconstructed_label, out_reconstructed_label_file);
    iftWriteImageByExt(basins, out_basins_file);

    /* Cleaning up! */
    free(out_seeds_img_file);
    free(out_reconstructed_label_file);
    free(out_basins_file);
    free(out_seeds_file_txt);
    free(resuming_type);
    free(stopping_criterion);
    free(secondary_stopping_criterion);

    iftDestroyLabeledSet(&seeds);
    iftDestroyImage(&reconstructed_label);
    iftDestroyImage(&seed_image);
    iftDestroyIGraph(&igraph);
    iftDestroyImageForest(&fst);
    iftDestroyMImage(&mimg);
}



void iftIGraphResumingDiffOrientedWatershed(iftIGraph *igraph, iftImage *input_label, iftLabeledSet *seeds,
                                            iftSet *trees_for_removal, iftDHeap *Q, double gamma)
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

    fprintf(stderr,"Is first oriented DIFT (IGraph Resuming Diff) iteration: %d\n", is_first_iteration);

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

                            double bndr = 1.0 - iftDiracDeltaFunction(input_label->val[p] - input_label->val[q]);

                            // If s is trying to conquer t, and the presegmentation labels of s and t differ,
                            // then bndr will be 1.0 and the weight will be multiplied by gamma. Otherwise,
                            // the regular weight is considered
                            tmp = iftMax(pvalue[s], (int)igraph->node[t].weight*pow(gamma, bndr));

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

                        double bndr = 1.0 - iftDiracDeltaFunction(input_label->val[p] - input_label->val[q]);

                        // If s is trying to conque t, and the presegmentation labels of s and t differ,
                        // then bndr will be 1.0 and the weight will be multiplied by gamma. Otherwise,
                        // the regular weight is considered
                        tmp = iftMax(pvalue[s], (int)igraph->node[t].weight*pow(gamma, bndr));

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

//    /* Fixing the label map after the first iteration */
//    if(!is_first_iteration)
//        iftIGraphFixLabelRootMap(igraph, &T);

    iftDestroySet(&adj);
    iftDestroySet(&T);
    iftDestroySet(&frontier_nodes);
}


/**
@brief Computes the watershed transform in a forest structure.

Complexity: O(|A| * |V|).

@param  fst  iftImageForest.    Forest structure created with a gradient image.
@param  seed iftLabeledSet.     List of spels with image index and seed label.
@param  removal_markers iftSet. List of spels marked for removal. NULL if empty

@return void. All forest maps are updated by reference.
*/
void iftDiffOrientedWatershedResuming(iftImageForest *fst, iftImage *input_label, iftLabeledSet *seed, iftSet * removal_markers, double gamma)
{
    iftAdjRel *A = fst->A;
    iftGQueue *Q = fst->Q;
    iftVoxel   u, v;
    int        i, p, q, tmp;
    iftSet    *Frontier = NULL;
    iftLabeledSet *S;
    iftBMap   *processed = fst->processed;
    iftImage  *pathval = fst->pathval, *pred = fst->pred, *label = fst->label;
    iftImage  *root = fst->root, *basins = fst->img, *marker = fst->marker;

    iftFillBMap(processed, 0);

    if (removal_markers != NULL)
    {
        Frontier = iftTreeRemoval(fst, removal_markers);
        while (Frontier != NULL)
        {
            p = iftRemoveSet(&Frontier);
            iftInsertGQueue(&Q, p);
        }
    }

    S = seed;
    while (S != NULL)
    {
        p = S->elem;

        if (Q->L.elem[p].color == IFT_GRAY)
        {
            /* p is also a frontier voxel, but the priority is it as a seed. */
            iftRemoveGQueueElem(Q, p);
        }

        label->val[p]   = S->label;
        pathval->val[p] = 0;
        root->val[p]    = p;
        pred->val[p]    = IFT_NIL;
        marker->val[p]  = S->marker;
        iftInsertGQueue(&Q, p);
        S = S->next;
    }

    /* Image Foresting Transform */
    while (!iftEmptyGQueue(Q))
    {
        p = iftRemoveGQueue(Q);
        u = iftGetVoxelCoord(basins, p);
        iftBMapSet1(processed, p);

        for (i = 1; i < A->n; i++)
        {
            v = iftGetAdjacentVoxel(A, u, i);
            if (iftValidVoxel(basins, v))
            {
                q = iftGetVoxelIndex(basins, v);

                if (Q->L.elem[q].color != IFT_BLACK) {
                    double bndr = 1.0-iftDiracDeltaFunction(input_label->val[p]-input_label->val[q]);

                    tmp = iftMax(pathval->val[p], iftRound(basins->val[q] * pow(gamma, bndr)));

                    /* if pred[q]=p then p and q belong to a tie-zone */
                    if ((tmp < pathval->val[q]) || ((pred->val[q] == p))) {
                        if (Q->L.elem[q].color == IFT_GRAY) {
                            iftRemoveGQueueElem(Q, q);
                        }
                        pred->val[q] = p;
                        root->val[q] = root->val[p];
                        label->val[q] = label->val[p];
                        marker->val[q] = marker->val[p];
                        pathval->val[q] = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }
    iftResetGQueue(Q);

}


int iftLocalMaximum(const iftImage *weights, int p, const iftAdjRel *disk)
{
    int local_max = p;
    iftVoxel u = iftGetVoxelCoord(weights, local_max);
    for (int i = 1; i < disk->n; i++) {
        iftVoxel v = iftGetAdjacentVoxel(disk, u, i);
        if (iftValidVoxel(weights, v)) {
            int q = iftGetVoxelIndex(weights, v);
            if (weights->val[local_max] < weights->val[q]) {
                local_max = q;
                u = iftGetVoxelCoord(weights, local_max);
            }
        }
    }
    return local_max;
}


int iftFurthestInError(const iftImage *source, const iftImage *target, const iftImage *mask, int from)
{
    iftVerifyImageDomains(source, mask, "iftFurthestInError");
    iftVerifyImageDomains(target, mask, "iftFurthestInError");

    iftAdjRel *disk = iftIs3DImage(source) ? iftSpheric(1.0) : iftCircular(1.0);
    iftImage *err = iftXor(source, target);
    iftImage *tmp = iftDilate(err, disk, NULL);
    iftDestroyImage(&err);
    err = tmp; tmp = NULL;

    iftSet *s = NULL; iftInsertSet(&s, from);
    iftImage *dist = iftGeodesicDistTransFromSet(err, s, NULL);

    iftImage *forbidden = iftMorphGrad(target, disk);

    int furthest = IFT_NIL, f_dist = IFT_NIL;
    for (int p = 0; p < source->n; p++) {
        if (err->val[p] && !forbidden->val[p] && (!mask || mask->val[p])) {
            if (f_dist < dist->val[p] && dist->val[p] != IFT_INFINITY_INT) {
                f_dist = dist->val[p];
                furthest = p;
            }
        }
    }

    iftDestroyImage(&err);
    iftDestroySet(&s);
    iftDestroyImage(&dist);
    iftDestroyAdjRel(&disk);

    return furthest;
}


static int _get_direction(const iftImage *obj, const iftAdjRel *A, int src)
{
    iftVoxel u = iftGetVoxelCoord(obj, src);

    // find object pixel
    int j = 1;
    for (; j < A->n; j++) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, j);
        if (iftValidVoxel(obj, v)) {
            int q = iftGetVoxelIndex(obj, v);
            if (obj->val[q])
                break;
        }
    }

    for (int i = 0; i < A->n; i++, j = (j + 1) % A->n) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, j);
        if (iftValidVoxel(obj, v)) {
            int q = iftGetVoxelIndex(obj, v);
            if (!obj->val[q])
                break;
        } else break;
    }

    return j;
}


iftSet *_sort_along_contour(const iftImage *obj, const iftImage *contour, const iftSet *S)
{
    iftAdjRel *A = iftClockCircular(sqrtf(2.0f));

    iftImage *count = iftCreateImage(contour->xsize, contour->ysize, contour->zsize);
    iftImage *pred = iftCreateImage(contour->xsize, contour->ysize, contour->zsize);

    iftFIFO *Q = iftCreateFIFO(contour->n);

    for (int i = 0; i < pred->n; i++) {
        pred->val[i] = IFT_NIL;
        count->val[i] = IFT_NIL;
    }

    int src = S->elem;
    int ref_orientation = _get_direction(obj, A, src);
    iftInsertFIFO(Q, src);

    int c = 0;
    while (!iftEmptyFIFO(Q))
    {
        int p = iftRemoveFIFO(Q);
        iftVoxel u = iftGetVoxelCoord(contour, p);
        count->val[p] = c++;

        int j = ref_orientation; // because of the selection of the starting point (no predecessor)
        if (pred->val[p] != IFT_NIL) {
            iftVoxel t = iftGetVoxelCoord(pred, pred->val[p]);
            j = ift6NeighborsClockFromLeft(&t, &u);
        }

        for (int i = 0; i < A->n; i++, j = (j + 1) % A->n)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, j);
            if (!iftValidVoxel(contour, v))
                continue;
            int q = iftGetVoxelIndex(contour, v);
            if (contour->val[q] != 0 && Q->color[q] != IFT_BLACK)
            {
                pred->val[q] = p;
                if (Q->color[q] == IFT_WHITE)
                    iftInsertFIFO(Q, q);
                break;
            }
        }
    }

    // trick to avoid unsorting nodes out of contour
    for (const iftSet *s = S, *prev = NULL; s; prev = s, s = s->next)
    {
        if (count->val[s->elem] == IFT_NIL && prev)
            count->val[s->elem] = count->val[prev->elem] + 1;
    }


    iftSet *out = NULL;
    iftSet *copy = iftSetCopyOrdered(S);
    while (copy)
    {
        iftSet *prev_max = NULL, *max = copy;
        for (iftSet *s = copy, *prev = NULL; s; prev = s, s = s->next) {
            if (count->val[max->elem] < count->val[s->elem]) {
                prev_max = prev;
                max = s;
            }
        }
        if (!prev_max)
            copy = max->next;
        else
            prev_max->next = max->next;

        max->next = out;
        out = max;
    }


    iftDestroyFIFO(&Q);
    iftDestroyImage(&count);
    iftDestroyImage(&pred);
    iftDestroyAdjRel(&A);

    return out;
}


iftSet *iftQueryForAnchorsPosition(const iftImage *gt_contour, const iftImage *source_mask,
                                   const iftImage *gt_mask, const iftSet *anchors)
{
    iftVerifyImageDomains(gt_contour, gt_mask, "iftQueryForAnchorsPosition");
    iftVerifyImageDomains(gt_contour, source_mask, "iftQueryForAnchorsPosition");

    iftImage *dist = iftCreateImageFromImage(gt_contour);
    iftImage *roots = iftCreateImageFromImage(gt_contour);
    iftImage *found = iftCreateImageFromImage(gt_contour);

    iftGQueue *Q = iftCreateGQueue(IFT_QSIZE, gt_contour->n, dist->val);

    for (int p = 0; p < gt_contour->n; p++) {
        roots->val[p] = IFT_NIL;
        dist->val[p] = IFT_INFINITY_INT;
        found->val[p] = 0;
    }

    iftImage *error = iftXor(gt_mask, source_mask);

    iftAdjRel *A = iftIs3DImage(gt_contour) ? iftSpheric(1.75) : iftCircular(1.5);
    iftImage *tmp = iftDilate(error, A, NULL);
    iftDestroyImage(&error);
    error = tmp; tmp = NULL;

    for (const iftSet *s = anchors; s; s = s->next)
    {
        int p = s->elem;
        dist->val[p] = 0;
        roots->val[p] = p;
        if (error->val[p]) {
            iftInsertGQueue(&Q, p);
        }
    }

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        int r = roots->val[p];
        // already found closest
        if (found->val[r])
            continue;

        iftVoxel u = iftGetVoxelCoord(gt_contour, p);
        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (!iftValidVoxel(gt_contour, v))
                continue;

            int q = iftGetVoxelIndex(gt_contour, v);
            if (error->val[q] && dist->val[q] > (dist->val[p] + 1))
            {
                if (Q->L.elem[q].color == IFT_GRAY)
                    iftRemoveGQueueElem(Q, q);
                dist->val[q] = dist->val[p] + 1;
                roots->val[q] = r;
                if (gt_contour->val[q])
                    found->val[r] = 1;
                // when closest found stop propagating from tree
                if (!found->val[r])
                    iftInsertGQueue(&Q, q);
            }
        }
    }

    for (int p = 0; p < gt_contour->n; p++) {
        if (gt_contour->val[p] && roots->val[p] != IFT_NIL) {
            int r = roots->val[p];
            roots->val[r] = p;
        }
    }

    iftSet *new_anchors = iftSetCopyOrdered(anchors);
    for (iftSet *s = new_anchors; s; s = s->next) {
        s->elem = roots->val[s->elem];
    }

    iftSet *validated = _sort_along_contour(gt_mask, gt_contour, new_anchors);
    iftDestroySet(&new_anchors);

    iftDestroyImage(&error);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&dist);
    iftDestroyImage(&roots);
    iftDestroyImage(&found);
    iftDestroyGQueue(&Q);

    return validated;
}


iftSet *iftFurtherThanThreshold(const iftVoxelArray *anchors,
                                const iftImage *mask,
                                float threshold)
{
    iftAdjRel *A = iftCircular(threshold);
    iftSet *set = NULL;
    for (int i = 0; i < anchors->n; i++) {
        bool found = false;
        for (int j = 0; j < A->n && !found; j++) {
            iftVoxel v = iftGetAdjacentVoxel(A, anchors->val[i], j);
            if (iftValidVoxel(mask, v)) {
                int p = iftGetVoxelIndex(mask, v);
                if (mask->val[p])
                    found = true;
            }
        }
        if (!found)
            iftInsertSet(&set, i);
    }
    iftDestroyAdjRel(&A);

    return set;
}
