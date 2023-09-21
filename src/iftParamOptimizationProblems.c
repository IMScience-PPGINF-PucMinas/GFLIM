#include "iftParamOptimizationProblems.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Stream.h"
#include "iftSegmentation.h"


/********************** PRIVATE FUNCTIONS *************************/
iftAdjRel *_iftTranslateModelSeedDisps(const iftAdjRel *seed_disps, iftVector trans_vec) {
    iftAdjRel *trans_seed_disps = iftCreateAdjRel(seed_disps->n);

    // translates the object model's seeds
    for (int i = 0; i < seed_disps->n; i++) {
        trans_seed_disps->dx[i] = seed_disps->dx[i] + trans_vec.x;
        trans_seed_disps->dy[i] = seed_disps->dy[i] + trans_vec.y;
        trans_seed_disps->dz[i] = seed_disps->dz[i] + trans_vec.z;
    }

    return trans_seed_disps;
}

/******************************************************************/


/********************** PUBLIC FUNCTIONS *************************/
iftMatrix *iftRegularMSPSDeltaMatrix(int n_params, int n_scales, int stride) {
    iftMatrix *delta = iftCreateMatrix(n_params, n_scales);

    for (int s = 0; s < n_scales; s++) {
        for (int i = 0; i < n_params; i++) {
            iftMatrixElem(delta, i, s) = (1 + s*stride);
        }
    }

    return delta;
}


iftGradMatchProb *iftCreateGradMatchProb(const iftImage *fixed_grad_img, const iftImage *mov_grad_img, iftVector max_disp) {
    iftVerifyImages(fixed_grad_img, mov_grad_img, "iftCreateGradMatchProb");

    iftGradMatchProb *prob = (iftGradMatchProb*) iftAlloc(1, sizeof(iftGradMatchProb));

    prob->fixed_grad_img = fixed_grad_img;
    prob->mov_grad_img   = mov_grad_img;
    prob->mov_grad_elems = NULL;

    iftList *L = iftCreateList();

    for (int p = 0; p < mov_grad_img->n; p++) {
        if (mov_grad_img->val[p] > 0) {
            iftInsertListIntoTail(L, p);
        }
    }

    prob->mov_grad_elems = iftCreateIntArray(L->n);
    size_t i             = 0;
    while (L->n != 0) {
        prob->mov_grad_elems->val[i++] = iftRemoveListHead(L);
    }
    iftDestroyList(&L);

    prob->max_disp.x = max_disp.x;
    prob->max_disp.y = max_disp.y;
    prob->max_disp.z = max_disp.z;


    return prob;
}


void iftDestroyGradMatchProb(iftGradMatchProb **prob) {
    if (prob != NULL) {
        iftGradMatchProb *aux = *prob;

        if (aux != NULL) {
            iftDestroyIntArray(&aux->mov_grad_elems);
            iftFree(aux);
            *prob = NULL;
        }
    }
}


iftPatchLocalizationNMIProb *iftCreatePatchLocalizationNMIProb(const iftImage *test_img,
                                    const iftImage *template_img, iftBoundingBox true_template_patch,
                                    iftBoundingBox test_patch, iftBoundingBox search_region) {
    iftPatchLocalizationNMIProb *prob = iftAlloc(1, sizeof(iftPatchLocalizationNMIProb));
    prob->test_img                          = iftCopyImage(test_img);
    prob->template_img                      = iftCopyImage(template_img);
    prob->true_template_patch               = true_template_patch;
    prob->test_patch                        = test_patch;
    prob->search_region                     = search_region;

    return prob;
}


void iftDestroyPatchLocalizationNMIProb(iftPatchLocalizationNMIProb **prob) {
    iftPatchLocalizationNMIProb *aux = *prob;
    if (aux != NULL) {
        iftDestroyImage(&aux->test_img);
        iftDestroyImage(&aux->template_img);
        iftFree(aux);
        *prob = NULL;
    }
}


float iftPatchNMI(void *void_prob, float *theta) {
    iftPatchLocalizationNMIProb *prob = (iftPatchLocalizationNMIProb *) void_prob;
    puts("patch");
    iftPrintBoundingBox(prob->test_patch);
    puts("template_patch");
    iftPrintBoundingBox(prob->true_template_patch);
    puts("search_region");
    iftPrintBoundingBox(prob->search_region);

    iftVoxel test_patch_center = iftBoundingBoxCenterVoxel(prob->test_patch);
    printf("disp: (%.0f, %.0f, %.0f)\n", theta[0], theta[1], theta[2]);
    printf("- test_patch_center: (%d, %d, %d)\n", test_patch_center.x, test_patch_center.y, test_patch_center.z);
    test_patch_center.x += iftRound(theta[0]);
    test_patch_center.y += iftRound(theta[1]);
    test_patch_center.z += iftRound(theta[2]);
    printf("- test_patch_center: (%d, %d, %d)\n", test_patch_center.x, test_patch_center.y, test_patch_center.z);


    if ((test_patch_center.x < prob->search_region.begin.x || test_patch_center.x > prob->search_region.end.x) ||
        (test_patch_center.y < prob->search_region.begin.y || test_patch_center.y > prob->search_region.end.y) ||
        (test_patch_center.z < prob->search_region.begin.z || test_patch_center.z > prob->search_region.end.z)) {
        printf("score = 0.0 (out of search region)\n\n");
        return 0.0;
    }

    iftBoundingBox new_test_patch = iftCentralizeBoundingBox(prob->test_patch, test_patch_center);

    iftImage *test_img_patch = iftExtractROI(prob->test_img, new_test_patch);
    iftImage *template_img_patch = iftExtractROI(prob->template_img, prob->true_template_patch);
    
    float score = iftNormalizedMutualInformation(test_img_patch, template_img_patch);

    iftDestroyImage(&test_img_patch);
    iftDestroyImage(&template_img_patch);

    printf("score = %f\n\n", score);

    return score;
}


iftBoundingBox iftMSPSMaxPatchNMI(const iftImage *img, const iftImage *template_img, iftBoundingBox true_template_patch,
                                        iftBoundingBox patch, iftBoundingBox search_region) {
    int n_params = 3; // disp vector: x, y, z
    int n_scales = 4;
    int stride   = 1;

    iftPatchLocalizationNMIProb *prob = iftCreatePatchLocalizationNMIProb(img, template_img, true_template_patch, patch, search_region);
    puts("true_template_patch");
    iftPrintBoundingBox(true_template_patch);
    puts("patch");
    iftPrintBoundingBox(patch);
    // iftPatchLocalizationNMIProb *prob = iftCreatePatchLocalizationNMIProb(iftMedianFilter(img, NULL), iftMedianFilter(template_img, NULL), true_template_patch, patch, search_region);
    iftVoxel test_patch_center = iftBoundingBoxCenterVoxel(prob->test_patch);
    iftMSPS *msps = iftCreateMSPS(n_params, n_scales, iftPatchNMI, prob);

    // computes a regular translation (displacement) of the moving image's coordinates for each scale.
    iftDestroyMatrix(&msps->delta);
    msps->delta = iftRegularMSPSDeltaMatrix(n_params, n_scales, stride);
    puts("Delta: MSPS");
    iftPrintMatrix(msps->delta);

    // run optimization
    iftMSPSMax(msps);
    printf("* best displacement: (%d, %d, %d)\n--------------------------\n\n",
            iftRound(msps->theta[0]), iftRound(msps->theta[1]), iftRound(msps->theta[2]));

    // move the input subcortical bb to the best position
    test_patch_center.x += iftRound(msps->theta[0]);
    test_patch_center.y += iftRound(msps->theta[1]);
    test_patch_center.z += iftRound(msps->theta[2]);

    iftBoundingBox best_patch = iftCentralizeBoundingBox(patch, test_patch_center);

    iftDestroyPatchLocalizationNMIProb(&prob);
    iftDestroyMSPS(&msps);

    return best_patch;
}


float iftMatchGradProbMean(void *problem, float *theta) {
    iftGradMatchProb *prob = (iftGradMatchProb*) problem;
    int dx = iftRound(theta[0]);
    int dy = iftRound(theta[1]);
    int dz = (iftIs3DImage(prob->fixed_grad_img)) ? iftRound(theta[2]) : 0;

    printf("disp: (%d, %d, %d)\n", dx, dy, dz);
    printf("max_disp: (%f, %f, %f)\n", prob->max_disp.x, prob->max_disp.y, prob->max_disp.z);

    // displacement vector out of the search region - it will be analyzed
    if ((abs(dx) > prob->max_disp.x) || (abs(dy) > prob->max_disp.y) || (abs(dz) > prob->max_disp.z)) {
        printf("* out trans_vec: (%d, %d, %d)\n", dx, dy, dz);
        printf("score = 0\n\n");
        return 0.0f;
    }

    float w_mean_grad  = 0.0;
    int n_valid_voxels = 0;

    // reduction operations may not be associative for real numbers, then I only use integer datatypes
    // #pragma omp parallel for reduction(+:w_mean_grad) reduction(+:n_valid_voxels)
    for (size_t i = 0; i < prob->mov_grad_elems->n; i++) {
        int p      = prob->mov_grad_elems->val[i]; // voxel index from the moving gradient
        iftVoxel u = iftGetVoxelCoord(prob->mov_grad_img, p); // voxel from the moving gradient

        // translated voxel position on the Fixed Image where Translated Moving Gradient's Voxel is
        iftVoxel v;
        v.x = u.x + dx;
        v.y = u.y + dy;
        v.z = u.z + dz;

        if (iftValidVoxel(prob->fixed_grad_img, v)) {
            int q = iftGetVoxelIndex(prob->fixed_grad_img, v);

            w_mean_grad += labs(prob->fixed_grad_img->val[q] - prob->mov_grad_img->val[p]);
            // w_mean_grad += prob->fixed_grad_img->val[q];
            n_valid_voxels++;
        }
    }

    if (n_valid_voxels == 0) {
        return IFT_INFINITY_INT;
    }

    float score = w_mean_grad / n_valid_voxels;
    // printf("* disp (%d, %d, %d) = %f\n", dx, dy, dz, score);

    // iftVector disp = {.x = dx, .y = dy, .z = dz};
    // iftImage *model_grad_trans = iftTranslateImageContent(prob->mov_grad_img, disp);
    // iftImage *merged = iftAdd(model_grad_trans, prob->fixed_grad_img);
    // iftWriteImageByExt(merged, "tmp/merged_%d_%d_%d.scn", dx, dy, dz);
    // iftDestroyImage(&model_grad_trans);
    // iftDestroyImage(&merged);

    printf("score = %f\n\n", score);

    return score;
}



iftVector iftMSPSMatchGradients(const iftImage *fixed_grad_img, const iftImage *mov_grad_img, int n_scales, int stride, iftVector max_disp) {
    // the coordinates x, y, and z (if 3D image) are the parameters to be optmized
    int n_params           = (iftIs3DImage(fixed_grad_img)) ? 3 : 2;
    iftGradMatchProb *prob = iftCreateGradMatchProb(fixed_grad_img, mov_grad_img, max_disp);
    iftMSPS *msps          = iftCreateMSPS(n_params, n_scales, iftMatchGradProbMean, prob);

    // computes a regular translation (displacement) of the moving image's coordinates for each scale.
    iftDestroyMatrix(&msps->delta);
    msps->delta = iftRegularMSPSDeltaMatrix(n_params, n_scales, stride);

    iftMSPSMax(msps);

    iftVector best_disp;
    best_disp.x = iftRound(msps->theta[0]);
    best_disp.y = iftRound(msps->theta[1]);
    best_disp.z = 0;
    if (n_params == 3)
        best_disp.z = iftRound(msps->theta[2]);

    iftDestroyGradMatchProb(&prob);
    iftDestroyMSPS(&msps);

    return best_disp;
}


iftObjModelMeanArcWeight *iftCreateObjModelMeanArcWeightProb(const iftImage *grad_test_img, int label,
                                                             const iftLabeledSet *seeds, const iftSet *certain_obj_region,
                                                             const iftSet *forbidden, iftVector max_disp) {
    iftObjModelMeanArcWeight *prob = (iftObjModelMeanArcWeight*) iftAlloc(1, sizeof(iftObjModelMeanArcWeight));

    prob->grad_test_img      = grad_test_img;
    prob->label              = label;
    prob->seeds              = seeds;
    prob->certain_obj_region = certain_obj_region;
    prob->forbidden          = forbidden;
    prob->max_disp           = max_disp;

    return prob;
}


void iftDestroyObjModelMeanArcWeightProb(iftObjModelMeanArcWeight **prob) {
    iftObjModelMeanArcWeight *aux = *prob;

    if (aux != NULL) {
        iftFree(aux);
        *prob = NULL;
    }
}


iftVector iftMSPSObjModel(const iftImage *grad_test_img, int label, const iftLabeledSet *seeds,
                          const iftSet *certain_obj_region, const iftSet *forbidden, 
                          int n_scales, int stride, iftVector max_disp) {
    // the coordinates x, y, and z (if 3D image) are the parameters to be optimized
    int n_params = (iftIs3DImage(grad_test_img)) ? 3 : 2;

    iftObjModelMeanArcWeight *prob = iftCreateObjModelMeanArcWeightProb(grad_test_img, label, seeds, certain_obj_region, forbidden, max_disp);
    iftMSPS *msps                  = iftCreateMSPS(n_params, n_scales, iftComputeObjModelMeanArcWeight, prob);

    // computes a regular translation (displacement) of the moving image's coordinates for each scale.
    iftDestroyMatrix(&msps->delta);
    msps->delta = iftRegularMSPSDeltaMatrix(n_params, n_scales, stride);

    iftMSPSMax(msps);

    iftVector best_disp;
    best_disp.x = iftRound(msps->theta[0]);
    best_disp.y = iftRound(msps->theta[1]);
    best_disp.z = 0;
    if (n_params == 3)
        best_disp.z = iftRound(msps->theta[2]);

    printf("- best disp (%0.f, %0.f, %0.f)\n\n", best_disp.x, best_disp.y, best_disp.z);

    // cleaning up
    iftDestroyObjModelMeanArcWeightProb(&prob);
    iftDestroyMSPS(&msps);

    return best_disp;
}


float iftComputeObjModelMeanArcWeight(void *problem, float *theta) {
    iftObjModelMeanArcWeight *prob = (iftObjModelMeanArcWeight*) problem;

    iftVector trans_vec;
    trans_vec.x = theta[0];
    trans_vec.y = theta[1];
    trans_vec.z = (iftIs3DImage(prob->grad_test_img)) ? theta[2] : 0;

    printf("trans_vec: (%f, %f, %f)\n", trans_vec.x, trans_vec.y, trans_vec.z);
    printf("max_disp: (%f, %f, %f)\n", prob->max_disp.x, prob->max_disp.y, prob->max_disp.z);

    // displacement vector out of the search region - it will be analyzed
    if ((fabs(trans_vec.x) > prob->max_disp.x) || (fabs(trans_vec.y) > prob->max_disp.y) || (fabs(trans_vec.z) > prob->max_disp.z)) {
        printf("* out trans_vec: (%f, %f, %f)\n", trans_vec.x, trans_vec.y, trans_vec.z);
        printf("score = 0\n\n");
        return 0.0f;
    }

    iftLabeledSet *trans_seeds       = iftTranslateLabeledSet(prob->seeds, prob->grad_test_img, trans_vec);
    iftSet *trans_certain_obj_region = iftTranslateSet(prob->certain_obj_region, prob->grad_test_img, trans_vec);
    iftSet *trans_forbidden          = iftTranslateSet(prob->forbidden, prob->grad_test_img, trans_vec);

    iftImage *seg_img = iftWatershed(prob->grad_test_img, NULL, trans_seeds, trans_forbidden);
    iftSetToImage(trans_certain_obj_region, seg_img, prob->label);
    // iftWriteImageByExt(seg_img, "tmp/obj%d_%.0f_%.0f_%.0f.hdr", prob->label, trans_vec.x, trans_vec.y, trans_vec.z);

    // iftImage *seg_img          = iftOrientedWatershed(prob->grad_test_img, prob->grad_test_img, NULL,
    //                                                   trans_seeds, 1, NULL);

    // char path[64];
    // sprintf(path, "tmp/%.0f_%.0f_%.0f_all_seeds.txt", trans_vec.x, trans_vec.y, trans_vec.z);
    // iftWriteSeeds(path, trans_seeds, prob->grad_test_img);


    // computes the score
    float score  = 0.0f;
    int n_voxels = 0;

    iftSet *borders = iftObjectBorderSet(seg_img, NULL);
    iftSet *S       = borders;
    while (S != NULL) {
        int p  = iftRemoveSet(&S);
        score += prob->grad_test_img->val[p];
        n_voxels++;
    }

    score /= (n_voxels*1.0);
    printf("* disp (%0.f, %0.f, %0.f) = %f\n", trans_vec.x, trans_vec.y, trans_vec.z, score);

    // DESTROYERS
    iftDestroyLabeledSet(&trans_seeds);
    iftDestroySet(&trans_certain_obj_region);
    iftDestroySet(&trans_forbidden);
    iftDestroyImage(&seg_img);

    return score;
}
/******************************************************************/























