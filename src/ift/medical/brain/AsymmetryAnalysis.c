#include "ift/medical/brain/AsymmetryAnalysis.h"

#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "ift/imgproc/basic/Histogram.h"
#include "iftGeometric.h"
#include "iftIGraph.h"
#include "iftImageMath.h"
#include "iftSeeds.h"


/************************** PRIVATE FUNCTIONS **************************/
/**
 * @brief Builds a Multi-Band Image for OISF algorithm by stacking the right and the flipped left brain sides, and the
 * (normalized) asymmetry map between them.
 *
 * @attention The function assumes that the image's Mid-Sagittal Plane is its central xslice.
 *
 * The resulting Multi-Band Image has 3 bands: the first is the input image, the second is the flipped input image,
 * the third is the normalized ([0, 1]) asymmetry map.
 * By assuming that the mid-sagittal plane is the central xslice of the input image, then each voxel
 * will have two values: one from the right brain side and the other from the flipped left brain side.
 *
 * Note that the input asymmetry map <asym_map> is not normalized.
 * We first get its maximum value (max asymmetry) and the used it to divide/normalize the asymmetry map.
 *
 * @param img Brain image.
 * @param asym_map Asymmetry Map.
 * @return Multi-Band Image with the stacked brain hemispheres.
 *
 * @author Samuel Martins
 * @date Sep 1, 2018
 */
iftMImage *_iftBuildBrainHemisphereMImageOISF(const iftImage *img, const iftImage *asym_map) {
    iftMImage *mimg = iftCreateMImage(img->xsize, img->ysize, img->zsize, 3);
    iftImage *flip_img = iftFlipImage(img, IFT_AXIS_X);
    float max_asym = iftMaximumValue(asym_map) * 1.0;

    #pragma omp parallel for
    for (int p = 0; p < mimg->n; p++) {
        mimg->val[p][0] = img->val[p];
        mimg->val[p][1] = flip_img->val[p];
        mimg->val[p][2] = asym_map->val[p] / max_asym;
    }
    iftDestroyImage(&flip_img);
    
    return mimg;
}


/************************** PUBLIC FUNCTIONS **************************/
iftImage *iftBrainAsymMap(const iftImage *img, const iftImage *bias) {
    iftImage *img_flip = iftFlipImage(img, IFT_AXIS_X);

    iftImage *diff = NULL;
    if (bias) {
        iftVerifyImages(img, bias, "iftBrainAsymMap");
        diff = iftCreateImageFromImage(img);

        #pragma omp parallel for
        for (int p = 0; p < diff->n; p++)
            diff->val[p] = iftMax(0, abs(img->val[p] - img_flip->val[p]) - bias->val[p]);
    }
    else diff = iftAbsSub(img, img_flip);

    iftDestroyImage(&img_flip);

    return diff;
}


iftImage *iftBrainDiffMap(const iftImage *img, const iftImage *template_img, const iftImage *bias) {
    iftVerifyImages(img, template_img, "iftBrainDiffMap");
    
    iftImage *diff = NULL;
    
    if (bias) {
        iftVerifyImages(img, bias, "iftBrainDiffMap");
        diff = iftCreateImageFromImage(img);

        #pragma omp parallel for
        for (int p = 0; p < diff->n; p++)
            diff->val[p] = iftMax(0, abs(img->val[p] - template_img->val[p]) - bias->val[p]);
    }
    else diff = iftAbsSub(img, template_img);
    
    return diff;
}



iftImage *iftMeanBrainAsymMap(const iftFileSet *img_set, bool add_stdev_asymmetries) {
    iftImage *img0 = iftReadImageByExt(img_set->files[0]->path);
    iftFImage *mean_asym = iftCreateFImage(img0->xsize, img0->ysize, img0->zsize);
    iftDestroyImage(&img0);
    
    for (int i = 0; i < img_set->n; i++) {
        iftImage *img = iftReadImageByExt(img_set->files[i]->path);
        iftImage *img_asym = iftBrainAsymMap(img, NULL);

        #pragma omp parallel for
        for (int p = 0; p < img_asym->n; p++)
            mean_asym->val[p] += (img_asym->val[p] / ((float) img_set->n));

        iftDestroyImage(&img);
        iftDestroyImage(&img_asym);
    }
    
    
    if (add_stdev_asymmetries) {
        iftFImage *stdev_asym = iftCreateFImage(mean_asym->xsize, mean_asym->ysize, mean_asym->zsize);
        
        for (int i = 0; i < img_set->n; i++) {
            iftImage *img = iftReadImageByExt(img_set->files[i]->path);
            iftImage *img_asym = iftBrainAsymMap(img, NULL);
            
            #pragma omp parallel for
            for (int p = 0; p < img_asym->n; p++)
                stdev_asym->val[p] += powf(img_asym->val[p] - mean_asym->val[p], 2);
            
            iftDestroyImage(&img);
            iftDestroyImage(&img_asym);
        }

        // adding stdev asymmetries to the mean asymmetries
        #pragma omp parallel for
        for (int p = 0; p < stdev_asym->n; p++)
            mean_asym->val[p] = mean_asym->val[p] + (sqrtf(stdev_asym->val[p] / img_set->n));
        
        iftDestroyFImage(&stdev_asym);
    }
    
    iftImage *mean_asym_map = iftRoundFImage(mean_asym);
    
    iftDestroyFImage(&mean_asym);
    
    return mean_asym_map;
}


iftImage *iftMeanBrainDiffMap(const iftFileSet *img_set, const iftImage *template_img, bool add_stdev_asymmetries) {
    iftFImage *mean_diff = iftCreateFImage(template_img->xsize, template_img->ysize, template_img->zsize);
    
    for (int i = 0; i < img_set->n; i++) {
        iftImage *img = iftReadImageByExt(img_set->files[i]->path);
        iftImage *img_diff = iftAbsSub(img, template_img);

        #pragma omp parallel for
        for (int p = 0; p < img_diff->n; p++)
            mean_diff->val[p] += (img_diff->val[p] / ((float) img_set->n));
        
        iftDestroyImage(&img);
        iftDestroyImage(&img_diff);
    }
    
    
    if (add_stdev_asymmetries) {
        iftFImage *stdev_asym = iftCreateFImage(mean_diff->xsize, mean_diff->ysize, mean_diff->zsize);
        
        for (int i = 0; i < img_set->n; i++) {
            iftImage *img = iftReadImageByExt(img_set->files[i]->path);
            iftImage *img_diff = iftAbsSub(img, template_img);

        #pragma omp parallel for
            for (int p = 0; p < img_diff->n; p++)
                stdev_asym->val[p] += powf(img_diff->val[p] - mean_diff->val[p], 2);
            
            iftDestroyImage(&img);
            iftDestroyImage(&img_diff);
        }
        
        // adding stdev asymmetries to the mean asymmetries
        #pragma omp parallel for
        for (int p = 0; p < stdev_asym->n; p++)
            mean_diff->val[p] = mean_diff->val[p] + (sqrtf(stdev_asym->val[p] / img_set->n));
        
        iftDestroyFImage(&stdev_asym);
    }
    
    iftImage *mean_diff_map = iftRoundFImage(mean_diff);
    
    iftDestroyFImage(&mean_diff);
    
    return mean_diff_map;
}


iftIntArray * iftGridSamplingByBrainAsymmetries(const iftImage *asym_map, const iftImage *bin_mask,
                                                int min_samples_on_symmetries, float thres_factor) {
    iftImage *asym_map_obj = iftMask(asym_map, bin_mask);
//    iftWriteImageByExt(asym_map_obj, "tmp/asym_map_obj.nii.gz");
    iftBoundingBox bb = iftMinBoundingBox(bin_mask, NULL);
    iftImage *asym_map_obj_roi = iftExtractROI(asym_map_obj, bb);
//    iftWriteImageByExt(asym_map_obj_roi, "tmp/asym_map_obj_roi.nii.gz");
    float otsu = iftOtsu(asym_map_obj_roi);
    float thres = thres_factor * otsu;
    printf("### otsu = %f\n", otsu);
    printf("### thres = %f\n\n", thres);
    iftDestroyImage(&asym_map_obj_roi);

    // mask with the asymmetries as binary componentes
    iftImage *asym_comps = iftThreshold(asym_map_obj, thres, IFT_INFINITY_INT, 1);
//    iftWriteImageByExt(asym_comps, "tmp/asym_comps.nii.gz");
    iftDestroyImage(&asym_map_obj);

    iftImage *sym_comps = iftComplement(asym_comps);
    iftImage *sym_comps_obj = iftMask(sym_comps, bin_mask);
//    iftWriteImageByExt(sym_comps_obj, "tmp/sym_comps_obj.nii.gz");
    iftDestroyImage(&sym_comps);
    
    iftImage *asym_comps_open = iftOpenBin(asym_comps, 1.0);
//    iftWriteImageByExt(asym_comps_open, "tmp/asym_comps_2xotsu_open.nii.gz");
    iftDestroyImage(&asym_comps);

    iftAdjRel *B = (iftIs3DImage(asym_map)) ? iftSpheric(1.74) : iftCircular(1.4);
    iftImage *label_comps = iftFastLabelComp(asym_comps_open, B);
//    iftWriteImageByExt(label_comps, "tmp/label_comps.nii.gz");
    
    iftDestroyAdjRel(&B);
    iftDestroyImage(&asym_comps_open);

    // grid points on the asymmetric regions
    iftIntArray *grid_asym = NULL;
    iftIntArray *max_obj_vals = iftMaximumObjectValues(asym_map, label_comps, &grid_asym);
    iftDestroyIntArray(&max_obj_vals);
    iftDestroyImage(&label_comps);

    int n_samples_on_asymmetries = grid_asym->n - 1; // ignore the value for the background [0]
//    int n_samples_on_symmetries = iftMax(n_samples - n_samples_on_asymmetries, min_samples_on_symmetries);
    int n_samples_on_symmetries = min_samples_on_symmetries;
     printf("n_samples_on_asymmetries: %d\n", n_samples_on_asymmetries);
     printf("n_samples_on_symmetries: %d\n", n_samples_on_symmetries);

    // grid points on the symmetric regions
    float radius = iftEstimateGridOnMaskSamplingRadius(sym_comps_obj, -1, n_samples_on_symmetries);
    iftIntArray *grid_sym = iftGridSamplingOnMask(sym_comps_obj, radius, -1, n_samples_on_symmetries);
    n_samples_on_symmetries = grid_sym->n;
    iftDestroyImage(&sym_comps_obj);
    printf("n_samples_on_symmetries: %d\n", n_samples_on_symmetries);

    iftIntArray *grid = iftCreateIntArray(n_samples_on_asymmetries + n_samples_on_symmetries);
    for (int i = 0; i < n_samples_on_asymmetries; i++) {
        int label = i + 1;
        grid->val[i] = grid_asym->val[label];
    }
    for (int i = 0; i < n_samples_on_symmetries; i++)
        grid->val[n_samples_on_asymmetries + i] = grid_sym->val[i];
    printf("** grid->n: %ld\n", grid->n);

    iftDestroyIntArray(&grid_asym);
    iftDestroyIntArray(&grid_sym);

    return grid;
}


iftMImage *iftBuildBrainHemisphereMImage(const iftImage *img) {
    iftMImage *mimg = iftCreateMImage(img->xsize, img->ysize, img->zsize, 2);
    iftImage *flip_img = iftFlipImage(img, IFT_AXIS_X);

    #pragma omp parallel for
    for (int p = 0; p < mimg->n; p++) {
        mimg->val[p][0] = img->val[p];
        mimg->val[p][1] = flip_img->val[p];
    }
    iftDestroyImage(&flip_img);

    return mimg;
}


iftMImage *iftBuildBrainHemisphereMImageAsym(const iftImage *img, const iftImage *asym_map) {
    iftMImage *mimg = iftCreateMImage(img->xsize, img->ysize, img->zsize, 3);
    iftImage *flip_img = iftFlipImage(img, IFT_AXIS_X);

    #pragma omp parallel for
    for (int p = 0; p < mimg->n; p++) {
        mimg->val[p][0] = img->val[p];
        mimg->val[p][1] = flip_img->val[p];
        mimg->val[p][2] = asym_map->val[p];
    }
    iftDestroyImage(&flip_img);

    return mimg;
}


/**
 * @brief Extracts the HAA feats (Histogram of Absolute Asymmetries) from an asymmetry map.
 *
 * This functions extracts HAA feats from an asymmetry map for each of its supervoxels.
 * The extracted supervoxels' feature vectors are stored in a set of datasets, one per supervoxel.
 * By assuming a set of <n_supervoxels> supervoxels, an array <Zarr> of <n_supervoxels> + 1 (for convenience), and
 * that the sample index from the asymmetry map is <s>, the extracted feature vector for each supervoxel i is stored
 * in the dataset Zarr[i] at the sample <s>.
 *
 * @warning The number of features of all supervoxel datasets must be equal to the number of histogram bins.
 *
 * @param asym_map Map with the absolute asymmetries between the brain sides.
 * @param svoxels_brain_side Label image with the symmetric supervoxels (from 1 to n_supervoxels) only from one brain side.
 * @param Zarr Array of datasets, one per supervoxel.
 * @param s Sample index in the datasets for the input asymmetry map.
 * @param n_bins Number of bins of the histograms.
 * @param n_supervoxels Number of supervoxels from the Supervoxel label image.
 *
 * @author Samuel Martins
 * @date Dec 11, 2018
 */
void iftSupervoxelHAAFeatsInAsymMap(const iftImage *asym_map, const iftImage *svoxels_brain_side,
                                    iftDataSet **Zarr, int s, int n_bins, int max_val, int n_supervoxels) {
    iftHist **hists = iftCalcGrayImageHistForLabels(asym_map, svoxels_brain_side, n_bins, max_val, true, NULL); // size n_supervoxels + 1

    #pragma omp parallel for
    for (int label = 1; label <= n_supervoxels; label++) {
        if (Zarr[label]->nfeats != n_bins)
            iftError("Number of Feats for Supervoxel DataSet [%d] is != from the number of histogram bins: %d != %d",
                        "iftSupervoxelHAAFeatsInAsymMap", label, Zarr[label]->nfeats, n_bins);
        for (int f = 0; f < Zarr[label]->nfeats; f++) {
            
            Zarr[label]->sample[s].feat[f] = hists[label]->val[f];
        }
        iftDestroyHist(&hists[label]);
    }
    iftDestroyHist(&hists[0]);
    iftFree(hists);
}

void iftSupervoxelBICFeatsInAsymMap(const iftImage *asym_map, const iftImage *svoxels_brain_side,
                                    iftDataSet **Zarr, int s, int n_bins_per_channel, int n_supervoxels) {
    iftMatrix *feats_mat = iftExtractBICForLabels(asym_map, svoxels_brain_side, n_bins_per_channel);

    #pragma omp parallel for
        for (int label = 1; label <= n_supervoxels; label++) {
            if (Zarr[label]->nfeats != feats_mat->ncols)
                iftError("Number of Feats for Supervoxel DataSet [%d] is != from the number of BIC feats: %d != %d",
                         "iftSupervoxelBICFeatsInAsymMap", label, Zarr[label]->nfeats, feats_mat->ncols);
            for (int f = 0; f < Zarr[label]->nfeats; f++) {
                Zarr[label]->sample[s].feat[f] = iftMatrixElem(feats_mat, f, label);
        }
    }

    iftDestroyMatrix(&feats_mat);
}

iftImage *iftSymmISF(const iftImage *img, const iftImage *bin_mask, float alpha, float beta, float thres_factor,
                     float min_dist_to_border, int n_seeds_on_symmetric_regions, const iftImage *normal_asymmap) {
    iftImage *asym_map = iftBrainAsymMap(img, normal_asymmap);
//    iftWriteImageByExt(asym_map, "tmp/asym_map.nii.gz");
    
    iftIntArray *grid = NULL;
    if (min_dist_to_border > 0.0) {
        puts("- eroding");
        iftSet *S = NULL;
        iftImage *bin_mask_eroded = iftErodeBin(bin_mask, &S, min_dist_to_border);
        
        grid = iftGridSamplingByBrainAsymmetries(asym_map, bin_mask_eroded, n_seeds_on_symmetric_regions, thres_factor);
        
        iftDestroySet(&S);
        iftDestroyImage(&bin_mask_eroded);
    }
    else grid = iftGridSamplingByBrainAsymmetries(asym_map, bin_mask, n_seeds_on_symmetric_regions, thres_factor);
    
    iftImage *seeds = iftCreateImageFromImage(img);
    iftIntArrayToImage(grid, seeds, 1);
    iftMImage *mimg = iftBuildBrainHemisphereMImage(img);
//    iftMImage *mimg = iftBuildBrainHemisphereMImageAsym(img, asym_map);
    
    iftAdjRel *A = iftSpheric(1.0);

    iftIGraph *igraph = iftExplicitIGraph(mimg, bin_mask, NULL, A);
    iftIGraphISF_Root(igraph, seeds, alpha, beta, 5);

    iftImage *svoxels_img = iftIGraphLabel(igraph);
    iftImage *final_svoxels_img = iftCreateImageFromImage(svoxels_img);
    
    int sagittal_midplane = img->xsize / 2;

    #pragma omp parallel for
    for (int p = 0; p < svoxels_img->n; p++) {
        if (svoxels_img->val[p]) {
            iftVoxel u = iftGetVoxelCoord(svoxels_img, p);
            
            int disp_x = sagittal_midplane - u.x;
            iftVoxel v = {.x = sagittal_midplane + disp_x, .y = u.y, .z = u.z};

            iftImgVoxelVal(final_svoxels_img, u) = iftImgVoxelVal(svoxels_img, u);
            iftImgVoxelVal(final_svoxels_img, v) = iftImgVoxelVal(svoxels_img, u);
        }
    }
    
    iftDestroyImage(&asym_map);
    iftDestroyIntArray(&grid);
    iftDestroyImage(&seeds);
    iftDestroyMImage(&mimg);
    iftDestroyAdjRel(&A);
    iftDestroyIGraph(&igraph);
    iftDestroyImage(&svoxels_img);

    return final_svoxels_img;
}



iftImage *iftSymmOISF(const iftImage *img, const iftImage *bin_mask, int n_supervoxels, float alpha, float beta,
                      float gamma, const iftImage *normal_asymmap, float thres_factor) {
    iftImage *asym_map = iftBrainAsymMap(img, normal_asymmap);
    
    iftIntArray *grid = iftGridSamplingByBrainAsymmetries(asym_map, bin_mask, 50, 2);
    iftImage *seeds = iftCreateImageFromImage(img);
    iftIntArrayToImage(grid, seeds, 1);
//    iftWriteImageByExt(seeds, "tmp/seeds.nii.gz");
    
    
    iftMImage *mimg = _iftBuildBrainHemisphereMImageOISF(img, asym_map);
    
    iftAdjRel *A = iftSpheric(1.0);
    iftIGraph *igraph = iftImplicitIGraph(mimg, bin_mask, A);
    iftIGraph *igraph_oisf = iftIGraphOISF(igraph, seeds, alpha, beta, gamma, 10);
    
    iftImage *supervoxels_img = iftIGraphLabel(igraph_oisf);
    
    int sagittal_midplane = img->xsize / 2;

#pragma omp parallel for
    for (int p = 0; p < supervoxels_img->n; p++) {
        if (supervoxels_img->val[p]) {
            iftVoxel u = iftGetVoxelCoord(supervoxels_img, p);
            
            int disp_x = sagittal_midplane - u.x;
            iftVoxel v = {.x = sagittal_midplane + disp_x, .y = u.y, .z = u.z};
            
            iftImgVoxelVal(supervoxels_img, v) = iftImgVoxelVal(supervoxels_img, u);
        }
    }
    
    iftDestroyImage(&asym_map);
    iftDestroyIntArray(&grid);
    iftDestroyImage(&seeds);
    iftDestroyMImage(&mimg);
    iftDestroyAdjRel(&A);
    iftDestroyIGraph(&igraph);
    iftDestroyIGraph(&igraph_oisf);
    
    return supervoxels_img;
}


iftImage *iftExtractBrainSide(const iftImage *img, iftBrainSide side) {
    int msp = img->xsize / 2;
    
    iftBoundingBox bb_opposite; // used to reset/set value 0 to the opposite side
    bb_opposite.begin.y = bb_opposite.begin.z = 0;
    bb_opposite.end.y = img->ysize - 1;
    bb_opposite.end.z = img->zsize - 1;
    
    if (side == IFT_RIGHT_BRAIN_SIDE) {
        bb_opposite.begin.x = msp;
        bb_opposite.end.x = img->xsize - 1;
    }
    else {
        bb_opposite.begin.x = 0;
        bb_opposite.end.x = msp;
    }
    
    iftImage *brain_side_img = iftCopyImage(img);
    iftFillBoundingBoxInImage(brain_side_img, bb_opposite, 0);
    
    return brain_side_img;
}


iftDataSet **iftExtractSupervoxelHAAFeats(const iftImage *test_img, const iftImage *test_sym_svoxels,
                                          const iftFileSet *train_set, int n_bins, const iftImage *normal_asym_map,
                                          int *n_svoxels_out) {
    iftImage *svoxels_right_side = iftExtractBrainSide(test_sym_svoxels, IFT_RIGHT_BRAIN_SIDE);
    int norm_val = iftMaximumValue(test_img);
    
    int n_svoxels = iftMaximumValue(svoxels_right_side);
    int n_samples = train_set->n + 1; // training samples + testing sample
    
    iftDataSet **Zarr = iftAlloc(n_svoxels + 1, sizeof(iftDataSet *));
    
    for (int label = 1; label <= n_svoxels; label++)
        Zarr[label] = iftCreateDataSet(n_samples, n_bins);

    #pragma omp parallel for
    for (long s = 0; s < train_set->n; s++) {
        iftImage *train_img = iftReadImageByExt(train_set->files[s]->path);
        iftImage *train_asym_map = iftBrainAsymMap(train_img, NULL);
        iftSupervoxelHAAFeatsInAsymMap(train_asym_map, svoxels_right_side, Zarr, s, n_bins, norm_val, n_svoxels);
        iftDestroyImage(&train_img);
        iftDestroyImage(&train_asym_map);
    }
    
    // testing sample is the last sample (index [n_samples - 1]) in all supervoxel datasets
    iftImage *test_asym_map = iftBrainAsymMap(test_img, normal_asym_map);
    iftSupervoxelHAAFeatsInAsymMap(test_asym_map, svoxels_right_side, Zarr, n_samples - 1, n_bins, norm_val, n_svoxels);

    iftDestroyImage(&test_asym_map);
    iftDestroyImage(&svoxels_right_side);
    
    if (n_svoxels_out)
        *n_svoxels_out = n_svoxels;
    
    return Zarr;
}


iftDataSet **iftExtractSupervoxelBICAsymmFeats(const iftImage *test_img, const iftImage *test_sym_svoxels,
                                               const iftFileSet *train_set, int n_bins_per_channel, const iftImage *normal_asym_map,
                                               int *n_svoxels_out) {
    iftImage *svoxels_right_side = iftExtractBrainSide(test_sym_svoxels, IFT_RIGHT_BRAIN_SIDE);
    int n_svoxels = iftMaximumValue(svoxels_right_side);
    int n_samples = train_set->n + 1; // training samples + testing sample
    int n_feats = 2 * n_bins_per_channel;

    iftDataSet **Zarr = iftAlloc(n_svoxels + 1, sizeof(iftDataSet *));
    for (int label = 1; label <= n_svoxels; label++)
        Zarr[label] = iftCreateDataSet(n_samples, n_feats);

    #pragma omp parallel for
    for (long s = 0; s < train_set->n; s++) {
        printf("[%ld/%ld]\n", s, train_set->n - 1);
        iftImage *train_img = iftReadImageByExt(train_set->files[s]->path);
        iftImage *train_asym_map = iftBrainAsymMap(train_img, NULL);
        iftSupervoxelBICFeatsInAsymMap(train_asym_map, svoxels_right_side, Zarr, s, n_bins_per_channel, n_svoxels);
        iftDestroyImage(&train_img);
        iftDestroyImage(&train_asym_map);
    }

    // testing sample is the last sample (index [n_samples - 1]) in all supervoxel datasets
    iftImage *test_asym_map = iftBrainAsymMap(test_img, normal_asym_map);
    iftSupervoxelBICFeatsInAsymMap(test_asym_map, svoxels_right_side, Zarr, n_samples - 1, n_bins_per_channel, n_svoxels);
    iftDestroyImage(&test_asym_map);

    iftDestroyImage(&svoxels_right_side);

    if (n_svoxels_out)
        *n_svoxels_out = n_svoxels;

    return Zarr;
}


void iftBuildRefDataSupervoxelDataSets(iftDataSet **Zarr, int n_supervoxels, const char *test_img_path,
                                       const char *test_supervoxels_path, const iftFileSet *train_set) {
    #pragma omp parallel for
    for (int label = 1; label <= n_supervoxels; label++) {
        if (Zarr[label] != NULL) {
            iftCSV *ref_data = iftCreateCSV(Zarr[label]->nsamples, 3);
            iftSetCSVHeader(ref_data, "image_path,supervoxels_path,supervoxel_label");
            
            if (Zarr[label]->nsamples != (train_set->n + 1))
                iftError("Number of samples for Supervoxel DataSet [%d] is != from the number of training samples + 1 ==> %d ! %d",
                        "iftBuildRefDataSupervoxelDataSets", label, Zarr[label]->nsamples, train_set->n + 1);

            #pragma omp parallel for
            for (int s = 0; s < train_set->n; s++) {
                strcpy(ref_data->data[s][0], train_set->files[s]->path);
                strcpy(ref_data->data[s][1], test_supervoxels_path);
                sprintf(ref_data->data[s][2], "%d", label);
            }
        
            strcpy(ref_data->data[train_set->n][0], test_img_path);
            strcpy(ref_data->data[train_set->n][1], test_supervoxels_path);
            sprintf(ref_data->data[train_set->n][2], "%d", label);
        
            iftSetRefData(Zarr[label], ref_data, IFT_REF_DATA_CSV);
        
            iftDestroyCSV(&ref_data);
        }
    }
}















