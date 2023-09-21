#include "ift/segm/GrabCut.h"

#include "ift/core/io/Stream.h"


#define CHI_SQ_DF 3
#define CHI_SQ_STDDEV 0.2721655
#define CHI_SQ_AVG 0.9259259

void iftGMMAddSample(iftGMM *gmm, const float *colors, int grp)
{
    double *sums = &gmm->sums[grp * 3];
    sums[0] += colors[0];
    sums[1] += colors[1];
    sums[2] += colors[2];

    double *prods = &gmm->prods[grp * 9];
    prods[0] += colors[0] * colors[0]; // [0][0]
    prods[1] += colors[0] * colors[1]; // [0][1]
    prods[2] += colors[0] * colors[2]; // [0][2]
    prods[3] += colors[1] * colors[0]; // [1][0]
    prods[4] += colors[1] * colors[1]; // [1][1]
    prods[5] += colors[1] * colors[2]; // [1][2]
    prods[6] += colors[2] * colors[0]; // [2][0]
    prods[7] += colors[2] * colors[1]; // [2][1]
    prods[8] += colors[2] * colors[2]; // [2][2]

    gmm->grp_size[grp]++;
    gmm->total_size++;
}


iftGMM *iftCreateGMM(const iftDataSet *Z)
{
    if (Z->nfeats != 3) {
        iftError("iftGMM only supports 3 dimensional datasets", "iftCreateGMM");
    }

    iftGMM *gmm = iftAlloc(1, sizeof *gmm);

    gmm->coefs = iftAlloc(Z->ngroups, sizeof *gmm->coefs);
    gmm->mean = iftAlloc(Z->ngroups * 3, sizeof *gmm->mean);
    gmm->cov = iftAlloc(Z->ngroups * 9, sizeof *gmm->cov);

    gmm->invCov = iftAlloc(Z->ngroups * 9, sizeof *gmm->invCov);
    gmm->covDet = iftAlloc(Z->ngroups, sizeof *gmm->covDet);

    gmm->sums = iftAlloc(Z->ngroups * 3, sizeof *gmm->sums);
    gmm->prods = iftAlloc(Z->ngroups * 9, sizeof *gmm->prods);

    gmm->grp_size = iftAlloc(Z->ngroups, sizeof *gmm->grp_size);
    gmm->total_size = 0;

    gmm->n_grp = Z->ngroups;

    for (int i = 0; i < Z->nsamples; i++) {
        iftGMMAddSample(gmm, Z->sample[i].feat, Z->sample[i].group - 1);
    }

    return gmm;

}


void iftDestroyGMM(iftGMM **gmm)
{
    iftGMM *aux = *gmm;
    if (aux != NULL) {
        iftFree(aux->coefs);
        iftFree(aux->mean);
        iftFree(aux->cov);
        iftFree(aux->invCov);
        iftFree(aux->covDet);
        iftFree(aux->sums);
        iftFree(aux->prods);
        iftFree(aux->grp_size);
        iftFree(aux);
    }
    aux = NULL;
}


void iftGMMBuildModel(iftGMM *gmm)
{
    const double var = 0.01;
    for (int i = 0; i < gmm->n_grp; i++)
    {
        int n = gmm->grp_size[i];
        if (n == 0) {
            gmm->coefs[i] = 0.0;
        } else {
            gmm->coefs[i] = n / ((double) gmm->total_size);

            double *mean = &gmm->mean[i * 3];
            double *sums = &gmm->sums[i * 3];
            double *cov  = &gmm->cov[i * 9];
            double *prods= &gmm->prods[i * 9];
            double *inv  = &gmm->invCov[i * 9];

            mean[0] = sums[0] / n;
            mean[1] = sums[1] / n;
            mean[2] = sums[2] / n;

            cov[0] = prods[0] / n - mean[0] * mean[0]; // [0][0]
            cov[4] = prods[4] / n - mean[1] * mean[1]; // [1][1]
            cov[8] = prods[8] / n - mean[2] * mean[2]; // [2][2]

            cov[1] = cov[3] = prods[1] / n - mean[0] * mean[1]; // [0][1] & [1][0]
            cov[2] = cov[6] = prods[2] / n - mean[0] * mean[2]; // [0][2] & [2][0]
            cov[5] = cov[7] = prods[5] / n - mean[1] * mean[2]; // [1][2] & [2][1]

            double det = cov[0] * (cov[4] * cov[8] - cov[5] * cov[7]) -
                         cov[1] * (cov[3] * cov[8] - cov[5] * cov[6]) +
                         cov[2] * (cov[3] * cov[7] - cov[4] * cov[6]);

            if (det < IFT_EPSILON) // adding variance to avoid singular covariance
            {
                cov[0] += var;
                cov[4] += var;
                cov[8] += var;
                det = cov[0] * (cov[4] * cov[8] - cov[5] * cov[7]) -
                      cov[1] * (cov[3] * cov[8] - cov[5] * cov[6]) +
                      cov[2] * (cov[3] * cov[7] - cov[4] * cov[6]);
            }

            gmm->covDet[i] = det;

            inv[0] =  (cov[4] * cov[8] - cov[5] * cov[7]) / det; // [0][0]
            inv[1] = -(cov[3] * cov[8] - cov[5] * cov[6]) / det; // [0][1]
            inv[2] =  (cov[3] * cov[7] - cov[4] * cov[6]) / det; // [0][2]
            inv[3] = -(cov[1] * cov[8] - cov[2] * cov[7]) / det; // [1][0]
            inv[4] =  (cov[0] * cov[8] - cov[2] * cov[6]) / det; // [1][1]
            inv[5] = -(cov[0] * cov[7] - cov[1] * cov[6]) / det; // [1][2]
            inv[6] =  (cov[1] * cov[5] - cov[2] * cov[4]) / det; // [2][0]
            inv[7] = -(cov[0] * cov[5] - cov[2] * cov[3]) / det; // [2][1]
            inv[8] =  (cov[0] * cov[4] - cov[1] * cov[3]) / det; // [2][2]
        }
    }
}


double iftGMMPixPDF(iftGMM *gmm, const double *colors, int comp)
{
    double res = 0;
    if (gmm->coefs[comp] > 0)
    {
        if (gmm->covDet[comp] < IFT_EPSILON) {
            printf("%lf\n", gmm->covDet[comp]);
            iftError("Computation with singular matrix", "iftGMMProb");
        }

        double *mean = &gmm->mean[3 * comp];
        double *inv = &gmm->invCov[9 * comp];

        double *diff = iftAlloc(3, sizeof *diff);
        diff[0] = colors[0] - mean[0];
        diff[1] = colors[1] - mean[1];
        diff[2] = colors[2] - mean[2];

        double mult = diff[0] * (diff[0] * inv[0] + diff[1] * inv[3] + diff[2] * inv[6]) +
                      diff[1] * (diff[0] * inv[1] + diff[1] * inv[4] + diff[2] * inv[7]) +
                      diff[2] * (diff[0] * inv[2] + diff[1] * inv[5] + diff[2] * inv[8]);
        res = (1.0 / sqrt(gmm->covDet[comp])) * exp(-0.5 * mult);

        iftFree(diff);
    }
    return res;
}


inline static double normalCDF(double value)
{
    return 0.5 * erfc(-value * M_SQRT1_2);
}


double iftGMMPixCDF(iftGMM *gmm, const double *colors, int comp, double jitter)
{
    double res = 0;
    if (gmm->coefs[comp] > 0)
    {
        if (gmm->covDet[comp] < IFT_EPSILON) {
            printf("%lf\n", gmm->covDet[comp]);
            iftError("Computation with singular matrix", "iftGMMProb");
        }

        double *mean = &gmm->mean[3 * comp];
        double *inv = &gmm->invCov[9 * comp];

        double *diff = iftAlloc(3, sizeof *diff);
        diff[0] = colors[0] - mean[0];
        diff[1] = colors[1] - mean[1];
        diff[2] = colors[2] - mean[2];

        double radii = diff[0] * (diff[0] * inv[0] + diff[1] * inv[3] + diff[2] * inv[6]) +
                       diff[1] * (diff[0] * inv[1] + diff[1] * inv[4] + diff[2] * inv[7]) +
                       diff[2] * (diff[0] * inv[2] + diff[1] * inv[5] + diff[2] * inv[8]);

        // https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution#Approximation

        radii += jitter;
        const double approx_rv = pow(radii / CHI_SQ_DF, 1.0 / 3.0);
        res = normalCDF( (approx_rv - CHI_SQ_AVG) / CHI_SQ_STDDEV);

        iftFree(diff);
    }

    return res;
}

int iftGMMGetComp(iftGMM *gmm, const iftMImage *mimg, int p)
{
    int idx = 0;
    double max = -1.0;

    double *cols = iftAlloc(3, sizeof *cols);
    for (int i = 0; i < gmm->n_grp; i++)
    {
        cols[0] = mimg->val[p][0];
        cols[1] = mimg->val[p][1];
        cols[2] = mimg->val[p][2];
        double prob = iftGMMPixPDF(gmm, cols, i);
        if (prob > max) {
            max = prob;
            idx = i;
        }
    }
    iftFree(cols);

    return idx;
}


iftImage *iftGMMAssignComp(iftGMM *gmm_obj, iftGMM *gmm_bkg, const iftMImage *mimg, iftImage *regions)
{
    iftImage *comp = iftCreateImage(regions->xsize, regions->ysize, regions->zsize);
    for (int p = 0; p < regions->n; p++)
    {
        if (regions->val[p] == SURE_OBJ || regions->val[p] == MAYBE_OBJ) {
            comp->val[p] = iftGMMGetComp(gmm_obj, mimg, p);
        } else {
            comp->val[p] = iftGMMGetComp(gmm_bkg, mimg, p);
        }
    }

    return comp;
}


double iftGMMPixProb(iftGMM *gmm, double *colors, int comp)
{
    double res = 0;
    const double jitter = 1e-2;
    if (gmm->coefs[comp] > 0)
        res = iftGMMPixCDF(gmm, colors, comp, jitter) - iftGMMPixCDF(gmm, colors, comp, -jitter);

    if (res != res) // is nan
        res = 0;

    return res;
}


iftFImage *iftGMMGetProb(iftGMM *gmm, const iftMImage *mimg)
{
    iftFImage *prob = iftCreateFImage(mimg->xsize, mimg->ysize, mimg->zsize);

    double *cols = iftAlloc(3, sizeof cols);
    for (int p = 0; p < prob->n; p++) {
        prob->val[p] = 0;
        for (int i = 0; i < gmm->n_grp; i++)
        {
            cols[0] = mimg->val[p][0];
            cols[1] = mimg->val[p][1];
            cols[2] = mimg->val[p][2];
            prob->val[p] += gmm->coefs[i] * iftGMMPixPDF(gmm, cols, i);
        }
    }
    iftFree(cols);

    return prob;
}


void iftGMMsUpdate(iftGMM *gmm_obj, iftGMM *gmm_bkg, const iftMImage *mimg, iftImage *regions, iftImage *comp)
{
    gmm_obj->total_size = 0;
    gmm_bkg->total_size = 0;

    for (int i = 0; i < gmm_obj->n_grp; i++) {
        gmm_obj->grp_size[i] = 0;
        gmm_bkg->grp_size[i] = 0;
    }

    for (int i = 0; i < gmm_obj->n_grp * 3; i++) {
        gmm_obj->sums[i] = 0;
        gmm_bkg->sums[i] = 0;

    }

    for (int i = 0; i < gmm_obj->n_grp * 9; i++) {
        gmm_obj->prods[i] = 0;
        gmm_bkg->prods[i] = 0;
    }

    float *colors = iftAlloc(3, sizeof *colors);
    for (int p = 0; p < mimg->n; p++) {
        colors[0] = mimg->val[p][0];
        colors[1] = mimg->val[p][1];
        colors[2] = mimg->val[p][2];
        if (regions->val[p] == SURE_OBJ || regions->val[p] == MAYBE_OBJ) {
            iftGMMAddSample(gmm_obj, colors, comp->val[p]);
        } else {
            iftGMMAddSample(gmm_bkg, colors, comp->val[p]);
        }
    }
    iftFree(colors);
}


iftImage *iftUncertMapToLabel(const iftImage *map)
{
    iftImage *label = iftCreateImage(map->xsize, map->ysize, map->zsize);

    for (int p = 0; p < map->n; p++) {
        if (map->val[p] == SURE_OBJ || map->val[p] == MAYBE_OBJ) {
            label->val[p] = 1;
        } else {
            label->val[p] = 0;
        }
    }

    return label;
}


iftImage *iftLabelToUncertMap(const iftImage *label, const iftAdjRel *A)
{
    iftImage *map = iftCreateImage(label->xsize, label->ysize, label->zsize);

    for (int p = 0; p < label->n; p++)
    {
        iftVoxel u = iftGetVoxelCoord(label, p);
        if (label->val[p] != 0) {
            map->val[p] = SURE_OBJ;
        } else {
            map->val[p] = SURE_BKG;
        }

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftValidVoxel(label, v)) {
                int q = iftGetVoxelIndex(label, v);
                if (label->val[q] != label->val[p]) {
                    if (label->val[p] != 0) map->val[p] = MAYBE_OBJ;
                    else map->val[p] = MAYBE_BKG;
                }
            }
        }
    }

    return map;
}


iftImage *iftMaybeForeground(const iftImage *label)
{
    iftImage *out = iftCreateImage(label->xsize, label->ysize, label->zsize);

    for (int p = 0; p < out->n; p++) {
        out->val[p] = ((label->val[p] == 0) ? SURE_BKG : MAYBE_OBJ);
    }

    return out;
}


iftImage *iftUpdateUncertMap(const iftImage *map, const iftImage *segm)
{
    iftImage *out = iftCopyImage(map);

    for (int p = 0; p < out->n; p++) {
        if (out->val[p] == MAYBE_OBJ || out->val[p] == MAYBE_BKG) {
            out->val[p] = ((segm->val[p] != 0) ? MAYBE_OBJ : MAYBE_BKG);
        }
    }

    return out;
}


iftImage *iftGrabCut(const iftMImage *mimg, const iftImage *regions_, double beta, int n_iters)
{
    iftDataSet *Z = iftMImageToDataSet(mimg, NULL, 0);

    iftImage *regions = iftCopyImage(regions_);

    for (int i = 0; i < Z->nsamples; i++) {
        if (regions->val[i] == SURE_OBJ || regions->val[i] == MAYBE_OBJ) {
            Z->sample[i].group = 1;
        } else {
            Z->sample[i].group = 0;
        }
    }
    Z->ngroups = 2;

    iftDataSet *Zbkg = iftExtractGroup(Z, 0);
    iftDataSet *Zobj = iftExtractGroup(Z, 1);

    iftDestroyDataSet(&Z);

    iftClusterDataSetByKMeans(Zbkg, 5, 10, 1e-4, 0, 0, 0);
    iftClusterDataSetByKMeans(Zobj, 5, 10, 1e-4, 0, 0, 0);

    iftGMM *gmm_obj = iftCreateGMM(Zobj);
    iftDestroyDataSet(&Zobj);
    iftGMMBuildModel(gmm_obj);

    iftGMM *gmm_bkg = iftCreateGMM(Zbkg);
    iftDestroyDataSet(&Zbkg);
    iftGMMBuildModel(gmm_bkg);

    for (int ite = 0; ite < n_iters; ite++)
    {
        iftImage *comp = iftGMMAssignComp(gmm_obj, gmm_bkg, mimg, regions);
        iftGMMsUpdate(gmm_obj, gmm_bkg, mimg, regions, comp);
        iftDestroyImage(&comp);

        iftGMMBuildModel(gmm_obj);
        iftGMMBuildModel(gmm_bkg);

        iftFImage *prob_obj = iftGMMGetProb(gmm_obj, mimg);
        iftFImage *prob_bkg = iftGMMGetProb(gmm_bkg, mimg);

        iftMaxflowGraph *graph = iftMaxflowGraphFromMImageProb(mimg, regions, prob_obj, prob_bkg, beta);
        iftMaxflow(graph);

        iftImage *segm = iftMMaxflowSegment(graph, mimg);
        iftImage *aux = regions;
        regions = iftUpdateUncertMap(aux, segm);
        iftDestroyImage(&aux);
        iftDestroyImage(&segm);

        iftDestroyFImage(&prob_obj);
        iftDestroyFImage(&prob_bkg);
        iftDestroyMaxflowGraph(&graph);
    }

    iftImage *out = iftUncertMapToLabel(regions);
    iftDestroyImage(&regions);
    iftDestroyGMM(&gmm_obj);
    iftDestroyGMM(&gmm_bkg);

    return out;
}


iftFImage *iftGMMDataSetDist(iftDataSet *train, const iftMImage *mimg, int n_comp)
{
    iftClusterDataSetByKMeans(train, n_comp, 10, 1e-4, 0, 0, 0);
    iftGMM *gmm = iftCreateGMM(train);
    iftGMMBuildModel(gmm);
    iftFImage *prob = iftGMMGetProb(gmm, mimg);
    iftDestroyGMM(&gmm);
    return prob;
}


#undef CHI_SQ_DF
#undef CHI_SQ_STDDEV
#undef CHI_SQ_AVG
