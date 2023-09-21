#include "iftDescriptors.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/io/Stream.h"
#include "ift/imgproc/basic/Histogram.h"
#include "iftClustering.h"
#include "iftDataSet.h"
#include "iftInterpolation.h"


/**
@file
@brief A description is missing here
*/

/**
 * @brief Returns the log2 for a value between [0, 255].
 * 
 * This function is a fast implementation for the log2 function for values in [0, 255].
 * 
 * @param val A value between [0, 255].
 * @return Log2 for the value.
 * 
 * @exception Raise an error if the value is out of [0, 255].
 * 
 * @author Samuel Martins
 * @date Oct 1, 2019
 */
unsigned char iftComputeLog2ForBIC(float val){
    if ((val < 0) || (val >= 256.0)) {
        iftError("Input value %f must be in [0, 255]",
                 "iftComputeLog2ForBIC", val);
    }

    if (iftAlmostZero(val)) return 0;
    else if (val < 1.0) return 1;
    else if (val < 2.0) return 2;
    else if (val < 4.0) return 3;
    else if (val < 8.0) return 4;
    else if (val < 16.0) return 5;
    else if (val < 32.0) return 6;
    else if (val < 64.0) return 7;
    else if (val < 128.0) return 8;
    else return 9;
}


iftImage* iftQuantize(const iftImage *img, int bins_per_band, int normalization_value){
  iftImage *quant = iftCreateImage(img->xsize,img->ysize,img->zsize);
  int p;

  if(iftIsColorImage(img)){
    for(p = 0; p < img->n; p++){
      iftColor color_ycbcr;
      color_ycbcr.val[0] = img->val[p];
      color_ycbcr.val[1] = img->Cb[p];
      color_ycbcr.val[2] = img->Cr[p];

      iftColor color_rgb = iftYCbCrtoRGB(color_ycbcr,normalization_value);

      int r = (color_rgb.val[0]*bins_per_band)/(normalization_value+1);
      int g = (color_rgb.val[1]*bins_per_band)/(normalization_value+1);
      int b = (color_rgb.val[2]*bins_per_band)/(normalization_value+1);

      quant->val[p] = r + g*bins_per_band + b*bins_per_band*bins_per_band;

    }
  }else{
    for(p = 0; p < img->n; p++){
      quant->val[p] = (img->val[p]*bins_per_band)/(normalization_value+1);
    }
  }

  return quant;
}

iftImage *iftQuantizeByClustering1D(iftImage *img, int bins_per_band) {
  iftImage *quant = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftDataSet *Z=NULL;
  iftKnnGraph *graph=NULL;

  iftRandomSeed(IFT_RANDOM_SEED);

  Z = iftImageToDataSet(img);
  iftSetStatus(Z,IFT_TRAIN);

  graph = iftCreateKnnGraphInt1D(Z, 15);

  iftUnsupTrainWithCClusters(graph, bins_per_band);

  quant = iftDataSetClustersToQuantizedImage(graph->Z, true);

  iftDestroyDataSet(&Z);
  iftDestroyKnnGraph(&graph);

  return quant;
}

/* ----------- Public functions --------------------*/



iftFeatures *iftCreateFeatures(int n)
{
  iftFeatures *feat=(iftFeatures *) iftAlloc(1,sizeof(iftFeatures));

  feat->val = iftAllocFloatArray(n);
  feat->n   = n;

  return(feat);
}

void iftDestroyFeatures(iftFeatures **feat)
{
  iftFeatures *aux=*feat;

  if (aux != NULL) {
    iftFree(aux->val);
    iftFree(aux);
    *feat = NULL;
  }

}


iftFeatures *iftReadFeatures(const char *filename) {
    if (filename == NULL)
        iftError("filename is NULL", "iftReadFeatures");

    FILE *fp = fopen(filename, "rb");

    struct stat st;
    size_t file_size = 0;

    if (stat(filename, &st) == 0)
        file_size = st.st_size;
    else iftError("Error when getting the Size of the File \"%s\"", "iftReadFeatures", filename);

    if (file_size % sizeof(float) != 0)
        iftError("Features stored seems not to be floats\n" \
                 "Size of the file \"%s\" is not multiple from float (%lu bytes)", "iftReadFeatures",
                 filename, sizeof(float));

    size_t n = file_size / sizeof(float);

    printf("file_size: %lu bytes\n", file_size);
    iftFeatures *feat = iftCreateFeatures(n);

    if (fread(feat->val, sizeof(float), feat->n, fp) != feat->n)
        iftError("Reading Error: feat->val", "iftReadFeatures");

    fclose(fp);

    return feat;
}


void iftWriteFeatures(const iftFeatures *feats, const char *filename) {
    if (filename == NULL)
        iftError("filename is NULL", "iftReadFeatures");

    FILE *fp = fopen(filename, "wb");

    if (fwrite(feats->val, sizeof(float), feats->n, fp) != feats->n) {
        iftError("Writing Error: feats->val", "iftReadFeatures");
    }

	fclose(fp);
}


iftFeatures *iftJoinFeatures(iftFeatures **array, int n)
{
    int nFeats = 0;
    for(int i = 0; i < n; i++)
        nFeats += array[i]->n;

    int s = 0;
    iftFeatures *feats = iftCreateFeatures(nFeats);
    for(int i = 0; i < n; i++) {
        for(int f = 0; f < array[i]->n; f++) {
            feats->val[s] = array[i]->val[f];
            s++;
        }
    }

    return feats;
}


iftImage *iftLocalBinaryPattern(const iftImage *img, const iftAdjRel *A) {
    iftImage *lbp = iftCreateImageFromImage(img);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++) {
        iftVoxel u = iftGetVoxelCoord(img, p);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(img, v)) {
                int q = iftGetVoxelIndex(img, v);

                if (img->val[q] > img->val[p]) {
                    lbp->val[p] |= 1 << (i - 1);
                }
            }
        }
    }

    return lbp;
}


iftFeatures *iftExtractLBP(iftImage *img, iftAdjRel *A)
{
  iftImage    *lbp=iftLocalBinaryPattern(img,A), *grad;
  iftFeatures *feat;

  grad = iftImageGradientMagnitude(lbp,A);
  feat = iftExtractBrightness(grad);

  iftDestroyImage(&lbp);
  iftDestroyImage(&grad);

  return(feat);
}


iftImageArray *ift3DLBPTOP(const iftImage *img) {
    if (!iftIs3DImage(img))
        iftError("Input image is no 3D.", "ift3DLBPTOP");
    if (iftIsColorImage(img))
        iftError("Input image must be in grayscale.", "ift3DLBPTOP");

    iftAdjRel *A = iftCircular(sqrtf(2));

    iftImageArray *lbps3D = iftCreateImageArray(3);

    iftImage *lbp3D_XY = iftCreateImageFromImage(img);
    iftImage *lbp3D_ZX = iftCreateImageFromImage(img);
    iftImage *lbp3D_YZ = iftCreateImageFromImage(img);


    #pragma omp parallel for
    for (int z = 0; z < img->zsize; z++) {
        iftImage *img2D_XY = iftGetXYSlice(img, z);

        iftImage *lbp2D_XY_slice = iftLocalBinaryPattern(img2D_XY, A);
        iftPutXYSlice(lbp3D_XY, lbp2D_XY_slice, z);

        iftDestroyImage(&img2D_XY);
        iftDestroyImage(&lbp2D_XY_slice);
    }

    #pragma omp parallel for
    for (int y = 0; y < img->ysize; y++) {
        iftImage *img2D_ZX = iftGetZXSlice(img, y);

        iftImage *lbp2D_ZX_slice = iftLocalBinaryPattern(img2D_ZX, A);
        iftPutZXSlice(lbp3D_ZX, lbp2D_ZX_slice, y);

        iftDestroyImage(&img2D_ZX);
        iftDestroyImage(&lbp2D_ZX_slice);
    }

    #pragma omp parallel for
    for (int x = 0; x < img->xsize; x++) {
        iftImage *img2D_YZ = iftGetYZSlice(img, x);

        iftImage *lbp2D_YZ_slice = iftLocalBinaryPattern(img2D_YZ, A);
        iftPutYZSlice(lbp3D_YZ, lbp2D_YZ_slice, x);

        iftDestroyImage(&img2D_YZ);
        iftDestroyImage(&lbp2D_YZ_slice);
    }

    lbps3D->val[0] = lbp3D_XY;
    lbps3D->val[1] = lbp3D_ZX;
    lbps3D->val[2] = lbp3D_YZ;

    iftDestroyAdjRel(&A);

    return lbps3D;
}


iftFeatures *iftExtract3DLBPTOPFeats(const iftImage *img, const iftImage *mask,
                                     int n_bins, bool normalize_histograms) {
    iftImageArray *lbps3D = ift3DLBPTOP(img);

    int max_val = 255; // the function iftVLBP builds the texture number/code with 8 bits (neighbors)
    if (n_bins > (max_val + 1)) {
        iftWarning("Number of Bins %d > max_val + 1: %d\n",
                   "iftExtract3DLBPTOPFeats", n_bins, max_val + 1);
    }

    iftFeatures *feats = iftCreateFeatures(3 * n_bins);

    iftHist *hist_XY = NULL;
    iftHist *hist_ZX = NULL;
    iftHist *hist_YZ = NULL;

    hist_XY = iftCalcGrayImageHist(lbps3D->val[0], mask, n_bins, max_val, normalize_histograms);
    hist_ZX = iftCalcGrayImageHist(lbps3D->val[1], mask, n_bins, max_val, normalize_histograms);
    hist_YZ = iftCalcGrayImageHist(lbps3D->val[2], mask, n_bins, max_val, normalize_histograms);

    #pragma omp parallel for
    for (int b = 0; b < hist_XY->nbins; b++)
        feats->val[b] = hist_XY->val[b];
    int offset = n_bins;

    #pragma omp parallel for
    for (int b = 0; b < hist_ZX->nbins; b++)
        feats->val[offset + b] = hist_ZX->val[b];
    offset = 2 * n_bins;

    #pragma omp parallel for
    for (int b = 0; b < hist_YZ->nbins; b++)
        feats->val[offset + b] = hist_YZ->val[b];

    iftDestroyImageArray(&lbps3D);
    iftDestroyHist(&hist_XY);
    iftDestroyHist(&hist_ZX);
    iftDestroyHist(&hist_YZ);

    return feats;
}


iftMatrix *iftExtract3DLBPTOPFeatsForLabels(const iftImage *img, const iftImage *label_img,
                                            int n_bins, bool normalize_histograms) {
    // puts("ift3DLBPTOP");
    iftImageArray *lbps3D = ift3DLBPTOP(img);
    // iftWriteImageByExt(lbps3D->val[0], "tmp/lbp_axial.nii.gz");
    // iftWriteImageByExt(lbps3D->val[1], "tmp/lbp_coronal.nii.gz");
    // iftWriteImageByExt(lbps3D->val[2], "tmp/lbp_sagittal.nii.gz");

    int max_val = 255; // the function iftVLBP builds the texture number/code with 8 bits (neighbors)
    if (n_bins > (max_val + 1)) {
        iftWarning("Number of Bins %d > max_val + 1: %d\n",
                   "iftExtract3DLBPTOPFeatsForLabels", n_bins, max_val + 1);
    }
    
    int max_label = iftMaximumValue(label_img);
    int n_feats = 3 * n_bins;

    iftMatrix *feats = iftCreateMatrix(n_feats, max_label + 1);

    // puts("Axial = iftCalcGrayImageHistForLabels");
    iftHist **hists_XY = iftCalcGrayImageHistForLabels(lbps3D->val[0], label_img,
                                                       n_bins, max_val, normalize_histograms, NULL);
    // puts("Coronal = iftCalcGrayImageHistForLabels");
    iftHist **hists_ZX = iftCalcGrayImageHistForLabels(lbps3D->val[1], label_img,
                                                       n_bins, max_val, normalize_histograms, NULL);
    // puts("Sagittal = iftCalcGrayImageHistForLabels");
    iftHist **hists_YZ = iftCalcGrayImageHistForLabels(lbps3D->val[2], label_img,
                                                       n_bins, max_val, normalize_histograms, NULL);
    // puts("end");

    int offset_ZX = n_bins;
    int offset_YZ = 2 * n_bins;

    #pragma omp parallel for
    for (int label = 0; label <= max_label; label++) {
        for (int b = 0; b < n_bins; b++) {
            iftMatrixElem(feats, b, label) = hists_XY[label]->val[b];
            iftMatrixElem(feats, b + offset_ZX, label) = hists_ZX[label]->val[b];
            iftMatrixElem(feats, b + offset_YZ, label) = hists_YZ[label]->val[b];
        }
        iftDestroyHist(&hists_XY[label]);
        iftDestroyHist(&hists_ZX[label]);
        iftDestroyHist(&hists_YZ[label]);
    }

    iftDestroyImageArray(&lbps3D);
    iftFree(hists_XY);
    iftFree(hists_ZX);
    iftFree(hists_YZ);

    return feats;
}


iftImage *iftVLBP(const iftImage *img) {
    iftAdjRel *A = iftCreateAdjRel(15);
    A->dx[0] = 0;  A->dy[0] = 0;  A->dz[0] = 0;
    
    A->dx[1] = 0;  A->dy[1] = 0;  A->dz[1] = -1;
    A->dx[2] = 1;  A->dy[2] = 0;  A->dz[2] = -1;
    A->dx[3] = 0;  A->dy[3] = -1; A->dz[3] = -1;
    A->dx[4] = -1; A->dy[4] = 0;  A->dz[4] = -1;
    A->dx[5] = 0;  A->dy[5] = 1;  A->dz[5] = -1;

    A->dx[6] = 1;  A->dy[6] = 0;  A->dz[6] = 0;
    A->dx[7] = 0;  A->dy[7] = -1; A->dz[7] = 0;
    A->dx[8] = -1; A->dy[8] = 0;  A->dz[8] = 0;
    A->dx[9] = 0;  A->dy[9] = 1;  A->dz[9] = 0;

    A->dx[10] = 1;   A->dy[10] = 0;   A->dz[10] = 1;
    A->dx[11] = 0;  A->dy[11] = -1; A->dz[11] = 1;
    A->dx[12] = -1; A->dy[12] = 0;  A->dz[12] = 1;
    A->dx[13] = 0;  A->dy[13] = 1;  A->dz[13] = 1;
    A->dx[14] = 0;  A->dy[14] = 0;  A->dz[14] = 1;

    iftImage *lbp = iftLocalBinaryPattern(img, A);

    iftDestroyAdjRel(&A);

    return lbp;
}


iftFeatures *iftExtractVLBPFeats(const iftImage *img, const iftImage *mask,
                                 int n_bins, bool normalize_hist) {
    if (n_bins <= 0) { n_bins = pow(2, 14); }

    iftImage *lbp = iftVLBP(img);
    int max_val = (int) pow(2, 14); // the function iftVLBP builds the texture number/code with 14 bits (neighbors)
    if (n_bins > (max_val + 1)) {
        iftWarning("Number of Bins %d > max_val + 1: %d\n",
                   "iftExtract3DLBPTOPFeatsForLabels", n_bins, max_val + 1);
    }

    iftHist *hist = iftCalcGrayImageHist(lbp, mask, n_bins, max_val, normalize_hist);

    iftFeatures *feats = iftCreateFeatures(hist->nbins);

    #pragma omp parallel for
    for (int b = 0; b < hist->nbins; b++) {
        feats->val[b] = hist->val[b];
    }

    iftDestroyImage(&lbp);
    iftDestroyHist(&hist);

    return feats;
}


iftMatrix *iftExtractVLBPFeatsForLabels(const iftImage *img, const iftImage *label_img,
                                        int n_bins, bool normalize_hist) {
    if (n_bins <= 0) { n_bins = pow(2, 14); }

    iftImage *lbp = iftVLBP(img);
    iftWriteImageByExt(lbp, "tmp/vlbp.nii.gz");
    int max_val = (int) pow(2, 14); // the function iftVLBP builds the texture number/code with 14 bits (neighbors)

    int max_label = iftMaximumValue(label_img);
    iftMatrix *feats = iftCreateMatrix(n_bins, max_label + 1);

    iftHist **hists = iftCalcGrayImageHistForLabels(lbp, label_img, n_bins, max_val, normalize_hist, NULL);

    #pragma omp parallel for
    for (int label = 0; label <= max_label; label++) {
        for (int b = 0; b < n_bins; b++) {
            iftMatrixElem(feats, b, label) = hists[label]->val[b];
        }
        iftDestroyHist(&hists[label]);
    }
    iftDestroyHist(&hists[0]);
    iftFree(hists);

    iftDestroyImage(&lbp);

    return feats;
}


iftFeatures *iftExtractBrightness(iftImage *img) {
    iftFeatures *feat = iftCreateFeatures(img->n);

    int img_max_val = iftMaximumValue(img);

    for (int p           = 0; p < img->n; p++)
        feat->val[p] = (float) img->val[p] / img_max_val;

    return (feat);
}

iftFeatures *iftExtractGradient(iftImage *img, iftAdjRel *A)
{
  iftImage *grad=iftImageGradientMagnitude(img,A);
  iftFeatures *feat=iftExtractBrightness(grad);

  iftDestroyImage(&grad);

  return(feat);
}

iftFeatures *iftMExtractBIC(iftMImage *mimg, int bins_per_band){
    iftFeatures *features = iftCreateFeatures(bins_per_band*mimg->m*2);

    #pragma omp parallel for
    for(int b = 0; b < mimg->m; b++) {
        iftImage *img = iftMImageToImage(mimg, 255, b);
        iftFeatures *feats = iftExtractBIC(img, NULL, bins_per_band);
        for(int f = 0; f < feats->n; f++)
            features->val[b*feats->n+f] = feats->val[f];
        iftDestroyFeatures(&feats);
        iftDestroyImage(&img);
    }

    return features;
}


iftFeatures *iftExtractBICClustering1D(iftImage *img, int bins_per_band) {
    int bins = bins_per_band, size=bins*bins*bins;
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));

    if(iftIsColorImage(img))
        bins = bins_per_band*bins_per_band*bins_per_band;

    iftImage* quant_img = iftQuantizeByClustering1D(img,bins_per_band);
    iftFeatures *features = iftCreateFeatures(size*2);

    iftAdjRel *A;
    if(iftIs3DImage(img))
        A = iftSpheric(1.0);
    else
        A = iftCircular(1.0);

    int p;
    for(p = 0; p < img->n; p++){
        iftVoxel u = iftGetVoxelCoord(img,p);

        int border = 0;
        int i;
        for(i = 1; i < A->n && !border; i++){
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);
            if(iftValidVoxel(img,v)){
                int q = iftGetVoxelIndex(img,v);

                if(quant_img->val[p] != quant_img->val[q]){
                    border = 1;
                }
            }
        }

        if(border){
            features->val[quant_img->val[p]]++;
        }else{
            features->val[quant_img->val[p] + size]++;
        }
    }

    int i;
    for(i = 0; i < features->n; i++)
      features->val[i] = 48+iftComputeLog2ForBIC(((float)features->val[i]/(float)img->n)*(float)normalization_value);

    iftDestroyImage(&quant_img);
    iftDestroyAdjRel(&A);

    return features;
}


iftFeatures *iftExtractBIC(const iftImage *img, const iftImage *bin_mask_in,
                           int n_bins_per_channel) {
    iftImage *bin_mask = NULL;
    if (bin_mask_in == NULL) {
        bin_mask = iftSelectImageDomain(img->xsize, img->ysize, img->zsize);
    }
    else {
        bin_mask = iftBinarize(bin_mask_in);
    }

    iftMatrix *feats_mat = iftExtractBICForLabels(img, bin_mask, n_bins_per_channel);
    iftFeatures *feats = iftCreateFeatures(feats_mat->ncols);

    #pragma omp parallel for
    for (int col = 0; col < feats_mat->ncols; col++) {
        feats->val[col] = iftMatrixElem(feats_mat, col, 1);
    }

    iftDestroyImage(&bin_mask);
    iftDestroyMatrix(&feats_mat);

    return feats;
}


iftMatrix *iftExtractBICForLabels(const iftImage *img, const iftImage *label_img,
                                  int n_bins_per_channel) {
    int n_bins = n_bins_per_channel;
    if (iftIsColorImage(img))
        n_bins = n_bins_per_channel * n_bins_per_channel * n_bins_per_channel;

    int norm_val = iftNormalizationValue(iftMaximumValue(img));
    iftImage *quant_img = iftQuantize(img, n_bins_per_channel, norm_val);
    iftAdjRel *A = (iftIs3DImage(img)) ? iftSpheric(1.0) : iftCircular(1.0);
    int n_objs = iftMaximumValue(label_img);

    iftMatrix *hist_borders = iftCreateMatrix(n_bins, n_objs + 1);
    iftIntArray *n_border_voxels_arr = iftCreateIntArray(n_objs + 1);

    iftMatrix *hist_interiors = iftCreateMatrix(n_bins, n_objs + 1);
    iftIntArray *n_interior_voxels_arr = iftCreateIntArray(n_objs + 1);

    for (int p = 0; p < quant_img->n; p++) {
        int label = label_img->val[p];
        
        iftVoxel u = iftGetVoxelCoord(quant_img, p);
        bool is_border = false;

        for (int i = 1; i < A->n && !is_border; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftValidVoxel(img, v) &&
                iftImgVoxelVal(label_img, u) == iftImgVoxelVal(label_img, v)) {
                is_border = (iftImgVoxelVal(quant_img, u) != iftImgVoxelVal(quant_img, v));
            }
            else is_border = true;  // u is on image's border or has a neighbor voxel with a different label
        }

        int bin = quant_img->val[p];
        if (is_border) {
            iftMatrixElem(hist_borders, bin, label)++;
            n_border_voxels_arr->val[label]++;
        }
        else {
            iftMatrixElem(hist_interiors, bin, label)++;
            n_interior_voxels_arr->val[label]++;
        }
    }

    iftMatrix *feats = iftCreateMatrix(2 * n_bins, n_objs + 1);

    // for each object, normalize the border and interior histograms to [0, 255],
    // apply log2 to them and merge them into a single feature vector
    #pragma omp parallel for
    for (int label = 0; label <= n_objs; label++) {
        // an object can have only border pixels so we must make sure
        // that the divisor will not be zero for interiors
        int n_border_voxels = iftMax(n_border_voxels_arr->val[label], 1);
        int n_interior_voxels = iftMax(n_interior_voxels_arr->val[label], 1);

        for (int bin = 0; bin < n_bins; bin++) {
            float norm_border_bin_val = (iftMatrixElem(hist_borders, bin, label) / n_border_voxels) * 255;
            iftMatrixElem(feats, bin, label) = iftComputeLog2ForBIC(norm_border_bin_val);

            float norm_interior_bin_val = (iftMatrixElem(hist_interiors, bin, label) / n_interior_voxels) * 255;
            iftMatrixElem(feats, n_bins + bin, label) = iftComputeLog2ForBIC(norm_interior_bin_val);
        }
    }

    iftDestroyImage(&quant_img);
    iftDestroyAdjRel(&A);
    iftDestroyMatrix(&hist_borders);
    iftDestroyIntArray(&n_border_voxels_arr);
    iftDestroyMatrix(&hist_interiors);
    iftDestroyIntArray(&n_interior_voxels_arr);

    return feats;
}


float iftFeatDistL1(iftFeatures *feat1, iftFeatures *feat2)
{
  int i;
  float dist=0.0f;

  if (feat1->n != feat2->n)
      iftError("Features are not compatible", "iftFeaturesDist");

  for (i=0; i < feat1->n; i++)
    dist += fabs(feat2->val[i]-feat1->val[i]);

  return(dist);
}


void iftWriteFeaturesInFile(iftFeatures *features, int truelabel, FILE *fp) {
	int i, distind = 0;

  /* Write features of a sample */
	fwrite(&distind, sizeof(int), 1, fp); /* Pre-computed distances are not being used. */
	fwrite(&truelabel, sizeof(int), 1, fp);

	for (i = 0; i < features->n; i++)
	  fwrite(&(features->val[i]), sizeof(float), 1, fp);
}


iftFeatures *iftMImageToFeatures(iftMImage *img) {
	iftFeatures *features = iftCreateFeatures(img->n * img->m);

//	int p,b;
//	for(b = 0; b < img->m; b++){
//		for(p = 0; p < img->n; p++){
//			features->val[ p + b*img->n] = img->val[p][b];
//		}
//	}
	int idx = 0;
	for (int p = 0; p < img->n; p++)
		for (int b = 0; b < img->m; b++)
			features->val[idx++] = img->val[p][b];

	return features;
}

iftFeatures *iftIntensityProfile(iftImage *img, iftPoint u1, iftPoint u2)
{
  iftFeatures *features=NULL;
  iftPoint     u;
  int          k, n;
  float        Dx=(u2.x-u1.x), Dy=(u2.y-u1.y), Dz=(u2.z-u1.z);
  float        dx=0, dy=0, dz=0;

  iftVoxel v1, v2;
  v1.x = iftRound(u1.x); v1.y = iftRound(u1.y); v1.z = iftRound(u1.z);
  v2.x = iftRound(u2.x); v2.y = iftRound(u2.y); v2.z = iftRound(u2.z);

  if ( (!iftValidVoxel(img,v1)) || (!iftValidVoxel(img,v2)) )
      iftError("Line has end point(s) outside the image domain", "iftIntensityProfile");


  /* DDA - Digital Differential Analyzer */

  if (iftVectorsAreEqual(u1, u2)) {
    n = 1;
  }else{ /* draw line from u1 to u2 */
    if ((fabs(Dx) >= fabs(Dy))&&(fabs(Dx) >= fabs(Dz))) { /* Dx is the maximum projection of
				     vector u1u2 */
      n  = (int)(fabs(Dx)+1);
      dx = iftSign(Dx);
      dy = dx*Dy/Dx;
      dz = dx*Dz/Dx;
    }else{
      if ((fabs(Dy) >= fabs(Dx))&&(fabs(Dy) >= fabs(Dz))) { /* Dy is the maximum projection of
				       vector u1u2 */
	n  = (int)(fabs(Dy)+1);
	dy = iftSign(Dy);
	dx = dy*Dx/Dy;
	dz = dy*Dz/Dy;
      } else { /* Dz is the maximum projection of vector u1u2 */
	n  = (int)(fabs(Dz)+1);
	dz = iftSign(Dz);
	dx = dz*Dx/Dz;
	dy = dz*Dy/Dz;
      }
    }
  }
  u.x = u1.x;  u.y = u1.y;  u.z = u1.z;
  features = iftCreateFeatures(n);

  for (k=0; k < n; k++) {
    features->val[k]=iftImageValueAtPoint(img,u);
    u.x = u.x + dx;
    u.y = u.y + dy;
    u.z = u.z + dz;
  }

  return(features);
}

iftHog* iftCreateHog2D(int cellSize, int blockSize, int blockStride, int nbins) {

    iftHog* hog = (iftHog*) iftAlloc(1, sizeof(iftHog));
    hog->blockSize.x = hog->blockSize.y = blockSize;
    hog->blockSize.z = 0;

    hog->cellSize.x = hog->cellSize.y = cellSize;
    hog->cellSize.z = 0;

    hog->blockStride.x = hog->blockStride.y = 1;
    hog->blockStride.z = 0;

    hog->nbins = nbins;

    return hog;
}

void iftDestroyHog(iftHog** hog) {
    iftFree(*hog);
    *hog=NULL;
}

//compute the 4 closest cell centers
void iftHogComputeClosestCellCenters(const iftHog* hog, const iftVoxel v, iftVoxel centers[4]) {

    float cX = round((float)v.x / hog->cellSize.x);
    float cY = round((float)v.y / hog->cellSize.y);

    centers[0].x = (cX-0.5) * hog->cellSize.x;
    centers[0].y = (cY-0.5) * hog->cellSize.y;
    centers[0].z = 0;

    centers[1].x = (cX + 0.5) * hog->cellSize.x;
    centers[1].y = (cY - 0.5) * hog->cellSize.y;
    centers[1].z = 0;

    centers[2].x = (cX - 0.5) * hog->cellSize.x;
    centers[2].y = (cY + 0.5) * hog->cellSize.y;
    centers[2].z = 0;

    centers[3].x = (cX + 0.5) * hog->cellSize.x;
    centers[3].y = (cY + 0.5) * hog->cellSize.y;
    centers[3].z = 0;
}

//2 closest orientation bins
void iftHogClosestOrientationBins(const iftHog *hog, float orientation, float *bins) {
    float binSize = 360.0/hog->nbins;
    int bin = round(orientation / binSize);
    bins[0] = (bin * binSize) - (binSize / 2);
    bins[1] = (bin * binSize) + (binSize / 2);

    if(bins[0]<0) {
        bins[0] = 360.0 - binSize/2;
    }
    if(bins[1]>360.0) {
        bins[1] = binSize/2;
    }
}

iftFeatures *iftExtractHOG2D(const iftHog *hog, iftImage *img) {

    int xcells = img->xsize/hog->cellSize.x;//number of cells in width
    int ycells = img->ysize/hog->cellSize.y;//number of cells in height

    int ncells = xcells*ycells;

    float **hist = (float**)iftAlloc(ncells, sizeof(float*));
    for (int i = 0; i < ncells; ++i) {
        hist[i] = iftAllocFloatArray(hog->nbins);
    }

    iftAdjRel* A = iftCircular(5.0);
    iftMImage* gradient = iftGradientVector(img, NULL, A);
    iftDestroyAdjRel(&A);

    iftFloatArray* orientation = iftCreateFloatArray(img->n);
    iftFloatArray* magnitude = iftCreateFloatArray(img->n);

    iftVector u;
    u.x = u.y = u.z = 0;
    for (int p = 0; p < magnitude->n; ++p) {
        u.x = gradient->val[p][0]; u.y = gradient->val[p][1];
        magnitude->val[p] = iftVectorMagnitude(u);
        orientation->val[p] = iftVectorXYAngle(u);
    }

    iftDestroyMImage(&gradient);

    int binSize = 360 / hog->nbins;

    float centerOrientation[2];
    iftVoxel centerCell[4];
    iftVoxel v;

    v.z = 0;
    for (v.y = 0; v.y < img->ysize; v.y+=hog->cellSize.y) {
            for (v.x = 0; v.x < img->xsize; v.x+=hog->cellSize.x) {//iterate over cells

            iftVoxel cell;
            cell.z = 0;
            for (cell.y = v.y; cell.y < v.y + hog->cellSize.y; ++cell.y) {
                for (cell.x = v.x; cell.x < v.x + hog->cellSize.x; ++cell.x) {//iterate inside cell

                        if(!iftValidVoxel(img, cell)) continue;

                        int p = iftGetVoxelIndex(img, cell);

                        float mag = magnitude->val[p];

                        if(mag<=0) continue;

                        iftHogComputeClosestCellCenters(hog, cell, centerCell);
                        iftHogClosestOrientationBins(hog, orientation->val[p], centerOrientation);

                        float orientationDist1 = fabsf(centerOrientation[0] - orientation->val[p])/binSize;
                        float orientationDist2 = fabsf(centerOrientation[1] - orientation->val[p])/binSize;

                        float w;
                        for(int i = 0; i < 4; i++)
                        {
                            int cellIndex = (centerCell[i].y / hog->cellSize.y) * xcells + (centerCell[i].x / hog->cellSize.x);
                            centerCell[i].z = 0;
                            if(cellIndex < ncells && cellIndex >= 0){
                                w = mag;
                                w = w * abs(centerCell[i].x - cell.x) / hog->cellSize.x;
                                w = w * abs(centerCell[i].y - cell.y) / hog->cellSize.y; // weight for the cell;

                                hist[cellIndex][(int)(centerOrientation[0]/binSize) ] += w * orientationDist1;//weight for the bin
                                hist[cellIndex][(int)(centerOrientation[1]/binSize) ] += w * orientationDist2;
                            }
                        }
                    }
                }
            }
    }

    iftFeatures* feats = iftCreateFeatures(ncells*hog->nbins);

    //number of times we pass each cell
    iftIntArray *normCount = iftCreateIntArray(ncells);

    //normalize feature vector
    for (int x = 0; x < xcells; x+=hog->blockStride.x) {
        for (int y = 0; y < ycells; y+=hog->blockStride.y) {
            double sum = 0.0;
            for (int x1 = x; x1 < (x + hog->blockSize.x) && x1 < xcells; ++x1) {
                for (int y1 = y; y1 < (y + hog->blockSize.y) && y1 < ycells; ++y1) {
                    int cellIndex = y1 * xcells + x1;
                    normCount->val[cellIndex]++;
                    for (int k = 0; k < hog->nbins; ++k)
                        sum+=hist[cellIndex][k]*hist[cellIndex][k];
                }
            }

            sum = sqrt(sum)+IFT_EPSILON;

            for (int x1 = x; x1 < (x + hog->blockSize.x) && x1 < xcells; ++x1) {
                for (int y1 = y; y1 < (y + hog->blockSize.y) && y1 < ycells; ++y1) {
                    int cellIndex = y1*xcells+x1;
                    for (int k = 0; k < hog->nbins; ++k)
                        feats->val[cellIndex*hog->nbins+k] += hist[cellIndex][k] / sum;
                }
            }
        }
    }

    for (int i = 0; i < ncells; ++i) {
        for (int j = 0; j < hog->nbins; ++j) {
            feats->val[i*hog->nbins+j]/=normCount->val[i];
        }
    }

    for (int i = 0; i < ncells; ++i) {
        iftFree(hist[i]);
    }
    iftFree(hist);

    iftDestroyIntArray(&normCount);
    iftDestroyFloatArray(&magnitude);
    iftDestroyFloatArray(&orientation);

    return feats;
}
