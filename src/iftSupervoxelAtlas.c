//
// Created by azaelmsousa on 12/11/20.
//

#include "iftMatrix.h"
#include "ift/core/io/File.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Stream.h"
#include "iftSupervoxelAtlas.h"
#include "iftImage.h"

iftMatrix *iftStandardizeTestMatrix(iftMatrix *M, char axis, iftMatrix *mean, iftMatrix *stdev)
{
    iftMatrix *standardized = iftCreateMatrix(M->ncols, M->nrows);

    if (axis == 'r') {

        for (int r = 0; r < M->nrows; r++) {
            for (int c = 0; c < M->ncols; c++) {
                iftMatrixElem(standardized,c,r) = (iftMatrixElem(M,c,r)-iftMatrixElem(mean,0,c))/iftMatrixElem(stdev,0,c);
            }
        }

    } else if (axis == 'c') {

        for (int r = 0; r < M->nrows; r++) {
            for (int c = 0; c < M->ncols; c++) {
                iftMatrixElem(standardized,c,r) = (iftMatrixElem(M,c,r)-iftMatrixElem(mean,c,0))/iftMatrixElem(stdev,c,0);
            }
        }

    } else
        iftError("Axis unavailable. Please, select either c or r as axis.","iftStandardizeTestMatrix");

    return standardized;

}

iftMatrix *iftStandardizeTrainMatrix(iftMatrix *M, char axis, iftMatrix **mean, iftMatrix **stdev)
{
    iftMatrix *standardized = iftCreateMatrix(M->ncols, M->nrows);
    iftDblArray *avg = NULL, *std = NULL;

    if (axis == 'r') {
        avg = iftComputeMeanVector(M,'r');
        std = iftCreateDblArray(M->ncols);

        for (int c = 0; c < M->ncols; c++){
            for (int r = 0; r < M->nrows; r++){
                std->val[c] += pow(iftMatrixElem(M,c,r)-avg->val[c],2);
            }
            std->val[c] /= M->nrows;
            std->val[c] = sqrtf(std->val[c]);
        }

        for (int r = 0; r < M->nrows; r++) {
            for (int c = 0; c < M->ncols; c++) {
                iftMatrixElem(standardized,c,r) = (iftMatrixElem(M,c,r)-avg->val[c])/std->val[c];
            }
        }

        *mean = iftCreateMatrix(1,avg->n);
        *stdev = iftCreateMatrix(1,std->n);
        for (int j = 0; j < avg->n; j++){
            iftMatrixElem(*mean,0,j) = avg->val[j];
            iftMatrixElem(*stdev,0,j) = std->val[j];
        }

    } else if (axis == 'c') {
        avg = iftComputeMeanVector(M,'c');
        std = iftCreateDblArray(M->nrows);

        for (int r = 0; r < M->nrows; r++){
            for (int c = 0; c < M->ncols; c++){
                std->val[r] += pow(iftMatrixElem(M,c,r)-avg->val[r],2);
            }
            std->val[r] /= M->ncols;
            std->val[r] = sqrtf(std->val[r]);
        }

        for (int r = 0; r < M->nrows; r++) {
            for (int c = 0; c < M->ncols; c++) {
                iftMatrixElem(standardized,c,r) = (iftMatrixElem(M,c,r)-avg->val[r])/std->val[r];
            }
        }

        *mean = iftCreateMatrix(avg->n,1);
        *stdev = iftCreateMatrix(std->n,1);
        for (int j = 0; j < avg->n; j++){
            iftMatrixElem(*mean,j,0) = avg->val[j];
            iftMatrixElem(*stdev,j,0) = std->val[j];
        }

    } else
        iftError("Axis unavailable. Please, select either c or r as axis.","iftStandardizeTrainMatrix");

    iftDestroyDblArray(&avg);
    iftDestroyDblArray(&std);

    return standardized;
}

iftMImage *iftAlignSupervoxelsHistogram(iftImage *svoxel, iftMImage *mimg, iftImage *mimg_label)
{
    iftMImage *aligned = iftCreateMImage(mimg->xsize,mimg->ysize,mimg->zsize,mimg->m);
    bool is_label_null = false;
    int n_svoxels = iftMaximumValue(svoxel);

    if (mimg_label == NULL) {
        mimg_label = iftThreshold(svoxel,1,n_svoxels,1);
        is_label_null = true;
    }

    iftMatrix *max = iftCreateMatrix(mimg->m,n_svoxels+1);
    iftMatrix *min = iftCreateMatrix(mimg->m,n_svoxels+1);
    iftSetMatrix(min,IFT_INFINITY_DBL);
    iftSetMatrix(max,IFT_INFINITY_DBL_NEG);

    // selecting min and max of every supervoxel
    for (int p = 0; p < svoxel->n; p++){
        if ((svoxel->val[p] > 0) && (mimg_label->val[p] > 0)){
            int svoxel_id = svoxel->val[p];
            for (int m = 0; m < mimg->m; m++){
                if (mimg->val[m][p] > iftMatrixElem(max,m,svoxel_id))
                    iftMatrixElem(max,m,svoxel_id) = mimg->val[m][p];
                if (mimg->val[m][p] < iftMatrixElem(min,m,svoxel_id))
                    iftMatrixElem(min,m,svoxel_id) = mimg->val[m][p];
            }
        }
    }

    // aligning svoxels histograms
    for (int p = 0; p < svoxel->n; p++){
        if ((svoxel->val[p] > 0) && (mimg_label->val[p] > 0)){
            int svoxel_id = svoxel->val[p];
            for (int m = 0; m < mimg->m; m++){
                aligned->val[m][p] = mimg->val[m][p] - (iftMatrixElem(max,m,svoxel_id) - iftMatrixElem(min,m,svoxel_id))/2;
            }
        }
    }
    iftDestroyMatrix(&max);
    iftDestroyMatrix(&min);

    if (is_label_null)
        iftDestroyImage(&mimg_label);

    return aligned;
}

iftHist *iftComputeSupervoxelHistogram(iftSupervoxelAtlas *atlas, int band, int svoxel_id)
{
    iftHist *h = NULL;

    iftImage *img = iftMImageToImage(atlas->mimg,4095,band);
    iftImage *mask = iftExtractObject(atlas->svoxel,svoxel_id);

    h = iftCalcGrayImageHist(img,mask,4095,iftMaximumValueInMask(img,mask),true);

    iftDestroyImage(&img);
    iftDestroyImage(&mask);

    return h;
}

iftSupervoxelAtlas *iftTrainSupervoxelAtlas(iftMImage *mimg, iftImage *svoxel)
{
    if (iftMaximumValue(svoxel) == 0)
        iftError("Supervoxel image appears to be empty.","iftCreateSupervoxelAtlas");

    iftSupervoxelAtlas *A = NULL;

    A = (iftSupervoxelAtlas *)calloc(1, sizeof(iftSupervoxelAtlas));

    A->svoxel = iftCopyImage(svoxel);
    A->mimg = iftCopyMImage(mimg);
    A->n_svoxels = iftMaximumValue(svoxel);
    A->covariance = (iftMatrix **)calloc(A->n_svoxels+1,sizeof(iftMatrix *));
    A->mean = (iftMatrix **)calloc(A->n_svoxels+1,sizeof(iftMatrix *));
    A->stdev = (iftMatrix **)calloc(A->n_svoxels+1,sizeof(iftMatrix *));

    iftIntArray *voxels_per_svoxel = iftCreateIntArray(A->n_svoxels+1);

    //Counting voxels per svoxel
    for (int p = 0; p < A->svoxel->n; p++){
        if (A->svoxel->val[p])
            voxels_per_svoxel->val[A->svoxel->val[p]]++;
    }

    //Turning every svoxel into a matrix data structure
    for (int i = 1; i < voxels_per_svoxel->n; i++){
        iftMatrix *data = iftCreateMatrix(mimg->m,voxels_per_svoxel->val[i]);
        int pos = 0;
        for (int p = 0; p < A->svoxel->n; p++){
            if (A->svoxel->val[p] == i){
                for (int m = 0; m < mimg->m; m++){
                    iftMatrixElem(data,m,pos) = mimg->val[p][m];
                }
                pos++;
            }
        }

        iftMatrix *standardized = iftStandardizeTrainMatrix(data,'r',&A->mean[i],&A->stdev[i]);
        iftDestroyMatrix(&data);

        A->covariance[i] = iftCovarianceMatrix(standardized);
        iftDestroyMatrix(&standardized);
    }
    iftDestroyIntArray(&voxels_per_svoxel);

    return A;
}

void iftDestroySupervoxelAtlas(iftSupervoxelAtlas **atlas)
{
    if (*atlas == NULL)
        return;

    for (int svoxel_id = 1; svoxel_id <= (*atlas)->n_svoxels; svoxel_id++) {
        iftDestroyMatrix(&(*atlas)->covariance[svoxel_id]);
        iftDestroyMatrix(&(*atlas)->mean[svoxel_id]);
        iftDestroyMatrix(&(*atlas)->stdev[svoxel_id]);
    }
    iftFree((*atlas)->covariance);
    iftFree((*atlas)->mean);

    iftDestroyImage(&(*atlas)->svoxel);
    iftDestroyMImage(&(*atlas)->mimg);
    iftFree(*atlas);
    *atlas = NULL;
}

void iftWriteSupervoxelAtlas(const iftSupervoxelAtlas *atlas, const char *filename, ...)
{
    if (!iftCompareStrings(iftFileExt(filename),".zip"))
    {
        iftError("File extension not supported (%s), please use a .zip file\n","iftWriteSupervoxelAtlas",iftFileExt(filename));
    }

    char filename_covariance[50]; //enough space to fit "tempdir_XXXXXX//covariance_%03.npy"
    char filename_mean[50]; //enough space to fit "tempdir_XXXXXX//mean_%03.npy"
    char filename_stdev[50]; //enough space to fit  "tempdir_XXXXXX//stdev_%03.npy"

    char *tmp_dir = iftMakeTempDir("tempdir_",NULL,NULL);

    for (int svoxel_id = 1; svoxel_id <= atlas->n_svoxels; svoxel_id++) {
        sprintf(filename_mean,"%s/mean_%03d.npy", tmp_dir, svoxel_id);
        sprintf(filename_stdev,"%s/stdev_%03d.npy", tmp_dir, svoxel_id);
        sprintf(filename_covariance,"%s/covariance_%03d.npy", tmp_dir, svoxel_id);

        iftWriteMatrix(atlas->mean[svoxel_id],filename_mean);
        iftWriteMatrix(atlas->stdev[svoxel_id],filename_stdev);
        iftWriteMatrix(atlas->covariance[svoxel_id],filename_covariance);
    }

    char filename_img[50];
    if (iftIs3DImage(atlas->svoxel))
        sprintf(filename_img,"%s/svoxel_model.nii.gz",tmp_dir);
    else
        sprintf(filename_img,"%s/svoxel_model.png",tmp_dir);
    iftWriteImageByExt(atlas->svoxel,filename_img);

    sprintf(filename_img,"%s/feats.mimg",tmp_dir);
    iftWriteMImage(atlas->mimg,filename_img);

    iftZipDirContent(tmp_dir,filename);
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
}

iftSupervoxelAtlas *iftReadSupervoxelAtlas(const char *filename)
{
    iftSupervoxelAtlas *atlas = NULL;
    atlas = (iftSupervoxelAtlas *)calloc(1, sizeof(iftSupervoxelAtlas));

    char *tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);
    iftUnzipFile(filename, tmp_dir);

    char filename_img[50];
    char filename_covariance[50];
    char filename_mean[50];
    char filename_stdev[50];

    sprintf(filename_img,"%s/svoxel_model.nii.gz",tmp_dir);
    if (!iftFileExists(filename_img))
        sprintf(filename_img,"%s/svoxel_model.png",tmp_dir);
    atlas->svoxel = iftReadImageByExt(filename_img);
    atlas->n_svoxels = iftMaximumValue(atlas->svoxel);

    sprintf(filename_img,"%s/feats.mimg",tmp_dir);
    atlas->mimg = iftReadMImage(filename_img);

    atlas->covariance = (iftMatrix **)calloc(atlas->n_svoxels+1,sizeof(iftMatrix *));
    atlas->mean = (iftMatrix **)calloc(atlas->n_svoxels+1,sizeof(iftMatrix *));
    atlas->stdev = (iftMatrix **)calloc(atlas->n_svoxels+1,sizeof(iftMatrix *));
    for (int svoxel_id = 1; svoxel_id <= atlas->n_svoxels; svoxel_id++) {
        sprintf(filename_mean,"%s/mean_%03d.npy", tmp_dir, svoxel_id);
        sprintf(filename_stdev,"%s/stdev_%03d.npy", tmp_dir, svoxel_id);
        sprintf(filename_covariance,"%s/covariance_%03d.npy", tmp_dir, svoxel_id);

        atlas->mean[svoxel_id] = iftReadMatrix(filename_mean);
        atlas->stdev[svoxel_id] = iftReadMatrix(filename_stdev);
        atlas->covariance[svoxel_id] = iftReadMatrix(filename_covariance);
    }

    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);

    return atlas;
}

iftImage *iftSegmentationBySupervoxelAtlas(iftMImage *test, iftSupervoxelAtlas *atlas, double bottom_threshold, double top_threshold)
{
    return NULL;
}

iftImage *iftProbMapBySupervoxelAtlas(iftMImage *test, iftImage *test_label, iftSupervoxelAtlas *atlas)
{
    if (test->m != atlas->covariance[1]->nrows){
        iftError("The number of channels on the test mimage (%d) must be the same as the atlas (%d)","iftProbMapBySupervoxelAtlas",test->m,atlas->covariance[0]->nrows);
    }

    iftImage *out = iftCreateImageFromImage(atlas->svoxel);

    for (int p = 0; p < atlas->svoxel->n; p++){

        if(atlas->svoxel->val[p] > 0) {

            //puts("------------");

            iftMatrix *data = iftCreateMatrix(test->m, 1);
            for (int m = 0; m < test->m; m++) {
                iftMatrixElem(data, m, 0) = test->val[p][m];
            }

            /*puts("\nSample bg");
            iftPrintMatrix(data);*/

            iftMatrix *standardized = iftStandardizeTestMatrix(data, 'r', atlas->mean[atlas->svoxel->val[p]],
                                                               atlas->stdev[atlas->svoxel->val[p]]);
            /*puts("\nMean bg");
            iftPrintMatrix(atlas->mean[atlas->svoxel->val[p]]);
            puts("\nstdev bg");
            iftPrintMatrix(atlas->stdev[atlas->svoxel->val[p]]);
            puts("\nStandardized bg");
            iftPrintMatrix(standardized);*/

            iftDestroyMatrix(&data);

            iftMatrix *mean = iftCreateMatrix(test->m, 1);
            double prob = iftMultivariateGaussPDF(standardized, atlas->covariance[atlas->svoxel->val[p]],mean);
            //printf("prob: %lf\n", prob);
            iftDestroyMatrix(&mean);
            iftDestroyMatrix(&standardized);

            out->val[p] = (int) 100 * (1 - prob);
        }
    }

    return out;
}

double iftMultivariateGaussPDF(iftMatrix *sample, iftMatrix *covariance, iftMatrix *mean)
{
    double normalization_factor = 1/(pow(2*PI,((double)covariance->nrows)/2)*sqrt(iftMatrixDeterminant(covariance)));

    iftMatrix *sample_minus_mean = iftSubtractMatrices(sample,mean);
    iftMatrix *sample_minus_meanT = iftTransposeMatrix(sample_minus_mean);
    iftMatrix *covariance_inv = iftInvertMatrix(covariance);

    /*puts("\nnormalization factor");
    printf("%lf\n",normalization_factor);
    puts("\nsample");
    iftPrintMatrix(sample);
    puts("\nsample_minus_mean");
    iftPrintMatrix(sample_minus_mean);
    puts("\ncovariance");
    iftPrintMatrix(covariance);
    puts("\ncovariance_inv");
    iftPrintMatrix(covariance_inv);*/

    iftMatrix *m1 = iftMultMatrices(sample_minus_meanT,covariance_inv);
    iftMatrix *m2 = iftMultMatrices(m1,sample_minus_mean);
    iftDestroyMatrix(&sample_minus_mean);
    iftDestroyMatrix(&sample_minus_meanT);
    iftDestroyMatrix(&covariance_inv);

    /*puts("m1");
    iftPrintMatrix(m1);
    puts("m2");
    iftPrintMatrix(m2);*/

    double e = exp(-0.5*iftMatrixElem(m2,0,0));
    iftDestroyMatrix(&m1);
    iftDestroyMatrix(&m2);

    /*puts("e");
    printf("%lf\n",e);*/

    return normalization_factor*e;
}

