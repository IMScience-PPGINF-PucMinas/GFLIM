#include "iftImageMath.h"

#include "ift/core/dtypes/Matrix.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/Numerical.h"
#include "iftMatrix.h"

/**
@file
@brief A description is missing here
*/

iftImage *iftAdd(const iftImage *img1, const iftImage *img2)
{
  iftImage *img3=NULL;
  int p;

  iftVerifyImageDomains(img1,img2,"iftAdd");

  if (iftIsColorImage(img1)){
    img3 = iftCreateColorImage(img1->xsize,img1->ysize,img1->zsize, iftImageDepth(img1));
    for (p=0; p < img1->n; p++){
      img3->val[p]=img1->val[p]+img2->val[p];
      img3->Cb[p] =img1->Cb[p];
      img3->Cr[p] =img1->Cr[p];
    }
  } else {
    img3 = iftCreateImage(img1->xsize,img1->ysize,img1->zsize);
    for (p=0; p < img1->n; p++){
      img3->val[p]=img1->val[p]+img2->val[p];
    }
  }
  iftCopyVoxelSize(img1,img3);

  return(img3);
}


void iftAddInPlace(const iftImage *src, const iftImage *dst) {
    iftVerifyImageDomains(src, dst, "iftAddInPlace");

    #pragma omp parallel for
    for (int p = 0; p < dst->n; p++) {
        dst->val[p] += src->val[p];
    }
}


iftImage *iftAddScalar(const iftImage *src, int scalar) {
    iftImage *dst = iftCopyImage(src);
    iftAddScalarInPlace(dst, scalar);

    return dst;
}


void iftAddScalarInPlace(iftImage *dst, int scalar) {
    for (int p = 0; p < dst->n; p++)
        dst->val[p] += scalar;
}


iftImage *iftSub(const iftImage *img1, const iftImage *img2) {
    iftVerifyImageDomains(img1, img2, "iftSub");
    
    iftImage *sub_img = iftCreateImageFromImage(img1);
    
    if (iftIsColorImage(img1)) {
        for (int p = 0; p < img1->n; p++) {
            sub_img->val[p] = img1->val[p] - img2->val[p];
            sub_img->Cb[p] = img1->Cb[p];
            sub_img->Cr[p] = img1->Cr[p];
        }
    } else {
        for (int p = 0; p < img1->n; p++) {
            sub_img->val[p] = img1->val[p] - img2->val[p];
        }
    }
    
    return sub_img;
}


iftImage *iftSubReLU(const iftImage *img1, const iftImage *img2) {
    iftVerifyImageDomains(img1, img2, "iftSubReLU");
    
    iftImage *sub_img = iftCreateImageFromImage(img1);
    
    if (iftIsColorImage(img1)) {
        #pragma omp parallel for
        for (int p = 0; p < img1->n; p++) {
            sub_img->val[p] = iftMax(0, img1->val[p] - img2->val[p]);
            sub_img->Cb[p] = img1->Cb[p];
            sub_img->Cr[p] = img1->Cr[p];
        }
    } else {
        #pragma omp parallel for
        for (int p = 0; p < img1->n; p++) {
            sub_img->val[p] = iftMax(0, img1->val[p] - img2->val[p]);
        }
    }
    
    return sub_img;
}



iftImage *iftAbsSub(const iftImage *img1, const iftImage *img2) {
    iftVerifyImageDomains(img1, img2, "iftAbsSub");
    
    iftImage *abs_sub_img = iftCreateImageFromImage(img1);
    
    if (iftIsColorImage(img1)) {
        #pragma omp parallel for
        for (int p = 0; p < img1->n; p++) {
            abs_sub_img->val[p] = abs(img1->val[p] - img2->val[p]);
            abs_sub_img->Cb[p] = img1->Cb[p];
            abs_sub_img->Cr[p] = img1->Cr[p];
        }
    }
    else {
        #pragma omp parallel for
        for (int p = 0; p < img1->n; p++)
            abs_sub_img->val[p] = abs(img1->val[p] - img2->val[p]);
    }
    
    return abs_sub_img;
}


iftImage *iftIntersec(const iftImage *label_img1, const iftImage *label_img2) {
    iftVerifyImageDomains(label_img1, label_img2, "iftIntersec");
    iftImage *intersec = iftCreateImageFromImage(label_img1);

    #pragma omp parallel for
    for (int p = 0; p < intersec->n; p++) {
        if (label_img1->val[p] == label_img2->val[p])
            intersec->val[p] = label_img1->val[p];
    }

    return intersec;
}


iftImage *iftAbs(const iftImage *img)
{
  iftImage *aimg=NULL;
  int p;

  if (iftIsColorImage(img)){
    aimg = iftCreateColorImage(img->xsize,img->ysize,img->zsize, iftImageDepth(img));
    for (p=0; p < img->n; p++){
      aimg->val[p]=abs(img->val[p]);
      aimg->Cb[p] =img->Cb[p];
      aimg->Cr[p] =img->Cr[p];
    }
  } else {
    aimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
    for (p=0; p < img->n; p++){
      aimg->val[p]=abs(img->val[p]);
    }
  }
  iftCopyVoxelSize(img,aimg);

  return(aimg);
}


iftImage *iftSquareSub(const iftImage *img1, const iftImage *img2) {
    iftVerifyImageDomains(img1, img2, "iftAbsSub");

    iftImage *square_sub_img = iftCreateImageFromImage(img1);

    if (iftIsColorImage(img1)) {
        #pragma omp parallel for
        for (int p = 0; p < img1->n; p++) {
            square_sub_img->val[p] = iftPowerOfTwo((img1->val[p] - img2->val[p]));
            square_sub_img->Cb[p] = img1->Cb[p];
            square_sub_img->Cr[p] = img1->Cr[p];
        }
    }
    else {
        #pragma omp parallel for
        for (int p = 0; p < img1->n; p++)
            square_sub_img->val[p] = iftPowerOfTwo((img1->val[p] - img2->val[p]));
    }

    return square_sub_img;
}


iftImage *iftAnd(const iftImage *img1, const iftImage *img2)
{
  iftImage *img3=NULL;
  int p;

  iftVerifyImageDomains(img1,img2,"iftAnd");

  if (iftIsColorImage(img1)){
    img3 = iftCreateColorImage(img1->xsize,img1->ysize,img1->zsize, iftImageDepth(img1));
    for (p=0; p < img1->n; p++){
      img3->val[p]=iftMin(img1->val[p],img2->val[p]);
      img3->Cb[p] =img1->Cb[p];
      img3->Cr[p] =img1->Cr[p];
    }
  } else {
    img3 = iftCreateImage(img1->xsize,img1->ysize,img1->zsize);
    for (p=0; p < img1->n; p++){
      img3->val[p]=iftMin(img1->val[p],img2->val[p]);
    }
  }
  iftCopyVoxelSize(img1,img3);

  return(img3);
}

iftImage *iftOr(const iftImage *img1, const iftImage *img2)
{
  iftImage *img3=NULL;
  int p;

  iftVerifyImageDomains(img1,img2,"iftOr");

  if (iftIsColorImage(img1)){
    img3 = iftCreateColorImage(img1->xsize,img1->ysize,img1->zsize, iftImageDepth(img1));
    for (p=0; p < img1->n; p++){
      img3->val[p]=iftMax(img1->val[p],img2->val[p]);
      img3->Cb[p] =img1->Cb[p];
      img3->Cr[p] =img1->Cr[p];
    }
  } else {
    img3 = iftCreateImage(img1->xsize,img1->ysize,img1->zsize);
    for (p=0; p < img1->n; p++){
      img3->val[p]=iftMax(img1->val[p],img2->val[p]);
    }
  }
  iftCopyVoxelSize(img1,img3);

  return(img3);
}

iftImage *iftMult(const iftImage *img1, const iftImage *img2)
{
  iftImage *img3=NULL;
  int p;

  iftVerifyImageDomains(img1,img2,"iftMult");

  if (iftIsColorImage(img1)){
    img3 = iftCreateColorImage(img1->xsize,img1->ysize,img1->zsize, iftImageDepth(img1));
    for (p=0; p < img1->n; p++){
      img3->val[p]=img1->val[p]*img2->val[p];
      img3->Cb[p] =img1->Cb[p];
      img3->Cr[p] =img1->Cr[p];
    }
  } else {
    img3 = iftCreateImage(img1->xsize,img1->ysize,img1->zsize);
    for (p=0; p < img1->n; p++){
      img3->val[p]=img1->val[p]*img2->val[p];
    }
  }
  iftCopyVoxelSize(img1,img3);

  return(img3);
}


iftFImage *iftMultByFImage(const iftImage *img, const iftFImage *fimg) {
    iftFImage *mult = iftCreateFImage(img->xsize, img->ysize, img->zsize);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++)
        mult->val[p] = img->val[p] * fimg->val[p];

    puts("oi");
    return mult;
}


iftImage *iftMultByIntScalar(const iftImage *img, int scalar) {
    iftImage *out_img = iftCreateImageFromImage(img);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++)
        out_img->val[p] = (img->val[p] * scalar);

    return out_img;
}


iftFImage *iftMultByScalar(const iftImage *img, float scalar) {
    iftFImage *out_fimg = iftCreateFImage(img->xsize, img->ysize, img->zsize);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++)
        out_fimg->val[p] = (img->val[p] * scalar);

    return out_fimg;
}


iftImage *iftXor(const iftImage *img1, const iftImage *img2) {
    iftVerifyImageDomains(img1, img2, "iftXor");

    iftImage *xor = iftCreateImageFromImage(img1);

    for (int p = 0; p != xor->n; p++)
        xor->val[p] = (img1->val[p] != img2->val[p]);

    return xor;
}




iftImage *iftComplement(const iftImage *img)
{
  iftImage *cimg=NULL;
  int p, maxval;

  if (iftIsColorImage(img)){
    cimg   = iftCreateColorImage(img->xsize,img->ysize,img->zsize, iftImageDepth(img));
    maxval = iftMaximumValue(img);
    for (p=0; p < img->n; p++){
      cimg->val[p]=maxval - img->val[p];
      cimg->Cb[p] =img->Cb[p];
      cimg->Cr[p] =img->Cr[p];
    }
  } else {
    cimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
    maxval = iftMaximumValue(img);
    for (p=0; p < img->n; p++){
      cimg->val[p]=maxval - img->val[p];
    }
  }
  iftCopyVoxelSize(img,cimg);

  return(cimg);
}


iftImage *iftReLUImage(const iftImage *img) {
    iftImage *relu = iftCreateImageFromImage(img);
    iftCopyCbCr(img, relu);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++)
        relu->val[p] = iftMax(0, img->val[p]);

    return relu;
}



iftImage *iftMask(const iftImage *img, const iftImage *mask) {
    iftVerifyImageDomains(img, mask, "iftMask");

    iftCopyImageVoxelFunc copyImageVoxelFunc = iftCopyGrayImageVoxel;
    if (iftIsColorImage(img))
        copyImageVoxelFunc = iftCopyColorImageVoxel;

    iftImage *mimg = iftCreateImageFromImage(img);
    
	#pragma omp parallel for
    for (int p = 0; p < img->n; p++)
        if (mask->val[p] != 0) {
            iftVoxel u = iftGetVoxelCoord(img, p);
            copyImageVoxelFunc(img, u, mimg, u);
        }

    return mimg;
}



iftImage *iftObjectMask(const iftImage *img, const iftImage *mask, int label) {
    iftVerifyImageDomains(img, mask, "iftObjectMask");

    iftCopyImageVoxelFunc copyImageVoxelFunc = iftCopyGrayImageVoxel;
    if (iftIsColorImage(img))
        copyImageVoxelFunc = iftCopyColorImageVoxel;

    iftImage *mimg = iftCreateImageFromImage(img);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++)
        if (mask->val[p] == label) {
            iftVoxel u = iftGetVoxelCoord(img, p);
            copyImageVoxelFunc(img, u, mimg, u);
        }

    return mimg;
}


iftImage *iftBinarize(const iftImage *label_img) {
    iftImage *bin_img = iftCreateImageFromImage(label_img);

    #pragma omp parallel for
    for (int p = 0; p < label_img->n; p++)
        bin_img->val[p] = (label_img->val[p] != 0);

    return bin_img;
}


iftImage *iftAddValue(const iftImage *img, int val)
{
  iftImage *nimg=NULL;
  int p;

  if (val < 0)
    iftWarning("Resulting image may have negative values","iftAddValue");

  if (iftIsColorImage(img)){
    nimg   = iftCreateColorImage(img->xsize,img->ysize,img->zsize, iftImageDepth(img));
    for (p=0; p < img->n; p++){
      nimg->val[p]=img->val[p] + val;
      nimg->Cb[p] =img->Cb[p];
      nimg->Cr[p] =img->Cr[p];
    }
  } else {
    nimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
    for (p=0; p < img->n; p++){
      nimg->val[p]=img->val[p] + val;
    }
  }
  iftCopyVoxelSize(img,nimg);

  return(nimg);
}

iftImage *iftBinaryFrom1To255(const iftImage *bin)
{
  iftImage *nbin=NULL;
  int p;

  nbin = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
    for (p=0; p < bin->n; p++){
      if (bin->val[p]!=0)
	nbin->val[p]=255;
  }
  iftCopyVoxelSize(bin,nbin);

  return(nbin);
}

iftImage *iftBinaryFrom255To1(const iftImage *bin)
{
  iftImage *nbin=NULL;
  int p;

  nbin = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
    for (p=0; p < bin->n; p++){
      if (bin->val[p]!=0)
	nbin->val[p]=1;
  }
  iftCopyVoxelSize(bin,nbin);

  return(nbin);
}

iftImage *iftAddValueInRegion(const iftImage *img, int value, const iftImage *mask)
{
  iftImage *nimg=NULL;
  int p;

  iftVerifyImageDomains(img,mask,"iftAddValueInRegion");

  if (value < 0)
    iftWarning("Resulting image may have negative values","iftAddValueInRegion");

  if (iftIsColorImage(img)){
    nimg   = iftCreateColorImage(img->xsize,img->ysize,img->zsize, iftImageDepth(img));
    for (p=0; p < img->n; p++){
      if (mask->val[p]!=0){
	nimg->val[p]=img->val[p] + value;
	nimg->Cb[p] =img->Cb[p];
	nimg->Cr[p] =img->Cr[p];
      }
    }
  } else {
    nimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
    for (p=0; p < img->n; p++){
      if (mask->val[p]!=0){
	nimg->val[p]=img->val[p] + value;
      }
    }
  }
  iftCopyVoxelSize(img,nimg);

  return(nimg);
}

iftFImage *iftSQRT(const iftImage *img)
{
  iftFImage *nimg=NULL;
  int p;

  if (iftIsColorImage(img))
    iftWarning("The color components are being lost","iftSQRT");

  nimg = iftCreateFImage(img->xsize,img->ysize,img->zsize);
  for (p=0; p < img->n; p++){
    nimg->val[p]=sqrt(img->val[p]);
  }
  
  nimg->dx = img->dx;
  nimg->dy = img->dy;
  nimg->dz = img->dz;

  return(nimg);
}

iftFImage *iftFSQRT(const iftFImage *img1) {
  iftFImage *img2 = iftCreateFImage(img1->xsize, img1->ysize, img1->zsize);
  int       p;
  
  for (p = 0; p < img2->n; p++) {
    img2->val[p] = sqrtf(img1->val[p]);
  }

  iftFCopyVoxelSize(img1, img2);

  return (img2);

}

iftImage* iftImageLinearCombination(const iftImage *img0, const iftImage *img1, double alpha) {
  double p = alpha;
  iftImage *combined = NULL;

  iftVerifyImageDomains(img0, img1, "iftLinearCombination");

  /* Differences in minimum and maximum values should not be a
     problem. If this is a problem, it should be treated outside the
     function. */

  if(iftMinimumValue(img0) != iftMinimumValue(img1))
    iftWarning("Minimum values of input images differ (%d vs %d)!\nNormalize the images first", "iftLinearCombination",
             iftMinimumValue(img0), iftMinimumValue(img1));

  if(iftMaximumValue(img0) != iftMaximumValue(img1))
    iftWarning("Maximum values of input images differ (%d vs %d)!\nNormalize the images first", "iftLinearCombination",
             iftMaximumValue(img0), iftMaximumValue(img1));


  if ((!iftIsColorImage(img0))&&(!iftIsColorImage(img1))){
    combined = iftCreateImage(img0->xsize, img0->ysize, img0->zsize);
#pragma omp parallel for
    for (int i = 0; i < img0->n; i++) {
      combined->val[i] = iftRound((img0->val[i]) * p + (img1->val[i]) * (1 - p));
    }
  } else {
    if ((iftIsColorImage(img0))&&(!iftIsColorImage(img1))){
      combined = iftCreateColorImage(img0->xsize, img0->ysize, img0->zsize, iftImageDepth(img0));
#pragma omp parallel for
      for (int i = 0; i < img0->n; i++) {
	combined->val[i] = iftRound((img0->val[i]) * p + (img1->val[i]) * (1 - p));
	combined->Cb[i] = img0->Cb[i];
	combined->Cr[i] = img0->Cr[i];
      }
    } else {
      if ((!iftIsColorImage(img0))&&(iftIsColorImage(img1))){
	combined = iftCreateColorImage(img0->xsize, img0->ysize, img0->zsize, iftImageDepth(img0));
#pragma omp parallel for
	for (int i = 0; i < img0->n; i++) {
	  combined->val[i] = iftRound((img0->val[i]) * p + (img1->val[i]) * (1 - p));
	  combined->Cb[i] = img1->Cb[i];
	  combined->Cr[i] = img1->Cr[i];
	}
      } else { /* both are color images */
	combined = iftCreateColorImage(img0->xsize, img0->ysize, img0->zsize, iftImageDepth(img0));
#pragma omp parallel for
	for (int i = 0; i < img0->n; i++) {
	  combined->val[i] = iftRound((img0->val[i]) * p + (img1->val[i]) * (1 - p));
	  combined->Cb[i]  = iftRound((img0->Cb[i]) * p + (img1->Cb[i]) * (1 - p));
	  combined->Cr[i]  = iftRound((img0->Cr[i]) * p + (img1->Cr[i]) * (1 - p));
	}
      }
    }
  }

  iftCopyVoxelSize(img0, combined);

  return(combined);
}


float iftMeanValue(const iftImage *img, iftImage *mask)
{

    char   no_mask=0;
    double mean=0.0, n=0;

    if (mask == NULL) {
        mask = iftCreateImage(img->xsize,img->ysize,img->zsize);
        iftSetImage(mask,1);
        no_mask = 1;
    } else {
        iftVerifyImageDomains(img,mask,"iftMeanValue");
    }

    for (int p = 0; p < img->n; p++) {
        if (mask->val[p]!=0){
            mean += img->val[p];
            n++;
        }
    }

    if (no_mask)
        iftDestroyImage(&mask);

    if (n != 0)
        return (mean / n);
    else
        return(0);
}


int iftModeValue(const iftImage *img, iftImage *mask)
{
    char  no_mask=0;
    int   mode, maxval, *hist = iftAllocIntArray((maxval=iftMaximumValue(img))+1);

    if (mask == NULL) {
        mask = iftCreateImage(img->xsize,img->ysize,img->zsize);
        iftSetImage(mask,1);
        no_mask = 1;
    } else {
        iftVerifyImageDomains(img,mask,"iftMeanValue");
    }

    for (int p = 0; p < img->n; p++) {
        if (mask->val[p]!=0){
            hist[img->val[p]]++;
        }
    }

    mode = 0;
    for (int i=1; i <= maxval; i++) {
        if (hist[i] > hist[mode]){
            mode = i;
        }
    }

    iftFree(hist);

    if (no_mask)
        iftDestroyImage(&mask);

    return(mode);
}


iftFloatArray *iftMeanValueInRegions(const iftImage *img, const iftImage *label_img) {
    iftVerifyImages(img, label_img, "iftMeanGrayImageLabels");

    int n_objects = iftMaximumValue(label_img);
    iftFloatArray *means = iftCreateFloatArray(n_objects + 1);
    iftIntArray *n_voxels_arr = iftCreateIntArray(n_objects + 1);

    for (int p = 0; p < img->n; p++) {
        int label = label_img->val[p];
        means->val[label] += img->val[p];
        n_voxels_arr->val[label]++;
    }

    for (int label = 1; label <= n_objects; label++) {
        if (n_voxels_arr->val[label] > 0)
            means->val[label] /= n_voxels_arr->val[label];
    }

    iftDestroyIntArray(&n_voxels_arr);

    return means;
}


iftMatrix *iftMeanStdevValuesInRegions(const iftImage *img, const iftImage *label_img, float norm_val) {
    iftVerifyImages(img, label_img, "iftMeanGrayImageLabels");
    iftFImage *img_norm = iftMultByScalar(img, 1.0 / norm_val);

    int n_objects = iftMaximumValue(label_img);

    // cols: [0] mean, [1] stdev
    iftMatrix *feats = iftCreateMatrix(2, n_objects + 1);
    iftIntArray *n_voxels_arr = iftCreateIntArray(n_objects + 1);

    // computes the mean
    for (int p = 0; p < img_norm->n; p++) {
        int label = label_img->val[p];
        iftMatrixElem(feats, 0, label) += img_norm->val[p];
        n_voxels_arr->val[label]++;
    }

    #pragma omp parallel for
    for (int label = 1; label <= n_objects; label++) {
        if (n_voxels_arr->val[label] > 0)
            iftMatrixElem(feats, 0, label) /= n_voxels_arr->val[label];
    }


    // computes the stdev
    for (int p = 0; p < img_norm->n; p++) {
        int label = label_img->val[p];
        float mean = iftMatrixElem(feats, 0, label);
        iftMatrixElem(feats, 1, label) += iftPowerOfTwo(img_norm->val[p] - mean);
    }

    #pragma omp parallel for
    for (int label = 1; label <= n_objects; label++) {
        if (n_voxels_arr->val[label] > 0)
            iftMatrixElem(feats, 1, label) = sqrtf(iftMatrixElem(feats, 1, label) / n_voxels_arr->val[label]);
    }

    iftDestroyFImage(&img_norm);
    iftDestroyIntArray(&n_voxels_arr);

    return feats;
}


iftImage *iftRoundFImage(const iftFImage *fimg) {
    iftImage *img = iftCreateImage(fimg->xsize, fimg->ysize, fimg->zsize);
    iftCopyVoxelSize(fimg, img);

    #pragma omp parallel for
    for (int p = 0; p < fimg->n; p++)
        img->val[p] = iftRound(fimg->val[p]);
    
    return img;
}

