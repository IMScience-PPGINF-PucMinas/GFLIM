#include "iftRadiometric.h"

#include "ift/core/dtypes/GQueue.h"
#include "ift/core/io/Stream.h"
#include "ift/imgproc/basic/Histogram.h"
#include "iftImageMath.h"

/**
@file
@brief A description is missing here
*/

/*------------------------------ Public Functions ------------------------- */

iftImage *iftLinearStretch(iftImage *img, double f1, double f2, double g1, double g2)
{
  iftImage *simg=iftCreateImage(img->xsize,img->ysize,img->zsize);
  double a;

  if (img->Cb != NULL) 
    iftCopyCbCr(img,simg);
  iftCopyVoxelSize(img,simg);

  if (f1 != f2) 
    a = (g2-g1)/(f2-f1);
  else
    a = IFT_INFINITY_DBL;

  #pragma omp parallel
  for (int p=0; p < img->n; p++){
    if (img->val[p] < f1)
      simg->val[p] = (int) iftRound(g1);
    else 
      if (img->val[p] > f2)
        simg->val[p] = (int) iftRound(g2);
      else {
	if (!iftAlmostZero(a - IFT_INFINITY_DBL))
      simg->val[p] = iftRound(a * (img->val[p] - f1) + g1);
	else{
	  simg->val[p] = (int) iftRound(g2);
	}   
      }
  }

  return(simg);
}

iftImage *iftGaussianStretch(iftImage *img, double mean, double stdev, int maxval)
{
  iftImage *simg=iftCreateImage(img->xsize,img->ysize,img->zsize);
  int p;
  float K = 2.0*stdev*stdev;

  if (img->Cb != NULL) 
    iftCopyCbCr(img,simg);
  iftCopyVoxelSize(img,simg);

  for (p=0; p < img->n; p++){
    simg->val[p] = (int)(maxval * exp(-(pow((img->val[p]-mean),2)/K)));
  }
  
  return(simg);
}

iftImage *iftExponenStretch(iftImage *img, double f1, double f2, double g1, double g2)
{
  iftImage *simg=iftCreateImage(img->xsize,img->ysize,img->zsize);
  int p;
  float a, b=exp(-1);

  if (img->Cb != NULL) 
    iftCopyCbCr(img,simg);
  iftCopyVoxelSize(img,simg);

  if (f1 != f2) 
    a = (f2-f1);
  else
    a = IFT_INFINITY_FLT;

  for (p=0; p < img->n; p++){
    if (img->val[p] < f1)
      simg->val[p] = (int)g1;
    else 
      if (img->val[p] > f2)
	simg->val[p] = (int)g2;
      else {
	if (a != IFT_INFINITY_FLT)	  
	  simg->val[p] = (int)((g2*(expf((img->val[p]-f2)/a)-b))+g1+b);
	else{
	  simg->val[p] = (int)g2;
	}   
      }
  }

  return(simg);
}


iftImage *iftNormalizeInRegion(const iftImage *img, const iftImage *region, double minval, double maxval) {
    iftImage *nimg = iftCopyImage(img);

    int img_min_val;
    int img_max_val;
    if (region) {
        iftVerifyImageDomains(img, region, "iftNormalizeInRegion");
        iftMinMaxValueInRegion(img, region, &img_min_val, &img_max_val);
    } else iftMinMaxValues(img, &img_min_val, &img_max_val);


    if (iftIsColorImage(img)){
        if (img_min_val < img_max_val) {
            iftCopyCbCr(img, nimg);

            int in_norm_value   = iftNormalizationValue(img_max_val);
            int out_norm_value  = iftNormalizationValue(maxval);

            for (int p = 0; p < img->n; p++) {
                iftColor rgb = iftGetRGB(img, p, in_norm_value);
                img_max_val = iftMax(img_max_val, iftMax(rgb.val[0], iftMax(rgb.val[1], rgb.val[2])));
            }

            #pragma omp parallel for
            for (int p = 0; p < img->n; p++){
                iftColor rgb = iftGetRGB(img, p, in_norm_value);

                rgb.val[0] = (int) ((maxval - minval) * ((double) rgb.val[0]) /
                                      ((double) img_max_val) + minval);
                rgb.val[1] = (int) ((maxval - minval) * ((double) rgb.val[1]) /
                                     ((double) img_max_val) + minval);
                rgb.val[2] = (int) ((maxval - minval) * ((double) rgb.val[2]) /
                                     ((double) img_max_val) + minval);

                iftSetRGB2(nimg, p, rgb, out_norm_value);
            }
        } else iftError("Cannot normalize empty image", "iftNormalizeInRegion");
    } else {
        if (img_min_val < img_max_val) {
            #pragma omp parallel for
            for (int p = 0; p < img->n; p++) {
                nimg->val[p] = (int) ((maxval - minval) * ((double) img->val[p] - (double) img_min_val) /
                                      ((double) img_max_val - (double) img_min_val) + minval);
            }
        } else iftError("Cannot normalize empty image", "iftNormalize");
    }


    // mask the normalized image only on the region
    // if (region) {
    //     iftImage *aux = nimg;
    //     nimg = iftMask(nimg, region);
    //     iftDestroyImage(&aux);
    // }

    return nimg;
}


iftImage *iftNormalize(const iftImage *img, double minval, double maxval) {
    return iftNormalizeInRegion(img, NULL, minval, maxval);
}


iftImage *iftImageDivisiveNormalization(iftImage *img, iftAdjRel *A, int minValNorm, int maxValNorm)
{
    int normValue = iftNormalizationValue(iftMaximumValue(img));
    iftImage *normImg = NULL;

    if(iftIsColorImage(img))
    {
        normImg = iftCreateColorImage(img->xsize, img->ysize, img->zsize, log2(maxValNorm));

        #pragma omp parallel for
        for (int p = 0; p < img->n; p++)
        {
            float sum = 0.0;
            iftVoxel u = iftGetVoxelCoord(img, p);
            for (int i = 1; i < A->n; i++)
            {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(img, v))
                {
                    int q = iftGetVoxelIndex(img, v);
                    iftColor rgb = iftGetRGB(img, q, normValue);
                    sum += rgb.val[0] * rgb.val[0];
                    sum += rgb.val[1] * rgb.val[1];
                    sum += rgb.val[2] * rgb.val[2];
                }
            }
            sum = sqrtf(sum);
            if (sum > IFT_EPSILON)
            {
                iftColor rgb = iftGetRGB(img, p, normValue);
                rgb.val[0] = (int) ((float)rgb.val[0] / sum * (float)maxValNorm);
                rgb.val[1] = (int) ((float)rgb.val[1] / sum * (float)maxValNorm);
                rgb.val[2] = (int) ((float)rgb.val[2] / sum * (float)maxValNorm);
                iftSetRGB(normImg, p, rgb.val[0], rgb.val[1], rgb.val[2], maxValNorm);
            }
        }
    } else
    {
        normImg = iftCreateImage(img->xsize, img->ysize, img->zsize);

        #pragma omp parallel for
        for (int p = 0; p < img->n; p++)
        {
            float sum = 0.0;
            iftVoxel u = iftGetVoxelCoord(img, p);
            for (int i = 1; i < A->n; i++)
            {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(img, v))
                {
                    int q = iftGetVoxelIndex(img, v);
                    sum += img->val[q] * img->val[q];
                }
            }
            sum = sqrtf(sum);
            if (sum > IFT_EPSILON)
            {
                normImg->val[p] = (int) ((float)img->val[p] / sum * (float)maxValNorm);
            }
        }
    }

    iftImage *normImg1 = iftNormalize(normImg, minValNorm, maxValNorm);
    iftDestroyImage(&normImg);

    return (normImg1);
}


//iftImage *iftNormalize(iftImage *img, double minval, double maxval) {
//  iftImage *nimg = iftCreateImage(img->xsize, img->ysize, img->zsize);
//  int      p;
//
//  int img_min_val = iftMinimumValue(img);
//  int img_max_val = iftMaximumValue(img);
//
//  iftCopyVoxelSize(img, nimg);
//
//  if (iftIsColorImage(img)){
//
//    iftCopyCbCr(img, nimg);
//
//    int img_max_Cb = iftMaximumCb(img);
//    int img_max_Cr = iftMaximumCr(img);
//
//    if (img_min_val < img_max_val) {
//      for (p = 0; p < img->n; p++){
//	nimg->val[p] = (int) ((maxval - minval) * ((double) img->val[p]) /
//                            ((double) img_max_val) + minval);
//	nimg->Cb[p] = (int) ((maxval - minval) * ((double) img->Cb[p]) /
//                            ((double) img_max_Cb) + minval);
//	nimg->Cr[p] = (int) ((maxval - minval) * ((double) img->Cr[p]) /
//			     ((double) img_max_Cr) + minval);
//      }
//    } else {
//      iftError("Cannot normalize empty image", "iftNormalize");
//    }
//  } else {
//
//    if (img_min_val < img_max_val) {
//      for (p = 0; p < img->n; p++){
//	nimg->val[p] = (int) ((maxval - minval) * ((double) img->val[p] - (double) img_min_val) /
//                            ((double) img_max_val - (double) img_min_val) + minval);
//      }
//    } else {
//      iftError("Cannot normalize empty image", "iftNormalize");
//    }
//  }
//
//  return (nimg);
//}

iftImage *iftWindowAndLevel(iftImage *img, int width, int level, int maxval) {
  double   f1, f2, g1, g2;
  iftImage *simg;

  f1 = (double) (level - width / 2.0);
  f2 = (double) (level + width / 2.0);
  int img_min_val = iftMinimumValue(img);
  int img_max_val = iftMaximumValue(img);

  if (f1 < img_min_val) f1 = img_min_val;
  if (f2 > img_max_val) f2 = img_max_val;

  g1   = 0.0;
  g2   = maxval;
  simg = iftLinearStretch(img, f1, f2, g1, g2);
  if (img->Cb != NULL)
    iftCopyCbCr(img, simg);
  iftCopyVoxelSize(img, simg);

  return (simg);
}


iftImage *iftEqualize(const iftImage *img, int max_val) {
  if (max_val <= 0) {
    max_val = iftMaxImageRange(iftImageDepth(img));    
  }

  iftImage *eimg = iftCreateImage(img->xsize,img->ysize, img->zsize);
  iftCopyVoxelSize(img,eimg);
  if (iftIsColorImage(img)) 
    iftCopyCbCr(img, eimg);
  
  iftGQueue *Q = iftCreateGQueue(max_val + 1, img->n, img->val);
  
  for (int p = 0; p < img->n; p++) 
    iftInsertGQueue(&Q, p);
  
  int nelems_per_bucket = (int) ((double) img->n / (double) (max_val+1));
  
  int nelems = 1;
  int i      = 0;
  
  while (!iftEmptyGQueue(Q)) {
    int p = iftRemoveGQueue(Q);
    
    if (nelems <= nelems_per_bucket) {
      eimg->val[p] = i;
      nelems++;
    } else {
      nelems = 1;
      i++;
    }
  }
  iftDestroyGQueue(&Q);
  
  /* Apply a median filter */
  
  iftImage  *filt = iftCreateImageFromImage(img);
  iftAdjRel *A;
  
  if (iftIs3DImage(eimg))
    A = iftSpheric(1.0);
  else
    A = iftCircular(1.0);
  
#pragma omp parallel
  for (int p = 0; p < eimg->n; p++) {
    iftVoxel u = iftGetVoxelCoord(eimg,p);
    int n = 0;
    int val[A->n];
    
    for (int i = 0; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, i);
      
      if (iftValidVoxel(eimg, v)) {
	int q = iftGetVoxelIndex(eimg, v);
	  val[i] = eimg->val[q];
	  n++;
	}
      }
      filt->val[p] = iftSelectTheKthElem(val, n, n/2);
    }

    if (iftIsColorImage(eimg))
      iftCopyCbCr(eimg, filt);

    iftDestroyAdjRel(&A);
    iftDestroyImage(&eimg);
    eimg = filt;
    
    return eimg;
}


iftImage *iftMatchHistogram(const iftImage *img, const iftImage *img_mask,
                            const iftImage *ref, const iftImage *ref_mask) {
    int norm_val = iftNormalizationValue(iftMaximumValue(ref));
    int nbins = norm_val + 1;

    iftHist *hist = iftCalcGrayImageHist(img, img_mask, nbins, norm_val, true);
    iftHist *acc_hist = iftCalcAccHist(hist);
    iftDestroyHist(&hist);

    iftHist *hist_ref = iftCalcGrayImageHist(ref, ref_mask, nbins, norm_val, true);
    iftHist *acc_hist_ref = iftCalcAccHist(hist_ref);
    iftDestroyHist(&hist_ref);

    iftImage *mimg = iftCreateImageFromImage(img);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++) {
        int start  = 0;
        int end    = acc_hist_ref->nbins - 1;
        int pos    = 0; // just to avoid warning in compilation
        bool found = false;

        // it tries to find the bucket in ref img's acc histogram
        // which contains the same value of the img's acc histogram
        // in bucket img->val[p]
        double acc_val = acc_hist->val[img->val[p]];

        // Binary search to find the corresponding accumulated histogram value
        // in the reference image
        while (start <= end) {
            pos = (start + end) / 2;

            double diff = acc_val - acc_hist_ref->val[pos]; 
            if (iftAlmostZero(diff)) {
                found = true;
                break;
            }
            else if (diff < 0.0) {
                end = pos - 1;
            } else {
                start = pos + 1;
            }
        }
        end = iftMax(pos, 0); // The end cannot be set lower than brightness 0        

        if (found) {
            mimg->val[p] = pos;
        }
        // IT'S ALMOST ALWAYS COMING HERE!
        else {
            // since it did not find a bucket in ref img's acc histogram with a value
            // exactly equal to the value in img's acc histogram, it gets the
            // bucket of ref img's hist with the closest value to the img's hist one
            if (fabs(acc_val - acc_hist_ref->val[start]) < fabs(acc_val - acc_hist_ref->val[end])) {
                mimg->val[p] = start;
            }
            else {
                mimg->val[p] = end;
            }
        }
    }
    
    iftDestroyHist(&acc_hist);
    iftDestroyHist(&acc_hist_ref);

    return mimg;
}


iftImage *iftNormalizeWithNoOutliersInRegion(const iftImage *img, const iftImage *region,
                                            int minval, int maxval, float perc) {
    int img_min_val;
    int img_max_val;

    if (region) {
        iftVerifyImageDomains(img, region, "iftNormalizeWithNoOutliersInRegion");
        iftMinMaxValueInRegion(img, region, &img_min_val, &img_max_val);
    }
    else iftMinMaxValues(img, &img_min_val, &img_max_val);

    int nbins             = (img_max_val - img_min_val + 1);
   
    iftHist *hist = iftCalcGrayImageHist(img, region, nbins, img_max_val, true);
    iftHist *acchist = iftCalcAccHist(hist);
    iftDestroyHist(&hist);


    iftImage *nimg = iftCopyImage(img);

    if (iftIsColorImage(img))
        iftCopyCbCr(img, nimg);

    // find the threshold/brightness to cut the outliers
    for (int i = nbins - 1; acchist->val[i] > perc; i--)
        img_max_val = i;

    if (img_min_val < img_max_val) {
        #pragma omp parallel for
        for (int p = 0; p < img->n; p++)
            if (img->val[p] <= img_max_val) {
                nimg->val[p] = (int) ((maxval - minval) * ((double) img->val[p] - (double) img_min_val) /
                                     ((double) img_max_val - (double) img_min_val) + minval);
            } else {
                nimg->val[p] = maxval;
            }
    } else iftError("Cannot normalize empty image", "iftNormalizeWithNoOutliers");
    
    iftDestroyHist(&acchist);

    // mask the normalized image only on the region
    // if (region) {
    //     iftImage *aux = nimg;
    //     nimg = iftMask(nimg, region);
    //     iftDestroyImage(&aux);
    // }

    return nimg;
}

iftImage *iftNormalizeWithNoOutliers(const iftImage *img, int minval, int maxval, float perc) {

    return iftNormalizeWithNoOutliersInRegion(img, NULL, minval, maxval, perc);
}


int     iftRadiometricResolution(iftImage *img)
{
    int maxval = iftMaximumValue(img);
    int bpp = 0;
    while (maxval > 0)
    {
        maxval >>= 1;
        bpp++;
    }

    return bpp;
}


iftFImage *iftFLinearStretch(iftFImage *img, double f1, double f2, double g1, double g2) {
  iftFImage *simg=iftCreateFImage(img->xsize,img->ysize,img->zsize);
  int p;
  double a;
  
  iftFCopyVoxelSize(img,simg);
  
  if (!iftAlmostZero(f1 - f2))
    a = (g2-g1)/(f2-f1);
  else
    a = IFT_INFINITY_FLT;
  
  for (p=0; p < img->n; p++) {
    if (img->val[p] < f1)
      simg->val[p] = g1;
    else
      if (img->val[p] > f2)
	simg->val[p] = g2;
      else {
	if (a != IFT_INFINITY_FLT)
	  simg->val[p] = a*(img->val[p]-f1)+g1;
	else{
	  simg->val[p] = g2;
	}
      }
  }
  
  return(simg);
}

iftImage *iftGaussianStandardization(iftImage *orig, float mean, float stdev)
{
  iftImage *simg = iftCopyImage(orig);
  int    max_val = iftMaxImageRange(iftImageDepth(orig));    

  float mean_aux=0.0, stdev_aux=0.0;

  for (int p=0; p < orig->n; p++) {
      mean_aux += orig->val[p];
  }
  mean_aux /= orig->n;

  for (int p=0; p < orig->n; p++) {
    stdev_aux += (orig->val[p]-mean_aux)*(orig->val[p]-mean_aux);
  }
  stdev_aux = sqrt(stdev_aux/orig->n);
  

  for (int p=0; p < orig->n; p++) {
    if (orig->val[p] > mean_aux){
      simg->val[p] = (int)(mean + sqrt(((orig->val[p] - mean_aux)*
					(orig->val[p] - mean_aux))*stdev/stdev_aux));
      
      if (simg->val[p] > max_val)
	simg->val[p] = max_val;
      
      } else {
      simg->val[p] = (int)(mean - sqrt(((orig->val[p] - mean_aux)*
					(orig->val[p] - mean_aux))*
				       stdev/stdev_aux));
      if (simg->val[p] < 0)
	simg->val[p] = 0;
    }
  }
  
  return(simg);
}

