#include "iftMImage.h"

#include "ift/metriclearn/MetricLearnCore.h"
#include "ift/metriclearn/NeighCompAnalysis.h"
#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"


/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Template for a function that copies the values from the image <img> at pixel <q> to the
 *        MImage <mimg> at pixel p starting at the band <initial_band>.
 * 
 * @param mimg                Target Multi-Band Image.
 * @param p                   Target pixel at the MImage.
 * @param img                 Source Image.
 * @param q                   Source pixel at the Source Image.
 * @param initial_band        Initial band where the values from the source image will copied at the MImage.
 * @param normalization_value Normalization value for the Copy. Use 1.0 for non-normalization.
 * 
 * @author Samuka Martins
 * @date July 21, 2018
 */
typedef void (*iftCopyVoxelToMImageFunc)(iftMImage *mimg, int p, const iftImage *img, int q, int initial_band,
                                         int normalization_value);



/**
 * @brief Function to copy a Gray Normalized value from an image <img> at pixel <q> to a MImage <mimg>
 *        at pixel <p> starting at band <initial_band>.
 * 
 * @param mimg                Target Multi-Band Image.
 * @param p                   Target pixel at the MImage.
 * @param img                 Source Image.
notebooks/deep-dynamic-unet.ipynb * @param q                   Source pixel at the Source Image.
 * @param initial_band        Initial band where the values from the source image will copied at the MImage.
 * @param normalization_value Normalization value for the Copy. Use 1.0 for non-normalization.
 *
 * @author Samuka Martins
 * @date July 21, 2018
 */
void iftCopyVoxelToMImageGRAYNorm(iftMImage *mimg, int p, const iftImage *img, int q, int initial_band,
                                  int normalization_value) {
      mimg->val[p][initial_band] = ((float) img->val[q]) / (float) normalization_value;
}


/**
 * @brief Function to copy a LAB2Norm Normalized value from a Color Image <img> at pixel <q> to a MImage <mimg>
 *        at pixel <p> starting at band <initial_band>.
 *
 * The image must be a COLOR Image.
 * Thus, the 3 LAB2Norm channels will be copied, respectively, at the bands: [initial_band, initial_band+1, initial_band+2].
 * 
 * 
 * @param mimg                Target Multi-Band Image.
 * @param p                   Target pixel at the MImage.
 * @param img                 Source Color Image.
 * @param q                   Source pixel at the Source Image.
 * @param initial_band        Initial band where the values from the source image will copied at the MImage.
 * @param normalization_value Normalization value for the Copy. Use 1.0 for non-normalization.
 *
 * @author Samuka Martins
 * @date July 21, 2018
 */
void iftCopyVoxelToMImageLAB2Norm(iftMImage *mimg, int p, const iftImage *img, int q, int initial_band,
                                         int normalization_value) {
    iftColor YCbCr;
    YCbCr.val[0] = img->val[q];
    YCbCr.val[1] = img->Cb[q];
    YCbCr.val[2] = img->Cr[q];

    iftFColor LabNorm = iftRGBtoLabNorm2(iftYCbCrtoRGB(YCbCr, normalization_value), normalization_value);

    mimg->val[p][initial_band]     = LabNorm.val[0];
    mimg->val[p][initial_band + 1] = LabNorm.val[1];
    mimg->val[p][initial_band + 2] = LabNorm.val[2];
}


/**
 * @brief Extend the feature vector from a voxel <p> from a source multidimensional image <mimg> by adding
 *        it voxel coordinates. The resulting feature vector is copied at the pixel <p> from the target multiband image <emimg>.
 *         
 * @param emimg    Target multiband image.
 * @param mimg     Source multiband image.
 * @param p        Source pixel for extension.
 * @param n_adj_voxels Number of times the coordinates must be
 * repeated must be the number of adjacent voxels.
 * @param norm_val Normalization factor. Use 1.0 to ignore the normalization.
 */
void iftExtendMVoxelFeatVector(iftMImage *emimg, const iftMImage *mimg, int p, int n_adj_voxels, float norm_val) {
    iftVoxel u = iftMGetVoxelCoord(mimg,p);

    if (iftIs3DMImage(mimg)) {
      int b;
      for (b = 0; b < mimg->m; b++) {
  emimg->val[p][b]   = mimg->val[p][b];
      }
      for (int i=0; i < n_adj_voxels; i += 3) {
  emimg->val[p][i+b]   = (float) u.x / norm_val;
  emimg->val[p][i+1+b] = (float) u.y / norm_val;
  emimg->val[p][i+2+b] = (float) u.z / norm_val;
      }
    } else {
      int b;
      for (b = 0; b < mimg->m; b++) {
  emimg->val[p][b]   = mimg->val[p][b];
      }
      for (int i=0; i < n_adj_voxels; i += 2) {
  emimg->val[p][i+b]   = (float) u.x / norm_val;
  emimg->val[p][i+1+b] = (float) u.y / norm_val;
      }
    }
}


/********************** PUBLIC FUNCTIONS *************************/

int iftMXSize(iftMImage *img)
{
  return(img->xsize);
}

int iftMYSize(iftMImage *img)
{
  return(img->ysize);
}

int iftMZSize(iftMImage *img)
{
  return(img->zsize);
}


inline iftVoxel iftMGetVoxelCoord(const iftMImage *img, int p)
{
    /* old
     * u.x = (((p) % (((img)->xsize)*((img)->ysize))) % (img)->xsize)
     * u.y = (((p) % (((img)->xsize)*((img)->ysize))) / (img)->xsize)
     * u.z = ((p) / (((img)->xsize)*((img)->ysize)))
     */
    iftVoxel u;
    div_t res1 = div(p, img->xsize * img->ysize);
    div_t res2 = div(res1.rem, img->xsize);

    u.x = res2.rem;
    u.y = res2.quot;
    u.z = res1.quot;

  return(u);
}


iftMImage  *iftCreateMImage(int xsize,int ysize,int zsize, int nbands)
{
  iftMImage *img=NULL;
  int        i,y,z,xysize;

  img = (iftMImage *) iftAlloc(1,sizeof(iftMImage));
  if (img == NULL){
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateMImage");
  }

  img->n       = xsize*ysize*zsize;
  img->m       = nbands;
  img->data    = iftCreateMatrix(img->m, img->n);
  img->val     = iftAlloc(img->n, sizeof *img->val);
  for (i = 0; i < img->n; i++)
      img->val[i] = iftMatrixRowPointer(img->data, i);
  img->xsize   = xsize;
  img->ysize   = ysize;
  img->zsize   = zsize;
  img->dx      = 1.0;
  img->dy      = 1.0;
  img->dz      = 1.0;
  img->tby     = iftAllocIntArray(ysize);
  img->tbz     = iftAllocIntArray(zsize);

  img->tby[0]=0;
  for (y=1; y < ysize; y++)
    img->tby[y]=img->tby[y-1] + xsize;

  img->tbz[0]=0; xysize = xsize*ysize;
  for (z=1; z < zsize; z++)
    img->tbz[z]=img->tbz[z-1] + xysize;

  return(img);
}

iftMImageArray *iftCreateMImageArray(long n)
{
    if(n <= 0)
        iftError("The number of mimages in the array must be >= 0 (n=%lu)", "iftCreateMImageArray", n);

    iftMImageArray *arr = (iftMImageArray*) iftAlloc(1, sizeof(iftMImageArray));

    arr->val = iftAlloc(n, sizeof(iftMImage*));
    arr->n = n;

    return arr;
}

void iftDestroyMImageArray(iftMImageArray **arr)
{
    iftMImageArray *aux = *arr;

    if (aux != NULL) {
        for(int i = 0; i < aux->n; i++) {
            if(aux->val[i] != NULL)
                iftDestroyMImage(&aux->val[i]);
        }
        iftFree(aux->val);
        iftFree(aux);
        *arr = NULL;
    }
}

iftMImageArray *iftCopyMImageArray(iftMImageArray *arr)
{
    iftMImageArray *copy = iftCreateMImageArray(arr->n);

    for(int i = 0; i < arr->n; i++)
        copy->val[i] = iftCopyMImage(arr->val[i]);

    return copy;
}

iftMImage *iftExtendMImageByLBP(iftMImage *img, iftAdjRel *A, char normalize)
{
  if (img->m!=1){
      iftError("The MImage must have only one band (MImage created from a grayscale image)", "iftExtendMImageByLBP");
  }

  iftMImage *eimg = iftCreateMImage(img->xsize, img->ysize, img->zsize, img->m+1);

  if (normalize){
    int max_lpb_val=INT_MIN;
#pragma omp parallel for shared(eimg,img,A)
    for(int p = 0; p < img->n; p++)
    {
      for(int b = 0; b < img->m ; b++)
        eimg->val[p][b] = img->val[p][b];
      iftVoxel u,v;
      u=iftMGetVoxelCoord(img,p);
      int lbp=0;
      for (int i=1; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A,u,i);
        if (iftMValidVoxel(img,v)){
          int q = iftMGetVoxelIndex(img,v);

          if (img->val[q][0] > img->val[p][0]) {
            lbp |= 1 << (i-1);
          }
        }
      }
      eimg->val[p][eimg->m-1]=lbp;
#pragma omp critical
      if (lbp > max_lpb_val)
        max_lpb_val=lbp;
    }
    /* normalize the lpb band*/
    for (int p=0;p<eimg->n;p++)
      eimg->val[p][eimg->m-1]/=max_lpb_val;
  }
  else {
#pragma omp parallel for shared(eimg,img,A)
    for (int p = 0; p < img->n; p++) {
      for (int b = 0; b < img->m; b++)
        eimg->val[p][b] = img->val[p][b];
      iftVoxel u, v;
      u = iftMGetVoxelCoord(img, p);
      int lbp = 0;
      for (int i = 1; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A, u, i);
        if (iftMValidVoxel(img, v)) {
          int q = iftMGetVoxelIndex(img, v);

          if (img->val[q][0] > img->val[p][0]) {
            lbp |= 1 << (i - 1);
          }
        }
      }
      eimg->val[p][eimg->m - 1] = lbp;
    }
  }

  return(eimg);
}


iftMImage *iftExtendMImageByAdjacencyAndVoxelCoord(const iftMImage *mimg, const iftAdjRel *A, char normalize_voxel_coord) {

  iftMImage *eimg = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, (mimg->m*A->n)+3);

  if (normalize_voxel_coord){
    float max_size = iftMax(iftMax(mimg->xsize,mimg->ysize),mimg->zsize);
#pragma omp parallel for shared(eimg,mimg,max_size,A)
    for(int p = 0; p < mimg->n; p++)
    {
      iftVoxel u = iftMGetVoxelCoord(mimg,p);
      for(int b = 0; b < mimg->m ; b++)
        for (int i = 0; i < A->n; i += 1) {
          iftVoxel v = iftGetAdjacentVoxel(A, u, i);
          if (iftMValidVoxel(mimg, v)) {
            int q = iftMGetVoxelIndex(mimg, v);
            eimg->val[p][i+b*A->n] = mimg->val[q][b];
          }
          else
            eimg->val[p][i+b*A->n] = mimg->val[p][b];
        }
      eimg->val[p][mimg->m*A->n] = (float)u.x/max_size;
      eimg->val[p][(mimg->m*A->n)+1] = (float)u.y/max_size;
      eimg->val[p][(mimg->m*A->n)+2] = (float)u.z/max_size;
    }
  }
  else {
#pragma omp parallel for shared(eimg,mimg,A)
    for(int p = 0; p < mimg->n; p++)
    {
      iftVoxel u = iftMGetVoxelCoord(mimg,p);
      for(int b = 0; b < mimg->m ; b++)
        for (int i = 0; i < A->n; i += 1) {
          iftVoxel v = iftGetAdjacentVoxel(A, u, i);
          if (iftMValidVoxel(mimg, v)) {
            int q = iftMGetVoxelIndex(mimg, v);
            eimg->val[p][i+b*A->n] = mimg->val[q][b];
          }
          else
            eimg->val[p][i+b*A->n] = mimg->val[p][b];
        }
      eimg->val[p][mimg->m*A->n] = (float)u.x;
      eimg->val[p][(mimg->m*A->n)+1] = (float)u.y;
      eimg->val[p][(mimg->m*A->n)+2] = (float)u.z;
    }
  }

  return(eimg);

}

iftMImage *iftExtendMImageByAdjacency(iftMImage *mimg, iftAdjRel *A){

  iftMImage *eimg = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m*A->n);

#pragma omp parallel for shared(eimg,mimg,A)
    for(int p = 0; p < mimg->n; p++)
    {
      iftVoxel u = iftMGetVoxelCoord(mimg,p);
      for(int b = 0; b < mimg->m ; b++)
        for (int i = 0; i < A->n; i += 1) {
          iftVoxel v = iftGetAdjacentVoxel(A, u, i);
          if (iftMValidVoxel(mimg, v)) {
            int q = iftMGetVoxelIndex(mimg, v);
            eimg->val[p][i+b*A->n] = mimg->val[q][b];
          }
          else
            eimg->val[p][i+b*A->n] = mimg->val[p][b];
        }
    }

  return(eimg);

}

iftMImage *iftExtendMImageByImage(iftMImage *mimg, iftImage *img)
{
    if (!iftIsDomainEqual(mimg,img))
        iftError("Input images must have the same domain","iftExtendMImageByImage");

    iftMImage *new = iftCreateMImage(mimg->xsize,mimg->ysize,mimg->zsize,mimg->m+1);

    for (int i = 0; i < mimg->m; i++)
        for (int j = 0; j < mimg->n; j++)
            new->val[j][i] = mimg->val[j][i];

    for (int j = 0; j < img->n; j++)
        new->val[j][mimg->m] = img->val[j];

    return new;
}

iftMImage *iftExtendMImageByMImage(iftMImage *mimg1, iftMImage *mimg2) {
    if (!iftIsDomainEqual(mimg1,mimg2))
        iftError("Input images must have the same domain","iftExtendMImageByMImage");

    iftMImage *new = iftCreateMImage(mimg1->xsize,mimg1->ysize,mimg1->zsize,mimg1->m+mimg2->m);
    iftMCopyVoxelSize(mimg1, new);

    for (int i = 0; i < mimg1->m; i++)
        for (int j = 0; j < mimg1->n; j++)
            new->val[j][i] = mimg1->val[j][i];

    for (int i = 0; i < mimg2->m; i++)
        for (int j = 0; j < mimg2->n; j++)
            new->val[j][mimg1->m + i] = mimg2->val[j][i];

    return new;
}


void  iftDestroyMImage(iftMImage **img)
{
    iftMImage *aux;

    aux = *img;
    if (aux != NULL) {
        if (aux->val != NULL)
            iftFree(aux->val);
        if (aux->data != NULL)
            iftDestroyMatrix(&aux->data);
        if (aux->tby != NULL)
            iftFree(aux->tby);
        if (aux->tbz != NULL)
            iftFree(aux->tbz);
        iftFree(aux);
        *img = NULL;
    }
}

inline char iftMValidVoxel(const iftMImage *img, iftVoxel v) {
  if ((v.x >= 0)&&(v.x < img->xsize)&&
      (v.y >= 0)&&(v.y < img->ysize)&&
      (v.z >= 0)&&(v.z < img->zsize))
    return(1);
  else
    return(0);
}



void iftMCopyVoxelSize(const iftMImage *img1, iftMImage *img2)
{
  img2->dx = img1->dx;
  img2->dy = img1->dy;
  img2->dz = img1->dz;
}

void iftMCopyVoxelSizeFromImage(const iftImage *img1, iftMImage *img2)
{
  img2->dx = img1->dx;
  img2->dy = img1->dy;
  img2->dz = img1->dz;
}

void iftMCopyVoxelSizeToImage(const iftMImage *img1, iftImage *img2)
{
  img2->dx = img1->dx;
  img2->dy = img1->dy;
  img2->dz = img1->dz;
}

iftMImage *iftImageToMImage(const iftImage *img1, char color_space)
{
  iftMImage *img2=NULL;
  iftImage  *input=NULL;
  int normalization_value = iftNormalizationValue(iftMaximumValue(img1));

  if ((iftIsColorImage(img1) == false) &&
      (color_space != GRAY_CSPACE) &&
      (color_space != GRAYNorm_CSPACE) ){ 
    iftColorTable *ctb = iftBlueToRedColorTable(normalization_value);
    input = iftGrayImageToColorImage(img1,ctb);
    iftDestroyColorTable(&ctb);
  } else {
    input = iftCopyImage(img1); 
  }
  
  switch (color_space) {

  case YCbCr_CSPACE:
    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,3);
#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)input->val[p]);
      img2->val[p][1]=((float)input->Cb[p]);
      img2->val[p][2]=((float)input->Cr[p]);
    }
    break;

  case YCbCrNorm_CSPACE:

    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,3);
#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)input->val[p])/(float)normalization_value;
      img2->val[p][1]=((float)input->Cb[p])/(float)normalization_value;
      img2->val[p][2]=((float)input->Cr[p])/(float)normalization_value;
    }
    break;

  case LAB_CSPACE:

    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,3);
#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor  YCbCr,RGB;
      iftFColor Lab;
      YCbCr.val[0] = input->val[p];
      YCbCr.val[1] = input->Cb[p];
      YCbCr.val[2] = input->Cr[p];
      RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
      Lab = iftRGBtoLab(RGB,normalization_value);
      img2->val[p][0]=Lab.val[0];
      img2->val[p][1]=Lab.val[1];
      img2->val[p][2]=Lab.val[2];
    }
    break;

    case LABNorm_CSPACE:

      img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,3);
#pragma omp parallel for shared(input, img2, normalization_value)
      for (int p=0; p < img2->n; p++) {
        iftColor  YCbCr,RGB;
        iftFColor Lab;
        YCbCr.val[0] = input->val[p];
        YCbCr.val[1] = input->Cb[p];
        YCbCr.val[2] = input->Cr[p];
        RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
        Lab = iftRGBtoLabNorm(RGB,normalization_value);
        img2->val[p][0]=Lab.val[0];
        img2->val[p][1]=Lab.val[1];
        img2->val[p][2]=Lab.val[2];
      }
      break;

  case LABNorm2_CSPACE:

    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,3);
#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor  YCbCr,RGB;
      iftFColor LabNorm;
      YCbCr.val[0] = input->val[p];
      YCbCr.val[1] = input->Cb[p];
      YCbCr.val[2] = input->Cr[p];
      RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
      LabNorm = iftRGBtoLabNorm2(RGB,normalization_value);
      img2->val[p][0]=LabNorm.val[0];
      img2->val[p][1]=LabNorm.val[1];
      img2->val[p][2]=LabNorm.val[2];
    }
    break;

  case RGB_CSPACE:

    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,3);
#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor YCbCr,RGB;
      YCbCr.val[0] = input->val[p];
      YCbCr.val[1] = input->Cb[p];
      YCbCr.val[2] = input->Cr[p];
      RGB          = iftYCbCrtoRGB(YCbCr,normalization_value);
      img2->val[p][0]=((float)RGB.val[0]);
      img2->val[p][1]=((float)RGB.val[1]);
      img2->val[p][2]=((float)RGB.val[2]);
    }
    break;

  case RGBNorm_CSPACE:

    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,3);
#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor YCbCr,RGB;
      YCbCr.val[0] = input->val[p];
      YCbCr.val[1] = input->Cb[p];
      YCbCr.val[2] = input->Cr[p];
      RGB          = iftYCbCrtoRGB(YCbCr,normalization_value);
      img2->val[p][0]=((float)RGB.val[0])/(float)normalization_value;
      img2->val[p][1]=((float)RGB.val[1])/(float)normalization_value;
      img2->val[p][2]=((float)RGB.val[2])/(float)normalization_value;
    }
    break;

  case GRAY_CSPACE:

    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,1);
#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)input->val[p]);
    }
    break;

  case GRAYNorm_CSPACE:

    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,1);
#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)input->val[p])/(float)normalization_value;
    }

    break;

  case WEIGHTED_YCbCr_CSPACE:

    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,3);

#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=(0.2/2.2)*((float)input->val[p]/(float)normalization_value);
      img2->val[p][1]=(1.0/2.2)*((float)input->Cb[p]/(float)normalization_value);
      img2->val[p][2]=(1.0/2.2)*((float)input->Cr[p]/(float)normalization_value);
    }
    break;
  
  case HSV_CSPACE:

    img2=iftCreateMImage(input->xsize,input->ysize,input->zsize,3);
#pragma omp parallel for shared(input, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor  YCbCr,RGB;
      iftColor HSV;
      YCbCr.val[0] = input->val[p];
      YCbCr.val[1] = input->Cb[p];
      YCbCr.val[2] = input->Cr[p];
      RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
      HSV = iftRGBtoHSV(RGB,normalization_value);
      img2->val[p][0]=HSV.val[0];
      img2->val[p][1]=HSV.val[1];
      img2->val[p][2]=HSV.val[2];
    }
    break;

  default:
      iftError("Invalid color space (see options in iftColor.h)", "iftImageToMImage");
  }

  img2->dx = input->dx;
  img2->dy = input->dy;
  img2->dz = input->dz;

  iftDestroyImage(&input);
  

  return(img2);
}

 iftImage   *iftMImageToImage(const iftMImage *img1, int Imax, int band)
{
  iftImage *img2=iftCreateImage(img1->xsize,img1->ysize,img1->zsize);
  int p,b=band;
  double min = IFT_INFINITY_FLT, max = IFT_INFINITY_FLT_NEG;

  if ((band < 0)||(band >= img1->m))
      iftError("Invalid band", "iftMImageToImage");

  for (p=0; p < img1->n; p++) {
    if (img1->val[p][b] < min)
      min = img1->val[p][b];
    if (img1->val[p][b] > max)
      max = img1->val[p][b];
  }

  //printf("min %lf max %lf\n",min,max);

  if (max > min){
    for (p=0; p < img2->n; p++) {
      img2->val[p]=(int)(Imax*(img1->val[p][b]-min)/(max-min));
    }
  }/*else{
    char msg[200];
    sprintf(msg,"Image is empty: max = %f and min = %f\n",max,min);
    puts(msg);
    // iftWarning(msg,"iftMImageToImage");
  }*/

  img2->dx = img1->dx;
  img2->dy = img1->dy;
  img2->dz = img1->dz;
  
  return(img2);
}

iftMImage  *iftCopyMImage(iftMImage *img)
{
  iftMImage *imgc=iftCreateMImage(img->xsize,img->ysize,img->zsize,img->m);
  int p, b;

  iftMCopyVoxelSize(img,imgc);

  for (b=0; b < img->m; b++)
    for (p=0; p < img->n; p++)
      imgc->val[p][b]=img->val[p][b];

  return(imgc);
}

iftMImage *iftReadMImage(const char *filename)
{

    iftMImage *mimg=NULL;
    int xsize,ysize,zsize;
    unsigned long m;
    char msg[200], type[10];
    FILE *fp = fopen(filename,"r");

    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadMImage", filename);

    if (fscanf(fp,"%s\n",type) != 1)
        iftError("Reading error", "iftReadMImage");

    if((strcmp(type,"MIG")==0)){
        // Write X Y Z nbands
        if(fscanf(fp,"%d %d %d %lu\n",&xsize,&ysize,&zsize,&m)!=4)
            iftError("Reading error", "iftReadMImage");
        mimg = iftCreateMImage(xsize,ysize,zsize,m);
        // Write dx dy dz
        if(fscanf(fp,"%f %f %f\n",&mimg->dx,&mimg->dy,&mimg->dz)!=3)
            iftError("Reading error", "iftReadMImage");

        fseek(fp,-(mimg->n*mimg->m*sizeof(float)),SEEK_END);

        if (fread(mimg->data->val, sizeof(float), mimg->n * mimg->m, fp) != mimg->n * mimg->m) {
            sprintf(msg, "Could not read pixels of the MImage");
            iftError(msg, "iftReadMImage");
        }
  } else {
      iftError("Input image must be a multi-band image", "iftReadMImage");
  }

  fclose(fp);
  return(mimg);
}

iftMImage *iftReadMImageFromNumPy(const char *format, ...) {
  va_list args;
  char npy_path[IFT_STR_DEFAULT_SIZE];

  va_start(args, format);
  vsprintf(npy_path, format, args);
  va_end(args);

  iftNumPyHeader *header = NULL;
  void *data = iftReadNumPy(npy_path, &header);
  iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);

  // (ysize, xsyze, nbands) or
  // (zsize, ysize, xsyze, nbands)
  int xsize = 0, ysize = 0, zsize = 0, n_bands = 0;

  if (header->n_dims == 3) {
    zsize = 1;
    ysize = header->shape[0];
    xsize = header->shape[1];
    n_bands = header->shape[2];
  }
  else if (header->n_dims == 4) {
    zsize = header->shape[0];
    ysize = header->shape[1];
    xsize = header->shape[2];
    n_bands = header->shape[3];
  }
  else {
    iftError("Number of dimensions %ld must be 3 or 4", "iftReadMImageNumPy", header->n_dims);
  }

  iftMImage *mimg = iftCreateMImage(xsize, ysize, zsize, n_bands);
  iftFloatArray *arr = NULL;

  if (cdtype == IFT_FLT_TYPE) {
    arr = iftCreateFloatArray(header->size);
    iftFree(arr->val);
    arr->val = data;
  }
  else {
    arr = iftConvertToFloatArray(data, header->size, cdtype);
    iftFree(data);
  }

  #pragma omp parallel for
  for (int p = 0; p < mimg->n; p++) {
    int start_voxel_idx = p * n_bands;

    for (int b = 0; b < n_bands; b++) {
      mimg->val[p][b] = arr->val[start_voxel_idx + b];
    }
  }

  iftDestroyFloatArray(&arr);
  iftDestroyNumPyHeader(&header);

  return mimg;
}


void  iftWriteMImage(iftMImage *mimg, const char *filename){
    FILE *fp = fopen(filename,"wb");

    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteMImage", filename);

    fprintf(fp,"MIG\n");
    // Write X Y Z nbands
    fprintf(fp,"%d %d %d %lu\n",mimg->xsize,mimg->ysize,mimg->zsize,mimg->m);
    // Write dx dy dz
    fprintf(fp,"%f %f %f\n",mimg->dx,mimg->dy,mimg->dz);
    // Write pixels
    if (mimg->n*mimg->m < 0)
        iftError("MImage size overflow, Please check if image size is correct.\nNumber of bands %d and image size %d.","iftWriteMImage",mimg->m,mimg->n);
    fwrite(mimg->data->val,sizeof(float),mimg->n*mimg->m,fp);
    fclose(fp);
}

void iftWriteMImageAsNumPy(const iftMImage *mimg, const char *format, ...) {
  va_list args;
  char npy_path[IFT_STR_DEFAULT_SIZE];

  va_start(args, format);
  vsprintf(npy_path, format, args);
  va_end(args);

  if (!iftEndsWith(npy_path, ".npy"))
    iftError("Invalid extension for numpy: %s\nTry *.npy",
             "iftWriteMImageAsNumPy", npy_path);

  int n_bands = mimg->m;

  int n_dims;
  long shape[4];

  if (iftIs3DMImage(mimg)) {
    n_dims = 4;
    shape[0] = mimg->zsize;
    shape[1] = mimg->ysize;
    shape[2] = mimg->xsize;
    shape[3] = n_bands;
  }
  else {
    n_dims = 3;
    shape[0] = mimg->ysize;
    shape[1] = mimg->xsize;
    shape[2] = n_bands;
  }

  iftFloatArray *arr = iftCreateFloatArray(mimg->n * n_bands);

  #pragma omp parallel for
  for (int p = 0; p < mimg->n; p++) {
    int start_voxel_idx = p * n_bands;

    for (int b = 0; b < n_bands; b++) {
      arr->val[start_voxel_idx + b] = mimg->val[p][b];
    }
  }

  iftNumPyHeader *header = iftCreateNumPyHeader(IFT_FLT_TYPE, shape, n_dims);
  iftWriteNumPy(header, arr->val, npy_path);
  iftDestroyNumPyHeader(&header);
  iftDestroyFloatArray(&arr);
}


void iftWriteMImageBands(iftMImage *mimg, char *base_filename){
  int i;
  char buffer[255];
  for(i = 0; i < mimg->m; i++){
    iftImage *img = iftMImageToImage(mimg,255,i);
    sprintf(buffer, "%s-%d.pgm", base_filename, i);
//		printf("%s\n",buffer);
    iftWriteImageP5(img,buffer);
    buffer[0] = 0;

    iftDestroyImage(&img);
  }
}


iftMImage *iftMAddFrame(iftMImage *img, int bx, int by, int bz, float value)
{
  iftMImage *fimg;


  if (iftIs3DMImage(img)) {
    fimg = iftCreateMImage(img->xsize+(2*bx),img->ysize+(2*by), img->zsize+(2*bz), img->m);
    iftMCopyVoxelSize(img,fimg);
    iftSetMImage(fimg,value);
#pragma omp parallel shared(fimg,img,bx,by,bz)
    for (int b=0; b < fimg->m; b++) {
      iftVoxel u; int p = 0;
      for (u.z=bz; u.z < fimg->zsize-bz; u.z++)
  for (u.y=by; u.y < fimg->ysize-by; u.y++)
    for (u.x=bx; u.x < fimg->xsize-bx; u.x++) {
      int q = iftMGetVoxelIndex(fimg,u);
      fimg->val[q][b] = img->val[p][b];
      p++;
    }
    }
  } else {
    fimg = iftCreateMImage(img->xsize+(2*bx),img->ysize+(2*by), 1, img->m);
    iftSetMImage(fimg,value);

#pragma omp parallel shared(fimg,img,bx,by)
    for (int b=0; b < fimg->m; b++) {
      iftVoxel u; u.z=0; int p=0;
      for (u.y=by; u.y < fimg->ysize-by; u.y++)
  for (u.x=bx; u.x < fimg->xsize-bx; u.x++) {
    int q = iftMGetVoxelIndex(fimg,u);
    fimg->val[q][b] = img->val[p][b];
    p++;
  }
    }
  }
  return(fimg);
}

iftMImage *iftMRemFrame(iftMImage *fimg, int bx, int by, int bz) {
  iftMImage *img;


  if (iftIs3DMImage(fimg)) {
    img = iftCreateMImage(fimg->xsize-(2*bx),fimg->ysize-(2*by),fimg->zsize-(2*bz), fimg->m);
    iftMCopyVoxelSize(fimg,img);
#pragma omp parallel shared(fimg,img,bx,by,bz)
    for (int b=0; b < fimg->m; b++) {
      iftVoxel u; int p = 0;
      for (u.z=bz; u.z < fimg->zsize-bz; u.z++)
  for (u.y=by; u.y < fimg->ysize-by; u.y++)
    for (u.x=bx; u.x < fimg->xsize-bx; u.x++) {
      int q = iftMGetVoxelIndex(fimg,u);
      img->val[p][b] = fimg->val[q][b];
      p++;
    }
    }
  } else {
    img = iftCreateMImage(fimg->xsize-(2*bx),fimg->ysize-(2*by),1, fimg->m);

#pragma omp parallel shared(fimg,img,bx,by)
    for (int b=0; b < fimg->m; b++) {
      iftVoxel u; u.z=0; int p = 0;
      for (u.y=by; u.y < fimg->ysize-by; u.y++)
  for (u.x=bx; u.x < fimg->xsize-bx; u.x++) {
    int q = iftMGetVoxelIndex(fimg,u);
    img->val[p][b] = fimg->val[q][b];
    p++;
  }
    }
  }
  return(img);
}

void iftSetMImage(iftMImage *img, float value)
{
  int p, b;

  for (b=0; b < img->m; b++)
    for (p=0; p < img->n; p++)
      img->val[p][b]=value;

}


iftImage *iftMImageBasins(const iftMImage *img, iftAdjRel *A)
{
   iftImage   *basins=iftCreateImage(img->xsize,img->ysize,img->zsize);
   double     *grad=iftAllocDoubleArray(img->n);
   float      *w, wt=0.0, K;
   int         dx, dy, dz;

   iftMaxAdjShifts(A, &dx, &dy, &dz);
   K = sqrtf(dx*dx + dy*dy + dz*dz);

   w = iftAllocFloatArray(A->n);
   for (int i=1; i < A->n; i++) {
     w[i] = K / sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i]);
     wt  += w[i];
   }
   for (int i=1; i < A->n; i++) {
     w[i] = w[i]/wt;
   }

   /* Compute basins image in the spatial domain */

#pragma omp parallel for shared(img,grad,A)
   for (int p=0; p < img->n; p++) {
     iftVoxel u   = iftMGetVoxelCoord(img,p);
     for (int i=1; i < A->n; i++) {
       iftVoxel v = iftGetAdjacentVoxel(A,u,i);
       double dist=0.0;
       if (iftMValidVoxel(img,v)){
   int q = iftMGetVoxelIndex(img,v);
   for (int b=0; b < img->m; b++) {
     dist += fabs(img->val[q][b]-img->val[p][b]);
   }
       }
       grad[p] += dist*w[i];
     }
   }
   
#pragma omp parallel for shared(grad,basins)
   for (int p=0; p < img->n; p++) {
     basins->val[p] = iftRound(grad[p]);
   }

   iftFree(grad);
   iftFree(w);
   
   return(basins);
 }


iftImage *iftMImageGradient(const iftMImage *img, iftAdjRel *A, int Imax) {
    iftFImage *gradI = iftCreateFImage(img->xsize,img->ysize,img->zsize);
    float Gmin=IFT_INFINITY_FLT, Gmax=IFT_INFINITY_FLT_NEG;

    for (ulong p=0; p < img->n; p++){
        iftVoxel u  = iftMGetVoxelCoord(img,p);
        float value = 0.0;
        for (int i=1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);
            if (iftMValidVoxel(img,v)){
            int q = iftMGetVoxelIndex(img,v);
                for (ulong b=0; b < img->m; b++)
                    value += (img->val[q][b]-img->val[p][b])*(img->val[q][b]-img->val[p][b]);
            }
        }
        value = sqrtf(value)/A->n;
        if (value > Gmax) Gmax = value;
        if (value < Gmin) Gmin = value;
        gradI->val[p] = value;
    }

    for (ulong p=0; p < img->n; p++){
        gradI->val[p] = (Imax*(gradI->val[p]-Gmin)/(Gmax-Gmin));
    }

    iftImage *out = iftCreateImage(gradI->xsize,gradI->ysize,gradI->zsize);
    for (ulong p=0; p < img->n; p++){
        out->val[p] = iftRound(gradI->val[p]);
    }
    iftDestroyFImage(&gradI);

    return(out);
}

iftImage *iftBorderProbImage(iftMImage *img) {
  iftAdjRel *A;
  iftImage  *prob;

  if (iftIs3DMImage(img))
    A = iftSpheric(sqrtf(3.0));
  else
    A = iftCircular(1.5);

  prob = iftMImageBasins(img, A);

  int      prob_max_val = iftMaximumValue(prob);
  for (int p            = 0; p < prob->n; p++)
    prob->val[p] = (int) ((float) prob->val[p] / prob_max_val * 100.0);

  iftDestroyAdjRel(&A);

  return (prob);
}

iftImage *iftUniformProbImage(iftMImage *img)
{
  iftImage  *prob = iftCreateImage(img->xsize, img->ysize, img->zsize);

  for (int p=0; p < prob->n; p++)
    prob->val[p] = 100;

  return(prob);
}

iftImage *iftRegionProbImage(iftMImage *img) {
  iftAdjRel *A;
  iftImage  *prob;

  if (iftIs3DMImage(img))
    A = iftSpheric(sqrt(3.0));
  else
    A = iftCircular(sqrt(2.0));

  prob = iftMImageBasins(img, A);
  int      max_val = iftMaximumValue(prob);
  for (int p       = 0; p < prob->n; p++)
    prob->val[p] = (int) ((float) (max_val - prob->val[p]) / max_val * 100.0);

  iftDestroyAdjRel(&A);

  return (prob);
}

void iftMultMImageByScalar(iftMImage *mimg, float scalar)
{
#pragma omp parallel for shared(mimg,scalar)
  for (int p=0; p < mimg->n; p++)
    for (int b=0; b < mimg->m; b++)
      mimg->val[p][b] *= scalar;
}

iftImage *iftSelectNonBorderVoxels(iftMImage *img, iftImage *mask1, int nsamples)
{
  iftImage  *prob  = iftBorderProbImage(img);
  iftImage  *mask2 = iftCreateImage(img->xsize,img->ysize,img->zsize);
  int        p, q, i, niters, xsize, ysize, zsize;
  iftAdjRel *A = NULL;
  iftVoxel   u, v, uo, uf;
  char       select;
  float      nsamples_per_axis, xspacing,yspacing,zspacing,spacing,radius;

  iftBoundingBox mbb = iftMinBoundingBox(mask1, NULL);
  uo = mbb.begin;
  uf = mbb.end;
  xsize = uf.x - uo.x + 1;
  ysize = uf.y - uo.y + 1;
  zsize = uf.z - uo.z + 1;


  if (iftIs3DMImage(img)){
    /* finds displacements along each axis in order to obtain nsamples_per_axis inside the mask */
    nsamples_per_axis = (float) pow((double)nsamples,1.0/3.0);
    xspacing          = (xsize/nsamples_per_axis);
    yspacing          = (ysize/nsamples_per_axis);
    zspacing          = (zsize/nsamples_per_axis);
    spacing           = (xspacing+yspacing+zspacing)/3.0;
    radius            = spacing/2.0;
    A                 = iftSpheric(radius);
  }else{
    nsamples_per_axis = (float) sqrtf((float)nsamples);
    xspacing          = (xsize/nsamples_per_axis);
    yspacing          = (ysize/nsamples_per_axis);
    spacing           = (xspacing+yspacing)/2.0;
    radius            = spacing/2.0;
    A                 = iftCircular(radius);
  }

  niters = 100*nsamples;
  while ((nsamples > 0)&&(niters > 0)){
    p = iftRandomInteger(0,mask1->n-1);
    if (mask1->val[p]>0){
      if (prob->val[p]>0){
  prob->val[p]--;
      }else{ /* select voxel */
  niters--;
  if (mask2->val[p]==0){
    select = 1;
    u = iftGetVoxelCoord(mask1,p);
    for (i=1; (i < A->n)&&(select); i++) {
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(mask1,v)){
        q = iftGetVoxelIndex(mask1,v);
        if (mask2->val[q]==1) {
    select=0;
        }
      }
    }
    if (select) {
      mask2->val[p]=1;
      nsamples--;
    }
  }
      }
    }
  }


  iftDestroyImage(&prob);
  iftDestroyAdjRel(&A);

  return(mask2);
}

iftImage *iftSelectNSamplesFromMask(iftMImage *img, iftImage *mask1, int nsamples)
{
  iftImage  *mask2 = iftCreateImage(img->xsize,img->ysize,img->zsize);
  int        p, q, i, niters, xsize, ysize, zsize;
  iftAdjRel *A;
  iftVoxel   u, v, uo, uf;
  char       select;
  float      nsamples_per_axis, xspacing,yspacing,zspacing,spacing,radius;

  /* compute space constraint */
  iftBoundingBox mbb = iftMinBoundingBox(mask1, NULL);
  uo = mbb.begin;
  uf = mbb.end;
  xsize = uf.x - uo.x + 1;
  ysize = uf.y - uo.y + 1;
  zsize = uf.z - uo.z + 1;

  if (iftIs3DMImage(img)){
    /* finds displacements along each axis in order to obtain nsamples_per_axis inside the mask */
    nsamples_per_axis = (float) pow((double)nsamples,1.0/3.0);
    xspacing          = (xsize/nsamples_per_axis);
    yspacing          = (ysize/nsamples_per_axis);
    zspacing          = (zsize/nsamples_per_axis);
    spacing           = (xspacing+yspacing+zspacing)/3.0;
    radius            = spacing/2.0;
    A = iftSpheric(radius);
  }else{
    nsamples_per_axis = (float) sqrtf((float)nsamples);
    xspacing          = (xsize/nsamples_per_axis);
    yspacing          = (ysize/nsamples_per_axis);
    spacing           = (xspacing+yspacing)/2.0;
    radius            = spacing/2.0;
    A = iftCircular(radius);
  }

  niters = 100*nsamples;
  while ((nsamples > 0)&&(niters > 0)){
    p = iftRandomInteger(0,mask1->n-1);
    if (mask1->val[p]>0){
      niters--;
      if (mask2->val[p]==0){
  select = 1;
  u = iftGetVoxelCoord(mask1,p);
  for (i=1; (i < A->n)&&(select); i++) {
    v = iftGetAdjacentVoxel(A,u,i);
    if (iftValidVoxel(mask1,v)){
      q = iftGetVoxelIndex(mask1,v);
      if (mask2->val[q]==1) {
        select=0;
      }
    }
  }
  if (select) {
    mask2->val[p]=1;
    nsamples--;
  }
      }
    }
  }

  //printf("nseeds %d\n",iftNumberOfElements(mask2));
  iftDestroyAdjRel(&A);

  return(mask2);
}

iftImage *iftGridSampling(iftMImage *img, iftImage *mask1, int nsamples)
{
  iftImage  *mask2 = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage  *prob  = iftBorderProbImage(img);
  int        p, q, i, qmin, xsize,ysize,zsize;
  float      xspacing,yspacing,zspacing,deltax,deltay,deltaz;
  iftAdjRel *A;
  iftVoxel   u, v, m, uo, uf;

  /* Compute the extreme voxels that define the region of interest in
     mask1, and then compute the xsize, ysize, and zsize of that
     ROI. */
  iftBoundingBox mbb = iftMinBoundingBox(mask1, NULL);
  uo = mbb.begin;
  uf = mbb.end;
  xsize = uf.x - uo.x + 1;
  ysize = uf.y - uo.y + 1;
  zsize = uf.z - uo.z + 1;

  if (iftIs3DMImage(img)){

    A  = iftSpheric(sqrtf(3.0));
    /* finds displacements along each axis */
    /* uncomment the next 4 lines to use same number o superpixels per axis */
    /*
    float nsamples_per_axis = (float) pow((double)nsamples,1.0/3.0);
    xspacing          = (xsize/nsamples_per_axis);
    yspacing          = (ysize/nsamples_per_axis);
    zspacing          = (zsize/nsamples_per_axis);
    */
    /* uncomment the next 5 lines to obtain equally spaced seeds in every axis */
    float superpixelsize = 0.5+(float)(xsize*ysize*zsize)/(float)nsamples;
    float step = (float) pow((double)superpixelsize,1.0/3.0)+0.5;
    xspacing = step;
    yspacing = step;
    zspacing = step;

    deltax            = xspacing/2.0;
    deltay            = yspacing/2.0;
    deltaz            = zspacing/2.0;

    if ((deltax < 1.0)||(deltay < 1.0)||(deltaz < 1.0))
        iftError("The number of samples is too high", "iftGridSampling");

    for (m.z=(int)deltaz; m.z < (zsize-deltaz); m.z = (int)(m.z + zspacing)) {
      for (m.y=(int)deltay; m.y < (ysize-deltay); m.y = (int)(m.y + yspacing)) {
  for (m.x=(int)deltax; m.x < (xsize-deltax); m.x = (int)(m.x + xspacing)) {
    u.z = uo.z + m.z; u.y = uo.y + m.y; u.x = uo.x + m.x;
    p = iftGetVoxelIndex(mask1,u);
    if (mask1->val[p]!=0){
      for (i=1, qmin=p; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A,u,i);
        q = iftGetVoxelIndex(mask1,v);
        if (iftValidVoxel(mask1,v) && (prob->val[q]<prob->val[qmin])&&
    (mask2->val[q]==0)&&
      (mask1->val[q]!=0))
    qmin=q;
      }
      mask2->val[qmin]=1;
    }
  }
      }
    }
  }else{
    A   = iftCircular(sqrtf(2.0));

    /* finds displacements along each axis  */
    /* uncomment the next 3 lines to use same number o superpixels per axis */
//    float nsamples_per_axis = (float) sqrt((double)nsamples);
//    xspacing          = (xsize/nsamples_per_axis);
//    yspacing          = (ysize/nsamples_per_axis);
    /* uncomment the next 4 lines to obtain equally spaced seeds in every axis */
    float superpixelsize = 0.5+(float)(xsize*ysize)/(float)nsamples;
    float step = sqrt(superpixelsize)+0.5;
    xspacing = step;
    yspacing = step;


    deltax            = xspacing/2.0;
    deltay            = yspacing/2.0;

    if ((deltax < 1.0)||(deltay < 1.0))
        iftError("The number of samples is too high", "iftGridSampling");

    u.z = 0;
    for (m.y=(int)deltay; m.y < ysize; m.y = (int)(m.y + yspacing)){
      for (m.x=(int)deltax; m.x < xsize; m.x = (int)(m.x + xspacing)) {
  u.y = uo.y + m.y; u.x = uo.x + m.x;
  p = iftGetVoxelIndex(mask1,u);
  if (mask1->val[p]!=0){
    for (i=1, qmin=p; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(mask1,v)){
        q = iftGetVoxelIndex(mask1,v);
        if ((prob->val[q]<prob->val[qmin])&&
            (mask2->val[q]==0)&&
            (mask1->val[q]!=0))
          qmin=q;
      }
    }
    mask2->val[qmin]=1;
  }
      }
    }
  }

  //printf("nseeds: %d\n",iftNumberOfElements(mask2));

  iftDestroyImage(&prob);
  iftDestroyAdjRel(&A);

  return(mask2);
}

float iftNormalizedShannonEntropy(float *arr, int size ) {
  int i, nbins, b;
  int *histogram;
  float im, entropy, binsize;
  float minVal, maxVal, range, factor;
  minVal = IFT_INFINITY_FLT;
  maxVal = IFT_INFINITY_FLT_NEG;
  binsize = 5.0;
  entropy = 0.0;
  factor = 0.0;
  // Quantize
  for (i=0; i < size; i++) {
    if (arr[i] < minVal)
      minVal = arr[i];
    if (arr[i] > maxVal)
      maxVal = arr[i];
  }
  range = maxVal - minVal;
  nbins = iftMax((int)(range / binsize), 1);
  histogram = iftAllocIntArray(nbins);

  if (range > 0)
    factor = ((float)nbins)/range;

  for (i=0; i < size; i++) {
    b = (int)( ( arr[i]- minVal )*factor );
    b = iftMin(b, nbins - 1);
    histogram[b]++;
  }

  // Compute entropy
  im = (float)size;
  for (i=0; i < nbins; i++) {
    if (histogram[i] > 0)
      entropy += -((float)histogram[i]/im) * log((float)histogram[i]/im);
  }
  iftFree(histogram);

  entropy /= log(size);
  return entropy;

}

iftImage *iftAltMixedSamplingSecondStage(iftMImage *img, iftImage *mask, int nsamples)
{
  iftImage  *out, *quadMask, *quadGrid;
  iftMImage *quadImg;
  int i,j,k, t, p, initX, initY, initZ, endX, endY, endZ, x, y, z, origp, quadValuesSize, indexQV, indexQuad;
  int *quadNSamples;
  int nquad, nquadz, quadImgSize, quadImgSizeX, quadImgSizeY, quadImgSizeXY;
  float *quadEntropy, totalEntropy, *quadValues;
  nquad = 2;
  nquadz = 2;
  indexQuad = 0;
  totalEntropy = 0.0;
  if (img->zsize == 1)
    nquadz = 1;
  quadEntropy = iftAllocFloatArray(nquad*nquad*nquadz);
  quadNSamples = iftAllocIntArray(nquad*nquad*nquadz);

  out = iftCreateImage(img->xsize, img->ysize, img->zsize);

  for(k=0; k<nquadz; k++) {
    for (i = 0; i < nquad; i++) {
      for (j = 0; j < nquad; j++) {
        // Compute init and end coordinates
        initZ = k * (img->zsize / nquadz);
        initY = i * (img->ysize / nquad);
        initX = j * (img->xsize / nquad);
        endZ = (k + 1) * (img->zsize / nquadz);
        endY = (i + 1) * (img->ysize / nquad);
        endX = (j + 1) * (img->xsize / nquad);
        if (k == (nquadz - 1))
          endZ = img->zsize;
        if (img->zsize == 1)
          endZ = img->zsize;
        if (i == (nquad - 1))
          endY = img->ysize;
        if (j == (nquad - 1))
          endX = img->xsize;

        // Divide image in quadrants
        quadImgSizeY = (endY - initY);
        quadImgSizeX = (endX - initX);
        quadImgSizeXY = quadImgSizeX * quadImgSizeY;
        quadImgSize = quadImgSizeXY * (endZ - initZ);
        quadValues = iftAllocFloatArray(quadImgSize);
        indexQV = 0;
        for (p = 0; p < quadImgSize; p++) {
          z = p / quadImgSizeXY;
          y = (p % quadImgSizeXY) / quadImgSizeX;
          x = (p % quadImgSizeXY) % quadImgSizeX;
          origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
          quadValues[indexQV] = img->val[origp][0];
          indexQV++;
        }
        quadValuesSize = indexQV;
        // Compute entropy
        quadEntropy[indexQuad] = iftNormalizedShannonEntropy(quadValues, quadValuesSize);
        totalEntropy += quadEntropy[indexQuad];

        indexQuad++;
        iftFree(quadValues);
      }
    }
  }

  indexQuad = 0;
  for(k=0; k<nquadz; k++) {
    for (i = 0; i < nquad; i++) {
      for (j = 0; j < nquad; j++) {
        if (totalEntropy == 0)
          quadNSamples[indexQuad] = iftRound(((float) nsamples / (float) (nquad * nquad * nquadz) ));
        else
          quadNSamples[indexQuad] = iftRound((quadEntropy[indexQuad] / totalEntropy) * ((float) nsamples));

        if (quadNSamples[indexQuad] == 0)
          quadNSamples[indexQuad] = 1;

        // Compute init and end coordinates
        initZ = k * (img->zsize / nquadz);
        initY = i * (img->ysize / nquad);
        initX = j * (img->xsize / nquad);
        endZ = (k + 1) * (img->zsize / nquadz);
        endY = (i + 1) * (img->ysize / nquad);
        endX = (j + 1) * (img->xsize / nquad);
        if (k == (nquadz - 1))
          endZ = img->zsize;
        if (img->zsize == 1)
          endZ = img->zsize;
        if (i == (nquad - 1))
          endY = img->ysize;
        if (j == (nquad - 1))
          endX = img->xsize;

        // Divide image in quadrants
        quadImgSizeY = (endY - initY);
        quadImgSizeX = (endX - initX);
        quadImgSizeXY = quadImgSizeX * quadImgSizeY;
        quadImg = iftCreateMImage(quadImgSizeX, quadImgSizeY, (endZ - initZ), img->m);
        quadMask = iftCreateImage(quadImgSizeX, quadImgSizeY, (endZ - initZ));

        for (p = 0; p < quadImg->n; p++) {
          z = p / quadImgSizeXY;
          y = (p % quadImgSizeXY) / quadImgSizeX;
          x = (p % quadImgSizeXY) % quadImgSizeX;
          origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
          for (t = 0; t < quadImg->m; t++) {
            quadImg->val[p][t] = img->val[origp][t];
          }
          if (mask->val[origp] != 0)
            quadMask->val[p] = 1;
        }
        quadGrid = iftGridSampling(quadImg, quadMask, quadNSamples[indexQuad]);
        // Write seeds in the final out image
        for (p = 0; p < quadGrid->n; p++) {
          if (quadGrid->val[p] > 0) {
            z = p / (quadGrid->xsize*quadGrid->ysize);
            y = (p % (quadGrid->xsize*quadGrid->ysize)) / quadGrid->xsize;
            x = (p % (quadGrid->xsize*quadGrid->ysize)) % quadGrid->xsize;
            origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
            out->val[origp] = 1;
          }
        }

        indexQuad++;
        iftDestroyMImage(&quadImg);
        iftDestroyImage(&quadMask);
        iftDestroyImage(&quadGrid);

      }
    }
  }
  iftFree(quadEntropy);
  iftFree(quadNSamples);

  return out;
}

iftImage *iftAltMixedSampling(iftMImage *img, iftImage *mask, int nsamples)
{
  iftImage  *out, *quadMask, *quadGrid;
  iftMImage *quadImg;
  int i,j,k, t, p, initX, initY, initZ, endX, endY, endZ, x, y, z, origp, quadValuesSize, indexQV, indexQuad;
  int *quadNSamples;
  int nquad, nquadz, quadImgSize, quadImgSizeX, quadImgSizeY, quadImgSizeXY;
  float *quadEntropy, totalEntropy, *quadValues, meanEntropy, thrEntropy, sdEntropy;
  nquad = 2;
  nquadz =2;
  indexQuad = 0;
  totalEntropy = 0.0;
  if (img->zsize == 1)
    nquadz = 1;
  quadEntropy = iftAllocFloatArray(nquad*nquad*nquadz);
  quadNSamples = iftAllocIntArray(nquad*nquad*nquadz);

  out = iftCreateImage(img->xsize, img->ysize, img->zsize);

  for(k=0; k<nquadz; k++) {
    for(i=0; i<nquad; i++) {
      for(j=0; j<nquad; j++) {
        // Compute init and end coordinates
        initZ = k * (img->zsize / nquadz);
        initY = i * (img->ysize / nquad);
        initX = j * (img->xsize / nquad);
        endZ = (k + 1) * (img->zsize / nquadz);
        endY = (i + 1) * (img->ysize / nquad);
        endX = (j + 1) * (img->xsize / nquad);
        if (k == (nquadz - 1))
          endZ = img->zsize;
        if (img->zsize == 1)
          endZ = img->zsize;
        if (i == (nquad - 1))
          endY = img->ysize;
        if (j == (nquad - 1))
          endX = img->xsize;

        // Divide image in quadrants
        quadImgSizeY = (endY - initY);
        quadImgSizeX = (endX - initX);
        quadImgSizeXY = quadImgSizeX * quadImgSizeY;
        quadImgSize = quadImgSizeXY * (endZ - initZ);
        quadValues = iftAllocFloatArray(quadImgSize);
        indexQV = 0;
        for (p = 0; p < quadImgSize; p++) {
          z = p / quadImgSizeXY;
          y = (p % quadImgSizeXY) / quadImgSizeX;
          x = (p % quadImgSizeXY) % quadImgSizeX;
          origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
          quadValues[indexQV] = img->val[origp][0];
          indexQV++;
        }
        quadValuesSize = indexQV;
        // Compute entropy
        quadEntropy[indexQuad] = iftNormalizedShannonEntropy(quadValues, quadValuesSize);
        totalEntropy += quadEntropy[indexQuad];

        indexQuad++;
        iftFree(quadValues);
      }
    }
  }

  indexQuad = 0;
  meanEntropy = totalEntropy / (float)(nquad*nquad*nquadz);
  sdEntropy = 0.0;
  for(k=0; k<nquadz; k++) {
    for (i = 0; i < nquad; i++) {
      for (j = 0; j < nquad; j++) {
        sdEntropy = (quadEntropy[indexQuad] - meanEntropy) * (quadEntropy[indexQuad] - meanEntropy);
        indexQuad++;
      }
    }
  }
  sdEntropy /= (float)(nquad*nquad*nquadz-1);
  sdEntropy = sqrtf(sdEntropy);
  thrEntropy = meanEntropy + sdEntropy;

  indexQuad = 0;
  for(k=0; k<nquadz; k++) {
    for (i = 0; i < nquad; i++) {
      for (j = 0; j < nquad; j++) {
        if (totalEntropy == 0)
          quadNSamples[indexQuad] = iftRound(((float) nsamples / (float) (nquad * nquad * nquadz) ));
        else
          quadNSamples[indexQuad] = iftRound((quadEntropy[indexQuad] / totalEntropy) * ((float) nsamples));
        if (quadNSamples[indexQuad] == 0)
          quadNSamples[indexQuad] = 1;

        // Compute init and end coordinates
        initZ = k * (img->zsize / nquadz);
        initY = i * (img->ysize / nquad);
        initX = j * (img->xsize / nquad);
        endZ = (k + 1) * (img->zsize / nquadz);
        endY = (i + 1) * (img->ysize / nquad);
        endX = (j + 1) * (img->xsize / nquad);
        if (k == (nquadz - 1))
          endZ = img->zsize;
        if (img->zsize == 1)
          endZ = img->zsize;
        if (i == (nquad - 1))
          endY = img->ysize;
        if (j == (nquad - 1))
          endX = img->xsize;

        // Divide image in quadrants
        quadImgSizeY = (endY - initY);
        quadImgSizeX = (endX - initX);
        quadImgSizeXY = quadImgSizeX * quadImgSizeY;
        quadImg = iftCreateMImage(quadImgSizeX, quadImgSizeY, (endZ - initZ), img->m);
        quadMask = iftCreateImage(quadImgSizeX, quadImgSizeY, (endZ - initZ));
        for (p = 0; p < quadImg->n; p++) {
          z = p / quadImgSizeXY;
          y = (p % quadImgSizeXY) / quadImgSizeX;
          x = (p % quadImgSizeXY) % quadImgSizeX;
          origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
          for (t = 0; t < quadImg->m; t++) {
            quadImg->val[p][t] = img->val[origp][t];
          }
          if (mask->val[origp] != 0)
            quadMask->val[p] = 1;
        }

        if (quadEntropy[indexQuad] > thrEntropy) {
          // Execute Second Stage of mix sampling
          quadGrid = iftAltMixedSamplingSecondStage(quadImg, quadMask, quadNSamples[indexQuad]);
        } else {
          quadGrid = iftGridSampling(quadImg, quadMask, quadNSamples[indexQuad]);
        }

        // Write seeds in the final out image
        for (p = 0; p < quadGrid->n; p++) {
          if (quadGrid->val[p] > 0) {
            z = p / (quadGrid->xsize*quadGrid->ysize);
            y = (p % (quadGrid->xsize*quadGrid->ysize)) / quadGrid->xsize;
            x = (p % (quadGrid->xsize*quadGrid->ysize)) % quadGrid->xsize;
            origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
            out->val[origp] = 1;
          }
        }

        indexQuad++;
        iftDestroyMImage(&quadImg);
        iftDestroyImage(&quadMask);
        iftDestroyImage(&quadGrid);

      }
    }
  }
  iftFree(quadEntropy);
  iftFree(quadNSamples);

  return out;
}

iftMImage *iftGradientVector(iftImage *img, iftImage *mask, iftAdjRel *A)
{
  iftMImage *grad=NULL;
  int        p, q, i, b;
  iftVoxel   u, v;
  float      mag[A->n], pq[3][A->n];

  for (i=0; i < A->n; i++)
    mag[i] = sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i]);

  for (i=0; i < A->n; i++) {
    pq[0][i] = A->dx[i]/mag[i];
    pq[1][i] = A->dy[i]/mag[i];
    pq[2][i] = A->dz[i]/mag[i];
  }

  if (iftIs3DImage(img)){
    grad=iftCreateMImage(img->xsize,img->ysize,img->zsize,3);
  }else{
    grad=iftCreateMImage(img->xsize,img->ysize,img->zsize,2);
  }

  for (p=0; p < img->n; p++) {
    if (mask==NULL || mask->val[p]) {//if no mask is passed, all the pixels are considered
      for (b=0; b < grad->m; b++)
  grad->val[p][b]=0;
      u = iftGetVoxelCoord(img,p);
      for (i=1; i < A->n; i++) {
  v = iftGetAdjacentVoxel(A,u,i);
  if (iftValidVoxel(img,v)){
    q = iftGetVoxelIndex(img,v);
    for (b=0; b < grad->m; b++)
      grad->val[p][b] += (img->val[q]-img->val[p])*pq[b][i];
  }
      }
    }
  }

  return(grad);
}


// If band is negative then the function will search for the maximum value among all bands
float iftMMaximumValue(const iftMImage *img, int band) {
  int b, i;
  float max_val = IFT_INFINITY_FLT_NEG;

  if(band < 0) {
    for(b = 0; b < img->m; b++) {
      for(i = 0; i < img->n; i++) {
  max_val = iftMax(img->val[i][b], max_val);
      }
    }
  } else {
    for(i = 0; i < img->n; i++) {
      max_val = iftMax(img->val[i][band], max_val);
    }
  }

  return max_val;
}

// If band is negative then the function will search for the minimum value among all bands
float iftMMinimumValue(iftMImage *img, int band) {
  int b, i;
  float min_val = IFT_INFINITY_FLT;

  if(band < 0) {
    for(b = 0; b < img->m; b++) {
      for(i = 0; i < img->n; i++) {
  min_val = iftMin(img->val[i][b], min_val);
      }
    }
  } else {
    for(i = 0; i < img->n; i++) {
      min_val = iftMin(img->val[i][band], min_val);
    }
  }

  return min_val;
}

iftImageTiles *iftMImageTilesByEquallyDividingAxes(const iftMImage *mimg, int ntiles_x, int ntiles_y, int ntiles_z) {
  iftBoundingBox bb;

  bb.begin.x = bb.begin.y = bb.begin.z = 0;
  bb.end.x = mimg->xsize - 1;
  bb.end.y = mimg->ysize - 1;
  bb.end.z = mimg->zsize - 1;

  return iftBoundingBoxImageTilesByEquallyDividingAxes(bb, ntiles_x, ntiles_y, ntiles_z);
}

iftImage *iftMImageTilesToLabelImage(iftMImage *mimg, int n_tiles){

  if (n_tiles <= 0)
      iftError("Number of tiles %d <= 0", "iftMImageTilesToLabelImage", n_tiles);

  if (mimg == NULL)
      iftError("Image is NULL", "iftMImageTilesToLabelImage");

  // tries to find regular tiles for the given number of tiles
  int block_vol = mimg->n / (1.0 * n_tiles);

  int n_blocks_x, n_blocks_y, n_blocks_z;
  iftNumberOfEquallyDimensionedTilesWithGivenMaximumSize(mimg->xsize, mimg->ysize, mimg->zsize,block_vol, &n_blocks_x, &n_blocks_y, &n_blocks_z);

  iftImageTiles *tiles = iftMImageTilesByEquallyDividingAxes(mimg, n_blocks_x, n_blocks_y, n_blocks_z);

  iftImage *out_img=iftCreateImage(mimg->xsize,mimg->ysize,mimg->zsize);

  #pragma omp parallel for
  for (int b = 0; b < tiles->ntiles; b++) {
    iftBoundingBox bb     = tiles->tile_coords[b];
    iftSetImageBoundingBoxValue(out_img,bb,b+1);
  }

  iftDestroyImageTiles(&tiles);

  return out_img;

}


float iftMImageSqDist(const iftMImage *mimg, int p, int q) {
    float dist = 0.0;

    // for each band from the MImage
    for (int b = 0; b < mimg->m; b++)
        dist += (mimg->val[q][b] - mimg->val[p][b]) * (mimg->val[q][b] - mimg->val[p][b]);

    return dist;
}


float iftMImageDist(const iftMImage *mimg, int p, int q) {
    float dist = 0.0;

    // for each band from the MImage
    for (int b = 0; b < mimg->m; b++)
        dist += (mimg->val[q][b] - mimg->val[p][b]) * (mimg->val[q][b] - mimg->val[p][b]);

    return sqrtf(dist);
}

int iftMImageDepth(iftMImage* mimg) {
    int max = IFT_INFINITY_INT_NEG;
    for(int b=0; b<mimg->m; b++) {
      int m = iftMMaximumValue(mimg, b);
      if(m > max)
        max = m;
    }
    max = iftNormalizationValue(max) + 1;
    double depth = iftLog(max, 2);
    return depth;
}



iftMImage *iftExtendMImageByVoxelCoord(const iftMImage *mimg, const iftImage *label_img, const int n_adj_voxels, bool normalize) {
    float norm_val = 1.0f;

    iftMImage *emimg;

    if (normalize)
        norm_val = iftMax(iftMax(mimg->xsize, mimg->ysize), mimg->zsize);

    if (iftIs3DMImage(mimg))
      emimg = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m + 3*n_adj_voxels);
    else
      emimg = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m + 2*n_adj_voxels);

    if (label_img == NULL) {
        #pragma omp parallel for
        for (int p = 0; p < mimg->n; p++)
    iftExtendMVoxelFeatVector(emimg, mimg, p, n_adj_voxels, norm_val);
    }
    else {
        #pragma omp parallel for
        for (int p = 0; p < mimg->n; p++)
            if (label_img->val[p] != 0)
        iftExtendMVoxelFeatVector(emimg, mimg, p, n_adj_voxels, norm_val);
    }

    return emimg;
}


iftMImage *iftExtractImageFeatures(const iftImage *img, const iftImage *label_img, const iftAdjRel *A, bool voxel_coord) {
    if ((label_img != NULL) && !iftIsDomainEqual(img, label_img))
        iftError("Image and Label Image have different domains: (%d, %d, %d) != (%d, %d, %d", "iftExtractImageFeatures",
                 img->xsize, img->ysize, img->zsize, label_img->xsize, label_img->ysize, label_img->zsize);

    int n_adj_voxels = 1;
    iftMImage *feat = NULL;
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));

    // only the own voxel is considered in feat. extraction
    if (iftIsColorImage(img)) {
        int n_channels = 3;
        int n_bands = 3;

        if (A == NULL) {
            feat = iftCreateMImage(img->xsize, img->ysize, img->zsize, n_bands);

            if (label_img == NULL) {
                #pragma omp parallel for
                for (int p = 0; p < img->n; p++)
                    iftCopyVoxelToMImageLAB2Norm(feat, p, img, p, 0, normalization_value);
            } else {
                #pragma omp parallel for
                for (int p = 0; p < img->n; p++)
                    if (label_img->val[p] != 0)
                        iftCopyVoxelToMImageLAB2Norm(feat, p, img, p, 0, normalization_value);
            }
        } else {
            // the own voxel and its neighbors defined by the adjacency are considered

            n_bands *= A->n; // one band for each channel from each adjacent pixel
            n_adj_voxels = A->n;
            feat = iftCreateMImage(img->xsize, img->ysize, img->zsize, n_bands);

            if (label_img == NULL) {
                #pragma omp parallel for
                for (int p = 0; p < img->n; p++) {
                    iftVoxel u = iftGetVoxelCoord(img, p);

                    for (int i = 0; i < A->n; i++) {
                        iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                        if (iftValidVoxel(img, v)) {
                            int q = iftGetVoxelIndex(img, v);
                            iftCopyVoxelToMImageLAB2Norm(feat, p, img, q, i * n_channels, normalization_value);
                        }
                    }
                }
            } else {
                #pragma omp parallel for
                for (int p = 0; p < img->n; p++) {
                    if (label_img->val[p] != 0) {
                        iftVoxel u = iftGetVoxelCoord(img, p);

                        for (int i = 0; i < A->n; i++) {
                            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                            if (iftValidVoxel(img, v)) {
                                int q = iftGetVoxelIndex(img, v);
                                iftCopyVoxelToMImageLAB2Norm(feat, p, img, q, i * n_channels, normalization_value);
                            }
                        }
                    }
                }
            }
        }
    } else {
        int n_channels = 1;
        int n_bands = 1;

        if (A == NULL) {
            feat = iftCreateMImage(img->xsize, img->ysize, img->zsize, n_bands);

            if (label_img == NULL) {
                #pragma omp parallel for
                for (int p = 0; p < img->n; p++)
                    iftCopyVoxelToMImageGRAYNorm(feat, p, img, p, 0, normalization_value);
            } else {
                #pragma omp parallel for
                for (int p = 0; p < img->n; p++)
                    if (label_img->val[p] != 0)
                        iftCopyVoxelToMImageGRAYNorm(feat, p, img, p, 0, normalization_value);
            }
        } else {
            // the own voxel and its neighbors defined by the adjacency are considered

            n_bands *= A->n; // one band for each channel from each adjacent pixel
            n_adj_voxels = A->n;
            feat = iftCreateMImage(img->xsize, img->ysize, img->zsize, n_bands);

            if (label_img == NULL) {
                #pragma omp parallel for
                for (int p = 0; p < img->n; p++) {
                    iftVoxel u = iftGetVoxelCoord(img, p);

                    for (int i = 0; i < A->n; i++) {
                        iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                        if (iftValidVoxel(img, v)) {
                            int q = iftGetVoxelIndex(img, v);
                            iftCopyVoxelToMImageGRAYNorm(feat, p, img, q, i * n_channels, normalization_value);
                        }
                    }
                }
            } else {
                #pragma omp parallel for
                for (int p = 0; p < img->n; p++) {
                    if (label_img->val[p] != 0) {
                        iftVoxel u = iftGetVoxelCoord(img, p);

                        for (int i = 0; i < A->n; i++) {
                            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                            if (iftValidVoxel(img, v)) {
                                int q = iftGetVoxelIndex(img, v);
                                iftCopyVoxelToMImageGRAYNorm(feat, p, img, q, i * n_channels, normalization_value);
                            }
                        }
                    }
                }
            }
        }
    }

    // adding the voxel coordinates
    if (voxel_coord) {
        iftMImage *aux = feat;
        feat = iftExtendMImageByVoxelCoord(feat, label_img, n_adj_voxels, true);
        iftDestroyMImage(&aux);
    }

    return feat;
}

iftMImage *iftTransformMImage(const iftMImage *input, const double *L, int L_nrow, int L_ncol)
{
    if (input->m != L_ncol)
        iftError("Cannot transform space, number of bands different from number os columns of L", "iftTransformMImage");

    iftMImage *output = iftCreateMImage(input->xsize, input->ysize, input->zsize, L_nrow);
    iftFree(output->data->val);

    output->data->val = iftSpaceTransformFloat(input->data->val, L, input->n, L_ncol, L_nrow);

    for (int i = 0; i < output->n; i++)
        output->val[i] = iftMatrixRowPointer(output->data, i);

    return output;
}


iftMImage *iftKernelTransformMImage(const iftMImage *input, const iftDoubleMatrix *train_data,
                                    const iftDoubleMatrix *omega, iftKernelFunction *K,
                                    double param1, double param2, double scale)
{
    if (train_data->ncols != input->m) {
        iftError("Input MImage must have the number of bands equal to train data number of features",
                "iftKernelTransformMImage");
    }

    if (omega->ncols != train_data->nrows) {
        iftError("Omega must have the number of columns equal to the number of train data samples",
                "iftKernelTransformMImage");
    }

    double *img_feat = iftAlloc(input->n * input->m, sizeof *img_feat);

    for (int i = 0; i < input->n; i++) {
        int row = i * input->m;
        for (int j = 0; j < input->m; j++) {
            img_feat[row + j] = input->val[i][j];
        }
    }

    double *space_feats = iftKernelFeatures(train_data->val, img_feat, train_data->nrows,
                                            input->n, train_data->ncols, K, param1, param2);

    iftFree(img_feat);

    double *Ldata = iftSpaceTransform(space_feats, omega->val, input->n, train_data->nrows, omega->nrows);

    iftFree(space_feats);

    iftMImage *out = iftCreateMImage(input->xsize, input->ysize, input->zsize, omega->nrows);

    for (int i = 0; i < out->n; i++) {
        int row = i * out->m;
        for (int j = 0; j < out->m; j++) {
            out->val[i][j] = (float) (Ldata[row + j] * scale);
        }
    }

    iftFree(Ldata);

    return out;
}


iftMImage *iftColorQuantization(const iftMImage *input, float lower_bound, float upper_bound, int bins_per_band)
{
    iftMImage *output = iftCreateMImage(input->xsize, input->ysize, input->zsize, input->m);

    float bin_width = (upper_bound - lower_bound) / bins_per_band;
    float *hist = iftAlloc(bins_per_band, sizeof *hist);

    hist[0] = lower_bound + (bin_width / 2);
    for (int i = 1; i < bins_per_band; i++) {
        hist[i] = hist[i - 1] + bin_width;
    }

    for (int b = 0; b < input->m; b++) {
        for (int p = 0; p < input->n; p++) {
            int idx = (int) ((input->val[p][b] - lower_bound) / bin_width);
            output->val[p][b] = hist[idx];
        }
    }

    iftFree(hist);

    return output;
}

iftMImage *iftMConvertColorSpace(iftMImage* mimg, iftColorSpace colSpaceIn, iftColorSpace colSpaceOut) {
    int normValue = iftNormalizationValue(iftMMaximumValue(mimg, -1));

    if(mimg->m != 1 && mimg->m != 3)
        iftError("Conversion between color spaces  can only be performed with mimages with 1 or 3 bands. The given mimage has %d bands",
            "iftMConvertColorSpace", mimg->m);

    /* validate the color spaces */
    bool isColorImgIn = (mimg->m == 3);

    if((isColorImgIn && (colSpaceIn == GRAY_CSPACE || colSpaceIn == GRAYNorm_CSPACE)) ||
        (!isColorImgIn && (colSpaceIn != GRAY_CSPACE && colSpaceIn != GRAYNorm_CSPACE))) {
        iftError("The number of image channels (%d) and the specified input color space does not match",
            "iftMConvertColorSpace", mimg->m);
    }

    if(!isColorImgIn && colSpaceOut != GRAY_CSPACE && colSpaceOut != GRAYNorm_CSPACE) {
        iftError("Color image representation can only be used with color images. The given image is grayscale.",
            "iftMConvertColorSpace");
    }

    if(colSpaceIn == colSpaceOut)
        return iftCopyMImage(mimg);

    /* create a new mimage whose number of bands depends on the chosen color space */
    iftMImage *mimgOut = NULL;
    if(colSpaceOut == GRAY_CSPACE || colSpaceOut == GRAYNorm_CSPACE)
        mimgOut = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, 1);
    else
        mimgOut = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, 3);
    
    /* convert the color space of each pixel */
    for(int p = 0; p < mimg->n; p++) {
        iftFColor colorIn, colorOut;
        for(int b = 0; b < mimg->n; b++)
            colorIn.val[b] = mimg->val[p][b];

        colorOut = iftConvertPixelColorSpace(colorIn, colSpaceIn, colSpaceOut, normValue);

        for(int b = 0; b < mimg->n; b++)
            mimgOut->val[p][b] = colorOut.val[b];
    }

    return mimgOut;
}

iftMImage *iftMImageCentralize(const iftMImage *mimg)
{
    double *mean = iftAlloc((size_t) mimg->m, sizeof *mean);

    for (int b = 0; b < mimg->m; b++) {
        for (int i = 0; i < mimg->n; i++) {
            mean[b] += mimg->val[i][b];
        }
        mean[b] /= mimg->n;
    }

    iftMImage *out = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m);

    for (int b = 0; b < mimg->m; b++) {
        for (int i = 0; i < mimg->n; i++) {
            out->val[i][b] = (float) (mimg->val[i][b] - mean[b]);
        }
    }

    return out;
}


iftMImage *iftMGaussianFilter(const iftMImage *mimg, float radius, float stddev)
{
    iftAdjRel *A = iftCircular(radius);
    float *w = iftAlloc(A->n, sizeof *w);
    float var = 2 * stddev * stddev, sum = 0.0f;

    for (int i = 0; i < A->n; i++) {
        float dist = A->dx[i] * A->dx[i] +
                     A->dy[i] * A->dy[i] + A
                     ->dz[i] * A->dz[i];
        w[i] = expf (- dist / var);
        sum += w[i];
    }

    for (int i = 0; i < A->n; i++) {
        w[i] /= sum;
    }

    iftMImage *out = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m);

    #pragma omp parallel for
    for (int b = 0; b < mimg->m; b++) {
        for (int p = 0; p < mimg->n; p++) {
            iftVoxel u = iftMGetVoxelCoord(mimg, p);
            out->val[p][b] = 0;
            for (int i = 0; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);
                if (!iftValidVoxel(mimg, v))
                    continue;
                int q = iftGetVoxelIndex(mimg, v);
                out->val[p][b] += w[i] * mimg->val[q][b];
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftFree(w);

    return out;
}


iftMImage *iftBinarySetsNCA(const iftMImage *mimg, const iftSet *obj, const iftSet *bkg,
                            int d_out, int n_iters, double learn_rate)
{
    int n = iftSetSize(obj) + iftSetSize(bkg);
    double *X = iftAlloc(n * mimg->m, sizeof *X);
    int *y = iftAlloc(n, sizeof *y);

    int i = 0;
    for (const iftSet *s = obj; s; s = s->next, i++) {
        int row = i * mimg->m;
        y[i] = 1;
        for (int j = 0; j < mimg->m; j++)
            X[row + j] = mimg->val[s->elem][j];
    }

    for (const iftSet *s = bkg; s; s = s->next, i++) {
        int row = i * mimg->m;
        y[i] = 0;
        for (int j = 0; j < mimg->m; j++)
            X[row + j] = mimg->val[s->elem][j];
    }

    double *L = iftNeighCompAnalysis(X, y, NULL, n, mimg->m, d_out, n_iters, learn_rate);

    iftFree(X);
    iftFree(y);

    iftMImage *Lmimg = iftTransformMImage(mimg, L, d_out, mimg->m);

    iftFree(L);

    return Lmimg;
}


void iftRescaleMImage(iftMImage *mimg, bool single_max)
{

    float *min = iftAlloc(mimg->m, sizeof *min);
    float *max = iftAlloc(mimg->m, sizeof *max);

    for (int i = 0; i < mimg->m; i++) {
        min[i] = IFT_INFINITY_FLT;
        max[i] = IFT_INFINITY_FLT_NEG;
    }

    for (int i = 0; i < mimg->n; i++) {
        for (int j = 0; j < mimg->m; j++) {
            if (mimg->val[i][j] < min[j]) {
                min[j] = mimg->val[i][j];
            }
            if (mimg->val[i][j] > max[j]) {
                max[j] = mimg->val[i][j];
            }
        }
    }

    if (single_max)
    {
        float min_min = min[0];
        float max_max = max[0];
        for (int i = 1; i < mimg->m; i++) {
            min_min = iftMin(min_min, min[i]);
            max_max = iftMax(max_max, max[i]);
        }

        for (int i = 0; i < mimg->m; i++) {
            max[i] = max_max;
            min[i] = min_min;
        }
    }


    for (int j = 0; j < mimg->m; j++)
        max[j] -= min[j];

    for (int i = 0; i < mimg->n; i++) {
        for (int j = 0; j < mimg->m; j++) {
            mimg->val[i][j] /= max[j];
        }
    }

    iftFree(min);
    iftFree(max);
}



iftMImage *iftStackGrayImages(int n, ...) {
    if (n <= 0)
        iftError("Invalid number of Gray Images %d <= 0", "iftConcatStrings", n);

    va_list img_list;
    va_start(img_list, n);
    const iftImage *img0 = va_arg(img_list, iftImage *);
    iftMImage *mimg = iftCreateMImage(img0->xsize, img0->ysize, img0->zsize, n);
    va_end(img_list);

    va_start(img_list, n);

    #pragma omp parallel for
    for (int b = 0; b < n; b++) {
        const iftImage *img = va_arg(img_list, iftImage *);

        if (!iftIsVoxelSizeEqual(img0, img)) {
            iftError("Different voxel size (dx, dy, dz)\nimg0: %.2f, %.2f, %.2f\nimg%d: %.2f, %.2f, %.2f\n",
                     "iftStackGrayImages", img0->dx, img0->dy, img0->dz, b, img->dx, img->dy, img->dz);
        }
        if (!iftIsDomainEqual(img0, img)) {
            iftError("Different domains (xsize, ysize, zsize)\nimg0: %d, %d, %d\nimg%d: %d, %d, %d\n",
                     "iftStackGrayImages", img0->xsize, img0->ysize, img0->zsize, b, img->xsize, img->ysize, img->zsize);
        }


        for (int p = 0; p < mimg->n; p++) {
            mimg->val[p][b] = img->val[p];
        }
    }
    va_end(img_list);
    
    return mimg;
}

int iftMGetVoxelIndex_pyift(iftMImage *s, iftVoxel v) {
  return ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)]);
}
