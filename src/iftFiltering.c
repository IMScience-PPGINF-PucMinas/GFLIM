#include "iftFiltering.h"

#include "ift/core/io/Stream.h"
#include "iftInterpolation.h"



iftImage  *iftLinearFilter(const iftImage *img, iftKernel *K)
{
  iftImage *fimg=iftCopyImage(img);
  iftAdjRel *A=K->A;

  if (iftIsColorImage(img)){
#pragma omp parallel for shared(img,fimg,K,A)
    for (int p=0; p < img->n; p++) {
      iftVoxel u   = iftGetVoxelCoord(img,p);
      float val1 = 0.0, val2 = 0.0, val3 = 0.0;
      for (int i=0; i < A->n; i++){
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(img,v)){
	  int q = iftGetVoxelIndex(img,v);
	  val1 += ((float) img->val[q] * K->weight[i]);
	  val2 += ((float) img->Cb[q]  * K->weight[i]);
	  val3 += ((float) img->Cr[q]  * K->weight[i]);
	}
      }
      fimg->val[p] = iftRound(val1);
      fimg->Cb[p]  = iftRound(val2);
      fimg->Cr[p]  = iftRound(val3);
    }
  }else{
#pragma omp parallel for shared(img,fimg,K,A)
    for (int p=0; p < img->n; p++) {
      iftVoxel u   = iftGetVoxelCoord(img,p);
      float val = 0.0;
      for (int i=0; i < A->n; i++){
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(img,v)){
	  int q = iftGetVoxelIndex(img,v);
	  val += ((float) img->val[q] * K->weight[i]);
	}
      }
      fimg->val[p] = iftRound(val);
    }
  }

  return(fimg);
}


iftImage  *iftLinearFilterInRegion(iftImage *img, iftImage *mask, iftKernel *K)
{
  iftImage *fimg=iftCopyImage(img);
  iftAdjRel *A=K->A;

  if (iftIsColorImage(img)){
#pragma omp parallel for shared(img,fimg,K,A)
    for (int p=0; p < img->n; p++) {
      if (mask->val[p]!=0) {
	iftVoxel u   = iftGetVoxelCoord(img,p);
	float val1 = 0.0, val2 = 0.0, val3 = 0.0;
	for (int i=0; i < A->n; i++){
	  iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(img,v)){
	    int q = iftGetVoxelIndex(img,v);
	    val1 += ((float) img->val[q] * K->weight[i]);
	    val2 += ((float) img->Cb[q]  * K->weight[i]);
	    val3 += ((float) img->Cr[q]  * K->weight[i]);
	  }
	}
	fimg->val[p] = iftRound(val1);
	fimg->Cb[p]  = iftRound(val2);
	fimg->Cr[p]  = iftRound(val3);
      }
    }
  }else{
#pragma omp parallel for shared(img,fimg,K,A)
    for (int p=0; p < img->n; p++) {
      if (mask->val[p]!=0) {
	iftVoxel u   = iftGetVoxelCoord(img,p);
	float val = 0.0;
	for (int i=0; i < A->n; i++){
	  iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(img,v)){
	    int q = iftGetVoxelIndex(img,v);
	    val += ((float) img->val[q] * K->weight[i]);
	  }
	}
	fimg->val[p] = iftRound(val);
      }
    }
  }
  return(fimg);
}

iftImage *iftFastLinearFilter(iftImage *img, iftKernel *K, char crop)
{
  iftMatrix     *A, *B, *M;
  iftImage      *fimg;
  iftFastAdjRel *F = iftCreateFastAdjRel(K->A,img->tby,img->tbz);
  int            xsize, ysize, zsize;
  iftImage      *gray;    

  if (crop){
    xsize = img->xsize-2*F->bx;
    ysize = img->ysize-2*F->by;
    zsize = img->zsize-2*F->bz;
  } else {
    xsize = img->xsize;
    ysize = img->ysize;
    zsize = img->zsize;
  }

  A    = iftKernelToMatrix(K);

  gray = iftImageGray(img);
  if (crop)
    B    = iftImageToMatrix(gray,F,crop);
  else
    B    = iftImageAdjToMatrix(gray,K->A);
  M    = iftMultMatrices(A,B);
  iftDestroyImage(&gray);
  gray = iftMatrixToImage(M,xsize,ysize,zsize);
  iftDestroyMatrix(&M);  
  iftDestroyMatrix(&B);


  if (iftIsColorImage(img)) {
    fimg = iftCreateColorImage(xsize,ysize,zsize,iftImageDepth(img));
    if (crop) {
      iftImage *aux = iftCropImage(img,F->bx,F->by,F->bz);

     for (int p=0; p < fimg->n; p++) {
	fimg->val[p] = gray->val[p];
	fimg->Cb[p]  = aux->Cb[p];
	fimg->Cr[p]  = aux->Cr[p];	
      }
      iftDestroyImage(&aux);
    } else {

      for (int p=0; p < fimg->n; p++) {
	fimg->val[p] = gray->val[p];
	fimg->Cb[p]  = img->Cb[p];
	fimg->Cr[p]  = img->Cr[p];	
      }
    }
  } else {
    fimg = iftCopyImage(gray);
  }

  iftDestroyMatrix(&A);
  iftDestroyFastAdjRel(&F);
  iftDestroyImage(&gray);

  if (iftIs3DImage(img)){
    iftCopyVoxelSize(img,fimg);
  }
  return(fimg);
}

iftFImage  *iftFLinearFilter(iftFImage *img, iftKernel *K)
{
  iftFImage *fimg=iftCreateFImage(img->xsize,img->ysize,img->zsize);
  iftAdjRel *A=K->A;

#pragma omp parallel for shared(img,fimg,K,A)
  for (int p=0; p < img->n; p++) {
    iftVoxel u   = iftFGetVoxelCoord(img,p);
    fimg->val[p] = 0.0;
    for (int i=0; i < A->n; i++){
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftFValidVoxel(img,v)){
	int q = iftFGetVoxelIndex(img,v);
	fimg->val[p] += ((float) img->val[q] * K->weight[i]);
      }
    }
  }

  return(fimg);
}

iftFImage  *iftFLinearFilterInRegion(iftFImage *img, iftImage *mask, iftKernel *K)
{
  iftFImage *fimg=iftFCopyImage(img);
  iftAdjRel *A=K->A;

#pragma omp parallel for shared(img,mask,fimg,K,A)
  for (int p=0; p < img->n; p++) {
    if (mask->val[p]!=0){
      iftVoxel u   = iftFGetVoxelCoord(img,p);
      fimg->val[p] = 0.0;
      for (int i=0; i < A->n; i++){
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftFValidVoxel(img,v)){
	  int q = iftFGetVoxelIndex(img,v);
	  fimg->val[p] += ((float) img->val[q] * K->weight[i]);
	}
      }
    }
  }

  return(fimg);
}

iftFImage  *iftFastFLinearFilter(iftFImage *img, iftKernel *K, char crop)
{
  iftMatrix *A, *B, *M;
  iftFImage   *fimg;
  iftFastAdjRel *F = iftCreateFastAdjRel(K->A,img->tby,img->tbz);


  A    = iftKernelToMatrix(K);
  B    = iftFImageToMatrix(img,F,crop);
  M    = iftMultMatrices(A,B);

  if (crop)
    fimg = iftMatrixToFImage(M,img->xsize-2*F->bx,img->ysize-2*F->by,img->zsize-2*F->bz);
  else
    fimg = iftMatrixToFImage(M,img->xsize,img->ysize,img->zsize);

  iftDestroyMatrix(&A);
  iftDestroyMatrix(&B);
  iftDestroyMatrix(&M);
  iftDestroyFastAdjRel(&F);

  if (iftIs3DFImage(img))
    iftFCopyVoxelSize(img,fimg);

  return(fimg);
}


iftMImage *iftMLinearFilter(iftMImage *img, iftMKernel *K)
{
  iftMatrix  *A, *B, *M;
  iftMImage   *fimg;
  iftFastAdjRel *F = iftCreateFastAdjRel(K->A,img->tby,img->tbz);

  if (img->m != K->nbands)
      iftError("Incompatible input", "iftMLinearFilter");

  A    = iftMKernelToMatrix(K);
  B    = iftMImageToMatrix(img,F);
  M    = iftMultMatrices(A,B);

  fimg = iftMatrixToMImage(M,img->xsize-2*F->bx,img->ysize-2*F->by,img->zsize-2*F->bz,1,'r');

  iftDestroyMatrix(&A);
  iftDestroyMatrix(&B);
  iftDestroyMatrix(&M);
  iftDestroyFastAdjRel(&F);

  if (iftIs3DMImage(img))
    iftMCopyVoxelSize(img,fimg);

  return(fimg);
}

iftMImage *iftMMLinearFilter(iftMImage *img, iftMMKernel *k_bank)
{
  iftMatrix  *A, *B, *M;
  iftMImage   *fimg;
  iftFastAdjRel *F = iftCreateFastAdjRel(k_bank->A,img->tby,img->tbz);

  if (img->m != k_bank->nbands)
      iftError("Incompatible input", "iftMMLinearFilter");

  A    = iftMMKernelToMatrix(k_bank);
  B    = iftMImageToMatrix(img,F);
  M    = iftMultMatrices(A,B);

  fimg = iftMatrixToMImage(M,img->xsize-2*F->bx,img->ysize-2*F->by,img->zsize-2*F->bz,k_bank->nkernels,'r');

  iftDestroyMatrix(&A);
  iftDestroyMatrix(&B);
  iftDestroyMatrix(&M);
  iftDestroyFastAdjRel(&F);

  if (iftIs3DMImage(img))
    iftMCopyVoxelSize(img,fimg);

  return(fimg);
}

iftMatrix *iftImageToMatrix(iftImage *img, iftFastAdjRel *F, char crop)
{
  iftMatrix *M;

  if (crop){
    M = iftCreateMatrix((img->xsize-2*F->bx)*(img->ysize-2*F->by)*(img->zsize-2*F->bz),F->n);
#pragma omp parallel for shared(img,M,F)
    for (int i=0; i < F->n; i++) {
      iftVoxel u;
      int c=0;
      for (u.z = F->bz; u.z < img->zsize-F->bz; u.z++)
	for (u.y = F->by; u.y < img->ysize-F->by; u.y++)
	  for (u.x = F->bx; u.x < img->xsize-F->bx; u.x++){
	    int p = iftGetVoxelIndex(img,u);
	    int q = p + F->dq[i];
	    M->val[iftGetMatrixIndex(M,c,i)] = (double) img->val[q];
	    c++;
	  }
    }
  }else{
    M = iftCreateMatrix(img->n,F->n);
#pragma omp parallel for shared(img,M,F)
    for (int i=0; i < F->n; i++) {
      iftVoxel u;
      for (u.z = F->bz; u.z < img->zsize-F->bz; u.z++)
	for (u.y = F->by; u.y < img->ysize-F->by; u.y++)
	  for (u.x = F->bx; u.x < img->xsize-F->bx; u.x++){
	    int p = iftGetVoxelIndex(img,u);
	    int q = p + F->dq[i];
	    M->val[iftGetMatrixIndex(M,p,i)] = (double) img->val[q];
	  }
    }
  }

  return(M);
}

iftMatrix *iftImageAdjToMatrix(iftImage *img, iftAdjRel *A)
{
  iftMatrix *M;

  M = iftCreateMatrix(img->n,A->n);
#pragma omp parallel for shared(img,M,A)
    for (int i=0; i < A->n; i++) {
      iftVoxel u;
      for (u.z = 0; u.z < img->zsize; u.z++)
	for (u.y = 0; u.y < img->ysize; u.y++)
	  for (u.x = 0; u.x < img->xsize; u.x++){
	    int p = iftGetVoxelIndex(img,u);
	    iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	    if (iftValidVoxel(img,v)){
	      int q = iftGetVoxelIndex(img,v);
	      M->val[iftGetMatrixIndex(M,p,i)] = (double) img->val[q];
	    }
	  }
    }

  return(M);
}

iftImage   *iftMatrixToImage(iftMatrix *M, int xsize, int ysize, int zsize)
{
  iftImage *img = iftCreateImage(xsize,ysize,zsize);

  if ( (xsize*ysize*zsize) != M->n )
      iftError("Incompatible image dimensions", "iftMatrixToImage");

  if (M->nrows != 1)
      iftError("Invalid matrix for convertion", "iftMatrixToImage");

  for (int c=0; c < M->ncols; c++) {
    img->val[c] = (int)(M->val[c] + 0.5);
  }

  return(img);
}

iftMatrix *iftFImageToMatrix(iftFImage *img, iftFastAdjRel *F, char crop)
{
  iftMatrix *M;

  if (crop){
    M = iftCreateMatrix((img->xsize-2*F->bx)*(img->ysize-2*F->by)*(img->zsize-2*F->bz),F->n);
#pragma omp parallel for shared(img,M,F)
    for (int i=0; i < F->n; i++) {
      iftVoxel u;
      int c=0;
      for (u.z = F->bz; u.z < img->zsize-F->bz; u.z++)
	for (u.y = F->by; u.y < img->ysize-F->by; u.y++)
	  for (u.x = F->bx; u.x < img->xsize-F->bx; u.x++){
	    int p = iftFGetVoxelIndex(img,u);
	    int q = p + F->dq[i];
	    M->val[iftGetMatrixIndex(M,c,i)] = (double) img->val[q];
	    c++;
	  }
    }
  }else{
    M = iftCreateMatrix(img->n,F->n);
#pragma omp parallel for shared(img,M,F)
    for (int i=0; i < F->n; i++) {
      iftVoxel u;
      for (u.z = F->bz; u.z < img->zsize-F->bz; u.z++)
	for (u.y = F->by; u.y < img->ysize-F->by; u.y++)
	  for (u.x = F->bx; u.x < img->xsize-F->bx; u.x++){
	    int p = iftFGetVoxelIndex(img,u);
	    int q = p + F->dq[i];
	    M->val[iftGetMatrixIndex(M,p,i)] = (double) img->val[q];
	  }
    }
  }

  return(M);
}

iftFImage   *iftMatrixToFImage(iftMatrix *M, int xsize, int ysize, int zsize)
{
  iftFImage *img = iftCreateFImage(xsize,ysize,zsize);

  if ( (xsize*ysize*zsize) != M->n )
      iftError("Incompatible image dimensions", "iftMatrixToFImage");

  if (M->nrows != 1)
      iftError("Invalid matrix for convertion", "iftMatrixToFImage");

  for (int c=0; c < M->ncols; c++) {
    img->val[c] = (float)M->val[c];
  }

  return(img);
}

iftMatrix *iftMImageToMatrix(iftMImage *img, iftFastAdjRel *F)
{
  iftMatrix *M;

  M = iftCreateMatrix((img->xsize-2*F->bx)*(img->ysize-2*F->by)*(img->zsize-2*F->bz),F->n*img->m);
#pragma omp parallel for shared(img,M,F)
  for (int i=0; i < F->n; i++) {
    iftVoxel u;
    for (int m=0; m < img->m; m++) {
      int r = i+m*F->n;
      int c = 0;
      for (u.z = F->bz; u.z < img->zsize-F->bz; u.z++)
	for (u.y = F->by; u.y < img->ysize-F->by; u.y++)
	  for (u.x = F->bx; u.x < img->xsize-F->bx; u.x++){
	    int p = iftGetVoxelIndex(img,u);
	    int q = p + F->dq[i];
	    M->val[iftGetMatrixIndex(M,c,r)] = (double) img->val[q][m];
	    c++;
	  }
    }
  }

  return(M);
}

iftMImage *iftMatrixToMImage(iftMatrix *M, int xsize, int ysize, int zsize, int nbands, char band_orientation)
{
  if(((long)xsize*(long)ysize*(long)zsize*(long)nbands) != M->n)
        iftError("Incompatible image dimensions", "iftMatrixToMImage");

    if(band_orientation == 'r') {
        if(M->nrows != nbands)
            iftError("Invalid matrix for convertion", "iftMatrixToMImage");
    }
    else if(band_orientation == 'c') {
        if(M->ncols != nbands)
            iftError("Invalid matrix for convertion", "iftMatrixToMImage");
    }

    iftMImage *mimg = iftCreateMImage(xsize,ysize,zsize,nbands);

    if(band_orientation == 'r') {
        for(int r=0; r < M->nrows; r++) {
            for(int c=0; c < M->ncols; c++) {
                mimg->val[c][r] = M->val[iftGetMatrixIndex(M,c,r)];
            }
        }
    }
    else if(band_orientation == 'c') {
        for(int r=0; r < M->nrows; r++) {
            for(int c=0; c < M->ncols; c++) {
                mimg->val[r][c] = M->val[iftGetMatrixIndex(M,c,r)];
            }
        }
    }

    return mimg;
}

iftMImageArray *iftMatrixToMImageArray(iftMatrix *M, int nImgs, int xsize, int ysize, int zsize, int nbands, char band_orientation)
{
    if((nImgs*xsize*ysize*zsize*nbands) != M->n)
        iftError("Incompatible image dimensions", "iftMatrixToMImageArray");

    if(band_orientation == 'r') {
        if(M->nrows != nbands)
            iftError("Invalid matrix for convertion", "iftMatrixToMImageArray");
    }
    else if(band_orientation == 'c') {
        if(M->ncols != nbands)
            iftError("Invalid matrix for convertion", "iftMatrixToMImageArray");
    }

    iftMImageArray *mimgArray = iftCreateMImageArray(nImgs);

    for(int i = 0; i < nImgs; i++)
        mimgArray->val[i] = iftCreateMImage(xsize,ysize,zsize,nbands);

    if(band_orientation == 'r') {
        for(int c=0; c < M->ncols; c++) {
            int i = c / mimgArray->val[0]->n;
            int idx = c % mimgArray->val[0]->n;
            for(int r=0; r < M->nrows; r++) {
                mimgArray->val[i]->val[idx][r] = M->val[iftGetMatrixIndex(M,c,r)];
            }
        }
    }
    else if(band_orientation == 'c') {
        for(int r=0; r < M->nrows; r++) {
            int i = r / mimgArray->val[0]->n;
            int idx = r % mimgArray->val[0]->n;
            for(int c=0; c < M->ncols; c++) {
                mimgArray->val[i]->val[idx][c] = M->val[iftGetMatrixIndex(M,c,r)];
            }
        }
    }

    return mimgArray;
}


iftMatrix *iftMImageToFeatureMatrix(iftMImage *mimg, iftAdjRel *A, float *fill_band_with)
{
    iftMatrix *matrix = iftCreateMatrix(A->n*mimg->m, mimg->n);

    #pragma omp parallel for
    for (int p = 0; p < mimg->n; p++) {
        iftVoxel u = iftMGetVoxelCoord(mimg, p);
        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v)) {
                int q = iftMGetVoxelIndex(mimg, v);
                for (int b = 0; b < mimg->m; b++) {
                    iftMatrixElem(matrix, b + i * mimg->m, p) = mimg->val[q][b];
                }
            }else if(fill_band_with != NULL){
              for (int b=0; b< mimg->m; b++){
                iftMatrixElem(matrix, b + i*mimg->m, p) = fill_band_with[b]; // uma linha para cada pixel, 
              }
            } 
        }
    }

    return matrix;
}

iftMatrix *iftMImageArrayToFeatureMatrix(iftMImageArray *mimgArray, iftAdjRel *A)
{
    iftMatrix *matrix = iftCreateMatrix(A->n*mimgArray->val[0]->m, mimgArray->n*mimgArray->val[0]->n);

    for (int m = 0; m < mimgArray->n; m++) {
        iftMImage *mimg = mimgArray->val[m];
        #pragma omp parallel for
        for (int p = 0; p < mimg->n; p++) {
            iftVoxel u = iftMGetVoxelCoord(mimg, p);
            for (int i = 0; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);
                if (iftMValidVoxel(mimg, v)) {
                    int q = iftMGetVoxelIndex(mimg, v);
                    for (int b = 0; b < mimg->m; b++) {
                        iftMatrixElem(matrix, b + i * mimg->m, p + m * mimg->n) = mimg->val[q][b];
                    }
                }
            }
        }
    }

    return matrix;
}

iftMatrix *iftKernelToMatrix(iftKernel *K)
{
  iftMatrix *M = iftCreateMatrix(K->A->n,1);

  for (int c=0; c < K->A->n; c++)
    M->val[c] = (double) K->weight[c];

  return(M);
}

iftMatrix *iftMKernelToMatrix(iftMKernel *K)
{
  iftMatrix *M = iftCreateMatrix(K->A->n*K->nbands,1);

  for (int m=0; m < K->nbands; m++) {
    int j=m*K->A->n;
    for (int i=0; i < K->A->n; i++)
      M->val[i+j] = (double) K->weight[m].val[i];
  }

  return(M);
}

iftMatrix *iftMMKernelToMatrix(iftMMKernel *k_bank)
{
  iftMatrix *M = iftCreateMatrix(k_bank->A->n*k_bank->nbands,k_bank->nkernels);

  for (int k=0; k < k_bank->nkernels; k++) {
    int Delta_k =  k*k_bank->A->n*k_bank->nbands;
    for (int b=0; b < k_bank->nbands; b++) {
      int Delta_bk = b*k_bank->A->n + Delta_k;
      for (int i=0; i < k_bank->A->n; i++)
	M->val[i+Delta_bk] = (double) k_bank->weight[k][b].val[i];
    }
  }

  return(M);
}


iftKernel *iftDoGKernel(float stdev1, float stdev2)
{
  float      stdev,radius;
  iftKernel *K1, *K2, *K;
  iftAdjRel *A;
  int        i;

  if (stdev1 > stdev2){
    stdev  = stdev1; stdev1 = stdev2; stdev2 = stdev;
  }

  radius = 3*stdev2;

  K1 = iftGaussianKernel(radius,stdev1);
  K2 = iftGaussianKernel(radius,stdev2);

  A = iftSpheric(radius);
  K = iftCreateKernel(A);

  for (i=0; i < K->A->n; i++) {
    K->weight[i] = K1->weight[i] - K2->weight[i];
  }

  iftDestroyKernel(&K1);
  iftDestroyKernel(&K2);
  iftDestroyAdjRel(&A);

  return(K);
}

iftKernel *iftDoGKernel2D(float stdev1, float stdev2)
{
  float      stdev,radius;
  iftKernel *K1, *K2, *K;
  iftAdjRel *A;
  int        i;

  if (stdev1 > stdev2){
    stdev  = stdev1; stdev1 = stdev2; stdev2 = stdev;
  }

  radius = 3*stdev2;

  K1 = iftGaussianKernel2D(radius,stdev1);
  K2 = iftGaussianKernel2D(radius,stdev2);

  A = iftCircular(radius);
  K = iftCreateKernel(A);

  for (i=0; i < K->A->n; i++) {
    K->weight[i] = K1->weight[i] - K2->weight[i];
  }

  iftDestroyKernel(&K1);
  iftDestroyKernel(&K2);
  iftDestroyAdjRel(&A);

  return(K);
}


iftImage *iftSobelGradientMagnitude(iftImage *img)
{
  iftKernel *K[3];
  iftImage  *grad[3];
  int  i, p;

  if (iftIs3DImage(img)){
    K[0] = iftSobelXKernel();
    K[1] = iftSobelYKernel();
    K[2] = iftSobelZKernel();
    for (i=0; i < 3; i++) {
      grad[i] = iftLinearFilter(img,K[i]);
      iftDestroyKernel(&K[i]);
    }

    for (p=0; p < img->n; p++)
      grad[2]->val[p] = (int)sqrt((double)grad[0]->val[p]*(double)grad[0]->val[p] +
				  (double)grad[1]->val[p]*(double)grad[1]->val[p] +
				  (double)grad[2]->val[p]*(double)grad[2]->val[p]);

  }else{
    grad[2] = iftCreateImage(img->xsize,img->ysize,img->zsize);
    K[0]    = iftSobelXKernel2D();
    K[1]    = iftSobelYKernel2D();
    for (i=0; i < 2; i++){
      grad[i] = iftLinearFilter(img,K[i]);
      iftDestroyKernel(&K[i]);
    }

    for (p=0; p < img->n; p++)
      grad[2]->val[p] = (int)sqrt((double)grad[0]->val[p]*(double)grad[0]->val[p] +
				  (double)grad[1]->val[p]*(double)grad[1]->val[p]);

  }

  iftDestroyImage(&grad[0]);
  iftDestroyImage(&grad[1]);

  return(grad[2]);
}


iftImage *iftMeanFilter(const iftImage *img, const iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = (iftIs3DImage(img)) ? iftSpheric(1.0) : iftCircular(1.0);
    else A = iftCopyAdjacency(Ain);


    iftImage *filt = iftCreateImageFromImage(img);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++) {
        iftVoxel u = iftGetVoxelCoord(img, p);
        int n_valid_neighbors = 0;

        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(img, v)) {
                iftImgVoxelVal(filt, u) += iftImgVoxelVal(img, v);
                n_valid_neighbors++;
            }
        }
        iftImgVoxelVal(filt, u) = iftRound(iftImgVoxelVal(filt, u) / ((float) n_valid_neighbors));
    }

    iftDestroyAdjRel(&A);

    return filt;
}


iftImage *iftMedianFilter(const iftImage *img, iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = (iftIs3DImage(img)) ? iftSpheric(1.0) : iftCircular(1.0);
    else A = iftCopyAdjacency(Ain);

    iftImage *filt = iftCreateImageFromImage(img);

    #pragma omp parallel
    for (int p = 0; p < img->n; p++) {
        iftVoxel u = iftGetVoxelCoord(img,p);
        int n = 0;
        int val[A->n];

        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            
            if (iftValidVoxel(img, v)) {
                int q = iftGetVoxelIndex(img, v);
                val[i] = img->val[q];
                n++;
            }
        }
        filt->val[p] = iftSelectTheKthElem(val, n, n/2);
    }

    if (iftIsColorImage(img))
        iftCopyCbCr(img, filt);

    iftDestroyAdjRel(&A);

    return filt;
}


iftMImage *iftMMedianFilter(iftMImage *img, iftAdjRel *A)
{
	int       n,j,i,b,p,q;
	float 		*value;
	iftVoxel   u,v;
	iftMImage  *med;

	value = iftAllocFloatArray(A->n);
	med   = iftCreateMImage(img->xsize,img->ysize,img->zsize, img->m);

	for (p=0; p < img->n; p++)
	{
		u = iftMGetVoxelCoord(img,p);
		for(b = 0; b < img->m; b++)
		{
			n = 0;
			for (i=0; i < A->n; i++)
			{
				v = iftGetAdjacentVoxel(A,u,i);
				if (iftMValidVoxel(img,v))
				{
					q = iftMGetVoxelIndex(img,v);
					/* insertion sort */
					value[n] = img->val[q][b];
					j        = n;
					while ((j > 0)&&(value[j] < value[j-1]))
					{
						iftSwap(value[j], value[j-1]);
						j--;
					}
					n++;
				}
			}
			med->val[p][b] = value[n/2];
		}
	}

	iftFree(value);

	iftMCopyVoxelSize(img,med);

  return(med);
}

iftImage *iftModaFilter(iftImage *img, iftAdjRel *A)
{
  int      *hist=iftAllocIntArray(iftMaximumValue(img)+1),i,p,q,imax,fmax;
  iftVoxel  u,v;
  iftImage *moda;

  moda  = iftCreateImage(img->xsize,img->ysize,img->zsize);

  for (p=0; p < img->n; p++) {
    u = iftGetVoxelCoord(img,p);
    for (i=0; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){
	q = iftGetVoxelIndex(img,v);
	hist[img->val[q]]++;
      }
    }
    imax = IFT_NIL; fmax = IFT_INFINITY_INT_NEG;
    for (i=0; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){
	q = iftGetVoxelIndex(img,v);
	if (hist[img->val[q]] > fmax) {
	  imax = img->val[q];
	  fmax = hist[img->val[q]];
	}
	hist[img->val[q]]=0;
      }
    }
    moda->val[p] = imax;
  }

  if (img->Cb != NULL)
    iftCopyCbCr(img,moda);
  iftCopyVoxelSize(img,moda);

  iftFree(hist);

  return(moda);
}

iftImage  *iftSmoothImage(iftImage *img, iftAdjRel *A, float sigma)
{
  iftImage  *fimg=iftCopyImage(img);
  float      Sigma = 2.0*sigma*sigma, *ws=iftAllocFloatArray(A->n);

  float sum=0.0;
  for (int i=0; i < A->n; i++){ 
    ws[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];
    sum  += ws[i];
  }
  for (int i=0; i < A->n; i++){
    ws[i] = (sum - ws[i])/sum;
  }
  
  if (iftIsColorImage(img)){
#pragma omp parallel for shared(img,fimg,A,Sigma,ws)
    for (int p=0; p < img->n; p++) {
      iftVoxel u   = iftGetVoxelCoord(img,p);
      float Yp     = img->val[p], Cbp    = img->Cb[p], Crp    = img->Cr[p];
      float valY   = 0.0, valCb = 0.0, valCr = 0.0, sum = 0.0;
      float    *wr = iftAllocFloatArray(A->n);
      for (int i=0; i < A->n; i++){
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(img,v)){
	  int q     = iftGetVoxelIndex(img,v);
	  float Yq  = img->val[q], Cbq = img->Cb[q], Crq = img->Cr[q];
	  wr[i] = (Yq-Yp)*(Yq-Yp) +
	    (Cbq-Cbp)*(Cbq-Cbp) + (Crq-Crp)*(Crq-Crp);
	  wr[i] = ws[i]*expf(-wr[i]/Sigma);
	  sum         += wr[i];
	  valY        += ((float) img->val[q] * wr[i]);
	  valCb       += ((float) img->Cb[q]  * wr[i]);
	  valCr       += ((float) img->Cr[q]  * wr[i]);
	}
      }
      fimg->val[p] = iftRound(valY / sum);
      fimg->Cb[p]  = iftRound(valCb / sum);
      fimg->Cr[p]  = iftRound(valCr / sum);
      iftFree(wr);
    }
  }else{
#pragma omp parallel for shared(img,fimg,A,Sigma,ws)
    for (int p=0; p < img->n; p++) {
      iftVoxel u   = iftGetVoxelCoord(img,p);
      float Yp     = img->val[p];
      float valY   = 0.0, sum = 0.0;
      float    *wr = iftAllocFloatArray(A->n);
      for (int i=0; i < A->n; i++){
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(img,v)){
	  int q     = iftGetVoxelIndex(img,v);
	  float Yq  = img->val[q];
	  wr[i] = (Yq-Yp)*(Yq-Yp);
	  wr[i] = ws[i]*expf(-wr[i]/Sigma);
	  sum         += wr[i];
	  valY        += ((float) img->val[q] * wr[i]);
	}
      }
      fimg->val[p] = iftRound(valY / sum);
      iftFree(wr);
    }
  }

  iftFree(ws);
  iftCopyVoxelSize(img,fimg);

  return(fimg);
}

iftImage  *iftSmoothImageInRegion(iftImage *img, iftImage *mask, iftAdjRel *A, float sigma)
{
  iftImage  *fimg=iftCopyImage(img);
  float      Sigma = 2.0*sigma*sigma, *ws=iftAllocFloatArray(A->n);

  for (int i=0; i < A->n; i++)
    ws[i] = 1.0/(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i] + 1);

  if (iftIsColorImage(img)){
#pragma omp parallel for shared(img,mask,fimg,A,Sigma,ws)
    for (int p=0; p < img->n; p++) {
      if (mask->val[p]) {
	iftVoxel u   = iftGetVoxelCoord(img,p);
	float Yp     = img->val[p], Cbp    = img->Cb[p], Crp    = img->Cr[p];
	float valY   = 0.0, valCb = 0.0, valCr = 0.0, sum = 0.0;
	float    *wr = iftAllocFloatArray(A->n);
	for (int i=0; i < A->n; i++){
	  iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(img,v)){
	    int q     = iftGetVoxelIndex(img,v);
	    float Yq  = img->val[q], Cbq = img->Cb[q], Crq = img->Cr[q];
	    wr[i] = (Yq-Yp)*(Yq-Yp) +
	      (Cbq-Cbp)*(Cbq-Cbp) + (Crq-Crp)*(Crq-Crp);
	    wr[i] = ws[i]*expf(-wr[i]/Sigma);
	    sum         += wr[i];
	    valY        += ((float) img->val[q] * wr[i]);
	    valCb       += ((float) img->Cb[q]  * wr[i]);
	    valCr       += ((float) img->Cr[q]  * wr[i]);
	  }
	}
	fimg->val[p] = iftRound(valY / sum);
	fimg->Cb[p]  = iftRound(valCb / sum);
	fimg->Cr[p]  = iftRound(valCr / sum);
	iftFree(wr);
      }
    }
  }else{
#pragma omp parallel for shared(img,fimg,A,Sigma,ws)
    for (int p=0; p < img->n; p++) {
      if (mask->val[p]) {
	iftVoxel u   = iftGetVoxelCoord(img,p);
	float Yp     = img->val[p];
	float valY   = 0.0, sum = 0.0;
	float    *wr = iftAllocFloatArray(A->n);
	for (int i=0; i < A->n; i++){
	  iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(img,v)){
	    int q     = iftGetVoxelIndex(img,v);
	    float Yq  = img->val[q];
	    wr[i] = (Yq-Yp)*(Yq-Yp);
	    wr[i] = ws[i]*expf(-wr[i]/Sigma);
	    sum         += wr[i];
	    valY        += ((float) img->val[q] * wr[i]);
	  }
	}
	fimg->val[p] = iftRound(valY / sum);
	iftFree(wr);
      }
    }
  }
  iftFree(ws);
  iftCopyVoxelSize(img,fimg);

  return(fimg);
}


/*@brief This function performs the bilateral filtering of a grayscale
  or color image. s_sigma is the spatial radius used to define the
  spatial adjacency in pixels. r_sigma is the feature space radius and
  should be between [0.0 and 1.0] because the channels are scaled
  accordingly.
*/

iftImage *iftFastBilateralFilter2D(iftImage *img, int s_sigma, float r_sigma)
{
  iftVoxel v;
  float * __restrict__ input;
  float * __restrict__ output;
  int c,p, nchannels = 1;
  iftImage *feats = NULL;
  int normalization_value = 255;

  if(iftIs3DImage(img))
      iftError("Fast bilateral filtering is only defined for 2D images",
               "iftFastBilateralFilter2D");

  // Determining the radiometric resolution and data necessary to convert the image to [0.0, 1.0]
  if(iftIsColorImage(img)) {
    nchannels = 3;
  }

  if(nchannels > 1)
    feats = iftCreateColorImage(img->xsize, img->ysize, img->zsize, iftImageDepth(img));
  else
    feats = iftCreateImage(img->xsize, img->ysize, img->zsize);

  input  = iftAllocAlignedFloatArray(img->n * nchannels, 32);
  output = iftAllocAlignedFloatArray(img->n * nchannels, 32);

  v.z = 0; // 2D image
  if(iftIsColorImage(img)) {
    for(v.y = 0; v.y < img->ysize; v.y++)
      for(v.x = 0; v.x < img->xsize; v.x++) {
        iftColor rgb;

        p = v.x+img->tby[v.y];

        rgb = iftGetRGB(img, p, 255);

        for(c = 0; c < nchannels; c++) {
          // Converting each band to range [0.0, 1.0]
          input[c+nchannels*p] = rgb.val[c]/(float)normalization_value;
        }
      }
  } else {
    for(v.y = 0; v.y < img->ysize; v.y++)
      for(v.x = 0; v.x < img->xsize; v.x++) {
        p = v.x+img->tby[v.y];

        // Converting each band to range [0.0, 1.0]
        input[p] = img->val[p]/(float)normalization_value;
      }
  }

  iftFastBilateralFilter2DAux(input, output, img->xsize, img->ysize, nchannels, s_sigma, r_sigma);


  v.z = 0; // 2D image
  if(iftIsColorImage(img)) {
    for(v.y = 0; v.y < img->ysize; v.y++)
      for(v.x = 0; v.x < img->xsize; v.x++) {
        iftColor rgb;
        // Initializing color to avoid compilation warning
        rgb.val[0] = 0;
        rgb.val[1] = 0;
        rgb.val[2] = 0;

        p = v.x+img->tby[v.y];

        for(c = 0; c < nchannels; c++) {
          // Converting each band back to the original range
          rgb.val[c] = output[c+nchannels*p]*normalization_value;
        }
        iftSetRGB2(feats, p, rgb, 255);
      }
  } else {
    for(v.y = 0; v.y < img->ysize; v.y++)
      for(v.x = 0; v.x < img->xsize; v.x++) {
        p = v.x+img->tby[v.y];
        // Converting each band back to the original range
        feats->val[p] = output[p]*normalization_value;
      }
  }

  _mm_free(input);
  _mm_free(output);

  return feats;
}


/*
  @brief This function performs the bilateral filtering of an
  iftMImage. s_sigma is the spatial radius used to define the spatial
  adjacency in pixels. r_sigma is the feature space radius and should
  be between [0.0 and 1.0] because the channels are scaled
  accordingly.
*/
iftMImage *iftFastMBilateralFilter2D(iftMImage *img, int s_sigma, float r_sigma)
{
  iftVoxel v;
  float * __restrict__ input;
  float * __restrict__ output;
  int c,p, nchannels = img->m;
  float max_val[nchannels];
  float min_val[nchannels];
  float diff;
  iftMImage *feats = NULL;

  if(iftIs3DMImage(img))
      iftError("Fast bilateral filtering is only defined for 2D images",
               "iftFastMBilateralFilter2D");

  // Determining the range of each channel
  for(c = 0; c < nchannels; c++) {
    max_val[c] = IFT_INFINITY_FLT_NEG;
    min_val[c] = IFT_INFINITY_FLT;

    for(p = 0; p < img->n; p++) {
      max_val[c] = iftMax(max_val[c], img->val[p][c]);
      min_val[c] = iftMin(min_val[c], img->val[p][c]);
    }
  }


  feats = iftCreateMImage(img->xsize, img->ysize, img->zsize, nchannels);


  input = iftAllocAlignedFloatArray(img->n * nchannels, 32);
  output = iftAllocAlignedFloatArray(img->n * nchannels, 32);

  v.z = 0; // 2D image
  for(v.y = 0; v.y < img->ysize; v.y++)
    for(v.x = 0; v.x < img->xsize; v.x++) {
      p = v.x+img->tby[v.y];

      for(c = 0; c < nchannels; c++) {
        diff = max_val[c] - min_val[c];
        // Converting each band to range [0.0, 1.0]
        input[c+nchannels*p] = (img->val[p][c]-min_val[c])/diff;
      }
    }

  iftFastBilateralFilter2DAux(input, output, img->xsize, img->ysize, nchannels, s_sigma, r_sigma);

  v.z = 0; // 2D image
  for(v.y = 0; v.y < img->ysize; v.y++)
    for(v.x = 0; v.x < img->xsize; v.x++)
      for(c = 0; c < nchannels; c++) {
        p = v.x+img->tby[v.y];

        diff = max_val[c] - min_val[c];
        // Converting each band back to the original range
        feats->val[p][c] = output[c+nchannels*p]*diff+min_val[c];
      }

  _mm_free(input);
  _mm_free(output);

  return feats;
}


/* Private functions obtained from https://github.com/victormatheus/halide-casestudies/blob/master/bilateral-filter-fast.cpp*/

/* This function is an image processing operation for GEGL
 *
 * GEGL is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * GEGL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GEGL; if not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2012 Victor Oliveira <victormatheus@gmail.com>
 */

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#if __GNUC_MINOR__ > 6
#include <stdalign.h>
#endif


inline static float lerp(float a, float b, float v) {
  return (1.0f - v) * a + v * b;
}

// inline static int clamp(int x, int low, int high) {
//   return (x < low) ? low : ((x > high) ? high : x);
// }

void iftFastBilateralFilter2DAux(const float * __restrict__ _input,
                                 float       * __restrict__ _output,
                                 int   width,
                                 int   height,
                                 int 	channels,
                                 int   s_sigma,
                                 float r_sigma) {
  const int padding_xy = 2;
  const int padding_z  = 2;

  const int sw = (width -1) / s_sigma + 1 + (2 * padding_xy);
  const int sh = (height-1) / 	 + 1 + (2 * padding_xy);
  const int depth = (int)(1.0f / r_sigma) + 1 + (2 * padding_z);
  int err;

  /* down-sampling */

#if __GNUC_MINOR__ > 6
  // enabling AVX (SSE evolution)
  float * __restrict__ input  = (float * __restrict__) __builtin_assume_aligned(_input,  32);
  float * __restrict__ output = (float * __restrict__) __builtin_assume_aligned(_output, 32);
#else
  float * __restrict__ input  = (float * __restrict__) _input;
  float * __restrict__ output = (float * __restrict__) _output;
#endif

  float * __restrict__ grid ;
  float * __restrict__ blurx;
  float * __restrict__ blury;
  float * __restrict__ blurz;

  err = posix_memalign((void **) &grid , 32, sw * sh * depth * channels * 2 * sizeof(float));
  if(err != 0) iftError("Memory alignment failed", "iftFastBilateralFilter2D");

  err =posix_memalign((void **) &blurx, 32, sw * sh * depth * channels * 2 * sizeof(float));
  if(err != 0) iftError("Memory alignment failed", "iftFastBilateralFilter2D");

  err =posix_memalign((void **) &blury, 32, sw * sh * depth * channels * 2 * sizeof(float));
  if(err != 0) iftError("Memory alignment failed", "iftFastBilateralFilter2D");

  err =posix_memalign((void **) &blurz, 32, sw * sh * depth * channels * 2 * sizeof(float));
  if(err != 0) iftError("Memory alignment failed", "iftFastBilateralFilter2D");

#define INPUT(x,y,c) input[(c)+channels*((x) + width * (y))]

#define GRID(x,y,z,c,i) grid [(i)+2*((c)+channels*((x)+sw*((y)+(z)*sh)))]
#define BLURX(x,y,z,c,i) blurx[(i)+2*((c)+channels*((x)+sw*((y)+(z)*sh)))]
#define BLURY(x,y,z,c,i) blury[(i)+2*((c)+channels*((x)+sw*((y)+(z)*sh)))]
#define BLURZ(x,y,z,c,i) blurz[(i)+2*((c)+channels*((x)+sw*((y)+(z)*sh)))]

  for (int k=0; k < (sw * sh * depth * channels * 2); k++) {
    grid [k] = 0.0f;
    blurx[k] = 0.0f;
    blury[k] = 0.0f;
    blurz[k] = 0.0f;
  }

  /* downsampling */

  for(int y = 0; y < height; y++)
    for(int x = 0; x < width; x++)
      for(int c = 0; c < channels; c++) {
        const float z = INPUT(x,y,c);

        const int small_x = (int)((float)(x) / s_sigma + 0.5f) + padding_xy;
        const int small_y = (int)((float)(y) / s_sigma + 0.5f) + padding_xy;
        const int small_z = (int)((float)(z) / r_sigma + 0.5f) + padding_z;

        assert(small_x >= 0 && small_x < sw);
        assert(small_y >= 0 && small_y < sh);
        assert(small_z >= 0 && small_z < depth);

        GRID(small_x, small_y, small_z, c, 0) += INPUT(x, y, c);
        GRID(small_x, small_y, small_z, c, 1) += 1.0f;
      }

  /* blur in x, y and z */
  /* XXX: we could use less memory, but at expense of code readability */

#pragma omp parallel for
  for (int z = 1; z < depth-1; z++)
    for (int y = 1; y < sh-1; y++)
      for (int x = 1; x < sw-1; x++)
        for(int c = 0; c < channels; c++)
          for (int i=0; i<2; i++)
            BLURX(x, y, z, c, i) = (GRID (x-1, y, z, c, i) + 2.0f * GRID (x, y, z, c, i) + GRID (x+1, y, z, c, i)) / 4.0f;

#pragma omp parallel for
  for (int z = 1; z < depth-1; z++)
    for (int y = 1; y < sh-1; y++)
      for (int x = 1; x < sw-1; x++)
        for(int c = 0; c < channels; c++)
          for (int i=0; i<2; i++)
            BLURY(x, y, z, c, i) = (BLURX (x, y-1, z, c, i) + 2.0f * BLURX (x, y, z, c, i) + BLURX (x, y+1, z, c, i)) / 4.0f;

#pragma omp parallel for
  for (int z = 1; z < depth-1; z++)
    for (int y = 1; y < sh-1; y++)
      for (int x = 1; x < sw-1; x++)
        for(int c = 0; c < channels; c++)
          for (int i=0; i<2; i++)
            BLURZ(x, y, z, c, i) = (BLURY (x, y, z-1, c, i) + 2.0f * BLURY (x, y, z, c, i) + BLURY (x, y, z+1, c, i)) / 4.0f;

  /* trilinear filtering */

#pragma omp parallel for
  for (int y=0; y < height; y++)
    for (int x=0; x < width; x++)
      for(int c = 0; c < channels; c++) {
        float xf = (float)(x) / s_sigma + padding_xy;
        float yf = (float)(y) / s_sigma + padding_xy;
        float zf = INPUT(x,y,c) / r_sigma + padding_z;

        int x1 = (int)xf;
        int y1 = (int)yf;
        int z1 = (int)zf;

        int x2 = x1+1;
        int y2 = y1+1;
        int z2 = z1+1;

        float x_alpha = xf - x1;
        float y_alpha = yf - y1;
        float z_alpha = zf - z1;

        float interpolated[2];

        assert(xf >= 0 && xf < sw);
        assert(yf >= 0 && yf < sh);
        assert(zf >= 0 && zf < depth);

        for (int i=0; i<2; i++)
          interpolated[i] =
                  lerp(lerp(lerp(BLURZ(x1, y1, z1, c, i), BLURZ(x2, y1, z1, c, i), x_alpha),
                            lerp(BLURZ(x1, y2, z1, c, i), BLURZ(x2, y2, z1, c, i), x_alpha), y_alpha),
                       lerp(lerp(BLURZ(x1, y1, z2, c, i), BLURZ(x2, y1, z2, c, i), x_alpha),
                            lerp(BLURZ(x1, y2, z2, c, i), BLURZ(x2, y2, z2, c, i), x_alpha), y_alpha), z_alpha);
        output[channels*(y*width+x)+c] = interpolated[0] / interpolated[1];
      }

  free (grid);
  free (blurx);
  free (blury);
  free (blurz);
}



/**
@file @brief It normalizes a gray-scale image within [0,Imax] by
reducing noise and illumination variation, and by enhancing the fine
texture details and object borders. The normalization uses a
convolution between an unit kernel and the squared values of the
image. This convolution is fast computed by multiplication of the
corresponding matrices.
*/

iftImage *iftNormalizeImage(iftImage *img, iftAdjRel *adj, int Imax)
{
  iftKernel     *K     = iftUniformKernel(adj);
  double        *value = iftAllocDoubleArray(img->n), minval, maxval;
  iftMatrix     *A     = iftKernelToMatrix(K);
  iftMatrix     *B     = iftImageAdjToMatrix(img,K->A);
  for (int p=0; p < B->n; p++) 
    B->val[p] = B->val[p]*B->val[p];
  iftMatrix     *M     = iftMultMatrices(A,B);
  for (int p=0; p < M->n; p++) 
    M->val[p] = sqrt(M->val[p]);
  iftImage      *sum   = iftMatrixToImage(M,img->xsize,img->ysize,img->zsize);
  iftImage      *nimg  = iftCreateImage(img->xsize,img->ysize,img->zsize);

  iftDestroyMatrix(&M);  
  iftDestroyMatrix(&B);
  iftDestroyMatrix(&A);
  iftDestroyKernel(&K);

  minval = IFT_INFINITY_DBL; maxval = IFT_INFINITY_DBL_NEG;
  for (int p=0; p < sum->n; p++) {
    if (sum->val[p] > IFT_EPSILON){
      value[p] = (double)img->val[p]/(double)sum->val[p];
      if (value[p] < minval)
	minval=value[p];
      if (value[p] > maxval)
	maxval=value[p];
    }
  }

  iftDestroyImage(&sum);

  double delta = maxval - minval;
  if (delta > IFT_EPSILON){
    for (int p=0; p < img->n; p++) {
      nimg->val[p] = (int)(Imax*(value[p]-minval)/delta);
    }
  }
  iftFree(value);

  return(nimg);
}

/**
@file @brief It computes the alpha-pooling of an input image for a given stride value. The pooling requires convolution between an unit kernel and the image values to the power of alpha. The convolution is computed by fast matrix multiplication. 
*/

iftImage *iftAlphaPooling(iftImage *img, iftAdjRel *adj, int stride, float alpha)
{
  iftKernel     *K     = iftUniformKernel(adj);
  iftMatrix     *A     = iftKernelToMatrix(K), *B, *M;
  iftFastAdjRel *F     = iftCreateFastAdjRel(adj, img->tby, img->tbz);
  iftImage      *pool_img=NULL, *output_img=NULL;
  int            xsize, ysize, zsize, p, q;
  iftVoxel       u;
  float          sx, sy, sz;

  if (iftIsColorImage(img))
      iftError("It is only suitable for grayscale images", "iftAlphaPooling");
  
  if (iftIs3DImage(img)){

    /* Reduces the image size to disconsider borders and to take into
       account the stride */

    xsize = ceilf(((float)(img->xsize - 2*F->bx))/stride);
    ysize = ceilf(((float)(img->ysize - 2*F->by))/stride);
    zsize = ceilf(((float)(img->zsize - 2*F->bz))/stride);

    pool_img = iftCreateImage(xsize, ysize, zsize);
    
    q=0;
    for(u.z = F->bz; u.z < img->zsize-F->bz; u.z += stride){
      for(u.y = F->by; u.y < img->ysize-F->by; u.y += stride){
	for(u.x = F->bx; u.x < img->xsize-F->bx; u.x += stride){
	  p = iftGetVoxelIndex(img,u);
	  pool_img->val[q] = img->val[p]; q++;
	}
      }
    }

    /* Apply the convolution */

    B     = iftImageAdjToMatrix(pool_img,K->A);
    for (p=0; p < B->n; p++) 
      B->val[p] = pow(B->val[p],alpha);
    M     = iftMultMatrices(A,B);
    for (p=0; p < M->n; p++) 
      M->val[p] = pow(M->val[p],1.0/alpha);

    iftDestroyImage(&pool_img);

    pool_img = iftMatrixToImage(M,xsize,ysize,zsize);

    /* Return the resulting image to its original size */
  
    sx = img->xsize/xsize;
    sy = img->ysize/ysize;
    sz = img->zsize/zsize;
    
    output_img = iftInterp(pool_img,sx,sy,sz);
    iftCopyVoxelSize(img,output_img);

  } else {

    /* Reduces the image size to disconsider borders and to take into
       account the stride */

    xsize = ceilf(((float)(img->xsize - 2*F->bx))/stride);
    ysize = ceilf(((float)(img->ysize - 2*F->by))/stride);
    zsize = 1;

    pool_img = iftCreateImage(xsize, ysize, 1); u.z = 0;
    q=0;
    for(u.y = F->by; u.y < img->ysize-F->by; u.y += stride){
      for(u.x = F->bx; u.x < img->xsize-F->bx; u.x += stride){
	p = iftGetVoxelIndex(img,u);
	pool_img->val[q] = img->val[p]; q++;
      }
    }

    /* Apply the convolution */

    B     = iftImageAdjToMatrix(pool_img,K->A);
    for (p=0; p < B->n; p++) {
      B->val[p] = pow(B->val[p],alpha);
    }
    M     = iftMultMatrices(A,B);
    for (p=0; p < M->n; p++) {
      M->val[p] = pow(M->val[p],1.0/alpha);
    }
    iftDestroyImage(&pool_img);

    pool_img = iftMatrixToImage(M,xsize,ysize,zsize);
    
    /* Return the resulting image to its original size */

    sx = img->xsize/xsize;
    sy = img->ysize/ysize;

    output_img = iftInterp2D(pool_img,sx,sy);
  }

  iftDestroyImage(&pool_img);
  iftDestroyMatrix(&M);  
  iftDestroyMatrix(&B);
  iftDestroyMatrix(&A);
  iftDestroyKernel(&K);
  iftDestroyFastAdjRel(&F);

  return(output_img);
}

iftMatrix *iftRandomKernelBankAsMatrix(int size, int nbands, int nKernels)
{
    iftMatrix *kernels = iftCreateMatrix(nKernels, size*nbands); /* ncols=nKernels, nrows=size*nbands */

    iftRandomSeed(time(NULL));
    
    for (int col = 0; col < nKernels; col++)
    {
      /* for each kernel */
        for (int row = 0; row < kernels->nrows; row++)
        { /* create its coefficients randomly */
            iftMatrixElem(kernels, col, row) = ((float)rand()) / ((float)RAND_MAX + 0.5);
        }
    }

    for (int col = 0; col < nKernels; col++)
    { /* for each kernel */
        float mean = 0;
        for (int row = 0; row < kernels->nrows; row++)
        { /* compute the mean */
            mean += iftMatrixElem(kernels, col, row);
        }
        mean /= kernels->nrows;
        for (int row = 0; row < kernels->nrows; row++)
        { /* centralize the coefficients */
            iftMatrixElem(kernels, col, row) = iftMatrixElem(kernels, col, row) - mean;
        }
    }

    for (int col = 0; col < nKernels; col++)
    { /* for each kernel */
        float norm = 0;
        for (int row = 0; row < kernels->nrows; row++)
        { /* compute the norm of the coefficients */
            norm += iftMatrixElem(kernels, col, row) * iftMatrixElem(kernels, col, row);
        }
        norm = sqrt(norm);
        if (norm > IFT_EPSILON)
        {
            for (int row = 0; row < kernels->nrows; row++)
            {   /* normalize its coefficients */
                iftMatrixElem(kernels, col, row) = iftMatrixElem(kernels, col, row) / norm;
            }
        }
    }

    return kernels;
}

iftFileSet *iftLocalBatchCentralization(iftFileSet *fs, int kernelSize, char *outputDir, char colSpace)
{
    int nClasses = iftFileSetLabelsNumber(fs);
    int nImgs = fs->n;
    int *nSampPerClass = iftCountSamplesPerClassFileSet(fs);
    iftImage *img0 = iftReadImageByExt(fs->files[0]->path);
    int nBands = iftIsColorImage(img0) ? 3 : 1;
    int nPositions = img0->n;
    iftAdjRel *A = iftRectangular(kernelSize, kernelSize);

    /* compute the mean feature for each band and position */
    iftMImage **mean = (iftMImage **)calloc(nClasses, sizeof(iftMImage *));
    for (int c = 0; c < nClasses; c++)
        mean[c] = iftCreateMImage(img0->xsize, img0->ysize, img0->zsize, nBands);

    for (int i = 0; i < nImgs; i++)
    {
        int c = fs->files[i]->label - 1;
        iftImage *img = iftReadImageByExt(fs->files[i]->path);

        iftMImage *mimg = iftImageToMImage(img, colSpace);
        #pragma omp parallel for
        for (int p = 0; p < nPositions; p++) {
            iftVoxel u = iftMGetVoxelCoord(mimg, p);
            for (int a = 0; a < A->n; a++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, a);
                if (iftMValidVoxel(mimg, v)) {
                    int q = iftGetVoxelIndex(mimg, v);
                    for (int b = 0; b < nBands; b++) {
                        mean[c]->val[p][b] += mimg->val[q][b];
                    }
                }
            }
        }
        iftDestroyImage(&img);
        iftDestroyMImage(&mimg);
    }

    #pragma omp parallel for
    for (int p = 0; p < nPositions; p++) {
        for (int b = 0; b < nBands; b++) {
            for (int c = 0; c < nClasses; c++) {
                mean[c]->val[p][b] /= (A->n * nSampPerClass[c + 1]);
            }
        }
    }

    /* compute the central point per position (origin of the local
       feature space) as the mean of the mean features from all
       classes */
    iftMImage *origin = iftCreateMImage(img0->xsize, img0->ysize, img0->zsize, nBands);
    char filename[2048];

    #pragma omp parallel for
    for (int p = 0; p < nPositions; p++) {
        for (int b = 0; b < nBands; b++) {
            for (int c = 0; c < nClasses; c++) {
                origin->val[p][b] += mean[c]->val[p][b];
            }
            origin->val[p][b] /= nClasses;
        }
    }

    iftFileSet *fsOut = iftCreateFileSet(fs->n);

    /* centralize feature mimages and output the resulting mimages */
    for (int i = 0; i < nImgs; i++) {
        char *basename = iftFilename(fs->files[i]->path, iftFileExt(fs->files[i]->path));
        iftImage *img = iftReadImageByExt(fs->files[i]->path);
        iftMImage *mimg = iftImageToMImage(img, colSpace);
        #pragma omp parallel for
        for (int p = 0; p < nPositions; p++) {
            for (int b = 0; b < nBands; b++) {
                mimg->val[p][b] = mimg->val[p][b] - origin->val[p][b];
            }
        }

        sprintf(filename, "%s/%s.mimg", outputDir, basename);
        iftWriteMImage(mimg, filename);
        fsOut->files[i] = iftCreateFile(iftAbsPathname(filename));
        iftDestroyImage(&img);
        iftDestroyMImage(&mimg);
        iftFree(basename);
    }

    return fsOut;
}

iftMImage *iftMConvolution(iftMImage *mimg, iftMatrix *kernels, int stride)
{
    int nKernels = kernels->ncols;
    int kernelSize = sqrt(kernels->nrows/mimg->m);
    iftAdjRel *A = iftRectangular(kernelSize, kernelSize);
    iftMatrix *vectors = iftMImageToFeatureMatrix(mimg, A, NULL);
    iftMatrix *conv = iftMultMatrices(vectors, kernels);
    iftMImage *convMImg = iftMatrixToMImage(conv, mimg->xsize, mimg->ysize, mimg->zsize, nKernels, 'c');

    iftMImage *mimgRet = iftCreateMImage(ceilf(mimg->xsize/(float)stride), ceilf(mimg->ysize/(float)stride), iftMax(ceilf(mimg->zsize/(float)stride), 1), convMImg->m);
    int q = 0;
    iftVoxel u;
    for (u.z = 0; u.z < convMImg->zsize; u.z = u.z + stride) {
        for (u.y = 0; u.y < convMImg->ysize; u.y = u.y + stride) {
            for (u.x = 0; u.x < convMImg->xsize; u.x = u.x + stride) {
                int p = iftMGetVoxelIndex(convMImg, u);
                for (int b = 0; b < convMImg->m; b++)
                    mimgRet->val[q][b] = convMImg->val[p][b];
                q++;
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftDestroyMatrix(&vectors);
    iftDestroyMatrix(&conv);
    iftDestroyMImage(&convMImg);

    return mimgRet;
}

iftMImageArray *iftMConvolutionArray(iftMImageArray *mimgArray, iftMatrix *kernels, int stride)
{
    iftMImage *mimg0 = mimgArray->val[0];

    int nKernels = kernels->ncols;
    int kernelSize = sqrt(kernels->nrows/mimg0->m);
    iftAdjRel *A = iftRectangular(kernelSize, kernelSize);
    iftMatrix *vectors = iftMImageArrayToFeatureMatrix(mimgArray, A);
    iftMatrix *conv = iftMultMatrices(vectors, kernels);
    iftMImageArray *convMImgArray = iftMatrixToMImageArray(conv, mimgArray->n, mimg0->xsize, mimg0->ysize, mimg0->zsize, nKernels, 'c');

    iftMImageArray *mimgArrayOut = iftCreateMImageArray(mimgArray->n);

    for (int i = 0; i < mimgArray->n; i++) {
        mimgArrayOut->val[i] = iftCreateMImage(ceilf(mimg0->xsize/(float)stride), ceilf(mimg0->ysize/(float)stride), iftMax(ceilf(mimg0->zsize/(float)stride), 1), convMImgArray->val[i]->m);
        int q = 0;
        iftVoxel u;
        for (u.z = 0; u.z < convMImgArray->val[i]->zsize; u.z = u.z + stride) {
            for (u.y = 0; u.y < convMImgArray->val[i]->ysize; u.y = u.y + stride) {
                for (u.x = 0; u.x < convMImgArray->val[i]->xsize; u.x = u.x + stride) {
                    int p = iftMGetVoxelIndex(convMImgArray->val[i], u);
                    for (int b = 0; b < convMImgArray->val[i]->m; b++)
                        mimgArrayOut->val[i]->val[q][b] = convMImgArray->val[i]->val[p][b];
                    q++;
                }
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftDestroyMatrix(&vectors);
    iftDestroyMatrix(&conv);
    iftDestroyMImageArray(&convMImgArray);

    return mimgArrayOut;
}

iftMImage *iftMMaxPooling(iftMImage *mimg, float radius, int stride)
{
    iftAdjRel *A;
    if (iftIs3DMImage(mimg))
      A = iftSpheric(radius);
    else
      A = iftCircular(radius);

    iftMImage *pool = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m);
    iftSetMImage(pool, IFT_INFINITY_FLT_NEG);

#pragma omp parallel for
    for (int p = 0; p < mimg->n; p++) {
        iftVoxel u = iftMGetVoxelCoord(mimg, p);
        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v)) {
                int q = iftMGetVoxelIndex(mimg, v);
                for (int b = 0; b < mimg->m; b++)
                    if (mimg->val[q][b] > pool->val[p][b])
                        pool->val[p][b] = mimg->val[q][b];
            }            
        }
    }
    
    iftDestroyAdjRel(&A);
    
    if (stride > 1){
      iftMImage *mimgRet = iftCreateMImage(ceilf(mimg->xsize/(float)stride), ceilf(mimg->ysize/(float)stride), iftMax(ceilf(mimg->zsize/(float)stride), 1), mimg->m);
      int q = 0;
      iftVoxel u;
      for (u.z = 0; u.z < pool->zsize; u.z = u.z + stride) {
        for (u.y = 0; u.y < pool->ysize; u.y = u.y + stride) {
	  for (u.x = 0; u.x < pool->xsize; u.x = u.x + stride) {
	    int p = iftMGetVoxelIndex(pool, u);
	    for (int b = 0; b < pool->m; b++)
	      mimgRet->val[q][b] = pool->val[p][b];
	    q++;
	  }
        }
      }

      iftDestroyMImage(&pool);

      return mimgRet;
    }

    return pool;
}

iftMImageArray *iftMMaxPoolingArray(iftMImageArray *mimgArray, float radius, int stride)
{
    iftMImageArray *mimgArrayOut = iftCreateMImageArray(mimgArray->n);

    for(int i = 0; i < mimgArray->n; i++)
        mimgArrayOut->val[i] = iftMMaxPooling(mimgArray->val[i], radius, stride);
    
    return mimgArrayOut;
}

iftMImage *iftMMaxPoolingRoi(iftMImage *mimg, iftRoiArray *roiArray, int stride)
{
    iftMatrix *maxMat = iftCreateMatrix(mimg->m, roiArray->n);
    iftSetMatrix(maxMat, IFT_INFINITY_FLT_NEG);

    #pragma omp parallel for
    for (int i = 0; i < roiArray->n; i++) {
        iftRoi *roi = roiArray->val[i];
        for (int v = 0; v < roi->n; v++) {
            int q = iftMGetVoxelIndex(mimg, roi->val[v]);
            for (int b = 0; b < mimg->m; b++) {
                if (mimg->val[q][b] > iftMatrixElem(maxMat, b, i))
                    iftMatrixElem(maxMat, b, i) = mimg->val[q][b];
            }
        }
    }

    iftMImage *pool = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m);
    for (int i = 0; i < roiArray->n; i++) {
        iftRoi *roi = roiArray->val[i];
        for (int v = 0; v < roi->n; v++) {
            int q = iftMGetVoxelIndex(mimg, roi->val[v]);
            for (int b = 0; b < mimg->m; b++) {
                pool->val[q][b] = iftMatrixElem(maxMat, b, i);
            }
        }
    }

    iftMImage *mimgRet = iftCreateMImage(ceilf(mimg->xsize/(float)stride), ceilf(mimg->ysize/(float)stride), iftMax(ceilf(mimg->zsize/(float)stride), 1), mimg->m);
    int q = 0;
    iftVoxel u;
    for (u.z = 0; u.z < pool->zsize; u.z = u.z + stride) {
        for (u.y = 0; u.y < pool->ysize; u.y = u.y + stride) {
            for (u.x = 0; u.x < pool->xsize; u.x = u.x + stride) {
                int p = iftMGetVoxelIndex(pool, u);
                for (int b = 0; b < pool->m; b++)
                    mimgRet->val[q][b] = pool->val[p][b];
                q++;
            }
        }
    }

    iftDestroyMImage(&pool);
    iftDestroyMatrix(&maxMat);

    return mimgRet;
}

iftMImageArray *iftMMaxPoolingRoiArray(iftMImageArray *mimgArray, iftRoiArray *roiArray, int stride)
{
    iftMImageArray *mimgArrayOut = iftCreateMImageArray(mimgArray->n);

    for(int i = 0; i < mimgArray->n; i++)
        mimgArrayOut->val[i] = iftMMaxPoolingRoi(mimgArray->val[i], roiArray, stride);
    
    return mimgArrayOut;
}

iftMImage *iftMMinPooling(iftMImage *mimg, float radius, int stride)
{
    iftAdjRel *A;
    if (iftIs3DMImage(mimg))
      A = iftSpheric(radius);
    else
      A = iftCircular(radius);

    iftMImage *pool = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m);
    iftSetMImage(pool, IFT_INFINITY_FLT);

    #pragma omp parallel for
    for (int p = 0; p < mimg->n; p++) {
        iftVoxel u = iftMGetVoxelCoord(mimg, p);
        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftMValidVoxel(mimg, v)) {
                int q = iftMGetVoxelIndex(mimg, v);
                for (int b = 0; b < mimg->m; b++)
                    if (mimg->val[q][b] < pool->val[p][b])
                        pool->val[p][b] = mimg->val[q][b];
            }
        }
    }

    iftDestroyAdjRel(&A);
    
    if (stride > 1) {
      iftMImage *mimgRet = iftCreateMImage(ceilf(mimg->xsize/(float)stride), ceilf(mimg->ysize/(float)stride), iftMax(ceilf(mimg->zsize/(float)stride), 1), mimg->m);
      int q = 0;
      iftVoxel u;
      for (u.z = 0; u.z < pool->zsize; u.z = u.z + stride) {
        for (u.y = 0; u.y < pool->ysize; u.y = u.y + stride) {
	  for (u.x = 0; u.x < pool->xsize; u.x = u.x + stride) {
	    int p = iftMGetVoxelIndex(pool, u);
	    for (int b = 0; b < pool->m; b++)
	      mimgRet->val[q][b] = pool->val[p][b];
	    q++;
	  }
        }
      }
      iftDestroyMImage(&pool);
      return mimgRet;
    }
    return(pool);
}

iftMImageArray *iftMMinPoolingArray(iftMImageArray *mimgArray, float radius, int stride)
{
    iftMImageArray *mimgArrayOut = iftCreateMImageArray(mimgArray->n);

    for(int i = 0; i < mimgArray->n; i++)
        mimgArrayOut->val[i] = iftMMinPooling(mimgArray->val[i], radius, stride);
    
    return mimgArrayOut;
}

iftMImage *iftMReLU(iftMImage *mimg) {
    iftMImage *relu = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m);

    #pragma omp parallel for
    for (int p = 0; p < mimg->n; p++)
        for (int b = 0; b < mimg->m; b++)
            relu->val[p][b] = iftMax(0, mimg->val[p][b]);

    return relu;
}

iftMImageArray *iftMReLUArray(iftMImageArray *mimgArray)
{
    iftMImageArray *mimgArrayOut = iftCreateMImageArray(mimgArray->n);

    for(int i = 0; i < mimgArray->n; i++)
        mimgArrayOut->val[i] = iftMReLU(mimgArray->val[i]);
    
    return mimgArrayOut;
}
