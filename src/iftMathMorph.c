#include "iftMathMorph.h"

#include "ift/core/dtypes/GQueue.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/io/Stream.h"


/*------------- Private functions for FastAreaOpen and FastAreaClose -------*/

static __inline__ int iftFindRoot(int *parent, int x) {
    if (parent[x] >= 0) {
        parent[x] = iftFindRoot(parent, parent[x]);
        return parent[x];
    }
    return x;
}

/* Sort voxel vector by radixSort */

static __inline__ void iftRadixSort (int byte, int N, iftVoxel *source, iftVoxel *dest) {
  int count[256];
  int index[256];
  int i;

  for (i=0; i < 256; i++) 
    count[i]=0;

  for (i = 0; i < N; i++) {
    count[(source[i].x >> (byte * 8)) & 0xff]++;
  }

  index[0] = 0;
  for (i = 1; i < 256; i++ ) {
    index[i] = index[i - 1] + count[i - 1];
  }
  for (i = 0; i < N; i++ ) {
    dest[index[(source[i].x >> (byte * 8)) & 0xff]++] = source[i];
  }
}

/*-------------------- Public functions ----------------------- */


iftImage *iftDilateWithKernel(const iftImage *img, const iftKernel *K, const iftImage *mask)
{
  iftImage *dil=iftCreateImage(img->xsize,img->ysize,img->zsize);

  if (mask==NULL) {
#pragma omp parallel for shared(img, dil, K)
    for (int p = 0; p < img->n; p++) {
      iftVoxel u = iftGetVoxelCoord(img, p);
      dil->val[p] = img->val[p] + (int) K->weight[0];
      for (int i = 1; i < K->A->n; i++) {
        iftVoxel v = iftGetAdjacentVoxel(K->A, u, i);
        if (iftValidVoxel(img, v)) {
          int q = iftGetVoxelIndex(img, v);
          if (img->val[q] > dil->val[p])
            dil->val[p] = img->val[q] + (int) K->weight[i];
        }
      }
    }
  } else {
    iftVerifyImageDomains(img,mask,"iftDilateWithKernel");
#pragma omp parallel for shared(img, dil, K, mask)
    for (int p = 0; p < img->n; p++) {
      if (mask->val[p]!=0) {
        iftVoxel u = iftGetVoxelCoord(img, p);
        dil->val[p] = img->val[p] + (int) K->weight[0];
        for (int i = 1; i < K->A->n; i++) {
          iftVoxel v = iftGetAdjacentVoxel(K->A, u, i);
          if (iftValidVoxel(img, v)) {
            int q = iftGetVoxelIndex(img, v);
            if ((mask->val[q]!=0)&&(img->val[q] > dil->val[p]))
              dil->val[p] = img->val[q] + (int) K->weight[i];
          }
        }
      }
    }
  }

  if (iftIsColorImage(img))
    iftCopyCbCr(img,dil);
  iftCopyVoxelSize(img,dil);
  
  return(dil);
}


iftImage *iftErodeWithKernel(const iftImage *img, const iftKernel *K, const iftImage *mask)
{
  iftImage *ero=iftCreateImage(img->xsize,img->ysize,img->zsize);

  if (mask == NULL) {
#pragma omp parallel for shared(img, ero, K)
    for (int p = 0; p < img->n; p++) {
      iftVoxel u = iftGetVoxelCoord(img, p);
      ero->val[p] = img->val[p] - (int) K->weight[0];
      for (int i = 1; i < K->A->n; i++) {
	iftVoxel v = iftGetAdjacentVoxel(K->A, u, i);
	if (iftValidVoxel(img, v)) {
	  int q = iftGetVoxelIndex(img, v);
	  if (img->val[q] < ero->val[p])
	    ero->val[p] = img->val[q] - K->weight[i];
	}
      }
	  }
	} else {
      iftVerifyImageDomains(img,mask,"iftErodeWithKernel");
#pragma omp parallel for shared(img, ero, K, mask)
      for (int p = 0; p < img->n; p++) {
	if (mask->val[p] != 0) {
	  iftVoxel u = iftGetVoxelCoord(img, p);
	  ero->val[p] = img->val[p] - (int) K->weight[0];
	  for (int i = 1; i < K->A->n; i++) {
	    iftVoxel v = iftGetAdjacentVoxel(K->A, u, i);
	    if (iftValidVoxel(img, v)) {
	      int q = iftGetVoxelIndex(img, v);
	      if ((mask->val[q] != 0)&&(img->val[q] < ero->val[p]))
		ero->val[p] = img->val[q] - K->weight[i];
	    }
	  }
	}
      }
    }

  if (iftIsColorImage(img))
    iftCopyCbCr(img,ero);
  iftCopyVoxelSize(img,ero);

  return(ero);
}

iftImage *iftCloseWithKernel(const iftImage *img, const iftKernel *K, const iftImage *mask)
{
  iftImage *dil,*ero;

  dil = iftDilateWithKernel(img,K,mask);
  ero = iftErodeWithKernel(dil,K,mask);
  iftDestroyImage(&dil);
  return(ero);

}


iftImage *iftOpenWithKernel(const iftImage *img, const iftKernel *K, const iftImage *mask)
{
  iftImage *dil,*ero;

  ero = iftErodeWithKernel(img,K, mask);
  dil = iftDilateWithKernel(ero,K,mask);
  iftDestroyImage(&ero);
  return(dil);
}

iftImage *iftDilate(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage *dil=iftCreateImage(img->xsize,img->ysize,img->zsize);

  if (mask == NULL){
#pragma omp parallel for shared(img,dil,A)
    for (int p=0; p < img->n; p++) {
      iftVoxel u = iftGetVoxelCoord(img,p);
      dil->val[p] = img->val[p];
      for (int i=1; i < A->n; i++) {
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(img,v)){
	  int q = iftGetVoxelIndex(img,v);
	  if (img->val[q]>dil->val[p])
	    dil->val[p] = img->val[q];
	}
      }
    }
  } else {
    iftVerifyImageDomains(img,mask,"iftDilate");
#pragma omp parallel for shared(img,dil,A,mask)
    for (int p=0; p < img->n; p++) {
      if (mask->val[p]!=0){
	iftVoxel u = iftGetVoxelCoord(img,p);
	dil->val[p] = img->val[p];
	for (int i=1; i < A->n; i++) {
	  iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(img,v)){
	    int q = iftGetVoxelIndex(img,v);
	    if ((mask->val[q]!=0)&&(img->val[q]>dil->val[p]))
	      dil->val[p] = img->val[q];
	  }
	}
      }
    }
  }
 
  if (iftIsColorImage(img))
    iftCopyCbCr(img,dil);
  iftCopyVoxelSize(img,dil);
  
  return(dil);
}

iftImage *iftErode(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage *ero=iftCreateImage(img->xsize,img->ysize,img->zsize);

  if (mask == NULL) {
  #pragma omp parallel for shared(img,ero,A)
    for (int p=0; p < img->n; p++) {
      iftVoxel u = iftGetVoxelCoord(img,p);
      ero->val[p] = img->val[p];
      for (int i=1; i < A->n; i++) {
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(img,v)){
	  int q = iftGetVoxelIndex(img,v);
	  if (img->val[q]<ero->val[p])
	    ero->val[p] = img->val[q];
	}
      }
    }
  } else {
    iftVerifyImageDomains(img,mask,"iftErode");
#pragma omp parallel for shared(img,ero,A,mask)
    for (int p=0; p < img->n; p++) {
      if (mask->val[p]!=0){
	iftVoxel u = iftGetVoxelCoord(img,p);
	ero->val[p] = img->val[p];
	for (int i=1; i < A->n; i++) {
	  iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(img,v)){
	    int q = iftGetVoxelIndex(img,v);
	    if ((mask->val[q]!=0)&&(img->val[q]<ero->val[p]))
	      ero->val[p] = img->val[q];
	  }
	}
      }
    }
  }
  if (iftIsColorImage(img))
    iftCopyCbCr(img,ero);
  iftCopyVoxelSize(img,ero);
  
  return(ero);
}

iftImage *iftClose(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage *dil,*ero;

  dil = iftDilate(img,A,mask);
  ero = iftErode(dil,A,mask);
  iftDestroyImage(&dil);
  return(ero);

}

iftImage *iftOpen(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage *dil,*ero;

  ero = iftErode(img,A,mask);
  dil = iftDilate(ero,A,mask);
  iftDestroyImage(&ero);
  return(dil);
  
}

iftImage *iftAsfOC(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage *open=NULL,*close=NULL;

  open  = iftOpen(img,A,mask);
  close = iftClose(open,A,mask);
  iftDestroyImage(&open);

  return(close);
}

iftImage *iftAsfCO(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage *open=NULL,*close=NULL;

  close = iftClose(img,A,mask);
  open  = iftOpen(close,A,mask);
  iftDestroyImage(&close);

  return(open);
}


iftImage *iftDilateBin(const iftImage *bin, iftSet **seed, float radius)
{
  iftImage  *dist=NULL,*root=NULL,*dil=NULL;
  iftGQueue *Q=NULL;
  int        i,p,q,tmp;
  iftVoxel   u,v,r;
  iftAdjRel *A=NULL,*B=NULL;
  float      maxdist=radius*radius; // the EDT computes squared
				    // Euclidean distance values

  if (iftIs3DImage(bin)){
    A = iftSpheric(1.74);
    B = iftSpheric(1.0);
  } else {
    A = iftCircular(1.5);
    B = iftCircular(1.0);
  }

  if (*seed == NULL) {
    (*seed) = iftObjectBorderSet(bin,B);
  }
  
  iftDestroyAdjRel(&B);

  // Initialization

  dil   = iftCopyImage(bin);
  dist  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  root  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  Q     = iftCreateGQueue(IFT_QSIZE,bin->n,dist->val);

  for (p=0; p < bin->n; p++) {
    if (bin->val[p] == 0)
      dist->val[p]= IFT_INFINITY_INT;
    }

  while ((*seed) != NULL) {
    p = iftRemoveSet(seed);
    dist->val[p]=0;
    root->val[p]=p;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform: Distance transform
  
  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    if (dist->val[p] <= maxdist){

      dil->val[p] = bin->val[root->val[p]]; // dilation
      u = iftGetVoxelCoord(bin,p);
      r = iftGetVoxelCoord(bin,root->val[p]);

      for (i=1; i < A->n; i++){
	v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(bin,v)){
	  q = iftGetVoxelIndex(bin,v);
	  if (dist->val[q] > dist->val[p]){
	    tmp = (v.x-r.x)*(v.x-r.x) + (v.y-r.y)*(v.y-r.y) + (v.z-r.z)*(v.z-r.z);
	    if (tmp < dist->val[q]){
	      if (dist->val[q] != IFT_INFINITY_INT)
		iftRemoveGQueueElem(Q, q);
	      dist->val[q]  = tmp;
	      root->val[q]  = root->val[p];
	      iftInsertGQueue(&Q, q);
	    }
	  }
	}
      }
    }else{ /* seeds for a possible subsequent erosion */
      iftInsertSet(seed,p);      
    }
  }
  
  iftDestroyGQueue(&Q);
  iftDestroyImage(&root);
  iftDestroyImage(&dist);
  iftDestroyAdjRel(&A);

  
  return(dil);
}


iftImage *iftErodeBin(const iftImage *bin, iftSet **seed, float radius)
{
  iftImage  *dist=NULL,*root=NULL,*ero=NULL;
  iftGQueue *Q=NULL;
  int        i,p,q,tmp;
  iftVoxel   u,v,r;
  iftAdjRel *A=NULL,*B=NULL;
  float      maxdist=radius*radius; // the EDT computes squared
				    // Euclidean distance values

  if (iftIs3DImage(bin)){
    A = iftSpheric(1.74);
    B = iftSpheric(1.0);
  } else {
    A = iftCircular(1.5);
    B = iftCircular(1.0);
  }
 
  if ((*seed) == NULL) {
    (*seed) = iftBackgroundBorderSet(bin,B);
  }
  
  iftDestroyAdjRel(&B);

  // Initialization

  ero   = iftCopyImage(bin);
  dist  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  root  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  Q     = iftCreateGQueue(IFT_QSIZE,bin->n,dist->val);

  for (p=0; p < bin->n; p++) {
    if (bin->val[p] != 0)
      dist->val[p]= IFT_INFINITY_INT;
    }

  while ((*seed) != NULL) {
    p = iftRemoveSet(seed);
    dist->val[p]=0;
    root->val[p]=p;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform: Distance transform
  
  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    if (dist->val[p] <= maxdist){

      ero->val[p] = bin->val[root->val[p]]; // erosion
      u = iftGetVoxelCoord(bin,p);
      r = iftGetVoxelCoord(bin,root->val[p]);

      for (i=1; i < A->n; i++){
	v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(bin,v)){
	  q = iftGetVoxelIndex(bin,v);
	  if (dist->val[q] > dist->val[p]){
	    tmp = (v.x-r.x)*(v.x-r.x) + (v.y-r.y)*(v.y-r.y) + (v.z-r.z)*(v.z-r.z);
	    if (tmp < dist->val[q]){
	      if (dist->val[q] != IFT_INFINITY_INT)
		iftRemoveGQueueElem(Q, q);
	      dist->val[q]  = tmp;
	      root->val[q]  = root->val[p];
	      iftInsertGQueue(&Q, q);
	    }
	  }
	}
      }
    }else{ /* seeds for a possible subsequent dilation */
      iftInsertSet(seed,p);
    }
  }
  
  iftDestroyGQueue(&Q);
  iftDestroyImage(&root);
  iftDestroyImage(&dist);
  iftDestroyAdjRel(&A);
  
  return(ero);
}


iftImage *iftErodeLabelImage(const iftImage *label_img, float radius) {
    iftImage *erode_img = iftCreateImageFromImage(label_img);

    iftIntArray *labels = iftGetObjectLabels(label_img);

    // since there will not be intersection between the objects after erosion, 
    // we can ignore the race conditions during copying into final erode image
    #pragma omp parallel for
    for (int o = 0; o < labels->n; o++) {
        iftImage *obj_img = iftExtractObject(label_img, labels->val[o]);

        iftSet *S = NULL;
        iftImage *erode_obj_img = iftErodeBin(obj_img, &S, radius);

        for (int p = 0; p < erode_obj_img->n; p++)
            if (erode_obj_img->val[p] != 0)
                erode_img->val[p] = erode_obj_img->val[p];

        iftDestroyImage(&obj_img);
        iftDestroyImage(&erode_obj_img);
        iftDestroySet(&S);
    }

    iftDestroyIntArray(&labels);

    return erode_img;
}


iftImage *iftCloseBin(const iftImage *bin, float radius)
{
  iftImage *close=NULL,*dil=NULL;
  iftSet   *seed=NULL;

  iftImage *aux = iftAddFrame(bin,iftRound(radius)+5,0);
  dil   = iftDilateBin(aux,&seed,radius);
  iftDestroyImage(&aux);
  aux   = iftErodeBin(dil,&seed,radius);
  close = iftRemFrame(aux,iftRound(radius)+5);
  iftDestroyImage(&aux);
  iftDestroyImage(&dil);
  iftDestroySet(&seed);

  return(close);
}

iftImage *iftOpenBin(const iftImage *bin, float radius)
{
  iftImage *open=NULL,*ero=NULL;
  iftSet   *seed=NULL;

  iftImage *aux = iftAddFrame(bin,iftRound(radius)+5,0);
  ero   = iftErodeBin(aux,&seed,radius);
  iftDestroyImage(&aux);
  aux   = iftDilateBin(ero,&seed,radius);
  open = iftRemFrame(aux,iftRound(radius)+5);  
  iftDestroyImage(&aux);
  iftDestroyImage(&ero);
  iftDestroySet(&seed);

  return(open);
}

iftImage *iftAsfOCBin(const iftImage *bin, float radius)
{
  iftImage *open=NULL,*ero=NULL;
  iftSet   *seed=NULL;

  iftImage *aux = iftAddFrame(bin,iftRound(radius)+5,0);
  ero   = iftErodeBin(aux,&seed,radius);
  open  = iftDilateBin(ero,&seed,2.0*radius);
  iftDestroyImage(&ero);
  iftDestroyImage(&aux);
  aux   = iftErodeBin(open,&seed,radius);
  ero = iftRemFrame(aux,iftRound(radius)+5);
  iftDestroyImage(&open);
  iftDestroySet(&seed);
  iftDestroyImage(&aux);

  return(ero);
}

iftImage *iftAsfCOBin(const iftImage *bin, float radius)
{
  iftImage *close=NULL,*dil=NULL;
  iftSet   *seed=NULL;

  iftImage *aux = iftAddFrame(bin,iftRound(radius)+5,0);
  dil    = iftDilateBin(aux,&seed,radius);
  close  = iftErodeBin(dil,&seed,2.0f*radius);
  iftDestroyImage(&dil);
  iftDestroyImage(&aux);
  aux    = iftDilateBin(close,&seed,radius);
  dil = iftRemFrame(aux,iftRound(radius)+5);
  iftDestroyImage(&close);
  iftDestroySet(&seed);
  iftDestroyImage(&aux);

  return(dil);
}

iftImage *iftCloseRecBin(const iftImage *bin, float radius)
{
  iftImage *close=NULL,*crec=NULL;

  iftImage *aux1 = iftAddFrame(bin,iftRound(radius)+5,0);
  close = iftCloseBin(aux1,radius);
  iftImage *aux2 = iftSuperiorRec(aux1,close,NULL);
  iftDestroyImage(&close);
  crec = iftRemFrame(aux2,iftRound(radius)+5);
  iftDestroyImage(&aux1);
  iftDestroyImage(&aux2);
  
  return(crec);
}


iftImage *iftOpenRecBin(const iftImage *bin, float radius)
{
  iftImage *orec=NULL,*aux[2];

  aux[0] = iftComplement(bin);
  aux[1] = iftCloseRecBin(aux[0],radius);
  orec   = iftComplement(aux[1]);
  iftDestroyImage(&aux[0]);
  iftDestroyImage(&aux[1]);

  return(orec);
}

iftImage *iftMorphGrad(const iftImage *img, const iftAdjRel *A)
{
  iftImage *grad=iftCreateImage(img->xsize,img->ysize,img->zsize);

  if (iftIsColorImage(img)){
#pragma omp parallel for shared(img,grad,A)
    for (int p=0; p < img->n; p++) {
      iftVoxel u = iftGetVoxelCoord(img,p);
      int min1 = img->val[p]; int max1 = min1;
      int min2 = img->Cb[p];  int max2 = min2;
      int min3 = img->Cr[p];  int max3 = min3;
      for (int i=1; i < A->n; i++) {
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(img,v)){
	  int q = iftGetVoxelIndex(img,v);
	  if (img->val[q]>max1) max1 = img->val[q];
	  if (img->val[q]<min1) min1 = img->val[q];
	  if (img->Cb[q]>max2)  max2 = img->Cb[q];
	  if (img->Cb[q]<min2)  min2 = img->Cb[q];
	  if (img->Cr[q]>max3)  max3 = img->Cr[q];
	  if (img->Cr[q]<min3)  min3 = img->Cr[q];
	}
      }
      grad->val[p]=iftMax(iftMax(max1-min1,max2-min2),max3-min3);
    }
  } else {
#pragma omp parallel for shared(img,grad,A)
    for (int p=0; p < img->n; p++) {
      iftVoxel u = iftGetVoxelCoord(img,p);
      int min = img->val[p]; int max = min;
      for (int i=1; i < A->n; i++) {
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(img,v)){
	  int q = iftGetVoxelIndex(img,v);
	  if (img->val[q]>max) max = img->val[q];
	  if (img->val[q]<min) min = img->val[q];
	}
      }
      grad->val[p]=max-min;
    }
  }
  
  iftCopyVoxelSize(img,grad);

  return(grad);
}

iftImage *iftSuperiorRec(const iftImage *img, const iftImage *marker, const iftImage *mask)
{
  iftImage   *srec=iftCopyImage(marker);
  iftGQueue  *Q=NULL;
  int         i,p,q,tmp;
  iftVoxel    u,v;
  iftAdjRel  *A;

  iftVerifyImageDomains(img,marker,"iftSuperiorRec");
  if (mask !=NULL)
    iftVerifyImageDomains(img,mask,"iftSuperiorRec");
    
  // Initialization 

  if (iftIs3DImage(img))
    A = iftSpheric(1.0);
  else
    A = iftCircular(1.0);

  Q     = iftCreateGQueue(iftMaximumValue(srec)+1,img->n,srec->val);

  for (p=0; p < img->n; p++) {
    if ((mask==NULL)||(mask->val[p]!=0)){
      iftInsertGQueue(&Q,p);
    }
  }
  
  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);
    u=iftGetVoxelCoord(img,p);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){	
	q = iftGetVoxelIndex(img,v);
	if ((mask==NULL)||(mask->val[q]!=0)){
	  if (srec->val[q] > srec->val[p]){
	    tmp = iftMax(srec->val[p], img->val[q]);
	    if (tmp < srec->val[q]){           
	      iftRemoveGQueueElem(Q,q);
	      srec->val[q] = tmp;
	      iftInsertGQueue(&Q, q); 
	    }
	  }
	}
      }
    }
  }
  
  iftDestroyGQueue(&Q);
  iftDestroyAdjRel(&A);

  return(srec);
}

iftImage *iftInferiorRec(const iftImage *img, const iftImage *marker, const iftImage *mask)
{
  iftImage   *irec=iftCopyImage(marker);
  iftGQueue  *Q=NULL;
  int         i,p,q,tmp;
  iftVoxel    u,v;
  iftAdjRel  *A;

  iftVerifyImageDomains(img,marker,"iftInferiorRec");
  if (mask !=NULL)
    iftVerifyImageDomains(img,mask,"iftInferiorRec");
    
  // Initialization 

  if (iftIs3DImage(img))
    A = iftSpheric(1.0);
  else
    A = iftCircular(1.0);
    
  Q     = iftCreateGQueue(iftMaximumValue(img)+1,img->n,irec->val);

  for (p=0; p < img->n; p++) {
    if ((mask==NULL)||(mask->val[p]!=0)){
      iftInsertGQueue(&Q,p);
    }
  }
  
  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);
    u=iftGetVoxelCoord(img,p);
    
    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){	
	q = iftGetVoxelIndex(img,v);
	if ((mask==NULL)||(mask->val[q]!=0)){
	  if (irec->val[q] < irec->val[p]){
	    tmp = iftMin(irec->val[p], img->val[q]);
	    if (tmp > irec->val[q]){           
	      iftRemoveGQueueElem(Q,q);
	      irec->val[q] = tmp;
	      iftInsertGQueue(&Q, q); 
	    }
	  }
	}
      }
    }
  }
  
  iftDestroyGQueue(&Q);
  iftDestroyAdjRel(&A);

  return(irec);
}


iftImage *iftCloseRec(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage   *crec=NULL,*marker=NULL;

  marker = iftClose(img,A,mask);
  crec   = iftSuperiorRec(img,marker,mask);
  iftDestroyImage(&marker);

  return(crec);
}


iftImage *iftOpenRec(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage   *marker, *orec=NULL;

  marker = iftOpen(img,A,mask);
  orec   = iftInferiorRec(img,marker,mask);
  iftDestroyImage(&marker);
  
  return(orec);
}

iftImage *iftAsfOCRec(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage *open=NULL,*close=NULL;

  open  = iftOpenRec(img,A,mask);
  close = iftCloseRec(open,A,mask);
  iftDestroyImage(&open);

  return(close);
}

iftImage *iftAsfCORec(const iftImage *img, const iftAdjRel *A, const iftImage *mask)
{
  iftImage *open=NULL,*close=NULL;

  close = iftCloseRec(img,A,mask);
  open  = iftOpenRec(close,A,mask);
  iftDestroyImage(&close);

  return(open);
}

iftImage *iftHClose(const iftImage *img, int H, const iftImage *mask)
{
  iftImage *marker,*suprec;

  if (H < 0)
      iftError("Invalid depth for closing basins", "iftHBasins");

  marker   = iftAddValue(img,H); 
  suprec   = iftSuperiorRec(img,marker,mask);
  iftDestroyImage(&marker);

  return(suprec);
}

iftImage *iftHOpen(const iftImage *img, int H, const iftImage *mask)
{
  iftImage *infrec,*suprec;

  if (H < 0)
      iftError("Invalid height for opening domes", "iftHDomes");

  iftImage *aux = iftComplement(img);
  suprec = iftHClose(aux,H,mask);
  iftDestroyImage(&aux);
  infrec = iftComplement(suprec);
  iftDestroyImage(&suprec);
    
  return(infrec);
}


iftImage *iftCloseBasins(const iftImage *img, iftSet *seed, const iftImage *_mask)
{
  iftImage   *cbas=NULL;
  iftGQueue  *Q=NULL;
  int         i,p,q,tmp;
  iftSet     *S=NULL;
  iftVoxel    u,v;
  iftAdjRel  *A;
  char        no_mask=0, no_seed=0;


  /* Initialization */

  if (iftIs3DImage(img)) 
    A = iftSpheric(1.0);    
  else
    A = iftCircular(1.0);

  iftImage *mask = NULL;
  if (_mask == NULL){ /* default region */
    mask = iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftSetImage(mask,1);
    no_mask=1;
  } else {
      iftVerifyImageDomains(img,_mask,"iftCloseBasins");
      mask = iftCopyImage(_mask);
  }

  if (seed == NULL){ /* default seed set */
    if (no_mask){
      seed=iftImageBorderSet(img);
    }else{
      seed=iftImageBorderSet(mask);
    }
    no_seed = 1;
  } 
    
  cbas  = iftCreateImage(img->xsize,img->ysize,img->zsize);
  if (iftIsColorImage(img)) 
    iftCopyCbCr(img,cbas);
  iftCopyVoxelSize(img,cbas);

  Q      = iftCreateGQueue(iftMaximumValue(img)+1,img->n,cbas->val);

  for (p=0; p < cbas->n; p++) {
    if(mask->val[p]!=0)
      cbas->val[p]= IFT_INFINITY_INT;    
  }

  S = seed;
  while (S != NULL){
    p = S->elem;
    if (mask->val[p]==0)
        iftError("Seed set is outside mask", "iftCloseBasins");
    cbas->val[p]=img->val[p];
    iftInsertGQueue(&Q,p);
    S = S->next;
  }


  /* Image Foresting Transform */

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    u = iftGetVoxelCoord(img,p);

    for (i=1; i < A->n; i++){
      v  = iftGetAdjacentVoxel(A,u,i);

      if (iftValidVoxel(img,v)){	
	q = iftGetVoxelIndex(img,v);
	if (mask->val[q]!=0){
	  if (cbas->val[q] > cbas->val[p]){
	    tmp = iftMax(cbas->val[p], img->val[q]);	    
	    if (tmp < cbas->val[q]){ /* This is a particular case of
					connectivity function that
					implies q has never been
					inserted in Q */ 
	      cbas->val[q]  = tmp;
	      iftInsertGQueue(&Q, q); 
	    }
	  }
	}
      }
    }
  }
  iftDestroyGQueue(&Q);
  iftDestroyAdjRel(&A);

  iftDestroyImage(&mask);

  if (no_seed)
    iftDestroySet(&seed);

  return(cbas);
}

iftImage *iftOpenDomes(const iftImage *img, iftSet *seed, const iftImage *_mask)
{
  iftImage   *odom=NULL;
  iftGQueue  *Q=NULL;
  int         i,p,q,tmp;
  iftSet     *S=NULL;
  iftVoxel    u,v;
  iftAdjRel  *A;
  char        no_mask=0, no_seed=0;

  /* Initialization */

  if (iftIs3DImage(img)) 
    A = iftSpheric(1.0);    
  else
    A = iftCircular(1.0);

  iftImage *mask = NULL;
  if (_mask == NULL){ /* default region */
    mask = iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftSetImage(mask,1);
    no_mask=1;
  } else {
    iftVerifyImageDomains(img,_mask,"iftOpenDomes");
    mask = iftCopyImage(_mask);
  }

  if (seed == NULL){ /* default seed set */
    if (no_mask){
      seed=iftImageBorderSet(img);
    }else{
      seed=iftImageBorderSet(mask);
    }
    no_seed = 1;
  } 
    
  odom  = iftCreateImage(img->xsize,img->ysize,img->zsize);
  if (iftIsColorImage(img)) 
    iftCopyCbCr(img,odom);
  iftCopyVoxelSize(img,odom);

  Q      = iftCreateGQueue(iftMaximumValue(img)+1,img->n,odom->val);
  iftSetRemovalPolicy(Q,MAXVALUE);

  for (p=0; p < odom->n; p++) {
    if(mask->val[p]!=0)
      odom->val[p]= IFT_INFINITY_INT_NEG;    
  }

  S = seed;
  while (S != NULL){
    p = S->elem;
    if (mask->val[p]==0)
        iftError("Seed set is outside mask", "iftCloseBasins");
    odom->val[p]=img->val[p];
    iftInsertGQueue(&Q,p);
    S = S->next;
  }


  /* Image Foresting Transform */

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    u = iftGetVoxelCoord(img,p);

    for (i=1; i < A->n; i++){
      v  = iftGetAdjacentVoxel(A,u,i);

      if (iftValidVoxel(img,v)){	
	q = iftGetVoxelIndex(img,v);
	if (mask->val[q]!=0){
	  if (odom->val[q] < odom->val[p]){
	    tmp = iftMin(odom->val[p], img->val[q]);	    
	    if (tmp > odom->val[q]){ /* This is a particular case of
					connectivity function that
					implies q has never been
					inserted in Q */ 
	      odom->val[q]  = tmp;
	      iftInsertGQueue(&Q, q); 
	    }
	  }
	}
      }
    }
  }
  iftDestroyGQueue(&Q);
  iftDestroyAdjRel(&A);

  iftDestroyImage(&mask);

  if (no_seed)
    iftDestroySet(&seed);

  return(odom);
}


/* Area opening via union-find algorithm */

iftImage *iftFastAreaOpen(const iftImage *img, int thres)
{
  int i, p, q, ip, iq;
  int *sortedVoxel, *area, *parent;
  iftVoxel *voxArray1, *voxArray2, u, v;
  iftAdjRel *A;
  iftImage *areaOpen;

  if (img->zsize==1) 
    A = iftCircular(1.0);
  else
    A = iftSpheric(1.0);
    

  area   = iftAllocIntArray(img->n);
  parent = iftAllocIntArray(img->n);
  
  /* Sort values by RadixSort */

  voxArray1 = (iftVoxel *) iftAlloc(img->n, sizeof(iftVoxel));
  voxArray2 = (iftVoxel *) iftAlloc(img->n, sizeof(iftVoxel));

  for (p = 0; p < img->n; p++) {
    voxArray1[p].x = img->val[p]; // x = voxel brightness
    voxArray1[p].y = p; // y = voxel index
  }
  sortedVoxel = iftAllocIntArray(img->n);
  /* if image is 8 bits, then call radixSort only once. */
  if (iftMaximumValue(img) <= 255) {
    iftRadixSort(0, img->n, voxArray1, voxArray2);
    /* store in sortedVoxel the indices of the voxels in the
       decreasing order. */
    for (p = img->n - 1; p >= 0; p--) {
      sortedVoxel[p] = voxArray2[img->n - p - 1].y;
    }
  } else {
    iftRadixSort(0, img->n, voxArray1, voxArray2);
    iftRadixSort(1, img->n, voxArray2, voxArray1);
    /* store in sortedVoxel the indices of the voxels in the
       decreasing order. */
    for (p = img->n - 1; p >= 0; p--) {
      sortedVoxel[p] = voxArray1[img->n - p - 1].y;
    }
  }
  iftFree(voxArray1);
  iftFree(voxArray2);

  static const int active = -1;
  static const int inactive = -2;

  /* process voxels in the decreasing order of brightness. */

  parent[sortedVoxel[0]] = active;
  area[sortedVoxel[0]] = 1;
  for (p = 1; p < img->n; p++) {
    ip = sortedVoxel[p];
    if (img->val[ip] != img->val[sortedVoxel[p - 1]]) {
      q = p - 1;
      while ((q >= 0) && (img->val[sortedVoxel[q]] == img->val[sortedVoxel[p - 1]])) {
	iq = sortedVoxel[q];
	if ((parent[iq] == active) && (area[iq] > thres)) {
	  parent[iq] = inactive;
	  area[iq]   = 1;
	}
	q--;
      }
    }
    parent[ip] = active;
    area[ip]   = 1;
    

    /* Process the neighbors of ip already visited. */

    u = iftGetVoxelCoord(img,ip);

    for (i = 1; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img, v)) {
	q = iftGetVoxelIndex(img,v);
	/* if already visited, try union */
	if (area[q] != 0) {
	  
	  iq = iftFindRoot(parent, q);
	  if (iq != ip) {
	    if ((img->val[iq] == img->val[ip]) || (parent[iq] == active)) {

	      
	      if ((parent[ip] == active) && (parent[iq] == active)) {
		area[ip] += area[iq]; // merge de x e y em y
		area[iq] = 1;
	      } else if (parent[iq] == active) {
		area[iq] = 1;
	      } else {
		area[ip] = 1;
		parent[ip] = inactive;
	      }
	      parent[iq] = ip;
	      
	    } else if (parent[ip] == active) {
	      parent[ip] = inactive;
	      area[ip] = 1;
	    }
	  }

	}
      }
    }
  }

  areaOpen = iftCreateImage(img->xsize, img->ysize, img->zsize);
  if (iftIsColorImage(img)) 
    iftCopyCbCr(img,areaOpen);
  iftCopyVoxelSize(img,areaOpen);
  for (p = img->n - 1; p >= 0; p--) {
    q = sortedVoxel[p];
    if (parent[q] < 0) {
      areaOpen->val[q] = img->val[q];
    } else {
      areaOpen->val[q] = areaOpen->val[parent[q]];
    }
  }
  
  iftFree(area);
  iftFree(parent);
  iftFree(sortedVoxel);
  iftDestroyAdjRel(&A);

  return areaOpen;
}

/* Area closing via union-find algorithm */

iftImage *iftFastAreaClose(const iftImage *img, int thres) {

  int i, p, q, ip, iq;
  int *sortedVoxel, *area, *parent;
  iftVoxel *voxArray1, *voxArray2, u, v;
  iftAdjRel *A;
  iftImage *areaClose;

  if (img->zsize==1) 
    A = iftCircular(1.0);
  else
    A = iftSpheric(1.0);
    

  area   = iftAllocIntArray(img->n);
  parent = iftAllocIntArray(img->n);
  
  /* Sort values by RadixSort */

  voxArray1 = (iftVoxel *) iftAlloc(img->n, sizeof(iftVoxel));
  voxArray2 = (iftVoxel *) iftAlloc(img->n, sizeof(iftVoxel));

  for (p = 0; p < img->n; p++) {
    voxArray1[p].x = img->val[p]; // x = voxel brightness
    voxArray1[p].y = p; // y = voxel index
  }
  sortedVoxel = iftAllocIntArray(img->n);
  /* if image is 8 bits, then call radixSort only once. */
  if (iftMaximumValue(img) <= 255) {
    iftRadixSort(0, img->n, voxArray1, voxArray2);
    /* store in sortedVoxel the indices of the voxels in the
       increasing order. */
    for (p = img->n - 1; p >= 0; p--) {
      sortedVoxel[p] = voxArray2[p].y;
    }
  } else {
    iftRadixSort(0, img->n, voxArray1, voxArray2);
    iftRadixSort(1, img->n, voxArray2, voxArray1);
    /* store in sortedVoxel the indices of the voxels in the
       increasing order. */
    for (p = img->n - 1; p >= 0; p--) {
      sortedVoxel[p] = voxArray1[p].y;
    }
  }
  iftFree(voxArray1);
  iftFree(voxArray2);

  static const int active = -1;
  static const int inactive = -2;

  /* process voxels in the decreasing order of brightness. */

  parent[sortedVoxel[0]] = active;
  area[sortedVoxel[0]]   = 1;
  for (p = 1; p < img->n; p++) {
    ip = sortedVoxel[p];
    if (img->val[ip] != img->val[sortedVoxel[p - 1]]) {
      q = p - 1;
      while ((q >= 0) && (img->val[sortedVoxel[q]] == img->val[sortedVoxel[p - 1]])) {
	iq = sortedVoxel[q];
	if ((parent[iq] == active) && (area[iq] > thres)) {
	  parent[iq] = inactive;
	  area[iq]   = 1;
	}
	q--;
      }
    }
    parent[ip] = active;
    area[ip]   = 1;
    

    /* Process the neighbors of ip already visited. */

    u = iftGetVoxelCoord(img,ip);

    for (i = 1; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img, v)) {
	q = iftGetVoxelIndex(img,v);
	/* if already visited, try union */
	if (area[q] != 0) {
	  
	  iq = iftFindRoot(parent, q);
	  if (iq != ip) {
	    if ((img->val[iq] == img->val[ip]) || (parent[iq] == active)) {

	      
	      if ((parent[ip] == active) && (parent[iq] == active)) {
		area[ip] += area[iq]; // merge de x e y em y
		area[iq] = 1;
	      } else if (parent[iq] == active) {
		area[iq] = 1;
	      } else {
		area[ip] = 1;
		parent[ip] = inactive;
	      }
	      parent[iq] = ip;
	      
	    } else if (parent[ip] == active) {
	      parent[ip] = inactive;
	      area[ip] = 1;
	    }
	  }

	}
      }
    }
  }

  areaClose = iftCreateImage(img->xsize, img->ysize, img->zsize);
  if (iftIsColorImage(img)) 
    iftCopyCbCr(img,areaClose);
  iftCopyVoxelSize(img,areaClose);

  for (p = img->n - 1; p >= 0; p--) {
    q = sortedVoxel[p];
    if (parent[q] < 0) {
      areaClose->val[q] = img->val[q];
    } else {
      areaClose->val[q] = areaClose->val[parent[q]];
    }
  }
  
  iftFree(area);
  iftFree(parent);
  iftFree(sortedVoxel);
  iftDestroyAdjRel(&A);

  return areaClose;
}

iftImage *iftAreaClose(const iftImage *img, int area_thres, iftCompTree *ctree)
{
  int i,p;
  iftImage *fimg;
  int *level=NULL;

  if (ctree == NULL) /* build the component tree from the minima */
    ctree = iftCreateMinTree(img);
  
  iftCumSize(ctree,ctree->root);
  level = iftAllocIntArray(ctree->numnodes);
  for (i=0; i < ctree->numnodes; i++) 
    level[i]=ctree->node[i].level;

  for (i=0; i < ctree->numnodes; i++) 
    if (ctree->node[i].numsons==0)
      level[i]=iftAreaLevel(ctree,level,i,area_thres);
  fimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
  if (iftIsColorImage(img)) 
    iftCopyCbCr(img,fimg);
  iftCopyVoxelSize(img,fimg);

  for (p=0; p < img->n; p++) 
    fimg->val[p]=level[ctree->cmap->val[p]];
  iftDestroyCompTree(&ctree);
  iftFree(level);

  return(fimg);
}


iftImage *iftVolumeClose(const iftImage *img, int volume_thres, iftCompTree *ctree)
{
  int i,p;
  iftImage *fimg;
  int *level=NULL;

  if (ctree == NULL) /* Build component tree from the minima */
    ctree = iftCreateMinTree(img);
  
  iftCumSize(ctree,ctree->root);
  level = iftAllocIntArray(ctree->numnodes);
  for (i=0; i < ctree->numnodes; i++) 
    level[i]=ctree->node[i].level;
  for (i=0; i < ctree->numnodes; i++) 
    if (ctree->node[i].numsons==0)
      level[i]=iftVolumeLevel(ctree,level,i,volume_thres,0);    
  fimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
  if (iftIsColorImage(img)) 
    iftCopyCbCr(img,fimg);
  iftCopyVoxelSize(img,fimg);
  for (p=0; p < img->n; p++) 
    fimg->val[p]=level[ctree->cmap->val[p]];
  iftDestroyCompTree(&ctree);
  iftFree(level);

  return(fimg);
}

iftImage *iftAreaOpen(const iftImage *img, int area_thres, iftCompTree *ctree)
{
  int i,p;
  iftImage  *fimg;
  int *level=NULL;

  if (ctree == NULL) /* Build component tree from the maxima */
    ctree = iftCreateMaxTree(img);
  
  iftCumSize(ctree,ctree->root);
  level = iftAllocIntArray(ctree->numnodes);
  for (i=0; i < ctree->numnodes; i++) 
    level[i]=ctree->node[i].level;
  
  for (i=0; i < ctree->numnodes; i++) 
    if (ctree->node[i].numsons==0)
      level[i]=iftAreaLevel(ctree,level,i,area_thres);
  fimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
  if (iftIsColorImage(img)) 
    iftCopyCbCr(img,fimg);
  iftCopyVoxelSize(img,fimg);
  for (p=0; p < img->n; p++) 
    fimg->val[p]=level[ctree->cmap->val[p]];
  iftDestroyCompTree(&ctree);
  iftFree(level);

  return(fimg);
}

iftImage *iftVolumeOpen(const iftImage *img, int volume_thres, iftCompTree *ctree)
{
  int i,p;
  iftImage *fimg;
  int *level=NULL;

  if (ctree == NULL) /* Build component tree from the maxima */
    ctree = iftCreateMaxTree(img);

  iftCumSize(ctree,ctree->root);
  level = iftAllocIntArray(ctree->numnodes);
  for (i=0; i < ctree->numnodes; i++) 
    level[i]=ctree->node[i].level;
  for (i=0; i < ctree->numnodes; i++) 
    if (ctree->node[i].numsons==0)
      level[i]=iftVolumeLevel(ctree,level,i,volume_thres,0);
  fimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
  if (iftIsColorImage(img)) 
    iftCopyCbCr(img,fimg);
  iftCopyVoxelSize(img,fimg);
  for (p=0; p < img->n; p++) 
    fimg->val[p]=level[ctree->cmap->val[p]];
  iftDestroyCompTree(&ctree);
  iftFree(level);

  return(fimg);
}


