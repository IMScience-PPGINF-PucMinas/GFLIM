#include "iftGeometric.h"

#include "ift/core/io/Stream.h"

iftImage *iftTransformImageByMatrix(const iftImage *img, const iftMatrix *M) {
    if (img == NULL)
        iftError("Input Image is NULL", "iftTransformImageByMatrix");
    if (M == NULL)
        iftError("Transformation Matrix is NULL", "iftTransformImageByMatrix");

    iftImage  *rimg;
    iftImage *Cb = NULL, *Cr = NULL;
    iftMatrix *InvM, *E, *InvE;
    iftPoint P1,P2;
    iftVoxel v,u;
    int q;
    int size = iftDiagonalSize(img);

    E    = iftExtendMatrix(M, 4, 4);
    InvM = iftInvertMatrix(M);
    InvE = iftExtendMatrix(InvM, 4, 4);

    if (iftIsColorImage(img)){
        Cb = iftImageCb(img);
        Cr = iftImageCr(img);
    }

    if (img->zsize > 1) { // 3D
        rimg = iftCreateImage(size, size, size);
        if (iftIsColorImage(img))
            iftSetCbCr(rimg,iftMaxImageRange(iftImageDepth(img)+1)/2);
        for (v.z = 0; v.z < rimg->zsize; v.z++)
            for (v.y = 0; v.y < rimg->ysize; v.y++)
                for (v.x = 0; v.x < rimg->xsize; v.x++) {
                    q    = iftGetVoxelIndex(rimg, v);
                    P1.x = v.x - rimg->xsize/2.0;
                    P1.y = v.y - rimg->ysize/2.0;
                    P1.z = v.z - rimg->zsize/2.0;
                    P2   = iftTransformPoint(InvE, P1);
                    P2.x = P2.x + img->xsize/2.0;
                    P2.y = P2.y + img->ysize/2.0;
                    P2.z = P2.z + img->zsize/2.0;
                    u.x  = iftRound(P2.x);
                    u.y  = iftRound(P2.y);
                    u.z  = iftRound(P2.z);
                    if (iftValidVoxel(img, u)) {
                        rimg->val[q] = iftImageValueAtPoint(img, P2);
                        if (iftIsColorImage(img)) {
                            rimg->Cb[q] = (ushort) iftImageValueAtPoint(Cb, P2);
                            rimg->Cr[q] = (ushort) iftImageValueAtPoint(Cr, P2);
                        }
                    }
                }
    }
    else { // 2D
        rimg = iftCreateImage(size, size, 1);
        if (iftIsColorImage(img)){
	  iftSetCbCr(rimg,iftMaxImageRange(iftImageDepth(img))/2+1);
	}
        v.z = u.z = P1.z = P2.z = 0;
        for (v.y = 0; v.y < rimg->ysize; v.y++)
            for (v.x = 0; v.x < rimg->xsize; v.x++) {
                q    = iftGetVoxelIndex(rimg, v);
                P1.x = v.x - rimg->xsize/2.0;
                P1.y = v.y - rimg->ysize/2.0;
                P2   = iftTransformPoint(InvE, P1);
                P2.x = P2.x + img->xsize/2.0;
                P2.y = P2.y + img->ysize/2.0;
                u.x  = iftRound(P2.x);
                u.y  = iftRound(P2.y);
                if (iftValidVoxel(img, u)) {
                  rimg->val[q] = iftImageValueAtPoint2D(img, P2);
                    if (iftIsColorImage(img)) {
                        rimg->Cb[q] = (ushort) iftImageValueAtPoint2D(Cb, P2);
                        rimg->Cr[q] = (ushort) iftImageValueAtPoint2D(Cr, P2);
                    }
                }
            }
    }

    iftDestroyImage(&Cb);
    iftDestroyImage(&Cr);
    
    P1.x = img->dx;
    P1.y = img->dy;
    P1.z = img->dz;
    
    P2 = iftTransformPoint(E, P1);
    
    rimg->dx = P2.x;
    rimg->dy = P2.y;
    rimg->dz = P2.z;
    
    
    iftDestroyMatrix(&InvM);
    iftDestroyMatrix(&InvE);
    iftDestroyMatrix(&E);
    
    return(rimg);
}

iftImage *iftAlignObject(iftImage *bin)
{
  iftDataSet *Z   = iftObjectToDataSet(bin);
  iftSetStatus(Z,IFT_TRAIN);
  iftDataSet *Zc  = iftCentralizeDataSet(Z);
  iftMatrix  *A   = iftDatasetCovarianceMatrix(Zc);
  iftMatrix  *U, *S, *Vt;
  iftSingleValueDecomp(A,&U,&S,&Vt);
  iftImage *obj   = iftTransformImageByMatrix(bin,Vt);
  int Imax        = iftMaximumValue(bin);
  if (Imax > 1) {
    for (int p=0; p < obj->n; p++)
      if (obj->val[p] >= Imax/2)
	obj->val[p] = Imax;
      else
	obj->val[p] = 0;
  }
  
  iftDestroyDataSet(&Z);
  iftDestroyDataSet(&Zc);
  iftDestroyMatrix(&A);
  iftDestroyMatrix(&U);
  iftDestroyMatrix(&S);
  iftDestroyMatrix(&Vt);

  return(obj);
}

iftImage *iftRotateImage2D(iftImage *img, float theta)
{
  iftImage  *rimg;
  iftMatrix *R=iftRotationMatrix(IFT_AXIS_Z, theta),*InvR=iftTransposeMatrix(R);
  iftPoint  P1,P2;
  iftVoxel   v,u;
  int        q, xsize, ysize;

  if (img->zsize != 1)
      iftError("Image must be 2D", "iftRotateImage2D");

  xsize=(int)(fabs(img->xsize*cos(theta * IFT_PI / 180.0)) + fabs(img->ysize * sin(theta * IFT_PI / 180.0)));
  ysize=(int)(fabs(img->xsize*sin(theta * IFT_PI / 180.0)) + fabs(img->ysize * cos(theta * IFT_PI / 180.0)));
  
  rimg = iftCreateImage(xsize,ysize,1);
  P2.z = P1.z = u.z = v.z = 0;

  if (iftIsColorImage(img)) {
    iftSetCbCr(rimg,(iftMaxImageRange(iftImageDepth(img))+1)/2);

    iftImage *Cb = iftImageCb(img);
    iftImage *Cr = iftImageCr(img);

    for (v.y = 0; v.y < rimg->ysize; v.y++){ 
      for (v.x = 0; v.x < rimg->xsize; v.x++){
	q    = iftGetVoxelIndex(rimg,v);
	P1.x = v.x - rimg->xsize/2.0;
	P1.y = v.y - rimg->ysize/2.0;
	P2   = iftTransformPoint(InvR,P1);
	P2.x = P2.x + img->xsize/2.0;
	P2.y = P2.y + img->ysize/2.0;
	u.x  = iftRound(P2.x);
	u.y  = iftRound(P2.y);	
	if (iftValidVoxel(img,u)){	  
	  int value    = iftImageValueAtPoint2D(img,P2);
	  rimg->val[q] = value;
	  rimg->Cb[q]  = (ushort) iftImageValueAtPoint2D(Cb,P2);
	  rimg->Cr[q]  = (ushort) iftImageValueAtPoint2D(Cr,P2);
	}
      }
    }
    iftDestroyImage(&Cb);
    iftDestroyImage(&Cr);
  }else{
    for (v.y = 0; v.y < rimg->ysize; v.y++) 
      for (v.x = 0; v.x < rimg->xsize; v.x++){
	q    = iftGetVoxelIndex(rimg,v);
	P1.x = v.x - rimg->xsize/2.0;
	P1.y = v.y - rimg->ysize/2.0;
	P2   = iftTransformPoint(InvR,P1);
	P2.x = P2.x + img->xsize/2.0;
	P2.y = P2.y + img->ysize/2.0;
	u.x  = iftRound(P2.x);
	u.y  = iftRound(P2.y);
	if (iftValidVoxel(img,u)){
	  rimg->val[q] = iftImageValueAtPoint2D(img,P2);
	}
      }
  }

  P1.x = img->dx;
  P1.y = img->dy;
  P1.z = img->dz;

  
  P2 = iftTransformPoint(R,P1);

  rimg->dx = P2.x;
  rimg->dy = P2.y;
  rimg->dz = P2.z;

  iftDestroyMatrix(&R);
  iftDestroyMatrix(&InvR);

  return(rimg);
}


iftImage *iftRotateImage(const iftImage *img, float theta_x, float theta_y) {
    if (img->zsize == 1)
        iftError("Use iftRotateImage2D", "iftRotateImage");
    
    iftMatrix *Rx = iftRotationMatrix(IFT_AXIS_X, theta_x);
    iftMatrix *InvRx = iftTransposeMatrix(Rx);
    iftMatrix *Ry = iftRotationMatrix(IFT_AXIS_Y, theta_y);
    iftMatrix *InvRy = iftTransposeMatrix(Ry);
    iftMatrix *InvM = iftMultMatrices(InvRx, InvRy);
    
    int xsize, ysize, zsize;
    xsize = ysize = zsize = iftDiagonalSize(img);
    
    iftImage *rimg = iftCreateImage(xsize, ysize, zsize);

    #pragma omp parallel for
    for (int q = 0; q < rimg->n; q++) {
        iftVoxel v = iftGetVoxelCoord(rimg, q);
        
        iftPoint Pv;
        Pv.x = v.x - rimg->xsize / 2.0;
        Pv.y = v.y - rimg->ysize / 2.0;
        Pv.z = v.z - rimg->zsize / 2.0;
        
        iftPoint Pu = iftTransformPoint(InvM, Pv);
        Pu.x = Pu.x + img->xsize / 2.0;
        Pu.y = Pu.y + img->ysize / 2.0;
        Pu.z = Pu.z + img->zsize / 2.0;
        iftVoxel u;
        u.x = iftRound(Pu.x);
        u.y = iftRound(Pu.y);
        u.z = iftRound(Pu.z);
        
        if (iftValidVoxel(img, u))
            iftImgVoxelVal(rimg, v) = iftImageValueAtPoint(img, Pu);
    }
    
    iftPoint P1 = {.x = img->dx, .y = img->dy, .z = img->dz};
    
    iftMatrix *M = iftMultMatrices(Ry, Rx);
    iftPoint P2 = iftTransformPoint(M, P1);
    iftDestroyMatrix(&M);
    
    rimg->dx = P2.x;
    rimg->dy = P2.y;
    rimg->dz = P2.z;
    
    iftDestroyMatrix(&Rx);
    iftDestroyMatrix(&InvRx);
    iftDestroyMatrix(&Ry);
    iftDestroyMatrix(&InvRy);
    iftDestroyMatrix(&InvM);
    
    return rimg;
}


iftImage *iftScaleImage(const iftImage *img, float sx, float sy, float sz) {
    if (img->zsize == 1) {
        iftWarning("Input Image is 2D. Then, using function iftScaleImage2D instead", "iftScaleImage");
        return iftScaleImage2D(img, sx, sy);
    }
    if ((sx == 0.0) || (sy == 0.0) || (sz == 0.0))
        iftError("Invalid scale factor(s): (sx, sy, sz) = (%f, %f, %f)... Try > 0.0 for all scale factors",
                 "iftScaleImage2D", sx, sy);

    iftImage dom;
    dom.xsize = iftRound(fabs(sx * img->xsize));
    dom.ysize = iftRound(fabs(sy * img->ysize));
    dom.zsize = iftRound(fabs(sz * img->zsize));

    iftMatrix *InvS = iftScaleMatrix(1.0/sx, 1.0/sy, 1.0/sz);
    iftImage *simg  = iftCreateImage(dom.xsize, dom.ysize, dom.zsize);
    simg->dx        = img->dx/fabs(sx);
    simg->dy        = img->dy/fabs(sy);
    simg->dz        = img->dz/fabs(sz);

    iftVoxel v;
    for (v.z = 0; v.z < simg->zsize; v.z++)
        for (v.y = 0; v.y < simg->ysize; v.y++) 
            for (v.x = 0; v.x < simg->xsize; v.x++) {
                int q = iftGetVoxelIndex(simg, v);
                iftPoint P1, P2;
                P1.x = v.x;
                P1.y = v.y;
                P1.z = v.z;
                P2   = iftTransformPoint(InvS, P1);
                
                iftVoxel u;
                u.x  = iftRound(P2.x);
                u.y  = iftRound(P2.y);
                u.z  = iftRound(P2.z);
                
                if (iftValidVoxel(img, u)) {
                    simg->val[q] = iftImageValueAtPoint(img, P2);
                }
            }

    iftDestroyMatrix(&InvS);

  return simg;
}

iftImage *iftScaleImage2D(const iftImage *img, float sx, float sy) {
    if (img->zsize != 1)
        iftError("Input Image is 3. You must use the function iftScaleImage instead", "iftScaleImage2D");
    if ((sx == 0.0) || (sy == 0.0))
        iftError("Invalid scale factor(s): (sx, sy) = (%f, %f)... Try > 0.0 for all scale factors", "iftScaleImage2D",
                 sx, sy);

    int xsize = iftRound(fabs(sx * img->xsize));
    int ysize = iftRound(fabs(sy * img->ysize));
    
    iftMatrix *InvS = iftScaleMatrix(1.0/sx, 1.0/sy, 1.0);

    iftImage *simg = iftCreateImage(xsize, ysize, 1);
    simg->dx       = img->dx / fabs(sx);
    simg->dy       = img->dy / fabs(sy);
    simg->dz       = img->dz;

    iftVoxel u, v;
    iftPoint P1, P2;
    P2.z = P1.z = u.z = v.z = 0;
  
    if (img->Cb != NULL) {
        iftSetCbCr(simg, (iftMaxImageRange(iftImageDepth(img)) + 1) / 2);
        iftImage *Cb = iftImageCb(img);
        iftImage *Cr = iftImageCr(img);
        
        for (v.y = 0; v.y < simg->ysize; v.y++) {
            for (v.x = 0; v.x < simg->xsize; v.x++) {
                int q = iftGetVoxelIndex(simg, v);
                P1.x  = v.x;
                P1.y  = v.y;
                P2    = iftTransformPoint(InvS, P1);
                u.x   = iftRound(P2.x);
                u.y   = iftRound(P2.y);
                
                if (iftValidVoxel(img, u)) {
                    simg->val[q] = iftImageValueAtPoint2D(img, P2);
                    simg->Cb[q]  = (ushort) iftImageValueAtPoint2D(Cb, P2);
                    
                    if (simg->Cb[q] == 0)
                        simg->Cb[q] = iftImageDepth(img);
                    simg->Cr[q] = (ushort) iftImageValueAtPoint2D(Cr, P2);
                    
                    if (simg->Cr[q] == 0)
                        simg->Cr[q] = iftImageDepth(img);
                } else {
                    simg->Cb[q]=simg->Cr[q]=iftImageDepth(img);
                }
            }
        }
        iftDestroyImage(&Cb);
        iftDestroyImage(&Cr);
    } else {
        for (v.y = 0; v.y < simg->ysize; v.y++) 
            for (v.x = 0; v.x < simg->xsize; v.x++) {
                int q = iftGetVoxelIndex(simg, v);
                P1.x  = v.x;
                P1.y  = v.y;
                P2    = iftTransformPoint(InvS, P1);
                u.x   = iftRound(P2.x);
                u.y   = iftRound(P2.y);

                if (iftValidVoxel(img, u)) {
                    simg->val[q] = iftImageValueAtPoint2D(img, P2);
                }
            }
    }
    iftDestroyMatrix(&InvS);

    return(simg);
}


iftVoxel iftFlipVoxelOnXAxis(iftVoxel u, int xsize, int ysize, int zsize) {
    iftVoxel v;
    v.x = xsize - 1 - u.x;
    v.y = u.y; v.z = u.z;

    return v;
}


iftVoxel iftFlipVoxelOnYAxis(iftVoxel u, int xsize, int ysize, int zsize) {
    iftVoxel v;
    v.y = ysize - 1 - u.y;
    v.x = u.x; v.z = u.z;
    return v;
}


iftVoxel iftFlipVoxelOnZAxis(iftVoxel u, int xsize, int ysize, int zsize) {
    iftVoxel v;
    v.z = zsize - 1 - u.z;
    v.x = u.x; v.y = u.y;

    return v;
}


iftImage *iftFlipImage(const iftImage *img, iftAxis axis) {
    iftImage *flip_img = iftCreateImageFromImage(img);

    iftFlipVoxelFunc flipVoxelFunc = iftFlipVoxelOnXAxis;
    if (axis == IFT_AXIS_X)
        flipVoxelFunc = iftFlipVoxelOnXAxis;
    else if (axis == IFT_AXIS_Y)
        flipVoxelFunc = iftFlipVoxelOnYAxis;
    else if (axis == IFT_AXIS_Z)
        flipVoxelFunc = iftFlipVoxelOnZAxis;

    iftCopyImageVoxelFunc copyImageVoxelFunc = (iftIsColorImage(img)) ? iftCopyColorImageVoxel : iftCopyGrayImageVoxel;
    
    #pragma omp parallel for
    for (int p = 0; p < img->n; p++) {
        iftVoxel u = iftGetVoxelCoord(img, p);
        iftVoxel v = flipVoxelFunc(u, img->xsize, img->ysize, img->zsize);
        copyImageVoxelFunc(img, u, flip_img, v);
    }

    return flip_img;
}


iftMImage *iftFlipMImage(const iftMImage *mimg, iftAxis axis) {
    iftMImage *flip_mimg = iftCreateMImage(mimg->xsize, mimg->ysize, mimg->zsize, mimg->m);

    iftFlipVoxelFunc flipVoxelFunc = iftFlipVoxelOnXAxis;
    if (axis == IFT_AXIS_X)
        flipVoxelFunc = iftFlipVoxelOnXAxis;
    else if (axis == IFT_AXIS_Y)
        flipVoxelFunc = iftFlipVoxelOnYAxis;
    else if (axis == IFT_AXIS_Z)
        flipVoxelFunc = iftFlipVoxelOnZAxis;

    for (int b = 0; b < mimg->m; b++) {
        #pragma omp parallel for
        for (int p = 0; p < mimg->n; p++) {
            iftVoxel u = iftMGetVoxelCoord(mimg, p);
            iftVoxel v = flipVoxelFunc(u, mimg->xsize, mimg->ysize, mimg->zsize);

            int p = iftMGetVoxelIndex(mimg, u);
            int q = iftMGetVoxelIndex(mimg, v);

            flip_mimg->val[q][b] = mimg->val[p][b];
        }
    }

    return flip_mimg;
}




iftImage *iftTransformImageClipping(iftImage *img, iftMatrix *InvE, int xsize, int ysize, int zsize)
{
   iftImage  *rimg = NULL;
   iftPoint P1,P2;
   iftVoxel v,u;
   int q;

   if (img->zsize > 1) 
   { // 3D
     rimg = iftCreateImage(xsize, ysize, zsize);
 
     for (v.z = 0; v.z < rimg->zsize; v.z++) 
        for (v.y = 0; v.y < rimg->ysize; v.y++) 
 	   for (v.x = 0; v.x < rimg->xsize; v.x++)
           {
 	      q    = iftGetVoxelIndex(rimg,v);
 	      P1.x = v.x;
   	      P1.y = v.y;
	      P1.z = v.z;
	      P2   = iftTransformPoint(InvE,P1);
	      P2.x = P2.x;
	      P2.y = P2.y;
	      P2.z = P2.z;
	      u.x  = iftRound(P2.x);
	      u.y  = iftRound(P2.y);
	      u.z  = iftRound(P2.z);
	      if (iftValidVoxel(img,u)){
	         rimg->val[q] = iftImageValueAtPoint(img,P2);
	   }
      }

     iftCopyVoxelSize(img,rimg);
   } 

   return(rimg);
}

iftFImage *iftTransformFImageClipping(iftFImage *img, iftMatrix *InvE, int xsize, int ysize, int zsize)
{
   iftFImage  *rimg = NULL;
   iftPoint P1,P2;
   iftVoxel v,u;
   int q;

   if (img->zsize > 1) 
   { // 3D
     rimg = iftCreateFImage(xsize, ysize, zsize);
 
     for (v.z = 0; v.z < rimg->zsize; v.z++) 
        for (v.y = 0; v.y < rimg->ysize; v.y++) 
 	   for (v.x = 0; v.x < rimg->xsize; v.x++)
           {
 	      q    = iftFGetVoxelIndex(rimg,v);
 	      P1.x = v.x;
   	      P1.y = v.y;
	      P1.z = v.z;
	      P2   = iftTransformPoint(InvE,P1);
	      P2.x = P2.x;
	      P2.y = P2.y;
	      P2.z = P2.z;
	      u.x  = iftRound(P2.x);
	      u.y  = iftRound(P2.y);
	      u.z  = iftRound(P2.z);
	      if (iftFValidVoxel(img,u)){
	         rimg->val[q] = iftFImageValueAtPoint(img,P2);
	      }
          }

     iftFCopyVoxelSize(img,rimg);
   } 

   return(rimg);
}

iftImage *iftTransformImageClip(iftImage *img, iftMatrix *M)
{
  iftImage  *rimg;
  iftMatrix *InvM, *E, *InvE;
  iftPoint P1,P2;
  iftVoxel v,u;
  int q;

  E    = iftExtendMatrix(M,4,4);
  InvM = iftInvertMatrix(M);
  InvE = iftExtendMatrix(InvM,4,4);

  // Only for 3D

    rimg = iftCreateImage(img->xsize,img->ysize,img->zsize);

    for (v.z = 0; v.z < rimg->zsize; v.z++) 
      for (v.y = 0; v.y < rimg->ysize; v.y++) 
	for (v.x = 0; v.x < rimg->xsize; v.x++){
	  q    = iftGetVoxelIndex(rimg,v);
	  P1.x = v.x;
	  P1.y = v.y;
	  P1.z = v.z;
	  P2   = iftTransformPoint(InvE,P1);
	  u.x  = iftRound(P2.x);
	  u.y  = iftRound(P2.y);
	  u.z  = iftRound(P2.z);
	  if (iftValidVoxel(img,u)){
	    rimg->val[q] = iftImageValueAtPoint(img,P2);
	  }
	}

  P1.x = img->dx;
  P1.y = img->dy;
  P1.z = img->dz;

  P2 = iftTransformPoint(E,P1);

  rimg->dx = P2.x;
  rimg->dy = P2.y;
  rimg->dz = P2.z;

  iftDestroyMatrix(&InvM);
  iftDestroyMatrix(&InvE);
  iftDestroyMatrix(&E);

  return(rimg);
}

