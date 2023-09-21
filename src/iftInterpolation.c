#include "iftInterpolation.h"

#include "ift/core/io/Stream.h"


/*-------------- Private Data Structures and Methods --------------------*/

/* Data structure to compute the fitness value for the best cut plane
	 problem at a given position */

typedef struct ift_BestCutPlaneProblem {
	iftPlane *pl;
	iftImage *weight;
	iftPoint  pos;
	int       xviewsize;
	int       yviewsize;
} iftBestCutPlaneProblem;

/* Create the best cut plane problem data structure */

iftBestCutPlaneProblem *iftCreateBestCutPlaneProblem(iftImage *weight, iftPoint pos, int xviewsize, int yviewsize)
{
	iftBestCutPlaneProblem *prob;

	prob     = (iftBestCutPlaneProblem *)iftAlloc(1,sizeof(iftBestCutPlaneProblem));
	prob->pl = iftCreatePlane();
	prob->xviewsize=xviewsize; 
	prob->yviewsize=yviewsize;
	prob->weight   = weight;

	iftSetPlanePos(prob->pl,pos.x,pos.y,pos.z);
	prob->pos    = pos;
	iftSetPlaneOrient(prob->pl,0.0,0.0,1.0);

	return(prob);
}

/* Destroy the best cut plane problem data structure */

void iftDestroyBestCutPlaneProblem(iftBestCutPlaneProblem **prob)
{
	iftBestCutPlaneProblem *aux=*prob;

	if (aux != NULL) {
		iftFree(aux->pl);
		iftFree(aux);
		*prob=NULL;
	}
}

/* Compute displacements for MSPS optimization */

void iftCutPlaneMSDeltas(iftMSPS *msps, float gamma, float *delta_max)
{
	int i,j;

	for (i=0; i < msps->n; i++) {
		for (j=0; j < msps->m; j++) {
		  msps->delta->val[iftGetMatrixIndex(msps->delta,i,j)] = pow(gamma,(j-msps->m+1))*delta_max[i];
		  printf("%lf ",msps->delta->val[iftGetMatrixIndex(msps->delta,i,j)]);
		}
		printf("\n");
	}
			
}


/* MSPS fitness function for the best cut plane problem */

float iftCutPlaneFitness(void *problem, float *theta)
{
	iftImage *slice;
	int p,i;
	float value=0.0;
	iftBestCutPlaneProblem *prob = (iftBestCutPlaneProblem *) problem;
	iftMatrix *R[3];

	/* Check the limits: constraints of the problem */

	for (i=0; i < 2; i++) {
		if (theta[i]<-180.0) theta[i]=-180.0;
		if (theta[i]>+180.0) theta[i]=+180.0;
	}

	/* Update plane orientation */

	iftSetPlaneOrient(prob->pl,0.0,0.0,1.0);

	R[0] = iftRotationMatrix(IFT_AXIS_X, (float)theta[0]);
	R[1] = iftRotationMatrix(IFT_AXIS_Y, (float)theta[1]);
	R[2] = iftMultMatrices(R[1],R[0]);
	prob->pl->normal = iftTransformVector(R[2],prob->pl->normal);
	for (i=0; i < 3; i++) iftDestroyMatrix(&R[i]);

	/* Compute fitness value */

	slice = iftGetSlice(prob->weight,prob->pl,prob->xviewsize,prob->yviewsize);

	for (p=0; p < slice->n; p++){
		value += slice->val[p];
	}

	iftDestroyImage(&slice);

	return(value);
}

/* Data structure to compute the fitness value for the best object cut
	 plane problem */

typedef struct ift_BestObjectCutPlaneProblem {
	iftPlane *pl;
	iftImage *weight;
	iftPoint  object_center;
	int       size;
} iftBestObjectCutPlaneProblem;

/* Create the best object cut plane problem data structure */

iftBestObjectCutPlaneProblem *iftCreateBestObjectCutPlaneProblem(iftImage *obj, iftVoxel pos, iftImage *weight)
{
	iftBestObjectCutPlaneProblem *prob;

	prob     = (iftBestObjectCutPlaneProblem *)iftAlloc(1,sizeof(iftBestObjectCutPlaneProblem));
	prob->pl = iftCreatePlane();
	prob->size=iftObjectDiagonal(obj);
	iftSetPlanePos(prob->pl,pos.x,pos.y,pos.z);
	iftSetPlaneOrient(prob->pl,0.0,0.0,1.0);
	prob->weight = weight;
	prob->object_center = iftGeometricCenter(obj);

	return(prob);
}

/* Destroy the best plane problem data structure */

void iftDestroyBestObjectCutPlaneProblem(iftBestObjectCutPlaneProblem **prob)
{
	iftBestObjectCutPlaneProblem *aux=*prob;

	if (aux != NULL) {
		iftFree(aux->pl);
		iftFree(aux);
		*prob=NULL;
	}
}


/* Compute displacements for MSPS optimization */

void iftObjectCutPlaneMSDeltas(iftMSPS *msps, float gamma, float *delta_max)
{
	int i,j;

	for (i=0; i < msps->n; i++) {
		for (j=0; j < msps->m; j++) {
		  msps->delta->val[iftGetMatrixIndex(msps->delta,i,j)] = pow(gamma,j-msps->m+1)*delta_max[i];
		  printf("%lf ",msps->delta->val[iftGetMatrixIndex(msps->delta,i,j)]);
		}
		printf("\n");
	}
			
}


/* MSPS fitness function for the best plane problem */

float iftObjectCutPlaneFitness(void *problem, float *theta)
{
	iftImage *weight;
	int p,i;
	float value=0.0;
	iftBestObjectCutPlaneProblem *prob = (iftBestObjectCutPlaneProblem *) problem;
	iftMatrix *R[3];

	
	/* Check the limits: constraints of the problem */

	for (i=0; i < 2; i++) {
		if (theta[i]<-180.0) theta[i]=-180.0;
		if (theta[i]>+180.0) theta[i]=+180.0;
	}

	for (i=2; i < 5; i++) {
		if (theta[i] < -prob->size/2.0) theta[i] = -prob->size/2.0;
		if (theta[i] > +prob->size/2.0) theta[i] = +prob->size/2.0;
	}

	/* Update plane orientation */

	iftSetPlaneOrient(prob->pl,0.0,0.0,1.0);
	R[0] = iftRotationMatrix(IFT_AXIS_X, (float)theta[0]);
	R[1] = iftRotationMatrix(IFT_AXIS_Y, (float)theta[1]);
	R[2] = iftMultMatrices(R[1],R[0]);
	prob->pl->normal = iftTransformVector(R[2],prob->pl->normal);
	for (i=0; i < 3; i++) iftDestroyMatrix(&R[i]);

	/* Update plane position */

	prob->pl->pos.x = (float)theta[2]+prob->object_center.x;
	prob->pl->pos.y = (float)theta[3]+prob->object_center.y;
	prob->pl->pos.z = (float)theta[4]+prob->object_center.z;

	/* Compute fitness value */

	weight = iftGetSlice(prob->weight,prob->pl,prob->size,prob->size);
	for (p=0; p < weight->n; p++) 
		value += weight->val[p];

	iftDestroyImage(&weight);

	return(value);
}


/*------------------- Public methods -----------------------------------*/

iftImage *iftResliceImage(iftImage *img, iftPlane *pl, int xviewsize, int yviewsize, int zviewsize)
{
	iftImage *rimg,*aux;
	iftPoint  P1,P0;
	iftVoxel  u;
	int       i,xyviewsize,offset,p,*src,*dst;

	/* Initialize output scene */

	rimg   = iftCreateImage(xviewsize,yviewsize,zviewsize);
	iftCopyVoxelSize(img,rimg);

	/* Starting at zviewsize/2 slices before, reslice the scene until zviewsize/2 slices
		 after the center of the plane and along the direction of v. */

	xyviewsize = xviewsize*yviewsize;

	P0.x   = pl->pos.x - (zviewsize/2.0)*pl->normal.x;
	P0.y   = pl->pos.y - (zviewsize/2.0)*pl->normal.y;
	P0.z   = pl->pos.z - (zviewsize/2.0)*pl->normal.z; 

	for (i=0,offset=0; i < zviewsize; i++,offset+=xyviewsize) {
		P1.x = P0.x + i*pl->normal.x; 
		P1.y = P0.y + i*pl->normal.y;
		P1.z = P0.z + i*pl->normal.z;
		u.x  = iftRound(P1.x);
		u.y  = iftRound(P1.y);
		u.z  = iftRound(P1.z);
		if (iftValidVoxel(img,u)){
			iftSetPlanePos(pl,P1.x,P1.y,P1.z);
			aux  = iftGetSlice(img,pl,xviewsize,yviewsize);
			src  = aux->val;
			dst  = rimg->val+offset;
			for (p=0; p < xyviewsize; p++) 
	dst[p] = src[p];
			iftDestroyImage(&aux);
		}
	}
	

	return(rimg);
}

/* The idea is to align the visualization plane with the cut plane and
   get the value of each pixel in the visualization plane by
   interpolation of the voxel values in the neighborhood of the mapped
   coordinate of that pixel in the cut plane */

iftImage *iftGetSlice(iftImage *img, iftPlane *pl, int usize, int vsize) 
{ 
	iftImage *slice;
	int p;
	iftVoxel u,v;
	iftPoint P1,P2;
	iftMatrix *A;
	
	slice = iftCreateImage(usize,vsize,1);
	iftCopyVoxelSize(img,slice);
	A     = iftRotationMatrixToAlignZWithVector(pl->normal);

	u.z = 0;
	for (u.y=0; u.y < vsize; u.y++)
	  for (u.x=0; u.x < usize; u.x++){

	    /* Consider the visualization plane initially placed in
	       parallel to the XY plane with z=0 and origin at the
	       center of the xsize x ysize region. This region must be
	       translated to the origin of the coordinate system,
	       rotated to align the z axis with the normal vector of
	       the cut plane, and then translated to the position of
	       the cut plane in order to compute the pixel intensity
	       of the mapped pixels by interpolation of the
	       intensities of their neighbors in the scene.
	    */
	    
	    p = iftGetVoxelIndex(slice,u);
	    P1.x = u.x - usize/2.0;
	    P1.y = u.y - vsize/2.0;
	    P1.z = u.z;

	    /* Rotation of the region such that the Z axis is aligned
	       to the normal of the cut plane */

	    P2 = iftTransformPoint(A,P1);

	    /* Translation of P2 to the position of the cut plane */

	    P2.x = P2.x + pl->pos.x;
	    P2.y = P2.y + pl->pos.y;
	    P2.z = P2.z + pl->pos.z;

	    /* Get its voxel intensity in the scene by trilinear
	       interpolation */
			
	    v.x = iftRound(P2.x);
	    v.y = iftRound(P2.y);
	    v.z = iftRound(P2.z);
			
	    if (iftValidVoxel(img,v)){
	      slice->val[p] = iftImageValueAtPoint(img,P2);
	    }
	    
	  }

	iftDestroyMatrix(&A);
	
	return(slice);
}



iftImage *iftInterp2D(const iftImage *img, float sx, float sy)
{
  iftImage  *ximg, *yimg;  
  int        xsize,ysize;
  
  if (img->zsize != 1)
      iftError("Image must be 2D", "iftInterp2D");
  if ((sx <= 0.0)||(sy <= 0.0))
      iftError("Invalid scale factors", "iftInterp2D");

  xsize= iftRound(fabs(sx * img->xsize));
  ysize= iftRound(fabs(sy * img->ysize));

  if (img->Cb != NULL) {
			

    if (sx != 1.0) {

      /* Interpolate along x */
      
      ximg = iftCreateImage(xsize,img->ysize,1);
      ximg->dx = img->dx/sx;
      ximg->dy = img->dy;
      ximg->dz = img->dz;      
      iftSetCbCr(ximg,(iftMaxImageRange(iftImageDepth(img))+1)/2);
      
#pragma omp parallel for shared(ximg, img, sx)
      for (int y = 0;  y < ximg->ysize; y++){
	iftVoxel u,v,w;
	u.z=v.z=w.z=0; 
	u.y = w.y = v.y = y;
	for (v.x = 0; v.x < ximg->xsize; v.x++){  
	  int q    = iftGetVoxelIndex(ximg,v); 
	  u.x  = (int)(v.x/sx);
	  float dx   = (v.x/sx) - u.x; 
	  w.x  = ((u.x+1)==img->xsize)?u.x:u.x+1;
	  int p    = iftGetVoxelIndex(img,u); 
	  int r    = iftGetVoxelIndex(img,w); 
	  ximg->val[q] = (int)(img->val[p]*(1.0-dx)+img->val[r]*dx);
	  ximg->Cb[q]  = (int)(img->Cb[p]*(1.0-dx)+img->Cb[r]*dx);
	  ximg->Cr[q]  = (int)(img->Cr[p]*(1.0-dx)+img->Cr[r]*dx);
	}
      }
    } else {
      ximg = iftCopyImage(img);
    }

    if ( sy != 1.0) {
 
      /* Interpolate along y */
      
      yimg = iftCreateImage(xsize,ysize,1);
      iftSetCbCr(yimg,(iftMaxImageRange(iftImageDepth(img))+1)/2);
      yimg->dx = ximg->dx;
      yimg->dy = ximg->dy/sy;
      yimg->dz = ximg->dz;
      
#pragma omp parallel for shared(yimg, ximg, sy)
      for (int x = 0; x < yimg->xsize; x++){
	iftVoxel u, v, w;
	u.z=v.z=w.z=0; 
	u.x = w.x = v.x = x;
	for (v.y = 0; v.y < yimg->ysize; v.y++){  
	  int q    = iftGetVoxelIndex(yimg,v); 
	  u.y  = (int)(v.y/sy);
	  float dy   = (v.y/sy) - u.y; 
	  w.y  = ((u.y+1)==ximg->ysize)?u.y:u.y+1;
	  int p    = iftGetVoxelIndex(ximg,u); 
	  int r    = iftGetVoxelIndex(ximg,w); 
	  yimg->val[q] = (int)(ximg->val[p]*(1.0-dy)+ximg->val[r]*dy);
	  yimg->Cb[q]  = (int)(ximg->Cb[p]*(1.0-dy)+ximg->Cb[r]*dy);
	  yimg->Cr[q]  = (int)(ximg->Cr[p]*(1.0-dy)+ximg->Cr[r]*dy);
	}
      }
    } else {
      yimg = iftCopyImage(ximg);
    }
    
  }else{ /* Gray-Scale Image */


    if (sx != 1.0) {

      /* Interpolate along x */
      
      ximg = iftCreateImage(xsize,img->ysize,1);
      ximg->dx = img->dx/sx;
      ximg->dy = img->dy;
      ximg->dz = img->dz;      
      
#pragma omp parallel for shared(ximg, img, sx)
      for (int y = 0; y < ximg->ysize; y++){
	iftVoxel u, v, w;
	u.z=v.z=w.z=0; 
 	u.y = w.y = v.y = y;
	for (v.x = 0; v.x < ximg->xsize; v.x++){  
	  int q    = iftGetVoxelIndex(ximg,v); 
	  u.x  = (int)(v.x/sx);
	  float dx   = (v.x/sx) - u.x; 
	  w.x  = ((u.x+1)==img->xsize)?u.x:u.x+1;
	  int p    = iftGetVoxelIndex(img,u); 
	  int r    = iftGetVoxelIndex(img,w); 
	  ximg->val[q] = (int)(img->val[p]*(1.0-dx)+img->val[r]*dx);
	}
      }
    } else {
      ximg = iftCopyImage(img);
    }
    
    if (sy != 1.0) {
			
      /* Interpolate along y */
			
      yimg = iftCreateImage(xsize,ysize,1);
      yimg->dx = ximg->dx;
      yimg->dy = ximg->dy/sy;
      yimg->dz = ximg->dz;

#pragma omp parallel for shared(yimg, ximg, sy)
      for (int x = 0; x < yimg->xsize; x++){
	iftVoxel u,v,w;
	u.z=v.z=w.z=0; 
	u.x = w.x = v.x = x;
	for (v.y = 0; v.y < yimg->ysize; v.y++){  
	  int q    = iftGetVoxelIndex(yimg,v); 
	  u.y  = (int)(v.y/sy);
	  float dy   = (v.y/sy) - u.y; 
	  w.y  = ((u.y+1)==ximg->ysize)?u.y:u.y+1;
	  int p    = iftGetVoxelIndex(ximg,u); 
	  int r    = iftGetVoxelIndex(ximg,w); 
	  yimg->val[q] = (int)(ximg->val[p]*(1.0-dy)+ximg->val[r]*dy);
	}
      }
    } else { 
      yimg = iftCopyImage(ximg);
    }  
  }

  iftDestroyImage(&ximg);

  return(yimg);
}


iftImage *iftInterpByNearestNeighbor(const iftImage *img, float sx, float sy, float sz) {
    if ((sx <= 0.0) || (sy <= 0.0) || (sz <= 0.0))
        iftError("Invalid scale factors: sx = %f, sy = %f, sz = %f\nTry > 0.0",
                 "iftInterpByNearestNeighbor", sx, sy, sz);

    iftImage *interp = NULL;

    if (sx == 1.0 && sy == 1.0 && sz == 1.0)
        interp = iftCopyImage(img);
    else {
        int xsize = iftRound(sx * img->xsize);
        int ysize = iftRound(sy * img->ysize);
        int zsize = iftRound(sz * img->zsize);

        iftCopyImageVoxelFunc copyImageVoxelFunc;
        if (iftIsColorImage(img)) {
            interp = iftCreateColorImage(xsize, ysize, zsize, iftImageDepth(img));
            copyImageVoxelFunc = iftCopyColorImageVoxel;
        }
        else {
            interp = iftCreateImage(xsize, ysize, zsize);
            copyImageVoxelFunc = iftCopyGrayImageVoxel;
        }

        interp->dx = img->dx / sx;
        interp->dy = img->dy / sy;
        interp->dz = img->dz / sz;

        #pragma omp parallel for
        for (int p = 0; p < interp->n; p++) {
            iftVoxel v = iftGetVoxelCoord(interp, p);
            iftVoxel u;
            u.x = iftMin(iftRound(v.x / sx), img->xsize-1);
            u.y = iftMin(iftRound(v.y / sy), img->ysize-1);
            u.z = iftMin(iftRound(v.z / sz), img->zsize-1);

            copyImageVoxelFunc(img, u, interp, v);
        }
    }


    return interp;
}

iftImage *iftInterp2DByNearestNeighbor(const iftImage *img, float sx, float sy) {
	if ((sx <= 0.0) || (sy <= 0.0))
		iftError("Invalid scale factors: sx = %f, sy = %f, sz = %f\nTry > 0.0",
				 "iftInterpByNearestNeighbor", sx, sy);

	iftImage *interp = NULL;

	if (sx == 1.0 && sy == 1.0)
		interp = iftCopyImage(img);
	else {
		int xsize = iftRound(sx * img->xsize);
		int ysize = iftRound(sy * img->ysize);

		iftCopyImageVoxelFunc copyImageVoxelFunc;
		if (iftIsColorImage(img)) {
			interp = iftCreateColorImage(xsize, ysize, 1, iftImageDepth(img));
			copyImageVoxelFunc = iftCopyColorImageVoxel;
		}
		else {
			interp = iftCreateImage(xsize, ysize, 1);
			copyImageVoxelFunc = iftCopyGrayImageVoxel;
		}

		interp->dx = img->dx / sx;
		interp->dy = img->dy / sy;

#pragma omp parallel for
		for (int p = 0; p < interp->n; p++) {
			iftVoxel v = iftGetVoxelCoord(interp, p);
			iftVoxel u;
			u.x = iftMin(iftRound(v.x / sx), img->xsize-1);
			u.y = iftMin(iftRound(v.y / sy), img->ysize-1);

			copyImageVoxelFunc(img, u, interp, v);
		}
	}


	return interp;
}

iftMImage * iftMInterp2D(iftMImage *mimg, float sx, float sy)
{
  iftMImage  *ximg, *yimg;  
  int        xsize,ysize;
	
  if (mimg->zsize != 1)
      iftError("Image must be 2D", "iftMInterp2D");
  if ((sx <= 0.0)||(sy <= 0.0))
      iftError("Invalid scale factors", "iftMInterp2D");
  
  xsize= iftRound(fabs(sx * mimg->xsize));
  ysize= iftRound(fabs(sy * mimg->ysize));
  
  if (sx != 1.0) {
    
    /* Interpolate along x */

    ximg = iftCreateMImage(xsize, mimg->ysize, 1, mimg->m);
    ximg->dx = mimg->dx/sx;
    ximg->dy = mimg->dy;
    ximg->dz = mimg->dz;      

#pragma omp parallel for shared(ximg, mimg, sx)    
    for (int b = 0; b < mimg->m ; b++){
      iftVoxel u, v, w;
      u.z=v.z=w.z=0; 
      for (v.y = 0; v.y < ximg->ysize; v.y++){
	u.y = w.y = v.y;
	for (v.x = 0; v.x < ximg->xsize; v.x++)
	  {  
	    int q    = iftMGetVoxelIndex(ximg,v); 
	    u.x  = (int)(v.x/sx);
	    float dx   = (v.x/sx) - u.x; 
	    w.x  = ((u.x+1)==mimg->xsize)?u.x:u.x+1;
	    int p    = iftMGetVoxelIndex(mimg,u); 
	    int r    = iftMGetVoxelIndex(mimg,w); 
	    ximg->val[q][b] = mimg->val[p][b]*(1.0-dx) + mimg->val[r][b]*dx;
	  }
      }
    }
  } else {
    ximg = iftCopyMImage(mimg);
  }
		
  if (sy != 1.0) {
		
    /* Interpolate along y */
		
    yimg = iftCreateMImage(xsize, ysize, 1, mimg->m);
    yimg->dx = ximg->dx;
    yimg->dy = ximg->dy/sy;
    yimg->dz = ximg->dz;

#pragma omp parallel for shared(yimg, ximg, sy)    
    for (int b = 0; b < mimg->m ; b++){
      iftVoxel u,v,w;
      u.z=v.z=w.z=0; 
      for (v.x = 0; v.x < yimg->xsize; v.x++){
	u.x = w.x = v.x;
	for (v.y = 0; v.y < yimg->ysize; v.y++){  
	  int q    = iftMGetVoxelIndex(yimg,v); 
	  u.y  = (int)(v.y/sy);
	  float dy   = (v.y/sy) - u.y; 
	  w.y  = ((u.y+1)==ximg->ysize)?u.y:u.y+1;
	  int p    = iftMGetVoxelIndex(ximg,u); 
	  int r    = iftMGetVoxelIndex(ximg,w); 
	  yimg->val[q][b] = ximg->val[p][b]*(1.0-dy) + ximg->val[r][b]*dy;
	}
      }
    }
  } else { 
    yimg = iftCopyMImage(ximg);
  }

  iftDestroyMImage(&ximg);

  return(yimg);
}

iftImage *iftInterp(const iftImage *img, float sx, float sy, float sz)
{
  iftImage  *ximg, *yimg, *zimg;  
  int        xsize,ysize,zsize;
  
  if (img->zsize == 1)
      iftError("Use iftInterp2D for 2D images", "iftInterp");
  if ((sx <= 0.0)||(sy <= 0.0)||(sz <= 0.0))
      iftError("Invalid scale factors", "iftInterp");
  
  xsize= iftRound(fabs(sx * img->xsize));
  ysize= iftRound(fabs(sy * img->ysize));
  zsize= iftRound(fabs(sz * img->zsize));
  
  if (sx != 1.0 ) {
    
    /* Interpolate along x */
	  
    ximg = iftCreateImage(xsize,img->ysize, img->zsize);
    ximg->dx = img->dx/sx;
    ximg->dy = img->dy;
    ximg->dz = img->dz;
    
    
#pragma omp parallel for shared(ximg, img, sx)
    for (int z = 0; z < ximg->zsize; z++){		  
      iftVoxel u, v, w;
      u.z = w.z = v.z = z;
      for (v.y = 0; v.y < ximg->ysize; v.y++){
	u.y = w.y = v.y;
	for (v.x = 0; v.x < ximg->xsize; v.x++){  
	  int q    = iftGetVoxelIndex(ximg,v); 
	  u.x  = (int)(v.x/sx);
	  float dx   = (v.x/sx) - u.x; 
	  w.x  = ((u.x+1)==img->xsize)?u.x:u.x+1;
	  int p    = iftGetVoxelIndex(img,u); 
	  int r    = iftGetVoxelIndex(img,w); 
	  ximg->val[q] = (int)(img->val[p]*(1.0-dx)+img->val[r]*dx);
	}
      }
    }
  }else{
    ximg = iftCopyImage(img);
  }
  
  if (sy != 1.0) { 
    
    /* Interpolate along y */
    
    yimg = iftCreateImage(xsize,ysize, img->zsize);
    yimg->dx = ximg->dx;
    yimg->dy = ximg->dy/sy;
    yimg->dz = ximg->dz;
    
#pragma omp parallel for shared(yimg, ximg, sy)
    for (int z = 0; z < yimg->zsize; z++){
      iftVoxel u,w,v;
      u.z = w.z = v.z = z;
      for (v.x = 0; v.x < yimg->xsize; v.x++){
	u.x = w.x = v.x;
	for (v.y = 0; v.y < yimg->ysize; v.y++){  
	  int q    = iftGetVoxelIndex(yimg,v); 
	  u.y  = (int)(v.y/sy);
	  float dy   = (v.y/sy) - u.y; 
	  w.y  = ((u.y+1)==ximg->ysize)?u.y:u.y+1;
	  int p    = iftGetVoxelIndex(ximg,u); 
	  int r    = iftGetVoxelIndex(ximg,w); 
	  yimg->val[q] = (int)(ximg->val[p]*(1.0-dy)+ximg->val[r]*dy);
	}
      }
    }
  } else {
    yimg = iftCopyImage(ximg);
  }
  
  iftDestroyImage(&ximg);
  
  if (sz != 1.0) {
    
    /* Interpolate along z */
    
    zimg = iftCreateImage(xsize,ysize,zsize);
    zimg->dx = yimg->dx;
    zimg->dy = yimg->dy;
    zimg->dz = yimg->dz/sz;
    
#pragma omp parallel for shared(zimg, yimg, sz)
    for (int y = 0; y < zimg->ysize; y++){
      iftVoxel u,v,w;
      u.y = w.y = v.y = y;
      for (v.x = 0; v.x < zimg->xsize; v.x++){
	u.x = w.x = v.x;
	for (v.z = 0; v.z < zimg->zsize; v.z++){  
	  int q    = iftGetVoxelIndex(zimg,v); 
	  u.z  = (int)(v.z/sz);
	  float dz   = (v.z/sz) - u.z; 
	  w.z  = ((u.z+1)==yimg->zsize)?u.z:u.z+1;
	  int p    = iftGetVoxelIndex(yimg,u); 
	  int r    = iftGetVoxelIndex(yimg,w); 
	  zimg->val[q] = (int)(yimg->val[p]*(1.0-dz)+yimg->val[r]*dz);
	}
      }
    }
    
  } else { 
    zimg = iftCopyImage(yimg);
  }
  
  iftDestroyImage(&yimg);
  
  return(zimg);
}

iftImage* iftResizeImage(const iftImage *img, int xsize, int ysize, int zsize) {
	if(iftIs3DImage(img)) {
		return iftInterp(img, ((float) xsize) / img->xsize, ((float) ysize) / img->ysize, ((float) zsize) / img->zsize);
	}
	else {
		return iftInterp2D(img, ((float) xsize) / img->xsize, ((float) ysize) / img->ysize);
	}
}

iftMImage *iftMInterp(iftMImage *mimg, float sx, float sy, float sz)
{
  iftMImage  *ximg, *yimg, *zimg;  
  int        xsize,ysize,zsize;
  
  if (mimg->zsize == 1)
      iftError("Use iftInterp2D for 2D images", "iftInterp");
  if ((sx <= 0.0)||(sy <= 0.0)||(sz <= 0.0))
      iftError("Invalid scale factors", "iftInterp");
  
  xsize= iftRound(fabs(sx * mimg->xsize));
  ysize= iftRound(fabs(sy * mimg->ysize));
  zsize= iftRound(fabs(sz * mimg->zsize));
  
  if (sx != 1.0 ) {
    
    /* Interpolate along x */
		
    ximg = iftCreateMImage(xsize, mimg->ysize, mimg->zsize, mimg->m);
    ximg->dx = mimg->dx/sx;
    ximg->dy = mimg->dy;
    ximg->dz = mimg->dz;
    
#pragma omp parallel for shared(ximg, mimg, sx)
    for (int b = 0; b < mimg->m ; b++){
      iftVoxel u,v,w;
      for (v.z = 0; v.z < ximg->zsize; v.z++){
	u.z = w.z = v.z;
	for (v.y = 0; v.y < ximg->ysize; v.y++){
	  u.y = w.y = v.y;
	  for (v.x = 0; v.x < ximg->xsize; v.x++){  
	    int q    = iftMGetVoxelIndex(ximg,v); 
	    u.x  = (int)(v.x/sx);
	    float dx   = (v.x/sx) - u.x; 
	    w.x  = ((u.x+1)==mimg->xsize)?u.x:u.x+1;
	    int p    = iftMGetVoxelIndex(mimg,u); 
	    int r    = iftMGetVoxelIndex(mimg,w); 
	    ximg->val[q][b] = mimg->val[p][b]*(1.0-dx) + mimg->val[r][b]*dx;
	  }
	}
      }
    }
  }else{
    ximg = iftCopyMImage(mimg);
  }
  
  if (sy != 1.0) { 

    /* Interpolate along y */

    yimg = iftCreateMImage(xsize, ysize, mimg->zsize, mimg->m);
    yimg->dx = ximg->dx;
    yimg->dy = ximg->dy/sy;
    yimg->dz = ximg->dz;
    
#pragma omp parallel for shared(yimg, ximg, sy)
    for (int b = 0; b < mimg->m ; b++){
      iftVoxel u,v,w;
      for (v.z = 0; v.z < yimg->zsize; v.z++){
	u.z = w.z = v.z;
	for (v.x = 0; v.x < yimg->xsize; v.x++){
	  u.x = w.x = v.x;
	  for (v.y = 0; v.y < yimg->ysize; v.y++){  
	    int q    = iftMGetVoxelIndex(yimg,v); 
	    u.y  = (int)(v.y/sy);
	    float dy   = (v.y/sy) - u.y; 
	    w.y  = ((u.y+1)==ximg->ysize)?u.y:u.y+1;
	    int p    = iftMGetVoxelIndex(ximg,u); 
	    int r    = iftMGetVoxelIndex(ximg,w); 
	    yimg->val[q][b] = ximg->val[p][b]*(1.0-dy) + ximg->val[r][b]*dy;
	  }
	}
      }
    }
  } else {
    yimg = iftCopyMImage(ximg);
  }
  
  iftDestroyMImage(&ximg);
  
  if (sz != 1.0) {
    
    /* Interpolate along z */
    
    zimg = iftCreateMImage(xsize,ysize,zsize, mimg->m);
    zimg->dx = yimg->dx;
    zimg->dy = yimg->dy;
    zimg->dz = yimg->dz/sz;
    
#pragma omp parallel for shared(zimg, yimg, sz)
    for (int b = 0; b < mimg->m ; b++){
      iftVoxel u,v,w;
      for (v.y = 0; v.y < zimg->ysize; v.y++){
	u.y = w.y = v.y;
	for (v.x = 0; v.x < zimg->xsize; v.x++){
	  u.x = w.x = v.x;
	  for (v.z = 0; v.z < zimg->zsize; v.z++){  
	    int q    = iftGetVoxelIndex(zimg,v); 
	    u.z  = (int)(v.z/sz);
	    float dz   = (v.z/sz) - u.z; 
	    w.z  = ((u.z+1)==yimg->zsize)?u.z:u.z+1;
	    int p    = iftGetVoxelIndex(yimg,u); 
	    int r    = iftGetVoxelIndex(yimg,w); 
	    zimg->val[q][b] = yimg->val[p][b]*(1.0-dz) + yimg->val[r][b]*dz;
	  }
	}
      }
    }
  } else { 
    zimg = iftCopyMImage(yimg);
  }
  
  iftDestroyMImage(&yimg);

  return(zimg);
}


iftMImage* iftResizeMImage(iftMImage* img, int xsize, int ysize, int zsize) {
	if(iftIs3DMImage(img)) {
		return iftMInterp(img, ((float) xsize) / img->xsize, ((float) ysize) / img->ysize, ((float) zsize) / img->zsize);
	}
	else {
		return iftMInterp2D(img, ((float) xsize) / img->xsize, ((float) ysize) / img->ysize);
	}
}


iftFImage *iftFGetSlice(iftFImage *img, iftPlane *pl, int xsize, int ysize) 
{ 
	iftFImage *slice;
	int p;
	iftVoxel u,v;
	iftPoint P1,P2;
	iftMatrix *A;

	slice = iftCreateFImage(xsize,ysize,1);
	iftFCopyVoxelSize(img,slice);
	A     = iftRotationMatrixToAlignZWithVector(pl->normal);

	u.z = 0;
	for (u.y=0; u.y < ysize; u.y++)
		for (u.x=0; u.x < xsize; u.x++){

		  /* Consider the visualization plane initially placed in parallel
	 to the XY plane with z=0 and origin at the center of the
	 xsize x ysize region. This region must be first centered at
	 the origin of the XYZ coordinate system. After that, compute
	 P1 as the coordinates of pixel u in that region. */

			p = iftGetVoxelIndex(slice,u);
			P1.x = u.x - xsize/2.0;
			P1.y = u.y - ysize/2.0;
			P1.z = u.z;


			/* Rotate the region such that the Z axis (normal of the plane)
			   is aligned to the normal of the desired slicing plane */

			P2 = iftTransformPoint(A,P1);

			/* Translate P2 onto the position of the desired slicing plane */

			P2.x = P2.x + pl->pos.x;
			P2.y = P2.y + pl->pos.y;
			P2.z = P2.z + pl->pos.z;

			/* Get its voxel intensity in the scene by trilinear
			   interpolation */
			
			v.x = iftRound(P2.x);
			v.y = iftRound(P2.y);
			v.z = iftRound(P2.z);
			
			if (iftFValidVoxel(img,v)){
			  slice->val[p] = iftFImageValueAtPoint(img,P2);
			}
			
		}

	iftDestroyMatrix(&A);
	
	return(slice);
}

iftFImage *iftFInterp2D(iftFImage *img, float sx, float sy)
{
  iftFImage  *ximg, *yimg;  
  int        xsize,ysize;
  
  if (img->zsize != 1)
      iftError("Image must be 2D", "iftFInterp2D");
  if ((sx <= 0.0)||(sy <= 0.0))
      iftError("Invalid scale factors", "iftFInterp2D");
  
  xsize= iftRound(fabs(sx * img->xsize));
  ysize= iftRound(fabs(sy * img->ysize));
  
  if (sx != 1.0) {

    /* Interpolate along x */
    
    ximg = iftCreateFImage(xsize,img->ysize,1);
    ximg->dx = img->dx/sx;
    ximg->dy = img->dy;
    ximg->dz = img->dz;      
    
#pragma omp parallel for shared(ximg, img, sx)
    for (int y = 0; y < ximg->ysize; y++){
      iftVoxel u,v,w;
      u.z=v.z=w.z=0; 
      u.y = w.y = v.y = y;
      for (v.x = 0; v.x < ximg->xsize; v.x++){  
	int q    = iftFGetVoxelIndex(ximg,v); 
	u.x  = (int)(v.x/sx);
	float dx   = (v.x/sx) - u.x; 
	w.x  = ((u.x+1)==img->xsize)?u.x:u.x+1;
	int p    = iftFGetVoxelIndex(img,u); 
	int r    = iftFGetVoxelIndex(img,w); 
	ximg->val[q] = (img->val[p]*(1.0-dx)+img->val[r]*dx);
      }
    }
  } else {
    ximg = iftFCopyImage(img);
  }
		
  if (sy != 1.0) {
		
    /* Interpolate along y */
    
    yimg = iftCreateFImage(xsize,ysize,1);
    yimg->dx = ximg->dx;
    yimg->dy = ximg->dy/sy;
    yimg->dz = ximg->dz;

#pragma omp parallel for shared(yimg, ximg, sy)
    for (int x = 0; x < yimg->xsize; x++){
      iftVoxel u,v,w;
      u.z=v.z=w.z=0; 
      u.x = w.x = v.x = x;
      for (v.y = 0; v.y < yimg->ysize; v.y++){  
	int q    = iftFGetVoxelIndex(yimg,v); 
	u.y  = (int)(v.y/sy);
	float dy   = (v.y/sy) - u.y; 
	w.y  = ((u.y+1)==ximg->ysize)?u.y:u.y+1;
	int p    = iftFGetVoxelIndex(ximg,u); 
	int r    = iftFGetVoxelIndex(ximg,w); 
	yimg->val[q] = (ximg->val[p]*(1.0-dy)+ximg->val[r]*dy);
      }
    }
  } else { 
    yimg = iftFCopyImage(ximg);
  }

  iftDestroyFImage(&ximg);

  return(yimg);
}

iftFImage *iftFInterp(iftFImage *img, float sx, float sy, float sz)
{
  iftFImage  *ximg, *yimg, *zimg;  
  int        xsize,ysize,zsize;
	
  if (img->zsize == 1)
      iftError("Use iftFInterp2D for 2D images", "iftFInterp");
  if ((sx <= 0.0)||(sy <= 0.0)||(sz <= 0.0))
      iftError("Invalid scale factors", "iftFInterp");

  xsize= iftRound(fabs(sx * img->xsize));
  ysize= iftRound(fabs(sy * img->ysize));
  zsize= iftRound(fabs(sz * img->zsize));

  if (sx != 1.0 ) {
		
    /* Interpolate along x */
		
    ximg = iftCreateFImage(xsize,img->ysize, img->zsize);
    ximg->dx = img->dx/sx;
    ximg->dy = img->dy;
    ximg->dz = img->dz;

#pragma omp parallel for shared(ximg, img, sx)
    for (int z = 0; z < ximg->zsize; z++){
      iftVoxel u,w,v;
      u.z = w.z = v.z = z;
      for (v.y = 0; v.y < ximg->ysize; v.y++){
	u.y = w.y = v.y;
	for (v.x = 0; v.x < ximg->xsize; v.x++){  
	  int q    = iftFGetVoxelIndex(ximg,v); 
	  u.x  = (int)(v.x/sx);
	  float dx   = (v.x/sx) - u.x; 
	  w.x  = ((u.x+1)==img->xsize)?u.x:u.x+1;
	  int p    = iftFGetVoxelIndex(img,u); 
	  int r    = iftFGetVoxelIndex(img,w); 
	  ximg->val[q] = (img->val[p]*(1.0-dx)+img->val[r]*dx);
	}
      }
    }
  }else{
    ximg = iftFCopyImage(img);
  }

  if (sy != 1.0) { 

    /* Interpolate along y */

    yimg = iftCreateFImage(xsize,ysize, img->zsize);
    yimg->dx = ximg->dx;
    yimg->dy = ximg->dy/sy;
    yimg->dz = ximg->dz;
    
#pragma omp parallel for shared(yimg, ximg, sy)
    for (int z = 0; z < yimg->zsize; z++){
      iftVoxel u,w,v;
      u.z = w.z = v.z = z;
      for (v.x = 0; v.x < yimg->xsize; v.x++){
	u.x = w.x = v.x;
	for (v.y = 0; v.y < yimg->ysize; v.y++){  
	  int q    = iftFGetVoxelIndex(yimg,v); 
	  u.y  = (int)(v.y/sy);
	  float dy   = (v.y/sy) - u.y; 
	  w.y  = ((u.y+1)==ximg->ysize)?u.y:u.y+1;
	  int p    = iftFGetVoxelIndex(ximg,u); 
	  int r    = iftFGetVoxelIndex(ximg,w); 
	  yimg->val[q] = (ximg->val[p]*(1.0-dy)+ximg->val[r]*dy);
	}
      }
    }
  } else {
    yimg = iftFCopyImage(ximg);
  }

  iftDestroyFImage(&ximg);
  
  if (sz != 1.0) {
    
    /* Interpolate along z */
    
    zimg = iftCreateFImage(xsize,ysize,zsize);
    zimg->dx = yimg->dx;
    zimg->dy = yimg->dy;
    zimg->dz = yimg->dz/sz;
    
#pragma omp parallel for shared(zimg, yimg, sz)
    for (int y = 0; y < zimg->ysize; y++){
      iftVoxel u,v,w;
      u.y = w.y = v.y = y;
      for (v.x = 0; v.x < zimg->xsize; v.x++){
	u.x = w.x = v.x;
	for (v.z = 0; v.z < zimg->zsize; v.z++){  
	  int q    = iftFGetVoxelIndex(zimg,v); 
	  u.z  = (int)(v.z/sz);
	  float dz   = (v.z/sz) - u.z; 
	  w.z  = ((u.z+1)==yimg->zsize)?u.z:u.z+1;
	  int p    = iftFGetVoxelIndex(yimg,u); 
	  int r    = iftFGetVoxelIndex(yimg,w); 
	  zimg->val[q] = (yimg->val[p]*(1.0-dz)+yimg->val[r]*dz);
	}
      }
    }		
  } else { 
    zimg = iftFCopyImage(yimg);
  }
  
  iftDestroyFImage(&yimg);

  return(zimg);
}


iftImage *iftShapeBasedInterp(iftImage *label, float sx, float sy, float sz) {
    char      *obj_code     = iftAllocCharArray(iftMaximumValue(label) + 1);
    int       l, p, xsize, ysize, zsize;
    iftVoxel  pos;
    iftImage  *bin, *nlabel = NULL;
    iftFImage *sdist, *inter;
    iftAdjRel *A            = iftSpheric(sqrtf(3.0));

    if (label->zsize == 1)
        iftError("Use iftShapeBasedInterp2D for 2D images", "iftShapeBasedInterp");

    if ((sx <= 0.0) || (sy <= 0.0) || (sz <= 0.0))
        iftError("Invalid scale factors", "iftShapeBasedInterp");

    /* Mark existing labels */

    for (p = 0; p < label->n; p++)
        obj_code[label->val[p]] = 1;

    /* Compute shape-based interpolation for each existing label */

    xsize = iftRound(fabs(sx * label->xsize));
    ysize = iftRound(fabs(sy * label->ysize));
    zsize = iftRound(fabs(sz * label->zsize));

    nlabel = iftCreateImage(xsize, ysize, zsize);
    nlabel->dx = label->dx / sx;
    nlabel->dy = label->dy / sy;
    nlabel->dz = label->dz / sz;

    int label_max_val = iftMaximumValue(label);
    for (l = 1; l <= label_max_val; l++) {
        if (obj_code[l] == 1) {
            iftBoundingBox mbb = iftMinObjectBoundingBox(label, l, NULL);
            pos = mbb.begin;
            bin = iftExtractROI(label, mbb);
            sdist = iftSignedDistTrans(bin, A);
            iftDestroyImage(&bin);
            inter = iftFInterp(sdist, sx, sy, sz);
            iftDestroyFImage(&sdist);
            bin = iftFThreshold(inter, 0, IFT_INFINITY_FLT, 1);
            iftDestroyFImage(&inter);
            pos.x = iftRound(pos.x * sx);
            pos.y = iftRound(pos.y * sy);
            pos.z = iftRound(pos.z * sz);
            iftInsertObject(bin, nlabel, l, pos);
            iftDestroyImage(&bin);
        }
    }

    iftFree(obj_code);
    iftDestroyAdjRel(&A);

    return (nlabel);
}

iftImage *iftShapeBasedInterp2D(iftImage *label, float sx, float sy)
{
  char      *obj_code=iftAllocCharArray(iftMaximumValue(label)+1);
  int        l,p,xsize,ysize;
  iftVoxel   pos;
  iftImage  *bin, *nlabel=NULL;
  iftFImage *sdist,*inter;
  iftAdjRel *A=iftSpheric(sqrtf(2.0));

  if (label->zsize > 1)
      iftError("Use iftShapeBasedInterp for 3D images", "iftShapeBasedInterp2D");
  
  if ((sx <= 0.0)||(sy <= 0.0))
      iftError("Invalid scale factors", "iftShapeBasedInterp2D");
  
  /* Mark existing labels */
  
  for (p=0; p < label->n; p++) 
    obj_code[label->val[p]]=1;
  
  /* Compute shape-based interpolation for each existing label */
  
  xsize= iftRound(fabs(sx * label->xsize));
  ysize= iftRound(fabs(sy * label->ysize));
  
  nlabel = iftCreateImage(xsize,ysize,1);

    int label_max_val = iftMaximumValue(label);
    for (l=1; l <= label_max_val; l++){
    if (obj_code[l]==1){ 
      iftBoundingBox mbb = iftMinObjectBoundingBox(label, l, NULL);
      pos = mbb.begin;
      bin = iftExtractROI(label, mbb);
      sdist = iftSignedDistTrans(bin,A);
      iftDestroyImage(&bin);
      inter = iftFInterp2D(sdist,sx,sy);
      iftDestroyFImage(&sdist);
      bin   = iftFThreshold(inter, 0, IFT_INFINITY_FLT, 1);
      iftDestroyFImage(&inter);
      pos.x = iftRound(pos.x * sx); pos.y = iftRound(pos.y * sy);
      iftInsertObject(bin,nlabel,l,pos);
      iftDestroyImage(&bin);
    }
  }
  
  iftFree(obj_code);
  iftDestroyAdjRel(&A);
  
  return(nlabel);
}

iftPlane *iftFindBestObjectCutPlane(iftImage *obj, iftImage *weight) {
    iftBestObjectCutPlaneProblem *prob       = NULL;
    iftMSPS                      *msps       = NULL;
    iftPlane                     *cutplane   = iftCreatePlane();
    float                        delta_max[] = {90.0, 90.0, 10.0, 10.0, 10.0};
    int                          i, wmax     = IFT_INFINITY_INT_NEG, p;
    iftVoxel                     pos         = {0, 0, 0};

    /* Find a local maximum as initial position */

    for (p = 0; p < weight->n; p++)
        if (weight->val[p] > wmax) {
            pos  = iftGetVoxelCoord(weight, p);
            wmax = weight->val[p];
        }

    /* Initialize the best object cut plane problem */

    prob = iftCreateBestObjectCutPlaneProblem(obj, pos, weight);

    /* Initialize context for MSPS optimization */

    msps = iftCreateMSPS(5, 3, iftObjectCutPlaneFitness, prob);

    /* Initialize theta and error */

    for (i = 0; i < 5; i++)
        msps->theta[i] = 0.0;

    iftObjectCutPlaneMSDeltas(msps, 5.0, delta_max);

    /* Find the best object cut plane */

    printf("Optimum fitness value is %lf\n", iftMSPSMax(msps));

    printf("Best plane %f %f %f %f %f %f\n", prob->pl->normal.x, prob->pl->normal.y,
           prob->pl->normal.z, prob->pl->pos.x, prob->pl->pos.y, prob->pl->pos.z);

    /* Return plane for object cut */

    cutplane->normal = prob->pl->normal;
    cutplane->pos    = prob->pl->pos;

    iftDestroyBestObjectCutPlaneProblem(&prob);
    iftDestroyMSPS(&msps);

    return (cutplane);
}

iftPlane *iftFindBestCutPlane(iftImage *img, iftPoint pos, int xviewsize, int yviewsize)
{
	iftBestCutPlaneProblem *prob=NULL; 
	iftMSPS      *msps=NULL;
	iftPlane     *cutplane=iftCreatePlane();
	float delta_max[] = {20.0,20.0};
	int i;
	iftVoxel u;

	u.x = iftRound(pos.x);
	u.y = iftRound(pos.y);
	u.z = iftRound(pos.z);

	if (!iftValidVoxel(img,u))
        iftError("Invalid plane position", "iftFindBestCutPlane");

	/* Initialize the best cut plane problem */

	prob = iftCreateBestCutPlaneProblem(img,pos,xviewsize,yviewsize);

	/* Initialize context for MSPS optimization */

	msps   = iftCreateMSPS(2,3,iftCutPlaneFitness,prob);
	
	/* Initialize theta and error */

	for (i=0; i < 2; i++) 
		msps->theta[i] = 0.0;

	iftCutPlaneMSDeltas(msps,2.0,delta_max);

	/* Find the best cut plane */

	printf("Optimum fitness value is %lf\n",iftMSPSMax(msps));

	printf("Best plane %f %f %f %f %f %f\n",prob->pl->normal.x,prob->pl->normal.y,prob->pl->normal.z,prob->pl->pos.x,prob->pl->pos.y,prob->pl->pos.z);

	/* Return plane for cut */

	cutplane->normal = prob->pl->normal;
	cutplane->pos    = prob->pl->pos;

	iftDestroyBestCutPlaneProblem(&prob);
	iftDestroyMSPS(&msps);
	 
	return(cutplane);
}

iftImage *iftResliceOnPrincipalAxis(iftImage *img, iftImage *bin)
{
	iftImage *rimg;
	iftVector paxis=iftPrincipalAxis(bin);
	iftPoint  center=iftGeometricCenter(bin);
	iftPlane *pl=iftCreatePlane();
	int       size = iftObjectDiagonal(bin);

	iftSetPlanePos(pl,center.x,center.y,center.z);
	iftSetPlaneOrient(pl,paxis.x,paxis.y,paxis.z);  
	rimg = iftResliceImage(img, pl, size, size, size);
	iftDestroyPlane(&pl);
	
	return(rimg);
}



iftImage *iftResliceImageSimple(const iftImage *img, int dx, int dy, int dz, iftAxisOrder axis_order)
{
  iftImage *rimg=NULL;
  iftVoxel  u, ui,uf;
  int       p, q; 

  if (!iftIs3DImage(img))
      iftError("Input image must be 3D", "iftResliceImageSimple");
  
  if (abs(dx*dy*dz)!=1)
      iftError("Increments dx, dy, and dz must be in {-1,1}", "iftResliceImageSimple");
    
  ui.x = ui.y = ui.z = 0;
  uf.x = img->xsize;
  uf.y = img->ysize;
  uf.z = img->zsize;

  if (dx == -1) {ui.x = img->xsize-1; uf.x = -1;}
  if (dy == -1) {ui.y = img->ysize-1; uf.y = -1;}
  if (dz == -1) {ui.z = img->zsize-1; uf.z = -1;}

  switch (axis_order) {
  case XYZ:

    rimg = iftCreateImage(img->xsize, img->ysize, img->zsize);
    q    = 0;
    for (u.z=ui.z; u.z != uf.z; u.z=u.z+dz){
      for (u.y=ui.y; u.y != uf.y; u.y=u.y+dy){
	for (u.x=ui.x; u.x != uf.x; u.x=u.x+dx){
	  p = iftGetVoxelIndex(img,u);
	  rimg->val[q]=img->val[p];
	  q++;
	}
      }
    }	  
    rimg->dx = img->dx;
    rimg->dy = img->dy;
    rimg->dz = img->dz;

    break;

  case XZY:

    rimg = iftCreateImage(img->xsize, img->zsize, img->ysize);
    q    = 0;
    for (u.y=ui.y; u.y != uf.y; u.y=u.y+dy){
      for (u.z=ui.z; u.z != uf.z; u.z=u.z+dz){
	for (u.x=ui.x; u.x != uf.x; u.x=u.x+dx){
	  p = iftGetVoxelIndex(img,u);
	  rimg->val[q]=img->val[p];
	  q++;
	}
      }
    }	  

    rimg->dx = img->dx;
    rimg->dy = img->dz;
    rimg->dz = img->dy;

    break;

  case YXZ:

    rimg = iftCreateImage(img->ysize, img->xsize, img->zsize);
    q    = 0;
    for (u.z=ui.z; u.z != uf.z; u.z=u.z+dz){
      for (u.x=ui.x; u.x != uf.x; u.x=u.x+dx){
	for (u.y=ui.y; u.y != uf.y; u.y=u.y+dy){
	  p = iftGetVoxelIndex(img,u);
	  rimg->val[q]=img->val[p];
	  q++;
	}
      }
    }	  

    rimg->dx = img->dy;
    rimg->dy = img->dx;
    rimg->dz = img->dz;

    break;

  case YZX:

    rimg = iftCreateImage(img->ysize, img->zsize, img->xsize);
    q    = 0;
    for (u.x=ui.x; u.x != uf.x; u.x=u.x+dx){
      for (u.z=ui.z; u.z != uf.z; u.z=u.z+dz){
	for (u.y=ui.y; u.y != uf.y; u.y=u.y+dy){
	  p = iftGetVoxelIndex(img,u);
	  rimg->val[q]=img->val[p];
	  q++;
	}
      }
    }	  

    rimg->dx = img->dy;
    rimg->dy = img->dz;
    rimg->dz = img->dx;

    break;

  case ZXY:

    rimg = iftCreateImage(img->zsize, img->xsize, img->ysize);
    q    = 0;
    for (u.y=ui.y; u.y != uf.y; u.y=u.y+dy){
      for (u.x=ui.x; u.x != uf.x; u.x=u.x+dx){
	for (u.z=ui.z; u.z != uf.z; u.z=u.z+dz){
	  p = iftGetVoxelIndex(img,u);
	  rimg->val[q]=img->val[p];
	  q++;
	}
      }
    }	  

    rimg->dx = img->dz;
    rimg->dy = img->dx;
    rimg->dz = img->dy;

    break;

  case ZYX:

    rimg = iftCreateImage(img->zsize, img->ysize, img->xsize);
    q    = 0;
    for (u.x=ui.x; u.x != uf.x; u.x=u.x+dx){
      for (u.y=ui.y; u.y != uf.y; u.y=u.y+dy){
	for (u.z=ui.z; u.z != uf.z; u.z=u.z+dz){
	  p = iftGetVoxelIndex(img,u);
	  rimg->val[q]=img->val[p];
	  q++;
	}
      }
    }	  
    rimg->dx = img->dz;
    rimg->dy = img->dy;
    rimg->dz = img->dx;

    break;

  default:
      iftError("Axis order must be in {XYZ, XZY, YXZ, YZX, ZXY, ZYX}", "iftResliceImageSimple");
      }
  
  iftCopyVoxelSize(img, rimg);

  
  return(rimg);
}



iftImage *iftResliceImageByTransMatrix(const iftImage *img, const iftMatrix *R) {
    if ((R->nrows != 4) || (R->ncols != 4))
        iftError("Transf. Matrix must be 4x4. Insted, it is %dx%d", "iftResliceImageByTransMatrix",
                  R->nrows, R->ncols);

    int sizes[3];
    sizes[0] = img->xsize;
    sizes[1] = img->ysize;
    sizes[2] = img->zsize;

    float voxel_sizes[3];
    voxel_sizes[0] = img->dx;
    voxel_sizes[1] = img->dy;
    voxel_sizes[2] = img->dz;

    iftImageDomain new_dom = {0, 0, 0};
    iftVoxelSize new_voxel_sizes = {1.0, 1.0, 1.0};

    // get the domain of the resliced image
    // assume that there only value != 0 (-1 or 1) in each row
    for (int j = 0; j <= 2; j++) {
        if ((iftMatrixElem(R, j, 0) == -1) || (iftMatrixElem(R, j, 0) == 1)) {
            new_dom.xsize = sizes[j];
            new_voxel_sizes.dx = voxel_sizes[j];
            continue;
        }
        if ((iftMatrixElem(R, j, 1) == -1) || (iftMatrixElem(R, j, 1) == 1)) {
            new_dom.ysize = sizes[j];
            new_voxel_sizes.dy = voxel_sizes[j];
            continue;
        }
        if ((iftMatrixElem(R, j, 2) == -1) || (iftMatrixElem(R, j, 2) == 1)) {
            new_dom.zsize = sizes[j];
            new_voxel_sizes.dz = voxel_sizes[j];
            continue;
        }
    }


    iftImage *out_img = iftCreateImage(new_dom.xsize, new_dom.ysize, new_dom.zsize);
    out_img->dx = new_voxel_sizes.dx;
    out_img->dy = new_voxel_sizes.dy;
    out_img->dz = new_voxel_sizes.dz;

    iftCopyImageVoxelFunc copy_func = iftCopyGrayImageVoxel;
    if (iftIsColorImage(img))
        copy_func = iftCopyColorImageVoxel;


    // after transformation, the points can have negative coordinates
    // we then use an offset to put it on image domain
    // offset is 0 or abs(size-1)
    iftVoxel offset = iftTransformVoxel(R, (iftVoxel) {img->xsize-1, img->ysize-1, img->zsize-1});
    offset.x = (offset.x < 0) ? abs(offset.x) : 0;
    offset.y = (offset.y < 0) ? abs(offset.y) : 0;
    offset.z = (offset.z < 0) ? abs(offset.z) : 0;

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++) {
        iftVoxel u = iftGetVoxelCoord(img, p);
        iftVoxel v = iftTransformVoxel(R, u);
        v = (iftVoxel) iftVectorSum(v, offset);

        copy_func(img, u, out_img, v);
    }

    return out_img;
}


int iftImageValueAtPoint(const iftImage *img, iftPoint P) {
	/*
     * Reference points considered in this function.
     *           u4-----------------u5
     *          /|                 /|
     *        /  |               /  |
     *      u0---|-------------u1   |
     *       |   |     P       |    |
     *       |   u6------------|---u7
     *       |  /              |   /
     *       |/                | /
     *      u2----------------u3
     *
     *      u0--------------u1
     *      |    ^          |
     *      |    | dy       |
     *      |    |          |   dx = P.x - u0.x
     *      |    V          |   dy = P.y - u0.y
     *      |<-->P          |   dz = P.z - u0.z
     *      | dx            |
     *      |               |
     *      u2--------------u3
     */
	if (!iftValidPoint(img, P))
	    return 0;
	
	iftVoxel u[8];
	u[0] = (iftVoxel) {(int) P.x, (int) P.y, (int) P.z};  // top-left most voxel to the point P
	
	u[1] = (iftVoxel) {u[0].x + 1, u[0].y,     u[0].z};
	u[2] = (iftVoxel) {u[0].x,     u[0].y + 1, u[0].z};
	u[3] = (iftVoxel) {u[0].x + 1, u[0].y + 1, u[0].z};
	u[4] = (iftVoxel) {u[0].x,     u[0].y,     u[0].z + 1};
	u[5] = (iftVoxel) {u[0].x + 1, u[0].y,     u[0].z + 1};
	u[6] = (iftVoxel) {u[0].x,     u[0].y + 1, u[0].z + 1};
	u[7] = (iftVoxel) {u[0].x + 1, u[0].y + 1, u[0].z + 1};
	
	for (int i = 0; i <= 7; i++)
		if (!iftValidVoxel(img, u[i]))
			return iftImgVoxelVal(img, u[0]);
	
	float dx = P.x - u[0].x;
	float dy = P.y - u[0].y;
	float dz = P.z - u[0].z;
	
	// after some algebric manipulations in the interpolation formulation, we have these equations
	float I01 = iftImgVoxelVal(img, u[0]) + (dx * (iftImgVoxelVal(img, u[1]) - iftImgVoxelVal(img, u[0])));
	float I23 = iftImgVoxelVal(img, u[2]) + (dx * (iftImgVoxelVal(img, u[3]) - iftImgVoxelVal(img, u[2])));
	float I45 = iftImgVoxelVal(img, u[4]) + (dx * (iftImgVoxelVal(img, u[5]) - iftImgVoxelVal(img, u[4])));
	float I67 = iftImgVoxelVal(img, u[6]) + (dx * (iftImgVoxelVal(img, u[7]) - iftImgVoxelVal(img, u[6])));
	float I0123 = I01 + (dy * (I23 - I01));
	float I4567 = I45 + (dy * (I67 - I45));
	
	return iftRound(I0123 + (dz * (I4567 - I0123)));
}


int iftImageValueAtPointNearestNeighbor(const iftImage *img, iftPoint P) {
    if (!iftValidPoint(img, P))
        return 0;
    
    iftVoxel u = {iftRound(P.x), iftRound(P.y), iftRound(P.z)};
    
    u.x = iftMin(iftMax(u.x, 0), img->xsize - 1);
    u.y = iftMin(iftMax(u.y, 0), img->ysize - 1);
    u.z = iftMin(iftMax(u.z, 0), img->zsize - 1);
    
    return iftImgVoxelVal(img, u);
}


int iftImageValueAtPoint2D(const iftImage *img, iftPoint P) {
    /*
     * Reference points considered in this function.
     *      u0--------------u1
     *      |    ^          |
     *      |    | dy       |
     *      |    |          |   dx = P.x - u0.x
     *      |    V          |   dy = P.y - u0.y
     *      |<-->P          |   
     *      | dx            |
     *      |               |
     *      u2--------------u3
     */
    if (!iftValidPoint(img, P))
        return 0;
    
    iftVoxel u[4];
    u[0] = (iftVoxel) {(int) P.x, (int) P.y};  // top-left most voxel to the point P
    
    u[1] = (iftVoxel) {u[0].x + 1, u[0].y};
    u[2] = (iftVoxel) {u[0].x,     u[0].y + 1};
    u[3] = (iftVoxel) {u[0].x + 1, u[0].y + 1};
    
    for (int i = 0; i <= 3; i++)
        if (!iftValidVoxel(img, u[i]))
            return iftImgVoxelVal(img, u[0]);
    
    float dx = P.x - u[0].x;
    float dy = P.y - u[0].y;
    
    // after some algebric manipulations in the interpolation formulation, we have these equations
    float I01 = iftImgVoxelVal(img, u[0]) + (dx * (iftImgVoxelVal(img, u[1]) - iftImgVoxelVal(img, u[0])));
    float I23 = iftImgVoxelVal(img, u[2]) + (dx * (iftImgVoxelVal(img, u[3]) - iftImgVoxelVal(img, u[2])));
    
    return iftRound(I01 + (dy * (I23 - I01)));
}


int iftImageValueAtPoint2DNearestNeighbor(const iftImage *img, iftPoint P) {
    iftVoxel u = {iftRound(P.x), iftRound(P.y), 0};
    
    u.x = iftMin(iftMax(u.x, 0), img->xsize-1);
    u.y = iftMin(iftMax(u.y, 0), img->ysize-1);
    
    return iftImgVoxelVal(img, u);
}




