#include "iftAdjacency.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/List.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/io/Stream.h"
#include "iftImage.h"
#include "iftSort.h"

/**
@file
@brief A description is missing here
*/

/*! \file
 * Author: Alexandre Falcao
 * Date: Aug 22, 2011
 * Last Update: Aug 22, 2011
 */

/* 3D adjacency relations */

iftAdjRel  *iftCreateAdjRel(int n) /*! \brief Allocates memory for a
				      3D adjacency relation */
{
  iftAdjRel *A=(iftAdjRel *)iftAlloc(1,sizeof(iftAdjRel ));

  A->dx = (int *)iftAllocIntArray(n);
  A->dy = (int *)iftAllocIntArray(n);
  A->dz = (int *)iftAllocIntArray(n);
  A->dt = (int *)iftAllocIntArray(n);
  A->n  = n;

  return(A);
}

void     iftDestroyAdjRel(iftAdjRel **A) /*! \brief Deallocates memory for a
					     3D adjacency relation */
{
  iftAdjRel *aux = *A;

  if (aux != NULL){
    if (aux->dx != NULL) iftFree(aux->dx);
    if (aux->dy != NULL) iftFree(aux->dy);
    if (aux->dz != NULL) iftFree(aux->dz);
    if (aux->dt != NULL) iftFree(aux->dt);    
    iftFree(aux);
    *A = NULL;
  }
}


iftAdjRel *iftSpheric(float r) {
    int r0 = (int) r;
    float r2 = (int) (r*r + 0.5);

    int n = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx =- r0; dx <= r0; dx++)
                if ( ((dx*dx) + (dy*dy) + (dz*dz)) <= r2)
                    n++;


    iftAdjRel *A = iftCreateAdjRel(n);

    int i = 0;
    int i0 = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx = -r0; dx <= r0; dx++)
                if ( ((dx*dx) + (dy*dy) + (dz*dz)) <= r2) {
                    A->dx[i] = dx;
                    A->dy[i] = dy;
                    A->dz[i] = dz;
                
                if ((dx == 0) && (dy == 0) && (dz == 0))
                    i0 = i;
                i++;
            }

    // shift to right and place central voxel at first
    for (int i = i0; i > 0; i--) {
        int dx = A->dx[i];
        int dy = A->dy[i];
        int dz = A->dz[i];
        A->dx[i] = A->dx[i-1];
        A->dy[i] = A->dy[i-1];
        A->dz[i] = A->dz[i-1];
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
        A->dz[i-1] = dz;
    }


    // sort by radius, so the 6 closest neighbors will come first
    float *dr = iftAllocFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dr[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];

    iftIntArray *idxs = iftIntRange(0, A->n-1, 1);
    iftFQuickSort(dr, idxs->val, 0, A->n-1, IFT_INCREASING);
    iftAdjRel *Asort = iftCreateAdjRel(A->n);

    for (int i = 0; i < A->n; i++) {
        int idx = idxs->val[i];
        Asort->dx[i] = A->dx[idx];
        Asort->dy[i] = A->dy[idx];
        Asort->dz[i] = A->dz[idx];
    }

    iftFree(dr);
    iftDestroyIntArray(&idxs);
    iftDestroyAdjRel(&A);

    return Asort;
}

iftAdjRel *iftHyperSpheric(float r) {
    int r0 = (int) r;
    float r2 = (int) (r*r + 0.5);

    int n = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx =- r0; dx <= r0; dx++)
               for (int dt =- r0; dt <= r0; dt++)
                	if ( ((dx*dx) + (dy*dy) + (dz*dz) + (dt*dt)) <= r2)
                    	n++;


    iftAdjRel *A = iftCreateAdjRel(n);

    int i = 0;
    int i0 = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx = -r0; dx <= r0; dx++)               
            	for (int dt =- r0; dt <= r0; dt++)
                	if ( ((dx*dx) + (dy*dy) + (dz*dz) + (dt*dt)) <= r2) {
		                 A->dx[i] = dx;
		                 A->dy[i] = dy;
		                 A->dz[i] = dz;
		                 A->dt[i] = dt;
                
				          if ((dx == 0) && (dy == 0) && (dz == 0) && (dt == 0))
				              i0 = i;
				          i++;
            		}

    // shift to right and place central voxel at first
    for (int i = i0; i > 0; i--) {
        int dx = A->dx[i];
        int dy = A->dy[i];
        int dz = A->dz[i];
        int dt = A->dt[i];        
        A->dx[i] = A->dx[i-1];
        A->dy[i] = A->dy[i-1];
        A->dz[i] = A->dz[i-1];
        A->dt[i] = A->dt[i-1];        
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
        A->dz[i-1] = dz;
        A->dt[i-1] = dt;        
    }


    // sort by radius, so the 6 closest neighbors will come first
    float *dr = iftAllocFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dr[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i] + A->dt[i]*A->dt[i];

    iftIntArray *idxs = iftIntRange(0, A->n-1, 1);
    iftFQuickSort(dr, idxs->val, 0, A->n-1, IFT_INCREASING);
    iftAdjRel *Asort = iftCreateAdjRel(A->n);

    for (int i = 0; i < A->n; i++) {
        int idx = idxs->val[i];
        Asort->dx[i] = A->dx[idx];
        Asort->dy[i] = A->dy[idx];
        Asort->dz[i] = A->dz[idx];
        Asort->dt[i] = A->dt[idx];        
    }

    iftFree(dr);
    iftDestroyIntArray(&idxs);
    iftDestroyAdjRel(&A);

    return Asort;
}

iftAdjRel *iftSemiHyperSpheric(float r, bool positive)
{
    int r0 = (int) r;
    float r2 = (int) (r*r + 0.5);

    int dir_factor = positive? 1:-1;

    int n = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx =- r0; dx <= r0; dx++)
                for (int dt =- r0; dt <= r0; dt++)
                    if ( ((dx*dx) + (dy*dy) + (dz*dz) + (dt*dt)) <= r2 && dt*dir_factor >= 0)
                        n++;


    iftAdjRel *A = iftCreateAdjRel(n);

    int i = 0;
    int i0 = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx = -r0; dx <= r0; dx++)
                for (int dt =- r0; dt <= r0; dt++)
                    if ( ((dx*dx) + (dy*dy) + (dz*dz) + (dt*dt)) <= r2 && dt*dir_factor >= 0) {
                        A->dx[i] = dx;
                        A->dy[i] = dy;
                        A->dz[i] = dz;
                        A->dt[i] = dt*dir_factor;

                        if ((dx == 0) && (dy == 0) && (dz == 0) && (dt == 0))
                            i0 = i;
                        i++;
                    }

    // shift to right and place central voxel at first
    for (int i = i0; i > 0; i--) {
        int dx = A->dx[i];
        int dy = A->dy[i];
        int dz = A->dz[i];
        int dt = A->dt[i];
        A->dx[i] = A->dx[i-1];
        A->dy[i] = A->dy[i-1];
        A->dz[i] = A->dz[i-1];
        A->dt[i] = A->dt[i-1];
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
        A->dz[i-1] = dz;
        A->dt[i-1] = dt;
    }


    // sort by radius, so the 6 closest neighbors will come first
    float *dr = iftAllocFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dr[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i] + A->dt[i]*A->dt[i];

    iftIntArray *idxs = iftIntRange(0, A->n-1, 1);
    iftFQuickSort(dr, idxs->val, 0, A->n-1, IFT_INCREASING);
    iftAdjRel *Asort = iftCreateAdjRel(A->n);

    for (int i = 0; i < A->n; i++) {
        int idx = idxs->val[i];
        Asort->dx[i] = A->dx[idx];
        Asort->dy[i] = A->dy[idx];
        Asort->dz[i] = A->dz[idx];
        Asort->dt[i] = A->dt[idx];
    }

    iftFree(dr);
    iftDestroyIntArray(&idxs);
    iftDestroyAdjRel(&A);

    return Asort;
}


iftAdjRel  *iftHemispheric(float r, char axis, int direction) /*! \brief Creates a 3D half-ball of radius r as
                                                            adjacency relation, in the corresponding axis and direction.
                                                            This adjacency is useful for segmenting a volume in a single
                                                            direction, e.g., a video-volume where z is time. */
{
  iftAdjRel *A=NULL;
  int i,j,k,n,r0,d,dx,dy,dz,dx0,dy0,dz0,xr1,yr1,zr1,i0=0;
  float *dr,r2,aux;

  n=0;
  r0 = (int)r;
  r2  = (int)(r*r + 0.5);

  dz0 = dy0 = dx0 = -r0;
  zr1 = yr1 = xr1 = r0;
  if(axis == 'z')
  {
    dz0 = (direction >= 0) ? 0 : -r0;
    zr1 = (direction >= 0) ? r0 : 0;
  }
  else if(axis == 'y')
  {
    dy0 = (direction >= 0) ? 0 : -r0;
    yr1 = (direction >= 0) ? r0 : 0;
  }
  else if(axis == 'x')
  {
    dx0 = (direction >= 0) ? 0 : -r0;
    xr1 = (direction >= 0) ? r0 : 0;
  }

  for(dz=dz0;dz<=zr1;dz++)
    for(dy=dy0;dy<=yr1;dy++)
      for(dx=dx0;dx<=xr1;dx++)
      if(((dx*dx)+(dy*dy)+(dz*dz)) <= r2)
	n++;

  A = iftCreateAdjRel(n);
  i=0;
  for(dz=dz0;dz<=zr1;dz++)
    for(dy=dy0;dy<=yr1;dy++)
      for(dx=dx0;dx<=xr1;dx++)
	if(((dx*dx)+(dy*dy)+(dz*dz)) <= r2){
	  A->dx[i]=dx;
	  A->dy[i]=dy;
	  A->dz[i]=dz;
	  if ((dx==0)&&(dy==0)&&(dz==0))
	    i0 = i;
	  i++;
	}

  /* shift to right and place central voxel at first */

  for (i=i0; i > 0; i--) {
    dx = A->dx[i];
    dy = A->dy[i];
    dz = A->dz[i];
    A->dx[i] = A->dx[i-1];
    A->dy[i] = A->dy[i-1];
    A->dz[i] = A->dz[i-1];
    A->dx[i-1] = dx;
    A->dy[i-1] = dy;
    A->dz[i-1] = dz;
  }

  /* sort by radius, so the 6 closest neighbors will come first */

  dr = iftAllocFloatArray(A->n);
  for (i=0; i < A->n; i++) {
    dr[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];
  }

  for (i=1; i < A->n-1; i++){
    k = i;
    for (j=i+1; j < A->n; j++)
      if (dr[j] < dr[k]){
	k = j;
      }
    aux   = dr[i];
    dr[i] = dr[k];
    dr[k] = aux;
    d        = A->dx[i];
    A->dx[i] = A->dx[k];
    A->dx[k] = d;
    d        = A->dy[i];
    A->dy[i] = A->dy[k];
    A->dy[k] = d;
    d        = A->dz[i];
    A->dz[i] = A->dz[k];
    A->dz[k] = d;
  }

  iftFree(dr);


  return(A);
}

iftAdjRel *iftSphericEdges(float r)
{
	int i,j,n=0;
	iftAdjRel* Ac = iftSpheric(r);

	for(i=0;i<Ac->n;i++)
		if (Ac->dz[i]*(2*r+1)*(2*r+1) + Ac->dy[i]*(2*r+1) + Ac->dx[i] >= 0)
			n++;

	iftAdjRel* A = iftCreateAdjRel(n);
	for(i=0,j=0;i<Ac->n;i++)
		if (Ac->dz[i]*(2*r+1)*(2*r+1) + Ac->dy[i]*(2*r+1) + Ac->dx[i] >= 0)
		{
			A->dx[j] = Ac->dx[i];
			A->dy[j] = Ac->dy[i];
			A->dz[j] = Ac->dz[i];
			j++;
		}

	iftDestroyAdjRel(&Ac);
	return A;
}


iftAdjRel *iftRectangular(int xsize, int ysize)
{
  iftAdjRel *A;
  int i,dx,dy,n,i0=0;

  n = 0;
  for(dy=-ysize/2;dy<=ysize/2;dy++)
    for(dx=-xsize/2;dx<=xsize/2;dx++)
      n++;

  A = iftCreateAdjRel(n);

  i = 0;
  for(dy=-ysize/2;dy<=ysize/2;dy++)
    for(dx=-xsize/2;dx<=xsize/2;dx++){
      A->dx[i] = dx;
      A->dy[i] = dy;
      A->dz[i] =  0;
      if ((dx==0)&&(dy==0))
	i0 = i;
      i++;
    }

  /* shift to right and place central point (origin) at first */

  for (i=i0; i > 0; i--) {
    dx = A->dx[i];
    dy = A->dy[i];
    A->dx[i] = A->dx[i-1];
    A->dy[i] = A->dy[i-1];
    A->dx[i-1] = dx;
    A->dy[i-1] = dy;
  }

  return(A);
}


iftAdjRel *iftRectangularWithDilation(int xsize, int ysize, int sx, int sy)
{
  /*
  Kernel_size: (int xsize, int ysize)
  dilation: (int sx, int sy)
  */

  iftAdjRel *A;
  int i,dx,dy,n,i0=0;

  if ((sx < 1)||(sy < 1))
    iftError("Dilation rates must be greater than 1","iftRectangularWithDilation");
  
  n = 0;
  for(dy=-ysize/2;dy<=ysize/2;dy++) // (-2,2) [-2,-1,0,1,2]
    for(dx=-xsize/2;dx<=xsize/2;dx++) // (-2,2) [-2,-1,0,1,2]
      n++; // 25

  A = iftCreateAdjRel(n);

  i = 0;
  for(dy=-ysize/2;dy<=ysize/2;dy++)
    for(dx=-xsize/2;dx<=xsize/2;dx++){
      A->dx[i] = dx*sx; // se sx=7, [-2*7,-1*7,0*7,1*7,2*7] = [-14,-7,0,7,14]
      A->dy[i] = dy*sy; // se sy=7, [-2*7,-1*7,0*7,1*7,2*7] = [-14,-7,0,7,14]
      A->dz[i] =  0;
      if ((dx==0)&&(dy==0))
	      i0 = i; // guarda a posição central
      i++;
    } 

  /* shift to right and place central point (origin) at first */

  for (i=i0; i > 0; i--) {
    dx = A->dx[i];
    dy = A->dy[i];
    A->dx[i] = A->dx[i-1];
    A->dy[i] = A->dy[i-1];
    A->dx[i-1] = dx;
    A->dy[i-1] = dy;
  }

  return(A);
}



iftAdjRel *iftCross(int xsize, int ysize)
{
  iftAdjRel *A;
  int i,dx,dy,n,i0=0;

  n = 0;
  for(dy=-ysize/2;dy<=ysize/2;dy++)
    n++;

  for(dx=-xsize/2;dx<=xsize/2;dx++)
    n++;

  n--;
  
  A = iftCreateAdjRel(n);

  i = 0;
  for(dx=-xsize/2;dx<=xsize/2;dx++){
    A->dx[i] = dx;
    A->dy[i] = 0;
    A->dz[i] = 0;
    if (dx==0)
      i0 = i;
    i++;
  }

  for(dy=-ysize/2;dy<=ysize/2;dy++){
    if (dy!=0){
      A->dx[i] = 0;
      A->dy[i] = dy;
      A->dz[i] = 0;
      i++;
    }
  }

  /* shift to right and place central point (origin) at first */

  for (i=i0; i > 0; i--) {
    dx = A->dx[i];
    dy = A->dy[i];
    A->dx[i] = A->dx[i-1];
    A->dy[i] = A->dy[i-1];
    A->dx[i-1] = dx;
    A->dy[i-1] = dy;
  }

  return(A);
}

iftAdjRel *iftCopyAdjacency(const iftAdjRel *A) {
  iftAdjRel *B = iftCreateAdjRel(A->n);
  int i;
  for (i=0; i < A->n; i++) {
    B->dx[i] = A->dx[i];
    B->dy[i] = A->dy[i];
    B->dz[i] = A->dz[i];
    B->dt[i] = A->dt[i];    
  }

  return(B);
}


iftAdjRel *iftCuboid(int xsize, int ysize, int zsize)
{
  iftAdjRel *A;
  int i,dx,dy,dz,n,i0=0;

  n = 0;
  for(dz=-zsize/2;dz<=zsize/2;dz++)
    for(dy=-ysize/2;dy<=ysize/2;dy++)
      for(dx=-xsize/2;dx<=xsize/2;dx++)
      n++;

  A = iftCreateAdjRel(n);

  i = 0;
  for(dz=-zsize/2;dz<=zsize/2;dz++)
    for(dy=-ysize/2;dy<=ysize/2;dy++)
      for(dx=-xsize/2;dx<=xsize/2;dx++){
	A->dx[i] = dx;
	A->dy[i] = dy;
	A->dz[i] = dz;
	if ((dx==0)&&(dy==0)&&(dz==0))
	  i0 = i;
	i++;
      }

  /* shift to right and place central point (origin) at first */

  for (i=i0; i > 0; i--) {
    dx = A->dx[i];
    dy = A->dy[i];
    dz = A->dz[i];
    A->dx[i] = A->dx[i-1];
    A->dy[i] = A->dy[i-1];
    A->dz[i] = A->dz[i-1];
    A->dx[i-1] = dx;
    A->dy[i-1] = dy;
    A->dz[i-1] = dz;
  }

  return(A);
}

iftAdjRel *iftHyperCuboid(int xsize, int ysize, int zsize, int tsize) {
    iftAdjRel *A;
    int i,dx,dy,dz,dt,n,i0=0;

    n = 0;
    for(dz=-zsize/2;dz<=zsize/2;dz++)
        for(dy=-ysize/2;dy<=ysize/2;dy++)
            for(dx=-xsize/2;dx<=xsize/2;dx++)
                for(dt=-tsize/2;dt<=tsize/2;dt++)
                    n++;

    A = iftCreateAdjRel(n);

    i = 0;
    for(dz=-zsize/2;dz<=zsize/2;dz++)
        for(dy=-ysize/2;dy<=ysize/2;dy++)
            for(dx=-xsize/2;dx<=xsize/2;dx++)
                for(dt=-tsize/2;dt<=tsize/2;dt++){
                    A->dx[i] = dx;
                    A->dy[i] = dy;
                    A->dz[i] = dz;
                    A->dt[i] = dt;
                    if ((dx==0)&&(dy==0)&&(dz==0)&&(dt==0))
                        i0 = i;
                    i++;
                }

    /* shift to right and place central point (origin) at first */

    for (i=i0; i > 0; i--) {
        dx = A->dx[i];
        dy = A->dy[i];
        dz = A->dz[i];
        dt = A->dt[i];
        A->dx[i] = A->dx[i-1];
        A->dy[i] = A->dy[i-1];
        A->dz[i] = A->dz[i-1];
        A->dt[i] = A->dt[i-1];
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
        A->dz[i-1] = dz;
        A->dt[i-1] = dt;
    }

    return(A);
}

iftAdjRel *iftCuboidWithDilation(int xsize, int ysize, int zsize, int sx, int sy, int sz)
{
    iftAdjRel *A;
    int i,dx,dy,dz,n,i0=0;

    if ((sx < 1)||(sy < 1)||(sz < 1))
        iftError("Dilation rates must be greater than 1","iftCuboidWithDilation");

    n = 0;
    for(dz=-zsize/2;dz<=zsize/2;dz++)
        for(dy=-ysize/2;dy<=ysize/2;dy++)
            for(dx=-xsize/2;dx<=xsize/2;dx++)
                n++;

    A = iftCreateAdjRel(n);

    i = 0;
    for(dz=-zsize/2;dz<=zsize/2;dz++)
        for(dy=-ysize/2;dy<=ysize/2;dy++)
            for(dx=-xsize/2;dx<=xsize/2;dx++){
                A->dx[i] = dx*sx;
                A->dy[i] = dy*sy;
                A->dz[i] = dz*sz;
                if ((dx==0)&&(dy==0)&&(dz==0))
                    i0 = i;
                i++;
            }

    /* shift to right and place central point (origin) at first */

    for (i=i0; i > 0; i--) {
        dx = A->dx[i];
        dy = A->dy[i];
        dz = A->dz[i];
        A->dx[i] = A->dx[i-1];
        A->dy[i] = A->dy[i-1];
        A->dz[i] = A->dz[i-1];
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
        A->dz[i-1] = dz;
    }

    return(A);
}

iftAdjRel *iftCircular(float r) {
    int r0 = (int) r;
    float r2 = (int) (r*r + 0.5);
    
    int n = 0;
    for (int dy = -r0; dy <= r0; dy++)
        for (int dx = -r0; dx <= r0; dx++)
            if (((dx*dx) + (dy*dy)) <= r2)
                n++;

    iftAdjRel *A = iftCreateAdjRel(n);
    int i = 0;
    int i0 = 0;
    for (int dy = -r0; dy <= r0; dy++)
        for (int dx = -r0; dx <= r0; dx++)
            if (((dx*dx) + (dy*dy)) <= r2) {
                A->dx[i] = dx;
                A->dy[i] = dy;
                A->dz[i] = 0;

                if ((dx==0) && (dy==0))
                    i0 = i;
                i++;
            }

    // shift to right and place central pixel at first
    for (int i = i0; i > 0; i--) {
        int dx     = A->dx[i];
        int dy     = A->dy[i];
        A->dx[i]   = A->dx[i-1];
        A->dy[i]   = A->dy[i-1];
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
    }


    // sort by radius, so the 4 closest neighbors will come first
    float *dr = iftAllocFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dr[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i];


    iftIntArray *idxs = iftIntRange(0, A->n-1, 1);
    iftFQuickSort(dr, idxs->val, 0, A->n-1, IFT_INCREASING);
    iftAdjRel *Asort = iftCreateAdjRel(A->n);

    for (int i = 0; i < A->n; i++) {
        int idx = idxs->val[i];
        Asort->dx[i] = A->dx[idx];
        Asort->dy[i] = A->dy[idx];
        Asort->dz[i] = A->dz[idx];
    }

    iftFree(dr);
    iftDestroyIntArray(&idxs);
    iftDestroyAdjRel(&A);

    return Asort;
}


iftAdjRel *iftCircularEdges(float r)
{
	int i,j,n=0;
	iftAdjRel* Ac = iftCircular(r);

	for(i=0;i<Ac->n;i++)
		if (Ac->dy[i]*(2*r+1) + Ac->dx[i] >= 0)
			n++;

	iftAdjRel* A = iftCreateAdjRel(n);
	for(i=0,j=0;i<Ac->n;i++)
		if (Ac->dy[i]*(2*r+1) + Ac->dx[i] >= 0)
		{
		  A->dx[j] = Ac->dx[i];
		  A->dy[j] = Ac->dy[i];
		  A->dz[j] = Ac->dz[i];
		  j++;
		}

	iftDestroyAdjRel(&Ac);
	return A;
}

iftAdjRel *iftClockCircular(float r)
{
  iftAdjRel *A=NULL;
  int i,j,k,n,dx,dy,r0,r2,d,i0=0;
  float *da,*dr,aux;

  n=0;

  r0  = (int)r;
  r2  = (int)(r*r + 0.5);
  for(dy=-r0;dy<=r0;dy++)
    for(dx=-r0;dx<=r0;dx++)
      if(((dx*dx)+(dy*dy)) <= r2)
	n++;

  A = iftCreateAdjRel(n);
  i=0;
  for(dy=-r0;dy<=r0;dy++)
    for(dx=-r0;dx<=r0;dx++)
      if(((dx*dx)+(dy*dy)) <= r2){
	A->dx[i]=dx;
	A->dy[i]=dy;
	A->dz[i]=0;
	if ((dx==0)&&(dy==0))
	  i0 = i;
	i++;
      }

  /* Set clockwise */

  da = iftAllocFloatArray(A->n);
  dr = iftAllocFloatArray(A->n);
  for (i=0; i < A->n; i++) {
    dx = A->dx[i];
    dy = A->dy[i];
    dr[i] = (float)sqrtf((dx*dx) + (dy*dy));
    if (i != i0){
      da[i] = atan2(-dy,-dx) * 180.0 / IFT_PI;
      if (da[i] < 0.0)
	da[i] += 360.0;
    }
  }
  da[i0] = 0.0;
  dr[i0] = 0.0;

  /* place central pixel at first */

  aux    = da[i0];
  da[i0] = da[0];
  da[0]  = aux;
  aux    = dr[i0];
  dr[i0] = dr[0];
  dr[0]  = aux;
  d         = A->dx[i0];
  A->dx[i0] = A->dx[0];
  A->dx[0]  = d;
  d         = A->dy[i0];
  A->dy[i0] = A->dy[0];
  A->dy[0]  = d;

  /* sort by angle */

  for (i=1; i < A->n-1; i++){
    k = i;
    for (j=i+1; j < A->n; j++)
      if (da[j] < da[k]){
	k = j;
      }
    aux   = da[i];
    da[i] = da[k];
    da[k] = aux;
    aux   = dr[i];
    dr[i] = dr[k];
    dr[k] = aux;
    d   = A->dx[i];
    A->dx[i] = A->dx[k];
    A->dx[k] = d;
    d        = A->dy[i];
    A->dy[i] = A->dy[k];
    A->dy[k] = d;
  }

  /* sort by radius for each angle */

  for (i=1; i < A->n-1; i++){
    k = i;
    for (j=i+1; j < A->n; j++)
      if ((dr[j] < dr[k])&&(da[j]==da[k])){
	k = j;
      }
    aux   = dr[i];
    dr[i] = dr[k];
    dr[k] = aux;
    d        = A->dx[i];
    A->dx[i] = A->dx[k];
    A->dx[k] = d;
    d        = A->dy[i];
    A->dy[i] = A->dy[k];
    A->dy[k] = d;
  }

  iftFree(dr);
  iftFree(da);

  return(A);
}

iftAdjRel *iftRightSide(iftAdjRel *A, float r){
  iftAdjRel *R=NULL;
  int i;
  float d;

  for (i=0; i < A->n; i++)
    if (A->dz[i]!=0)
        iftError("It must be a 2D adjacency relation", "iftRightSide2D");

  /* Let p -> q be an arc represented by the increments dx,dy,dz. Its
     right side at distance r is given by the increments Dx = (-dy/d +
     dx/2)*r and Dy = dx/d + dy/2, where d=sqrt(dx²+dy²). */

  R = iftCreateAdjRel(A->n);
  for (i=0; i < R->n; i++){
    d  = sqrt(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
    if (d != 0){
      R->dx[i] = (int) ( (float)A->dx[i] / 2.0 * r - (float)A->dy[i] / d   * r );
      R->dy[i] = (int) ( (float)A->dx[i] / d   * r + (float)A->dy[i] / 2.0 * r );
      R->dz[i] = 0.0;
    }
  }

  return(R);
}

iftAdjRel *iftLeftSide(iftAdjRel *A, float r)
{
  iftAdjRel *L=NULL;
  int i;
  float d;

  for (i=0; i < A->n; i++)
    if (A->dz[i]!=0)
        iftError("It must be a 2D adjacency relation", "iftLeftSide2D");

  /* Let p -> q be an arc represented by the increments dx,dy. Its
     left side is given by the increments Dx = dy/d + dx/2 and Dy =
     -dx/d + dy/2, where d=sqrt(dx²+dy²). */

  L = iftCreateAdjRel(A->n);
  for (i=0; i < L->n; i++){
    d  = sqrt(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
    if (d != 0){
      L->dx[i] = (int)( (float)A->dx[i] / 2.0 * r + (float)A->dy[i] / d * r );
      L->dy[i] = (int)( (float)A->dy[i] / 2.0 * r - (float)A->dx[i] / d * r );
      L->dz[i] = 0;
    }
  }

  return(L);
}

void iftMaxAdjShifts(const iftAdjRel *A, int *dx, int *dy, int *dz)
{
  int i, d[3];

  *dx = *dy = *dz = 0.0;

  for (i=0; i < A->n; i++) {
    d[0] = abs(A->dx[i]);
    d[1] = abs(A->dy[i]);
    d[2] = abs(A->dz[i]);
    if (*dx < d[0]) *dx = d[0];
    if (*dy < d[1]) *dy = d[1];
    if (*dz < d[2]) *dz = d[2];
  }

}

iftFastAdjRel *iftCreateFastAdjRel(iftAdjRel *A, int *tby, int *tbz)
{
  iftFastAdjRel *F = (iftFastAdjRel *) iftAlloc(1,sizeof(iftFastAdjRel ));
  int i, po, qo;
  iftVoxel v;

  F->n  = A->n;
  F->dq = iftAllocIntArray(F->n);

  /* compute maximum displacements */

  iftMaxAdjShifts(A,&F->bx,&F->by,&F->bz);

  /* compute displacements to adjacent voxels */

  po = F->bx + tby[F->by] + tbz[F->bz];

  for (i=0; i < F->n; i++) {
    v.x = F->bx + A->dx[i];
    v.y = F->by + A->dy[i];
    v.z = F->bz + A->dz[i];
    qo  = v.x + tby[v.y] + tbz[v.z];
    F->dq[i] = qo-po;
  }

  return(F);

}

void iftDestroyFastAdjRel(iftFastAdjRel **F)
{
  iftFastAdjRel *aux=*F;

  if (aux != NULL) {
    iftFree(aux->dq);
    iftFree(aux);
    *F = NULL;
  }
}


void iftWriteAdjRel(iftAdjRel *A, char* filename){
    if (A != NULL) {
        FILE *fp = fopen(filename, "w");
        if(!fp)
            iftError("Unable to open file", "iftWriteAdjRel");

        fprintf(fp, "%d\n", A->n);
        for (int i = 0; i < A->n; i++)
            fprintf(fp,"%d %d %d\n",A->dx[i],A->dy[i],A->dz[i]);

        fclose(fp);
    }
}

void iftPrintAdjRel(iftAdjRel *A){
    if (A != NULL) {
        printf("%d\n", A->n);
        for (int i = 0; i < A->n; i++)
	  printf("%d %d %d\n",A->dx[i],A->dy[i],A->dz[i]);
    }
}

iftAdjRel* iftReadAdjRel(char* filename){
  FILE *fp = fopen(filename, "r");
  if(!fp)
      iftError("Unable to open file", "iftReadAdjRel");

  int n;
  if (fscanf(fp,"%d",&n)!=1)
      iftError("Reading error", "iftReadAdjRel");

  iftAdjRel *A = iftCreateAdjRel(n);

  for (int i=0; i < A->n; i++)
    if (fscanf(fp,"%d %d %d",&A->dx[i],&A->dy[i],&A->dz[i])!=3)
        iftError("Reading error", "iftReadAdjRel");

  fclose(fp);

  return(A);
}


iftAdjRel *iftReadAdjRelBinFile(const char *path) {
    if (path == NULL) {
        iftWarning("Pathname of the Adjacency Relation is NULL. Returning a NULL pointer",
                   "iftReadAdjRelBinFile");
        return NULL;
    }
    else {
        FILE *fp = fopen(path, "rb");
        if (!fp)
            iftError("Unable to open file", "iftReadAdjRelBinFile");
        
        int n;
        if (fread(&n, sizeof(int), 1, fp) != 1)
            iftError("Error when reading the Adj. Relation size", "iftReadAdjRelBinFile");
        
        iftAdjRel *A = iftCreateAdjRel(n);

        if (fread(A->dx, sizeof(int), A->n, fp) != A->n)
            iftError("Error when reading the dx array", "iftReadAdjRelBinFile");
        if (fread(A->dy, sizeof(int), A->n, fp) != A->n)
            iftError("Error when reading the dy array", "iftReadAdjRelBinFile");
        if (fread(A->dz, sizeof(int), A->n, fp) != A->n)
            iftError("Error when reading the dz array", "iftReadAdjRelBinFile");

        fclose(fp);

        return A;
    }
}



void iftWriteAdjRelBinFile(const iftAdjRel *A, const char *filename) {
    if ((A != NULL) && (filename != NULL)) {
        FILE *fp = fopen(filename, "wb");
        if (!fp)
            iftError("Unable to open file", "iftWriteAdjRelBinFile");

        if (fwrite(&A->n, sizeof(int), 1, fp) != 1)
            iftError("Error when writing adjacency size", "iftWriteAdjRelBinFile");
        
        if (fwrite(A->dx, sizeof(int), A->n, fp) != A->n)
            iftError("Error when writing dx array", "iftWriteAdjRelBinFile");

        if (fwrite(A->dy, sizeof(int), A->n, fp) != A->n)
            iftError("Error when writing dy array", "iftWriteAdjRelBinFile");

        if (fwrite(A->dz, sizeof(int), A->n, fp) != A->n)
            iftError("Error when writing dz array", "iftWriteAdjRelBinFile");

        fclose(fp);
    }
}



inline iftVoxel iftGetAdjacentVoxel(const iftAdjRel *A, iftVoxel u, int adj)
{
  iftVoxel v;

  v.x = u.x + A->dx[adj];
  v.y = u.y + A->dy[adj];
  v.z = u.z + A->dz[adj];
  v.t = u.t + A->dt[adj];

  return(v);
}

int iftIsAdjRel3D(iftAdjRel *A)
{
  int dx,dy,dz;
  iftMaxAdjShifts(A,&dx,&dy,&dz);
  if (dz > 0)
    return 1;
  else
    return 0;
}



iftAdjRel *iftExtractAdjacencyShell(iftAdjRel *A, double radius_min, double radius_max, uchar include_center) {
  int i, i0, n, found, dx, dy, dz;
  double r;
  iftAdjRel *Aout = NULL;

  /* Counting the number of elements within the closed interval [radius_min, radius_max] */
  n = 0;
  for(i = 0; i < A->n; i++) {
    r = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];

    if(r >= radius_min*radius_min && r <= radius_max*radius_max) {
      n++;
    }
  }


  i0 = 0;
  found = IFT_NIL;

  /* If include_center is True then we add the displacement (0,0,0) to the adjacency relation even if its magnitude
   * is not within the closed interval [radius_min, radius_max]
   */
  if(include_center) {
    /* If radius_min is 0.0 then we know that the center will be added eventually and need not do that now */
    if(!iftAlmostZero(radius_min)) {
      /* For the sake of generality, we search for the displacement (0,0,0) and add it first if include_center
         is True. This way we need not assume that (0,0,0) is in the first position
       */
      for(i = 0; i < A->n && !(A->dx[i] == A->dy[i] && A->dy[i] == A->dz[i] && A->dx[i] == 0); i++);

      // If the displacement (0,0,0) is found, we add it in the first position of the resulting adjacency relation
      if(i < A->n) {
        found = i;
        n++;
        i0++;
      }
    }
  }

  Aout = iftCreateAdjRel(n);

  /* If found >= 0, then we know that include_center is True, radius_min is != 0.0, and the displacement (0,0,0)
    was encountered. Hence, we add it in the first position because it will not be added later.
   */
  if(found >= 0) {
    Aout->dx[i0] = A->dx[found];
    Aout->dy[i0] = A->dy[found];
    Aout->dz[i0] = A->dz[found];
    i0++;
  }

  /* Adding the voxel displacements within the specified interval [radius_min, radius_max] */
  for(i = 0; i < A->n; i++) {
    r = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];

    if(r >= radius_min*radius_min && r <= radius_max*radius_max) {
      Aout->dx[i0] = A->dx[i];
      Aout->dy[i0] = A->dy[i];
      Aout->dz[i0] = A->dz[i];

      i0++;
    }
  }

  /* Recycling <found> to ensure that (0,0,0) be placed in the first index */
  found = IFT_NIL;
  /* Searching for the displacement (0,0,0) in <Aout> to place it in the first index in case that was not done before */
  for(i0 = 0; i0 < Aout->n && found < 0; i0++) {
    if (Aout->dx[i0] == Aout->dy[i0] && Aout->dy[i0] == Aout->dz[i0] && Aout->dx[i0] == 0)
      found = i0;
  }

  /* Shift to right and place central voxel at first, just in case that was not done before*/
  for (i = found; i > 0; i--) {
    dx = Aout->dx[i];
    dy = Aout->dy[i];
    dz = Aout->dz[i];
    Aout->dx[i] = Aout->dx[i - 1];
    Aout->dy[i] = Aout->dy[i - 1];
    Aout->dz[i] = Aout->dz[i - 1];
    Aout->dx[i - 1] = dx;
    Aout->dy[i - 1] = dy;
    Aout->dz[i - 1] = dz;
  }

  return Aout;
}


iftAdjRel *iftAdjacencyBoundaries(const iftAdjRel *A, const iftAdjRel *B_in) {
    bool is_3D_adj = false;
    for (int i = 0; i < A->n; i++) {
        if (A->dz[i] != 0) {
            is_3D_adj = true;
            break;
        }
    }

    iftAdjRel *B = NULL;
    if (B_in)
        B = iftCopyAdjacency(B_in);
    else
        B = (is_3D_adj) ? iftSpheric(1.0) : iftCircular(1.0);

    // get the maximum distance between the center of A and the other displacement vectors
    float max_dist = sqrtf(A->dx[A->n-1]*A->dx[A->n-1] + A->dy[A->n-1]*A->dy[A->n-1] + A->dz[A->n-1]*A->dz[A->n-1]);

    int xsize, ysize, zsize;
    xsize = ysize = zsize = (3 * ceil(max_dist)); // find out a given size to create an image that fits A
    if (!is_3D_adj)
        zsize = 1;

    // create "an image" and centralize A on it
    int ***M = iftAlloc(zsize, sizeof(int**));
    for (int z = 0; z < zsize; z++) {
        M[z] = iftAlloc(ysize, sizeof(int*));
        
        for (int y = 0; y < ysize; y++)
            M[z][y] = iftAlloc(xsize, sizeof(int));
    }

    iftVoxel center = {xsize/2, ysize/2, zsize/2};

    // label the adjacent on the "image"
    for (int i = 0; i < A->n; i++) {
        iftVoxel u;
        u.x = center.x + A->dx[i];
        u.y = center.y + A->dy[i];
        u.z = center.z + A->dz[i];

        M[u.z][u.y][u.x] = 1;
    }
    

    iftList *boundaries = iftCreateList();

    for (int i = 0; i < A->n; i++) {
        iftVoxel u;
        u.x = center.x + A->dx[i];
        u.y = center.y + A->dy[i];
        u.z = center.z + A->dz[i];

        for (int j = 1; j < B->n; j++) {
            // adjacent voxel of v
            iftVoxel v;
            v.x = u.x + B->dx[j];
            v.y = u.y + B->dy[j];
            v.z = u.z + B->dz[j];

            // if v is out of domain or it has label != of u
            if ((v.x < 0 || v.x >= xsize) || (v.y < 0 || v.y >= ysize) || (v.z < 0 || v.z >= zsize) ||
                (M[u.z][u.y][u.x] != M[v.z][v.y][v.x])) {
                iftInsertListIntoTail(boundaries, i);
                break;
            }
        }
    }

    iftAdjRel *Abound = iftCreateAdjRel(boundaries->n);

    int i = 0;
    while (!iftIsEmptyList(boundaries)) {
        int j = iftRemoveListTail(boundaries);
        Abound->dx[i] = A->dx[j];
        Abound->dy[i] = A->dy[j];
        Abound->dz[i] = A->dz[j];
        i++;
    }
    iftDestroyList(&boundaries);


    for (int z = 0; z < zsize; z++) {
        for (int y = 0; y < ysize; y++)
            iftFree(M[z][y]);
        iftFree(M[z]);
    }        
    iftFree(M);

    if (B_in)
        iftDestroyAdjRel(&B);

    return Abound;
}


int ift6NeighborsClockFromLeft(iftVoxel *pred, iftVoxel *cur)
{
    int map[3][3];
    map[0][0] = 7; // -1, -1
    map[1][0] = 8; //  0, -1
    map[2][0] = 1; //  1, -1
    map[2][1] = 2; //  1,  0
    map[2][2] = 3; //  1,  1
    map[1][2] = 4; //  0,  1
    map[0][2] = 5; // -1,  1
    map[0][1] = 6; // -1,  0
    map[1][1] = 0; //  0,  0
    int dx = cur->x - pred->x;
    int dy = cur->y - pred->y;
    return map[dx + 1][dy + 1];
}


void iftSidePixels(const iftVoxel *u, const iftVoxel *v, float radii, iftVoxel *left, iftVoxel *right)
{
    int dx = v->x - u->x;
    int dy = v->y - u->y;
    float d = sqrt(dx * dx + dy * dy);
    if (d > IFT_EPSILON) {
        right->x = u->x + (int) roundf(dx / 2.0f - radii * dy / d);
        right->y = u->y + (int) roundf(dx / d + radii * dy / 2.0f);
        right->z = 0;
        left->x = u->x + (int) roundf(dx / 2.0f + radii * dy / d);
        left->y = u->y + (int) roundf(dy / 2.0f - radii * dx / d);
        left->z = 0;
    } else {
        right->x = left->x = u->x;
        right->y = left->y = u->y;
        right->z = left->z = 0;
    }
}

static inline int m_jitter(int value, int mod)
{
    if (value < -IFT_EPSILON) {
        return -1 * (random() % mod);
    } else if (value > IFT_EPSILON) {
        return random() % mod;
    }
    return 0;
}

void iftSidePixelsWithJitter(const iftVoxel *u, const iftVoxel *v, float radii,
                             int mod, iftVoxel *left, iftVoxel *right)
{
    int dx = v->x - u->x;
    int dy = v->y - u->y;
    float d = sqrt(dx * dx + dy * dy);
    if (d > IFT_EPSILON) {
        right->x = u->x + m_jitter((int) roundf(dx / 2.0f - radii * dy / d), mod);
        right->y = u->y + m_jitter((int) roundf(dx / d + radii * dy / 2.0f), mod);
        right->z = 0;
        left->x = u->x + m_jitter((int) roundf(dx / 2.0f + radii * dy / d), mod);
        left->y = u->y + m_jitter((int) roundf(dy / 2.0f - radii * dx / d), mod);
        left->z = 0;
    } else {
        right->x = left->x = u->x;
        right->y = left->y = u->y;
        right->z = left->z = 0;
    }
}





