#include "iftSpectrum.h"



iftMImage *iftFFTImageToMImage(iftImage *img)
{
  iftMImage *mimg = iftCreateMImage(img->xsize, img->ysize, img->zsize, 2);

  for (int p=0; p < img->n; p++){
    mimg->val[p][0] = img->val[p];
    mimg->val[p][1] = 0.0;
  }

  return(mimg);
}

iftImage  *iftFFTMImageToImage(iftMImage *mimg, iftImage *orig)
{
  iftImage *img = iftCreateImageFromImage(orig);
  iftVoxel  u;
  float     min=IFT_INFINITY_FLT;
  float     max=IFT_INFINITY_FLT_NEG;
  
  for (u.z=0; u.z < mimg->zsize; u.z++) 
    for (u.y=0; u.y < mimg->ysize; u.y++) 
      for (u.x=0; u.x < mimg->xsize; u.x++) {
	int p = iftMGetVoxelIndex(mimg,u);
	float mag = sqrtf(mimg->val[p][0]*mimg->val[p][0] + 
			  mimg->val[p][1]*mimg->val[p][1]);
	if (mag < min)
	  min = mag;
	if (mag > max)
	  max = mag;
      }
  
  int Imax = iftNormalizationValue(iftImageDepth(orig));
  for (u.z=0; u.z < orig->zsize; u.z++) 
    for (u.y=0; u.y < orig->ysize; u.y++) 
      for (u.x=0; u.x < orig->xsize; u.x++) {
	int p = iftMGetVoxelIndex(mimg,u);
	float mag = sqrtf(mimg->val[p][0]*mimg->val[p][0] + 
			  mimg->val[p][1]*mimg->val[p][1]);
	int q       = iftGetVoxelIndex(orig,u);
	img->val[q] = Imax * (mag - min) / (max - min);
      }

  return(img);
}

void iftFFT1D(int dir, long nn, double *x, double *y)
{
  
/*
   This computes an in-place complex-to-complex FFT x and y are the
   real and imaginary arrays of nn points.  dir = 1 gives forward
   transform dir = -1 gives reverse transform
*/

   int m;
   long i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   m = (int)(log(nn)/log(2)+.00001);

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 =  0.0;
   l2 =  1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for reverse transform */
   if (dir == -1) {
      for (i=0;i<nn;i++) {
         x[i] /= (double)nn;
         y[i] /= (double)nn;
      }
   }

}

iftMImage *iftFFT2D(iftMImage *mimg)
{
  iftMImage *spec;
  int        p, xsize, ysize;
  double    *real, *imag;
  iftVoxel   u;
  
  /* Enlarge xsize and ysize to the nearest power of 2 */

  xsize = (int)ceil(log(mimg->xsize)/log(2));
  ysize = (int)ceil(log(mimg->ysize)/log(2));
  xsize = 1 << xsize;
  ysize = 1 << ysize;
  
  spec = iftCreateMImage(xsize, ysize, 1, 2);

  real = (double*)iftAllocDoubleArray(spec->xsize);
  imag = (double*)iftAllocDoubleArray(spec->xsize);

  u.z = 0;
  
  for (u.y=0; u.y<spec->ysize; u.y++)
    {
      for (u.x=0; u.x<spec->xsize; u.x++)
	{
	  if (u.y<mimg->ysize && u.x<mimg->xsize){
	    p         = iftMGetVoxelIndex(mimg,u);
	    real[u.x] = (double)mimg->val[p][0];
	    imag[u.x] = (double)mimg->val[p][1];
	  }else{
	    real[u.x] = 0.0;
	    imag[u.x] = 0.0;
	  }
	}
      iftFFT1D(1, (long)spec->xsize, real, imag);
      for (u.x=0; u.x<spec->xsize; u.x++)
	{
	  p   = iftMGetVoxelIndex(spec,u);
	  spec->val[p][0] = real[u.x];
	  spec->val[p][1] = imag[u.x];
	}
    }
  
  iftFree(real);
  iftFree(imag);
  
  real = iftAllocDoubleArray(spec->ysize);
  imag = iftAllocDoubleArray(spec->ysize);
  
  for (u.x=0; u.x<spec->xsize; u.x++)
    {
      for (u.y=0; u.y<spec->ysize; u.y++)
	{
 	  p         = iftMGetVoxelIndex(spec,u); 
	  real[u.y] = spec->val[p][0];
	  imag[u.y] = spec->val[p][1];
	}
      iftFFT1D(1, (long)spec->ysize, real, imag);
      for (u.y=0; u.y<spec->ysize; u.y++)
	{
	  p         = iftMGetVoxelIndex(spec,u); 
	  spec->val[p][0] = real[u.y];
	  spec->val[p][1] = imag[u.y];
	}
    }
  
  iftFree(real);
  iftFree(imag);
  
  return(spec);
}

iftMImage *iftFFT3D(iftMImage *mimg)
{
  iftMImage *spec;
  int        xsize, ysize, zsize, p;
  double    *real, *imag;
  iftVoxel   u;

  /* Enlarge xsize, ysize and zsize to the nearest power of 2 */

  xsize   = (int)ceil(log(mimg->xsize)/log(2));
  ysize   = (int)ceil(log(mimg->ysize)/log(2));
  zsize   = (int)ceil(log(mimg->zsize)/log(2));
  xsize   =  1 << xsize;
  ysize   =  1 << ysize;
  zsize   =  1 << zsize;

  spec = iftCreateMImage(xsize, ysize, zsize, 2);

  for(u.z=0; u.z < spec->zsize; u.z++){
    
    real = iftAllocDoubleArray(spec->xsize);
    imag = iftAllocDoubleArray(spec->xsize);
    
    for (u.y=0; u.y<spec->ysize; u.y++)
      {
	for (u.x=0; u.x<spec->xsize; u.x++)
	  {
	    if (u.y<mimg->ysize && u.x<mimg->xsize && u.z<mimg->zsize){
	      p = iftMGetVoxelIndex(mimg,u);
	      real[u.x] = (double)mimg->val[p][0];
	      imag[u.x] = (double)mimg->val[p][1];
	    }else{
	      real[u.x] = 0.0;
	      imag[u.x] = 0.0;
	    }
	  }
	iftFFT1D(1, (long)spec->xsize, real, imag);
	for (u.x=0; u.x<spec->xsize; u.x++)
	  {
	    p = iftMGetVoxelIndex(spec,u);
	    spec->val[p][0] = real[u.x];
	    spec->val[p][1] = imag[u.x];
	  }
      }
      
    iftFree(real);
    iftFree(imag);

    real = iftAllocDoubleArray(spec->ysize);
    imag = iftAllocDoubleArray(spec->ysize);
    
    for (u.x=0; u.x<spec->xsize; u.x++)
      {
	for (u.y=0; u.y<spec->ysize; u.y++)
	  {
	    p = iftMGetVoxelIndex(spec,u);
	    real[u.y] = spec->val[p][0];
	    imag[u.y] = spec->val[p][1];
	  }
	iftFFT1D(1, (long)spec->ysize, real, imag);
	for (u.y=0; u.y<spec->ysize; u.y++)
	  {
	    p = iftMGetVoxelIndex(spec,u);
	    spec->val[p][0] = real[u.y];
	    spec->val[p][1] = imag[u.y];
	  }
      }
    
    iftFree(real);
    iftFree(imag);    
  }
  
  real = iftAllocDoubleArray(spec->zsize);
  imag = iftAllocDoubleArray(spec->zsize);

  for (u.y=0; u.y<spec->ysize; u.y++)
    {
      for (u.x=0; u.x<spec->xsize; u.x++)
	{
	  for(u.z=0; u.z<spec->zsize;u.z++)
	    {
	      p = iftMGetVoxelIndex(spec,u);
	      real[u.z]= spec->val[p][0];
	      imag[u.z]= spec->val[p][1];
	    }
	  iftFFT1D(1, (long)spec->zsize, real, imag);
	  for (u.z=0; u.z<spec->zsize; u.z++)
	    {
	      p = iftMGetVoxelIndex(spec,u);
	      spec->val[p][0] = real[u.z];
	      spec->val[p][1] = imag[u.z];
	    }
	}
    }
  
  iftFree(real);
  iftFree(imag);
  
  return(spec);
}

iftMImage *iftInvFFT2D(iftMImage *spec)
{
  iftMImage  *mimg;
  double     *real, *imag;
  iftMImage  *specaux;
  iftVoxel    u;
  int         p;
  
  mimg     = iftCreateMImage(spec->xsize, spec->ysize, 1, 2);
  specaux  = iftCreateMImage(spec->xsize, spec->ysize, 1, 2);
  
  real    = iftAllocDoubleArray(spec->xsize);
  imag    = iftAllocDoubleArray(spec->xsize);

  u.z = 0;
  
  for (u.y=0; u.y<spec->ysize; u.y++)
    {
      for (u.x=0; u.x<spec->xsize; u.x++)
	{
	  p         = iftMGetVoxelIndex(spec,u);
	  real[u.x] = spec->val[p][0];
	  imag[u.x] = spec->val[p][1];
	}
      iftFFT1D(-1, (long)spec->xsize, real, imag);
      for (u.x=0; u.x<spec->xsize; u.x++)
	{
	  p = iftMGetVoxelIndex(spec,u);
	  specaux->val[p][0] = real[u.x];
	  specaux->val[p][1] = imag[u.x];
	}
      }
  
    iftFree(real);
    iftFree(imag);
    
    real = iftAllocDoubleArray(specaux->ysize);
    imag = iftAllocDoubleArray(specaux->ysize);

    for (u.x=0; u.x<mimg->xsize; u.x++){
      for (u.y=0; u.y<mimg->ysize; u.y++){
	p = iftMGetVoxelIndex(spec,u);
	real[u.y] = specaux->val[p][0];
	imag[u.y] = specaux->val[p][1];
      }
      iftFFT1D(-1, (long)specaux->ysize, real, imag);
      for (u.y=0; u.y<mimg->ysize; u.y++){
	p = iftGetVoxelIndex(specaux,u);
	mimg->val[p][0] = real[u.y];
	mimg->val[p][1] = imag[u.y];
      }
    }
 
    iftFree(real);
    iftFree(imag);
    
    iftDestroyMImage(&specaux);
    
    return(mimg);
}

iftMImage *iftInvFFT3D(iftMImage *spec)
{
  iftMImage  *mimg;
  int         p;
  double     *real, *imag;
  iftMImage  *specaux;
  iftVoxel    u;
  
  mimg     = iftCreateMImage(spec->xsize, spec->ysize, spec->zsize, 2);
  specaux  = iftCreateMImage(spec->xsize, spec->ysize, spec->zsize, 2); 
    
  for(u.z=0;u.z<spec->zsize;u.z++){

    real = iftAllocDoubleArray(spec->xsize);
    imag = iftAllocDoubleArray(spec->xsize);
    
    for (u.y=0; u.y<spec->ysize; u.y++)
      {
	for (u.x=0; u.x<spec->xsize; u.x++)
	  {
	    p         = iftMGetVoxelIndex(spec,u);
	    real[u.x] = spec->val[p][0];
	    imag[u.x] = spec->val[p][1];
	  }
	iftFFT1D(-1, (long)spec->xsize, real, imag);
	for (u.x=0; u.x<spec->xsize; u.x++)
	  {
	    p         = iftMGetVoxelIndex(spec,u);	    
	    specaux->val[p][0] = real[u.x];
	    specaux->val[p][1] = imag[u.x];
	  }
      }
    
    iftFree(real);
    iftFree(imag);
    
    real = iftAllocDoubleArray(specaux->ysize);
    imag = iftAllocDoubleArray(specaux->ysize);
    
    for (u.x=0; u.x<specaux->xsize; u.x++)
      {
	for (u.y=0; u.y<specaux->ysize; u.y++)
	  {
	    p         = iftMGetVoxelIndex(spec,u);
	    real[u.y] = specaux->val[p][0];
	    imag[u.y] = specaux->val[p][1];
	  }
	iftFFT1D(-1, (long)specaux->ysize, real, imag);
	for (u.y=0; u.y<specaux->ysize; u.y++)
	  {
	    p         = iftMGetVoxelIndex(spec,u);
	    specaux->val[p][0] = real[u.y];
	    specaux->val[p][1] = imag[u.y];
	  }
      }
    
    iftFree(real);
    iftFree(imag);
  }
  
  real = iftAllocDoubleArray(spec->zsize);
  imag = iftAllocDoubleArray(spec->zsize);
  
  for (u.y=0; u.y<spec->ysize; u.y++){
      for (u.x=0; u.x<spec->xsize; u.x++)
	{
	  for(u.z=0; u.z<spec->zsize;u.z++)
	    {
	      p = iftMGetVoxelIndex(spec,u);
	      real[u.z]= specaux->val[p][0];
	      imag[u.z]= specaux->val[p][1];
	    }
	  iftFFT1D(-1, (long)spec->zsize, real, imag);
	  for (u.z=0; u.z<spec->zsize; u.z++)
	    {
	      p = iftGetVoxelIndex(mimg,u);
	      mimg->val[p][0] = real[u.z];
	      mimg->val[p][1] = imag[u.z];
	    }
	}
  }

  iftFree(real);
  iftFree(imag);
  iftDestroyMImage(&specaux);
  
  return(mimg);
}
iftImage *iftViewMagnitude(iftMImage *spec)
{
  iftImage *img;
  double *aux; // could a aux iftImage be used?
  int p,q;
  iftVoxel u,v;
  double max;

    
  max = IFT_INFINITY_DBL_NEG;
  img = iftCreateImage(spec->xsize, spec->ysize, spec->zsize);
  aux = iftAllocDoubleArray(img->n);

  // Calculate de magnitude sqrt(real² + imag²) and centralize the image
  if (iftIs3DMImage(spec)){
    for (u.z=0; u.z < spec->zsize; u.z++){
      for (u.y=0; u.y < spec->ysize; u.y++){
        for (u.x=0; u.x < spec->xsize; u.x++){
          p = iftMGetVoxelIndex(spec,u); // original image index (centralized)
          // Getting the transformed image coordinates
          v.x = (u.x + img->xsize/2)%img->xsize;
          v.y = (u.y + img->ysize/2)%img->ysize;
          v.z = (u.z + img->zsize/2)%img->zsize;
          q = iftMGetVoxelIndex(spec, v); // trnasformed image index

          aux[q] = sqrt(spec->val[p][0]*spec->val[p][0] + spec->val[p][1]*spec->val[p][1]);
          if (aux[q] > max)
            max = aux[q];
        }
      }
    }
  }
  else { // 2D
    u.z = v.z = 0;
    for (u.y=0; u.y < spec->ysize; u.y++){
      for (u.x=0; u.x < spec->xsize; u.x++){

        p = iftMGetVoxelIndex(spec,u); // original image index (centralized)
        // Getting the transformed image coordinates
        v.x = (u.x + img->xsize/2)%img->xsize;
        v.y = (u.y + img->ysize/2)%img->ysize;
        v.z = 0;
        q = iftMGetVoxelIndex(spec, v); // trnasformed image index

        aux[q] = sqrt(spec->val[p][0]*spec->val[p][0] + spec->val[p][1]*spec->val[p][1]);
        if (aux[q] > max)
          max = aux[q];
      }
    }
  }

  // Normalize the image to 255 colors
  max = 255.0/max;
  for (p=0; p<img->n; p++)
    img->val[p] = (int)(max*aux[p]);

  iftFree(aux);

  return(img);
}

iftImage *iftViewLogMagnitude(iftMImage *spec)
{
  iftImage  *img;
  double    *aux; 
  int        p, q;
  iftVoxel   u, v;
  double     max;

  max = IFT_INFINITY_DBL_NEG;
  img = iftCreateImage(spec->xsize, spec->ysize, spec->zsize);
  aux = iftAllocDoubleArray(img->n);
  
  /* Calculate the magnitude log(1+sqrt(real^2+imag^2)) and centralize
     the image */

  if (iftIs3DMImage(spec)) {
    for (u.z=0; u.z < spec->zsize; u.z++){
      for (u.y=0; u.y < spec->ysize; u.y++){
	for (u.x=0; u.x < spec->xsize; u.x++){
	  p   = iftMGetVoxelIndex(spec,u);
	  v.x = (u.x + img->xsize/2)%img->xsize;
	  v.y = (u.y + img->ysize/2)%img->ysize;
	  v.z = (u.z + img->zsize/2)%img->zsize;
	  q   = iftMGetVoxelIndex(spec,v);
	  aux[q] = log(1+sqrt(spec->val[p][0]*spec->val[p][0] + spec->val[p][1]*spec->val[p][1]));
	  if (aux[q] > max)
	    max = aux[q];
	}
      }
    }
  } else { // 2D 
    u.z = 0;
    for (u.y=0; u.y < spec->ysize; u.y++){
      for (u.x=0; u.x < spec->xsize; u.x++){
	p   = iftMGetVoxelIndex(spec,u);
	v.x = (u.x + img->xsize/2)%img->xsize;
	v.y = (u.y + img->ysize/2)%img->ysize;
	v.z = 0;
	q   = iftMGetVoxelIndex(spec,v);
	aux[q] = log(1+sqrt(spec->val[p][0]*spec->val[p][0] + spec->val[p][1]*spec->val[p][1]));
	if (aux[q] > max)
	  max = aux[q];
      }
    }
  }
    
  /* Normalize the image to 255 colors */
  max = 255.0/max;
  for (p=0; p<img->n; p++)
    img->val[p] = (int)(max*aux[p]);

  iftFree(aux);

  return(img);
}

iftImage *iftViewPhase(iftMImage *spec)
{
  iftImage  *img;
  double    *aux; 
  int        p, q;
  iftVoxel   u, v;
  double     max, min;

  max = IFT_INFINITY_DBL_NEG;
  min = IFT_INFINITY_DBL;
  img = iftCreateImage(spec->xsize, spec->ysize, spec->zsize);
  aux = iftAllocDoubleArray(img->n);
  
  /* Calculate the phase atan(imag/real) and centralize
     the image */

  if (iftIs3DMImage(spec)) {
    for (u.z=0; u.z < spec->zsize; u.z++){
      for (u.y=0; u.y < spec->ysize; u.y++){
	for (u.x=0; u.x < spec->xsize; u.x++){
	  p   = iftMGetVoxelIndex(spec,u);
	  v.x = (u.x + img->xsize/2)%img->xsize;
	  v.y = (u.y + img->ysize/2)%img->ysize;
	  v.z = (u.z + img->zsize/2)%img->zsize;
	  q   = iftMGetVoxelIndex(spec,v);
	  aux[q] = atan(spec->val[p][1]/spec->val[p][0]);
	  if (aux[q] > max)
	    max = aux[q];
	  if (aux[q] < min)
	    min = aux[q];
	}
      }
    }
  } else { // 2D 
    u.z = 0;
    for (u.y=0; u.y < spec->ysize; u.y++){
      for (u.x=0; u.x < spec->xsize; u.x++){
	p   = iftMGetVoxelIndex(spec,u);
	v.x = (u.x + img->xsize/2)%img->xsize;
	v.y = (u.y + img->ysize/2)%img->ysize;
	v.z = 0;
	q   = iftMGetVoxelIndex(spec,v);
	aux[q] = atan(spec->val[p][1]/spec->val[p][0]);
	if (aux[q] > max)
	  max = aux[q];
	if (aux[q] < min)
	  min = aux[q];
      }
    }
  }
    
  /* Normalize the image to 255 colors */
  max = 255.0/(max-min);
  for (p=0; p<img->n; p++)
    img->val[p] = (int)(max*(aux[p]-min));

  iftFree(aux);

  return(img);
}

iftMImage *iftMultSpectra(iftMImage *spec, iftMImage *filt)
{
  iftMImage *fspec;
  
  fspec = iftCreateMImage(spec->xsize, spec->ysize, spec->zsize, 2);
    
  for (int p=0; p < spec->n; p++) 
    for (int b=0; b < spec->m; b += 2) {
      fspec->val[p][0] = (spec->val[p][0]*filt->val[p][0]) -
	(spec->val[p][1]*filt->val[p][1]);
      fspec->val[p][1] = (spec->val[p][0]*filt->val[p][1]) +
	(spec->val[p][1]*filt->val[p][0]);
    }
  
  return(fspec);
}

iftMImage *iftGaborFilter2D(iftMImage *spec, float u0, float v0, float P, float x0, float y0, float a, float b, float theta) {
    /*
    Gabor filter implementation in the space domain.
    Follow the reference for more: https://inc.ucsd.edu/mplab/75/media//gabor.pdf

    Parameters:
    float (u0,v0): sinusoid carrier space frequency along x and y, respectively;
    float P: phase of the sinusoid carrier
    float (x0,y0): location of the peak of the Gaussian envelope (should be filt->xsize/2 and filt->ysize/2 if you want to centralize it)
    float (a,b): scales the 2 axis of the Gaussian envelope
    float theta: rotation angle of the gaussian envelope

    gabor(x,y, parameters) = sinusoid(x,y).envelope(x,y)
    */
    // Create filter templates with same size and nbands as given spec image
    iftMImage *filt = iftCreateMImage(spec->xsize, spec->ysize, spec->zsize, 2); // filt to move centered filter to (0,0)
    iftVoxel u; // map coordinates to indexes
    int p;

    for (u.z = 0; u.z < filt->zsize; u.z++) {
        for (u.y = 0; u.y < filt->ysize; u.y++) {
            for (u.x = 0; u.x < filt->xsize; u.x++) {
                // calculating the gaussian envelope
                float x_rotation = ((u.x - spec->xsize/2) - x0) * cos(theta) +
                                   ((u.y - spec->ysize/2) - y0) * sin(theta); // (x - x0)r = (x - x0)cos(theta) + (y - y0)sen(theta)
                float y_rotation = -((u.x - spec->xsize/2) - x0) * sin(theta) +
                                   ((u.y - spec->ysize/2) - y0) * cos(theta); // (y - y0)r = -(x - x0)sin(theta) + (y - y0)cos(theta)
                float envelope = exp(-IFT_PI * (pow(a, 2.0) * pow(x_rotation, 2.0) + pow(b, 2.0) * pow(y_rotation, 2.0)));

                // calculating the sinusoid carrier (real and imaginary parts)
                //printf("Spectrum\n");
                //printf("u0:%f\t(u.x - spec->xsize/2):%d\tv0:%f\t(u.y - spec->xsize/2):%d\tResultado do produto:%f \n",u0, (u.x - spec->xsize/2), v0, (u.y - spec->xsize/2), u0 * (u.x - spec->xsize/2) + v0 * (u.y - spec->xsize/2));
                float real_carrier = cos(2.0 * IFT_PI * (u0 * (u.x - spec->xsize/2) + v0 * (u.y - spec->ysize/2)) + P);
                float imag_carrier = sin(2.0 * IFT_PI * (u0 * (u.x - spec->xsize/2) + v0 * (u.y - spec->ysize/2)) + P);

                // finally, calculating real and imag parts of gabor itself
                float real_gabor = envelope * real_carrier;
                float imag_gabor = envelope * imag_carrier;

                // getting coordinate index and filling the bands
                p = iftMGetVoxelIndex(filt, u);
                filt->val[p][0] = real_gabor;
                filt->val[p][1] = imag_gabor;
            }
        }
    }
    return (filt);
}



iftMImage *iftButterworthHighPass(iftMImage *spec, int n, float Dh)
{
  iftMImage *filt = iftCreateMImage(spec->xsize, spec->ysize, spec->zsize,2);
  iftMImage *aux  = iftCreateMImage(spec->xsize, spec->ysize, spec->zsize,2);
  iftVoxel   u, v;
  int        p, q;
  
  /* Create the filter with origin at the center of the image domain */

  Dh = Dh * sqrt((filt->xsize/2*filt->xsize/2) +
		 (filt->ysize/2*filt->ysize/2) +
		 (filt->zsize/2*filt->zsize/2));
  n  = 2*n;
  for (u.z=0; u.z < filt->zsize; u.z++) {
    for (u.y=0; u.y < filt->ysize; u.y++) {
      for (u.x=0; u.x < filt->xsize; u.x++) {
	float D = sqrt((u.x - filt->xsize/2)*(u.x - filt->xsize/2) + 
		       (u.y - filt->ysize/2)*(u.y - filt->ysize/2) +
		       (u.z - filt->zsize/2)*(u.z - filt->zsize/2));
	p = iftMGetVoxelIndex(filt, u);
	if (!iftAlmostZero(D)){
	  aux->val[p][0] = 1.0 / (1.0 + 0.414 * pow((Dh / D),n));
	} 	  
      }
    }
  }

  /* Move the origin to the (0,0) coordinate considering that the
     spectrum is periodic with period 2*pi */

  if (iftIs3DMImage(filt)){
    for (u.z=0; u.z < filt->zsize; u.z++) {
      for (u.y=0; u.y < filt->ysize; u.y++) {
	for (u.x=0; u.x < filt->xsize; u.x++) {
	  v.x = (filt->xsize + (u.x - filt->xsize/2))%filt->xsize;
	  v.y = (filt->ysize + (u.y - filt->ysize/2))%filt->ysize;
	  v.z = (filt->zsize + (u.z - filt->zsize/2))%filt->zsize;
	  p   = iftMGetVoxelIndex(filt,u);
	  q   = iftMGetVoxelIndex(filt,v);
	  filt->val[q][0] = aux->val[p][0];
	}
      }
    }
  } else {
    u.z = v.z = 0;
    for (u.y=0; u.y < filt->ysize; u.y++) {
      for (u.x=0; u.x < filt->xsize; u.x++) {
	v.x = (filt->xsize + (u.x - filt->xsize/2))%filt->xsize;
	v.y = (filt->ysize + (u.y - filt->ysize/2))%filt->ysize;
	p   = iftMGetVoxelIndex(filt,u);
	q   = iftMGetVoxelIndex(filt,v);
	filt->val[q][0] = aux->val[p][0];
      }
    }
  }

  iftDestroyMImage(&aux);
  return(filt);
}

iftMImage *iftButterworthLowPass(iftMImage *spec, int n, float Dl)
{
  iftMImage *filt = iftCreateMImage(spec->xsize, spec->ysize, spec->zsize,2);
  iftMImage *aux  = iftCreateMImage(spec->xsize, spec->ysize, spec->zsize,2);
  iftVoxel   u, v;
  int        p, q;
  
  /* Create the filter with origin at the center of the image domain */

  Dl = Dl * sqrt((filt->xsize/2*filt->xsize/2) +
		 (filt->ysize/2*filt->ysize/2) +
		 (filt->zsize/2*filt->zsize/2));
  n  = 2*n;
  for (u.z=0; u.z < filt->zsize; u.z++) {
    for (u.y=0; u.y < filt->ysize; u.y++) {
      for (u.x=0; u.x < filt->xsize; u.x++) {
	float D = sqrt((u.x - filt->xsize/2)*(u.x - filt->xsize/2) + 
		       (u.y - filt->ysize/2)*(u.y - filt->ysize/2) +
		       (u.z - filt->zsize/2)*(u.z - filt->zsize/2));
	p = iftMGetVoxelIndex(filt, u);
	if (!iftAlmostZero(D)){
	  aux->val[p][0] = 1.0 / (1.0 + 0.414 * pow((D / Dl),n));
	} 	  
      }
    }
  }

  /* Move the origin to the (0,0) coordinate considering that the
     spectrum is periodic with period 2*pi */

  if (iftIs3DMImage(filt)){
    for (u.z=0; u.z < filt->zsize; u.z++) {
      for (u.y=0; u.y < filt->ysize; u.y++) {
	for (u.x=0; u.x < filt->xsize; u.x++) {
	  v.x = (filt->xsize + (u.x - filt->xsize/2))%filt->xsize;
	  v.y = (filt->ysize + (u.y - filt->ysize/2))%filt->ysize;
	  v.z = (filt->zsize + (u.z - filt->zsize/2))%filt->zsize;
	  p   = iftMGetVoxelIndex(filt,u);
	  q   = iftMGetVoxelIndex(filt,v);
	  filt->val[q][0] = aux->val[p][0];
	}
      }
    }
  } else {
    u.z = v.z = 0;
    for (u.y=0; u.y < filt->ysize; u.y++) {
      for (u.x=0; u.x < filt->xsize; u.x++) {
	v.x = (filt->xsize + (u.x - filt->xsize/2))%filt->xsize;
	v.y = (filt->ysize + (u.y - filt->ysize/2))%filt->ysize;
	p   = iftMGetVoxelIndex(filt,u);
	q   = iftMGetVoxelIndex(filt,v);
	filt->val[q][0] = aux->val[p][0];
      }
    }
  }

  iftDestroyMImage(&aux);
  return(filt);
}

iftMImage *iftButterworthBandPass(iftMImage *spec, int n, float Dl, float Dh)
{
  if (Dl > Dh) {
    
    iftMImage *lowpass  = iftButterworthLowPass(spec,n,Dl);
    iftMImage *highpass = iftButterworthHighPass(spec,n,Dh);
    iftMImage *filt = iftCreateMImage(spec->xsize,spec->ysize,spec->zsize,2);
   
    for (int p=0; p < spec->n; p++) {
      filt->val[p][0] = lowpass->val[p][0]*highpass->val[p][0];
    }

    iftDestroyMImage(&lowpass);
    iftDestroyMImage(&highpass);
    return(filt);
    
  } else {
    iftError("Bandpass filter requires Dl > Dh","iftButterworthBandPass");
  }

  return(NULL);
}

iftMImage *iftButterworthBandReject(iftMImage *spec, int n, float Dl, float Dh)
{
  if (Dl < Dh) {
    
    iftMImage *lowpass  = iftButterworthLowPass(spec,n,Dl);
    iftMImage *highpass = iftButterworthHighPass(spec,n,Dh);
    iftMImage *filt = iftCreateMImage(spec->xsize,spec->ysize,spec->zsize,2);
   
    for (int p=0; p < spec->n; p++) {
      filt->val[p][0] = lowpass->val[p][0] + highpass->val[p][0];
    }

    iftDestroyMImage(&lowpass);
    iftDestroyMImage(&highpass);
    return(filt);
    
  } else {
    iftError("Bandreject filter requires Dl < Dh","iftButterworthBandReject");
  }
    
  return(NULL);
}
