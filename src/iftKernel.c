#include "iftKernel.h"


#include "ift/core/io/Stream.h"


iftKernel *iftReadKernel(char *filename)
{
  FILE      *fp=fopen(filename,"r");
  iftKernel *K=NULL;
  iftAdjRel *A=NULL;
  int        i,n;

  if(fscanf(fp,"%d",&n)!=1) iftError("Reading error", "iftReadKernel");

  A = iftCreateAdjRel(n);
  K = iftCreateKernel(A);

  /* Read displacements and weights of the adjacent voxels */

  for (i=0; i < K->A->n; i++) {
    if(fscanf(fp,"%d %d %d %f",&K->A->dx[i],&K->A->dy[i],&K->A->dz[i],&K->weight[i])!=4)
        iftError("Reading error", "iftReadKernel");
  }
  
  fclose(fp);

  iftDestroyAdjRel(&A); 

  return(K);
}

void iftWriteKernel(iftKernel *K, char *filename)
{
  FILE   *fp=fopen(filename,"w");
  int     i;

  fprintf(fp,"%d\n",K->A->n);
  
  /* Write displacements and weights of the adjacent voxels */

  for (i=0; i < K->A->n; i++) {
    fprintf(fp,"%d %d %d %f\n",K->A->dx[i],K->A->dy[i],K->A->dz[i],K->weight[i]); 
  }
  
  fclose(fp);

}

iftMKernel *iftCreateMKernel(iftAdjRel *A, int nbands){
	iftMKernel* kernel = (iftMKernel*)iftAlloc(1,sizeof(iftMKernel));

	kernel->A = iftCopyAdjacency(A);
	kernel->nbands = nbands;

	kernel->weight = (iftBand*)iftAlloc(nbands,sizeof(iftBand));

	int i;
	for(i = 0; i < nbands; i++){
		kernel->weight[i].val = iftAllocFloatArray(A->n);
	}

	return kernel;
}


void iftDestroyMKernel(iftMKernel **K){
	iftMKernel *kernel = *K;

	int i;
	for(i = 0; i < kernel->nbands; i++)
	  iftFree(kernel->weight[i].val);
	iftFree(kernel->weight);
	
	iftDestroyAdjRel(&kernel->A);
	
	iftFree(kernel);
	*K = 0;
}

iftMKernel *iftCopyMKernel(const iftMKernel *K) {
    iftMKernel *cpy = NULL;

    cpy = iftCreateMKernel(K->A, K->nbands);

    for(int i = 0; i < K->nbands; i++) {
        for(int j = 0; j < K->A->n; j++) {
            cpy->weight[i].val[j] = K->weight[i].val[j];
        }
    }

    return cpy;
}

iftMKernel *iftRandomMKernel(iftAdjRel *A, int nbands){
  iftMKernel *kernel = iftCreateMKernel(A,nbands);

  int p,b;
  float average = 0.0;
  for(b = 0; b < nbands; b++){
    for(p = 0; p < A->n; p++){
      float w = ((float)rand())/((float)RAND_MAX + 0.5);
      kernel->weight[b].val[p] = w;
      average += w;
    }
  }
  average = average/(A->n * nbands);
  
  //The norm cannot be computed before we have the average
  float norm = 0.0;
  for(b = 0; b < nbands; b++){
    for(p = 0; p < A->n; p++){
      kernel->weight[b].val[p] = kernel->weight[b].val[p] - average;
      norm += (kernel->weight[b].val[p])*(kernel->weight[b].val[p]);
    }
  }
  norm = sqrt(norm);
  if (norm > IFT_EPSILON) {
    //The final value cannot be computed before the norm
    for(b = 0; b < nbands; b++){
      for(p = 0; p < A->n; p++){
	kernel->weight[b].val[p] = kernel->weight[b].val[p]/norm;
      }
    }
  }
  
  return kernel;
}


iftMMKernel *iftCreateMMKernel(iftAdjRel *A, int nbands, int nkernels){
  iftMMKernel * k_bank = (iftMMKernel*) iftAlloc(1, sizeof(iftMMKernel));

  k_bank->A        = iftCopyAdjacency(A);
  k_bank->nkernels = nkernels;
  k_bank->nbands   = nbands;
  k_bank->W        = NULL;
  k_bank->mean     = NULL;

  k_bank->bias = iftAllocFloatArray(nkernels);

  k_bank->weight   = (iftBand **) iftAlloc(nkernels, sizeof(iftBand *));
  if (k_bank->weight == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateMMKernel");
  for (int k = 0 ; k < nkernels ; k++){
    k_bank->weight[k] = (iftBand *) iftAlloc(nbands, sizeof(iftBand ));
    if (k_bank->weight[k] == NULL)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateMMKernel");
    for (int b = 0 ; b < nbands ; b++){
      k_bank->weight[k][b].val = iftAllocFloatArray(A->n);
    }
  }
  
  return k_bank;

}

iftMMKernel *iftReadMMKernel(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadMMKernel", filename);

    iftMMKernel *k_bank = iftReadMMKernelFILE(fp);

    fclose(fp);

    return k_bank;
}



iftMMKernel *iftReadMMKernelFILE(FILE* fp)
{
  int whitening, dim,size,nkernels,nbands;

  if (fread(&dim     ,sizeof(int),1,fp) != 1) iftError("Reading error", "iftReadMMKernel");
  if (fread(&size    ,sizeof(int),1,fp) != 1) iftError("Reading error", "iftReadMMKernel");
  if (fread(&nbands  ,sizeof(int),1,fp) != 1) iftError("Reading error", "iftReadMMKernel");
  if (fread(&nkernels,sizeof(int),1,fp) != 1) iftError("Reading error", "iftReadMMKernel");

  /* always considering squared adjacencies */
  iftAdjRel* A;
  if (dim == 3) {
    size = powf(size,1./3);
    A = iftCuboid(size,size,size);
  }
  else {
    size = powf(size,0.5);
    A = iftRectangular(size,size);
  }
  //  fprintf(stderr,"dim: %d, size: %d, nbands: %d, nkernels: %d\n",dim,size,nbands,nkernels);
  iftMMKernel* K = iftCreateMMKernel(A,nbands,nkernels);

  /* Read ts and weights of the adjacent voxels */
  for (int k=0; k < K->nkernels; k++) {
    for (int b=0; b < K->nbands; b++) {	  
      if (fread(K->weight[k][b].val,sizeof(float),A->n,fp) != K->A->n) iftError("Reading error", "iftReadMMKernel");
    }
  }

  if (fread(K->bias, sizeof(float), K->nkernels, fp) != K->nkernels) iftError("Reading error", "iftReadMMKernel");

  if (fread(&whitening,sizeof(int),1,fp) != 1) iftError("Reading error", "iftReadMMKernel");

  if (whitening) {
    K->W     = iftReadRawMatrixFILE(fp);
    //    iftPrintMatrix(K->W);
    K->mean  = iftAllocFloatArray(A->n*nbands);
    if (fread(K->mean ,sizeof(float),A->n*nbands,fp) != A->n*nbands) iftError("Reading error", "iftReadMMKernel");
    //    iftPrintFloatArray(K->mean,A->n*nbands);
  }
  
  iftDestroyAdjRel(&A); 

  return(K);
}

void iftWriteMMKernel(const iftMMKernel *k_bank, const char *filename) {
    FILE *fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteMMKernel", filename);

    iftWriteMMKernelFILE(k_bank, fp);

    fclose(fp);
}



void iftWriteMMKernelFILE(const iftMMKernel *K, FILE* fp)
{
  int dim;
  if ( iftIsAdjRel3D(K->A) ) 
    dim = 3;
  else
    dim = 2;

  if (fwrite(&dim          ,sizeof(int),1,fp) != 1) iftError("Writing error", "iftWriteMMKernel");
  if (fwrite(&(K->A->n)    ,sizeof(int),1,fp) != 1) iftError("Writing error", "ittWriteMMKernel");
  if (fwrite(&(K->nbands  ),sizeof(int),1,fp) != 1) iftError("Writing error", "iftWriteMMKernel");
  if (fwrite(&(K->nkernels),sizeof(int),1,fp) != 1) iftError("Writing error", "iftWriteMMKernel");

  for (int k=0; k < K->nkernels; k++) {
    for (int b=0; b < K->nbands; b++) {
      if (fwrite(K->weight[k][b].val,sizeof(float),K->A->n,fp) != K->A->n) iftError("Writing error", "iftWriteMMKernel");
    }
  }

  if (fwrite(K->bias, sizeof(float), K->nkernels, fp) != K->nkernels) iftError("Writing error", "iftWriteMMKernel");

  if (K->W != NULL){
    int whitening = 1;
    if (fwrite(&whitening,sizeof(int),1,fp) != 1) iftError("Writening error", "iftWriteMMKernel");
    iftWriteRawMatrixFILE(K->W   ,fp);
    //iftPrintMatrix(K->W);
    if (fwrite(K->mean ,sizeof(float),K->A->n*K->nbands,fp) != K->A->n*K->nbands)
        iftError("Writening error", "iftWriteMMKernel");
    //iftPrintFloatArray(K->mean,K->A->n*K->nbands);
  }
  else {
    int whitening = 0;
    if (fwrite(&whitening,sizeof(int),1,fp) != 1) iftError("Writening error", "iftWriteMMKernel");
  }
}


iftMMKernel *iftCopyMMKernel(iftMMKernel *k)
{

  iftMMKernel* k_bank;
  k_bank = iftCreateMMKernel(k->A,k->nbands,k->nkernels);
  for (int u=0; u < k->nkernels; u++)  {
    for (int b=0; b < k->nbands; b++) {
      for (int i=0; i < k->A->n; i++) {
	k_bank->weight[u][b].val[i] = k->weight[u][b].val[i];
      }
    }
    k_bank->bias[u] = k->bias[u];
  }

  if (k->W != NULL) {
    k_bank->W = iftCopyMatrix(k->W);
    k_bank->mean  = iftAllocFloatArray(k->nbands*k->A->n);
    for (int i=0; i < k->nbands*k->A->n; i++) {
      k_bank->mean[i]  = k->mean[i];
    }
  }

  return k_bank;

}

iftMMKernel* iftUnionMMKernel(iftMMKernel *k1, iftMMKernel *k2)
{
  iftMMKernel* kU;

  if ( (k1->A->n != k2->A->n) || (k1->nbands != k2->nbands) ) {
    char msg[200];
    sprintf(msg,"Kernel dimensions are incompatible!\nNeighboors (%d/%d), Bands (%d/%d)",
	    k1->A->n, k2->A->n, k1->nbands, k2->nbands);
      iftError(msg, "iftUnionMMKernel");
  }    

  kU = iftCreateMMKernel(k1->A,k1->nbands,k1->nkernels+k2->nkernels);

  for(int u=0; u < k1->nkernels ; u++) {
    for (int b=0; b < k1->nbands; b++) {
      for (int i=0; i < k1->A->n; i++) {
	kU->weight[u][b].val[i] = k1->weight[u][b].val[i];
      }
    }
    kU->bias[u] = k1->bias[u];
  }

  for(int u=0; u < k2->nkernels ; u++) {
    for (int b=0; b < k2->nbands; b++) {
      for (int i=0; i < k2->A->n; i++) {
	kU->weight[k1->nkernels+u][b].val[i] = k2->weight[u][b].val[i];
      }
    }
    kU->bias[k1->nkernels+u] = k2->bias[u];
  }

  if ( ((k1->W != NULL) && (k2->W == NULL)) || ((k2->W != NULL) && (k1->W == NULL)) ) {
    char msg[200];
    sprintf(msg,"Kernels are from different nature (Whitening/NonWhitening)");
      iftError(msg, "iftUnioMMKernel");
  }

  if ( (k1->W != NULL) && (k2->W != NULL) ) {
    // Call - k1 parameters are chosen
    kU->W = iftCopyMatrix(k1->W);
    kU->mean  = iftAllocFloatArray(k1->nbands*k1->A->n);
    for (int i=0; i < k1->nbands*k1->A->n; i++) {
      kU->mean[i]  = k1->mean[i];
    }
  }

  return kU;
}

iftMMKernel *iftRandomMMKernel(iftAdjRel *A, int nbands, int nkernels){
  iftMMKernel * k_bank = NULL;

  k_bank = iftCreateMMKernel(A, nbands, nkernels);

  for (int k=0; k < nkernels; k++) {
    float mean=0.0, norm=0.0;
    for (int b=0; b < nbands; b++) {
      for (int i=0; i < A->n; i++) {
	k_bank->weight[k][b].val[i] = 
	  ((float)rand())/((float)RAND_MAX + 0.5); // generate random numbers in 0..1 range
	mean += k_bank->weight[k][b].val[i];
      }
    }
    mean /= (A->n *nbands); 

    for (int b=0; b < nbands; b++) {
      for (int i=0; i < A->n; i++) {
	k_bank->weight[k][b].val[i] = (k_bank->weight[k][b].val[i]-mean); 
	norm += (k_bank->weight[k][b].val[i]*k_bank->weight[k][b].val[i]);
      }
    }
    norm = sqrtf(norm);
    if (norm > IFT_EPSILON) {
      for (int b=0; b < nbands; b++) {
	for (int i=0; i < A->n; i++) {
	  k_bank->weight[k][b].val[i] /= norm;
	}
      }
    }
  }

  return(k_bank);

}

iftMMKernel* iftV1likeMMKernel2D(int sizex, int sizey, int nbands, int norients,int nfreqs, float *freqs){
  iftMMKernel * k_bank = NULL;

  int nkernels = norients*nfreqs;

  iftAdjRel* A = iftRectangular(sizex,sizey);
  k_bank = iftCreateMMKernel(A, nbands, nkernels);

  float orients[norients];
  for (int o=0;o<norients;o++)
    orients[o] = (o+1) * IFT_PI / norients;

  int dx,dy,dz;
  float gx0,gy0,gw,gh;
  iftMaxAdjShifts(A,&dx,&dy,&dz);
  gx0 = dx;gy0 = dy;dz=0;gw = gx0; gh = gy0;

  for (int k=0; k < nkernels; k++) {

    iftKernel* gabor = iftGabor2D(gw,gh,gx0,gy0,freqs[k/norients],orients[k%norients],0.,A);

    // normalizing to zero mean and unit length the k-th gabor filter
    float mean=0.0, norm=0.0;
    for (int b=0; b < nbands; b++) {
      for (int i=0; i < A->n; i++) {
	k_bank->weight[k][b].val[i] = gabor->weight[i];
	mean += k_bank->weight[k][b].val[i];
      }
    }
    iftDestroyKernel(&gabor);
    mean /= ( A->n * nbands); 

    for (int b=0; b < nbands; b++) {
      for (int i=0; i < A->n; i++) {
	k_bank->weight[k][b].val[i] = (k_bank->weight[k][b].val[i]-mean); 
	norm += (k_bank->weight[k][b].val[i]*k_bank->weight[k][b].val[i]);
      }
    }
    norm = sqrtf(norm);
    if (norm > 1.) {
      for (int b=0; b < nbands; b++) {
	for (int i=0; i < A->n; i++) {
	  k_bank->weight[k][b].val[i] /= norm;
	}
      }
    }
  }

  iftDestroyAdjRel(&A);

  return(k_bank);
}


iftMMKernel *iftRandomZMMMKernel(iftAdjRel *A, int nbands, int nkernels){
  iftMMKernel * k_bank = NULL;

  k_bank = iftCreateMMKernel(A, nbands, nkernels);

  for (int k=0; k < nkernels; k++) {
    float mean=0.0;
    for (int b=0; b < nbands; b++) {
      for (int i=0; i < A->n; i++) {
	k_bank->weight[k][b].val[i] = 
	  ((float)rand())/((float)RAND_MAX + 0.5);
	mean += k_bank->weight[k][b].val[i];
      }
    }
    mean /= (A->n *nbands); 
  }

  return(k_bank);

}

iftMMKernel *iftOnesMMKernel(iftAdjRel *A, int nbands, int nkernels){
  iftMMKernel * k_bank = NULL;

  k_bank = iftCreateMMKernel(A, nbands, nkernels);

#pragma omp parallel for shared(k_bank,A,nbands,nkernels)
  for (int k=0; k < nkernels; k++) {
    for (int b=0; b < nbands; b++) {
      for (int i=0; i < A->n; i++) {
        if (b*(nkernels/nbands)+i==k)
          k_bank->weight[k][b].val[i] = 1.0;
        else
          k_bank->weight[k][b].val[i] = 0.0;
      }
    }
  }

  return(k_bank);

}

void iftDestroyMMKernel(iftMMKernel **k_bank){
  iftMMKernel *aux = *k_bank;

  if (aux->A != NULL) iftDestroyAdjRel(&aux->A);
  for (int k=0; k < aux->nkernels; k++) {
    for (int b=0; b < aux->nbands; b++) {
      iftFree(aux->weight[k][b].val);
    }
    iftFree(aux->weight[k]);
  }
  iftFree(aux->weight);
  iftFree(aux->bias);
  if (aux->W != NULL) iftDestroyMatrix(&aux->W);
  if (aux->mean  != NULL) iftFree(aux->mean);

  iftFree(aux);
  *k_bank = NULL;
  
}

iftKernel *iftCreateKernel(iftAdjRel *A)
{
  iftKernel *K=(iftKernel *) iftAlloc(1,sizeof(iftKernel));

  K->A      = iftCopyAdjacency(A); 
  K->weight = iftAllocFloatArray(K->A->n);

  return(K);
}

void iftDestroyKernel(iftKernel **K)
{
  iftKernel *aux=*K;

  if (aux != NULL) {
    iftDestroyAdjRel(&aux->A);
    iftFree(aux->weight);
    iftFree(aux);
    *K = NULL;
  }
}

 float *GaborEnvelope2D(iftAdjRel *A, int x0, int y0, float a, float b, float theta){
  /*
  Auxiliar function to create a GaborKernel. Calculates only the Gaussian envelope part.
  Parameters:
  (x0, y0): Location of the peak of the Gaussian envelope
  (a,b): Scales the 2 axis x and y respectively. This basically changes the covariance matrix of 2D Gaussian. 
  theta: Rotation angle of Gaussian envelope

  Returns a float pointer containing the weights of gaussian envelope
  Read: https://inc.ucsd.edu/mplab/75/media//gabor.pdf section 2.2 for further reference
  */ 

 float *envelope = (float *)(calloc(A->n, sizeof(float)));

 for (int i = 0; i < A->n; i++){
   float x_rotation = (A->dx[i] - x0)*cos(theta) + (A->dy[i] - y0)*sin(theta); // (x - x0)r = (x - x0)cos(theta) + (y - y0)sen(theta)
   float y_rotation = -(A->dx[i] - x0)*sin(theta) + (A->dy[i] - y0)*cos(theta); // (y - y0)r = -(x - x0)sin(theta) + (y - y0)cos(theta)
   
   envelope[i] = exp(-IFT_PI * (pow(a,2.0) * pow(x_rotation,2.0) + pow(b, 2.0) * pow(y_rotation, 2.0)));
 }
 
 return(envelope);
}

float *RealGabor2D(iftAdjRel *A, float *envelope, float u0, float v0, float P){
  /*
  Auxiliar function to create a GaborKernel. Calculates only the real part of the complex Gabor Kernel.
  Parameters:
  (u0, v0): Spatial frequencies of the Gabor Kernel in Cartesian coordinates
  P: Phase of the Gabor Kernel

  Returns a floit ponter containing the weights of the real part of gabor kernel.
  Read: https://inc.ucsd.edu/mplab/75/media//gabor.pdf section 2.3 for further reference
  */

 float *real_gabor = (float *)(calloc(A->n, sizeof(float)));
 
 for (int i = 0; i < A->n; i++){
    real_gabor[i] = envelope[i] * cos(2.0 * IFT_PI * (u0 * A->dx[i] + v0 * A->dy[i]) + P);
    
 }

 return(real_gabor);

}

float *ImagGabor2D(iftAdjRel *A, float *envelope, float u0, float v0, float P){
  /*
  Auxiliar function to create a GaborKernel. Calculates only the imaginary part of the complex Gabor kernel.
  Parameters:
  (u0, v0): Spatial frequencies of the Gabor kernel in Cartesian coordinates
  P: Phase of the Gabor kernel

  Returns a floit ponter containing the weights of the imag part of gabor kernel.
  Read: https://inc.ucsd.edu/mplab/75/media//gabor.pdf section 2.3 for further reference
  */

 float *imag_gabor = (float *)(calloc(A->n, sizeof(float)));
 
 for (int i = 0; i < A->n; i++){
    imag_gabor[i] = envelope[i] * sin(2.0 * IFT_PI * (u0 * A->dx[i] + v0 * A->dy[i]) + P);

 }

 return(imag_gabor);

}


iftMKernel *iftGaborMKernel2D(int xsize, int ysize, float u0, float v0, float P, int x0, int y0, float a, float b, float theta)
{
  iftAdjRel *A = iftRectangular(xsize,ysize);
  iftMKernel *K = iftCreateMKernel(A, 2);
  float *gabor_envelope = GaborEnvelope2D(A, x0, y0, a, b, theta);

  K->weight[0].val = RealGabor2D(A, gabor_envelope, u0, v0, P);
  K->weight[1].val = ImagGabor2D(A, gabor_envelope, u0, v0, P);

  
  iftDestroyAdjRel(&A);

  return(K);
  }



iftMMKernel *iftGaborMMKernelBank2D(int xsize, int ysize, int n_kernels, float **gabor_parameters, int nbands)
{
  /*
  Gabor filter implementation in the space domain.
  Follow the reference for more: https://inc.ucsd.edu/mplab/75/media//gabor.pdf

  Function parameters:
  int xsize, int ysize: size of gabor space kernel
  int n_kernels: number of gabor kernels in the bank
  float **gabor_parameters: list of parameters listed below

  Gabor Parameters:
  float (u0,v0): sinusoid carrier space frequency along x and y, respectively;
  float (a,b): scales the 2 axis of the Gaussian envelope
  float theta: rotation angle of the gaussian envelope
  float P: phase of the sinusoid carrier
  float (x0,y0): location of the peak of the Gaussian envelope
 
  gabor(x,y, parameters) = sinusoid(x,y).envelope(x,y)
  */

  iftAdjRel   *A = iftRectangular(xsize, ysize);
  iftMMKernel *K = iftCreateMMKernel(A, nbands, 2 * n_kernels);
  float *gabor_envelope=NULL;
  int    param_counter = 0; // since i doesn't increments by 1, we need another counter variable

  for (int i = 0; i < 2 * n_kernels; i = i + 2){
    gabor_envelope = GaborEnvelope2D(A, gabor_parameters[param_counter][6], gabor_parameters[param_counter][7], gabor_parameters[param_counter][2], gabor_parameters[param_counter][3], gabor_parameters[param_counter][4]);
    
    // propagate the real and imaginary parts through every band of filter
    for (int b = 0; b < nbands; b++){
      K->weight[i][b].val = RealGabor2D(A, gabor_envelope, gabor_parameters[param_counter][0], gabor_parameters[param_counter][1], gabor_parameters[param_counter][5]);
      K->weight[i + 1][b].val = ImagGabor2D(A, gabor_envelope, gabor_parameters[param_counter][0], gabor_parameters[param_counter][1], gabor_parameters[param_counter][5]);
    }
    param_counter++;
    iftFree(gabor_envelope);
  }

  iftDestroyAdjRel(&A);

  return(K);
  
}




iftKernel * iftGaussianKernel(float radius, float stdev)
{
  iftAdjRel *A = iftSpheric(radius);
  iftKernel *K = iftCreateKernel(A);
  int       i;
  float     dist, var = stdev*stdev, sum=0.0;

  for (i=0; i < A->n; i++) {
    dist = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];
    K->weight[i] = exp(-dist/(2*var));
    sum += K->weight[i];
  }

  for (i=0; i < A->n; i++) {
    K->weight[i] /= sum;
  }

  iftDestroyAdjRel(&A);

  return(K);
}

iftKernel *iftGaussianKernel2D(float radius, float stdev)
{
  iftAdjRel *A = iftCircular(radius);
  iftKernel *K = iftCreateKernel(A);
  int       i;
  float     dist, var = stdev*stdev, sum=0.0;

  for (i=0; i < A->n; i++) {
    dist = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i];
    K->weight[i] = exp(-dist/(2*var));
    sum += K->weight[i];
  }

  for (i=0; i < A->n; i++) {
    K->weight[i] /= sum;
  }

  iftDestroyAdjRel(&A);

  return(K);
}

iftKernel *iftSobelXKernel2D()
{
  iftAdjRel *A = iftRectangular(3,3);
  iftKernel *K = iftCreateKernel(A);
  int i;

  for (i=0; i < A->n; i++) {
    if (A->dx[i]==0)
      K->weight[i]=0;
    else{
      if (A->dx[i]==-1){
	if (A->dy[i]==0){
	  K->weight[i]=-2;
	}else{
	  K->weight[i]=-1;
	}
      }else{ /* A->dx[i]==1 */
	if (A->dy[i]==0){
	  K->weight[i]=2;
	}else{
	  K->weight[i]=1;
	}
      }
    }
  }
	
  iftDestroyAdjRel(&A);

  return(K);
}
iftKernel *iftSobelYKernel2D()
{
  iftAdjRel *A = iftRectangular(3,3);
  iftKernel *K = iftCreateKernel(A);
  int i;

  for (i=0; i < A->n; i++) {
    if (A->dy[i]==0)
      K->weight[i]=0;
    else{
      if (A->dy[i]==-1){
	if (A->dx[i]==0){
	  K->weight[i]=-2;
	}else{
	  K->weight[i]=-1;
	}
      }else{ /* A->dy[i]==1 */
	if (A->dx[i]==0){
	  K->weight[i]=2;
	}else{
	  K->weight[i]=1;
	}
      }
    }
  }
	
  iftDestroyAdjRel(&A);

  return(K);
}

iftKernel *iftSobelXKernel()
{
  iftAdjRel *A = iftCuboid(3,3,3);
  iftKernel *K = iftCreateKernel(A);
  int i;

  for (i=0; i < A->n; i++) {
    if (A->dx[i]==0)
      K->weight[i]=0;
    else{
      if (A->dx[i]==-1){
	if ((A->dy[i]==0)&&(A->dz[i]==0)){
	  K->weight[i]=-4;
	}else{
	  if ((abs(A->dy[i])+abs(A->dz[i]))==1){
	    K->weight[i]=-2;
	  }else{ /* (abs(A->dy[i])+abs(A->dz[i]))==2 */
	    K->weight[i]=-1;
	  }
	}
      }else{ /* A->dx[i]==1 */
	if ((A->dy[i]==0)&&(A->dz[i]==0)){
	  K->weight[i]=4;
	}else{
	  if ((abs(A->dy[i])+abs(A->dz[i]))==1){
	    K->weight[i]=2;
	  }else{ /* (abs(A->dy[i])+abs(A->dz[i]))==2 */
	    K->weight[i]=1;
	  }
	}
      }
    }
  }
	
  iftDestroyAdjRel(&A);

  return(K);
}

iftKernel *iftSobelYKernel()
{
  iftAdjRel *A = iftCuboid(3,3,3);
  iftKernel *K = iftCreateKernel(A);
  int i;

  for (i=0; i < A->n; i++) {
    if (A->dy[i]==0)
      K->weight[i]=0;
    else{
      if (A->dy[i]==-1){
	if ((A->dx[i]==0)&&(A->dz[i]==0)){
	  K->weight[i]=-4;
	}else{
	  if ((abs(A->dx[i])+abs(A->dz[i]))==1){
	    K->weight[i]=-2;
	  }else{ /* (abs(A->dx[i])+abs(A->dz[i]))==2 */
	    K->weight[i]=-1;
	  }
	}
      }else{ /* A->dy[i]==1 */
	if ((A->dx[i]==0)&&(A->dz[i]==0)){
	  K->weight[i]=4;
	}else{
	  if ((abs(A->dx[i])+abs(A->dz[i]))==1){
	    K->weight[i]=2;
	  }else{ /* (abs(A->dx[i])+abs(A->dz[i]))==2 */
	    K->weight[i]=1;
	  }
	}
      }
    }
  }
	
  iftDestroyAdjRel(&A);

  return(K);
}

iftKernel *iftSobelZKernel()
{
  iftAdjRel *A = iftCuboid(3,3,3);
  iftKernel *K = iftCreateKernel(A);
  int i;

  for (i=0; i < A->n; i++) {
    if (A->dz[i]==0)
      K->weight[i]=0;
    else{
      if (A->dz[i]==-1){
	if ((A->dx[i]==0)&&(A->dy[i]==0)){
	  K->weight[i]=-4;
	}else{
	  if ((abs(A->dx[i])+abs(A->dy[i]))==1){
	    K->weight[i]=-2;
	  }else{ /* (abs(A->dx[i])+abs(A->dy[i]))==2 */
	    K->weight[i]=-1;
	  }
	}
      }else{ /* A->dz[i]==1 */
	if ((A->dx[i]==0)&&(A->dy[i]==0)){
	  K->weight[i]=4;
	}else{
	  if ((abs(A->dx[i])+abs(A->dy[i]))==1){
	    K->weight[i]=2;
	  }else{ /* (abs(A->dx[i])+abs(A->dy[i]))==2 */
	    K->weight[i]=1;
	  }
	}
      }
    }
  }
	
  iftDestroyAdjRel(&A);

  return(K);
}

iftKernel* iftGabor2D(float gw,float gh,float gx0,float gy0,float wfreq,float worient,float wphase,iftAdjRel* A){
// Inputs
//     gw -- width of the gaussian envelope
//     gh -- height of the gaussian envelope
//     gx0 -- x indice of center of the gaussian envelope
//     gy0 -- y indice of center of the gaussian envelope
//     wfreq -- frequency of the 2d wave
//     worient -- orientation of the 2d wave
//     wphase -- phase of the 2d wave
//     shape -- shape tuple (height, width)
//
// Outputs:
//     gabor -- 2d gabor with zero mean and unit-variance
  int l;

  iftKernel* gabor = iftCreateKernel(A);

  // [x,y] = meshgrid(0:width-1,0:height-1);
  int *x = iftAllocIntArray(gabor->A->n);
  int *y = iftAllocIntArray(gabor->A->n);
  for(l=0;l<gabor->A->n;l++) {
    x[l] = gabor->A->dx[l]+gx0;
    y[l] = gabor->A->dy[l]+gy0;
  }

  // X = x .* (cos(worient) * wfreq);
  // Y = y .* (sin(worient) * wfreq);
  float* X = iftAllocFloatArray(gabor->A->n);
  float* Y = iftAllocFloatArray(gabor->A->n);
  float cosw = cos(worient) * wfreq;
  float sinw = sin(worient) * wfreq;
  for(l=0;l<gabor->A->n;l++) {
    X[l] = x[l] * cosw;
    Y[l] = y[l] * sinw;
  }

  // env  = exp( -pi * ( (x-gx0).^2 ./ (gw^2) + (y-gy0).^2 / (gh^2) ) );
  // wave = exp( j * (2.*pi*(X+Y) + wphase) )
  float* env   = iftAllocFloatArray(gabor->A->n);
  float* wave  = iftAllocFloatArray(gabor->A->n);
  for(l=0;l<gabor->A->n;l++) {
    env[l]  = exp(-IFT_PI * ((x[l] - gx0) * (x[l] - gx0) / (gw * gw) + (y[l] - gy0) * (y[l] - gy0) / (gh * gh) ) );
    // For considering only the real part of e^{j*ang} = cos(ang) + j sin(ang), we take only the cosine part
    wave[l] = cos(2 * IFT_PI * (X[l] + Y[l]) + wphase);
 
    // gabor = real(env .* wave);
    gabor->weight[l] = env[l] * wave[l];
  }
  iftFree(x);iftFree(y); iftFree(X);iftFree(Y);

  iftFree(env);iftFree(wave);

  return gabor;
}

iftKernel *iftUniformKernel(iftAdjRel *A)
{
  iftKernel *K=iftCreateKernel(A);
  int i;
  
  for (i=0; i < A->n; i++) 
    K->weight[i]=1.0;

  return(K);
}

iftHessianKernel *iftCreateHessianKernel()
{
    iftHessianKernel *Dk = (iftHessianKernel*)iftAlloc(1,sizeof(iftHessianKernel));
    return Dk;
}

void iftDestroyHessianKernel(iftHessianKernel **K)
{
    iftDestroyKernel(&(*K)->Kxx);
    iftDestroyKernel(&(*K)->Kxy);
    iftDestroyKernel(&(*K)->Kxz);
    iftDestroyKernel(&(*K)->Kyy);
    iftDestroyKernel(&(*K)->Kyz);
    iftDestroyKernel(&(*K)->Kzz);
    iftFree(*K);
    *K=NULL;
}


iftHessianKernel *iftGaussianHessianKernels(float radius, float stdev)
{
  iftHessianKernel *Dg = iftCreateHessianKernel();
  iftAdjRel *A = iftSpheric(radius);
  float     dist, var = stdev*stdev;

  Dg->Kxx = iftCreateKernel(A);
  Dg->Kyy = iftCreateKernel(A);
  Dg->Kzz = iftCreateKernel(A);
  Dg->Kxy = iftCreateKernel(A);
  Dg->Kxz = iftCreateKernel(A);
  Dg->Kyz = iftCreateKernel(A);
  for (int i = 0; i < A->n; i++){
    dist = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];

    Dg->Kxx->weight[i] = 1/sqrt((2*IFT_PI*var)*(2*IFT_PI*var)*(2*IFT_PI*var))*exp(-dist/(2*var))*((A->dx[i]*A->dx[i])/(var*var)-1/var);
    Dg->Kyy->weight[i] = 1/sqrt((2*IFT_PI*var)*(2*IFT_PI*var)*(2*IFT_PI*var))*exp(-dist/(2*var))*((A->dy[i]*A->dy[i])/(var*var)-1/var);
    Dg->Kzz->weight[i] = 1/sqrt((2*IFT_PI*var)*(2*IFT_PI*var)*(2*IFT_PI*var))*exp(-dist/(2*var))*((A->dz[i]*A->dz[i])/(var*var)-1/var);
    Dg->Kxy->weight[i] = 1/sqrt((2*IFT_PI*var)*(2*IFT_PI*var)*(2*IFT_PI*var))*exp(-dist/(2*var))*(A->dx[i]*A->dy[i]/var);
    Dg->Kxz->weight[i] = 1/sqrt((2*IFT_PI*var)*(2*IFT_PI*var)*(2*IFT_PI*var))*exp(-dist/(2*var))*(A->dx[i]*A->dz[i]/var);
    Dg->Kyz->weight[i] = 1/sqrt((2*IFT_PI*var)*(2*IFT_PI*var)*(2*IFT_PI*var))*exp(-dist/(2*var))*(A->dy[i]*A->dz[i]/var);

    /*
     * Normalization of the gaussian second order derivative for cross-scale comparison.
     * For more information, please refer to the paper:
     *  Marcin Rudzki, "Vessel Detection Method Based on Eigenvalues of the Hessian Matrix
     *  and its Applicability to Airway Tree Segmentation". 2009.
     */
    Dg->Kxx->weight[i] *= var;
    Dg->Kyy->weight[i] *= var;
    Dg->Kzz->weight[i] *= var;
    Dg->Kxy->weight[i] *= var;
    Dg->Kxz->weight[i] *= var;
    Dg->Kyz->weight[i] *= var;

  }

  return Dg;
}
