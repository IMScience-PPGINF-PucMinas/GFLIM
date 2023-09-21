#ifndef IFT_KERNEL_H_
#define IFT_KERNEL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftAdjacency.h"
#include "iftMatrix.h"
#include "iftSpectrum.h"

/* Multiband kernel: iftKernel *K; K->A->n; K->A->dx[i],
   K->A->dy[i],K->A->dz[i]; K->weight[i]  */
//! swig(destroyer = iftDestroyKernel)
typedef struct ift_kernel {
  iftAdjRel *A;
  float     *weight;
} iftKernel; 

/* Multiband kernel: iftMKernel *K; K->A->n; K->A->dx[i],
   K->A->dy[i],K->A->dz[i]; K->weight[b].val[i]  */

//! swig(destroyer = iftDestroyMKernel)
typedef struct ift_mkernel {
  iftAdjRel *A;
  iftBand   *weight;
  int        nbands;
} iftMKernel;

/* Kernel bank of multiband kernels: iftMMKernel *K; K->A->n;
   K->A->dx[i], K->A->dy[i],K->A->dz[i]; K->nkernels; K->nbands;
   K->weight[k][b].val[i]  */

//! swig(destroyer = iftDestroyMMKernel)
typedef struct ift_mmkernel {
  iftAdjRel *A;
  iftBand  **weight; /* one kernel per row and one band per column */
  float	    *bias; /* bias provided from SVM -> one bias per Kernel Multiband */
  int        nkernels;
  int        nbands;
  iftMatrix *W; /* Whitening transform */
  float     *mean; /* for centralization before whitening (size is W->ncols) */ 
} iftMMKernel;

typedef struct ift_hessian_kernel {
    /* Second order derivative kernel of a function F in relation to xx */
    iftKernel *Kxx;
    /* Second order derivative kernel of a function F in relation to yy */
    iftKernel *Kyy;
    /* Second order derivative kernel of a function F in relation to zz */
    iftKernel *Kzz;
    /* Second order derivative kernel of a function F in relation to xy */
    iftKernel *Kxy;
    /* Second order derivative kernel of a function F in relation to xz */
    iftKernel *Kxz;
    /* Second order derivative kernel of a function F in relation to yz */
    iftKernel *Kyz;
} iftHessianKernel;

iftKernel   *iftCreateKernel(iftAdjRel *A);
void         iftDestroyKernel(iftKernel **K);
iftKernel   *iftReadKernel(char *filename);
void         iftWriteKernel(iftKernel *K, char *filename); 
iftMKernel  *iftCreateMKernel(iftAdjRel *A, int nbands);
void         iftDestroyMKernel(iftMKernel **K);
iftMKernel  *iftRandomMKernel(iftAdjRel *A, int nbands);
iftMMKernel *iftCreateMMKernel(iftAdjRel *A, int nbands, int nkernels);
iftHessianKernel *iftCreateHessianKernels();
void iftDestroyHessianKernel(iftHessianKernel **K);
/**
 * @brief Copies an existing iftMKernel.
 *
 * @author Thiago Vallin Spina
 * @date Mar 01, 2016
 *
 */
iftMKernel  *iftCopyMKernel(const iftMKernel *K);
iftMMKernel *iftReadMMKernel(const char *filename);
iftMMKernel *iftReadMMKernelFILE(FILE* fp);
void         iftWriteMMKernel(const iftMMKernel *k_bank, const char *filename);
void         iftWriteMMKernelFILE(const iftMMKernel *k_bank,FILE* p);
iftMMKernel *iftCopyMMKernel(iftMMKernel *k);
iftMMKernel * iftUnionMMKernel(iftMMKernel *k1, iftMMKernel *k2);

//! swig(newobject, stable)
iftMMKernel *iftRandomMMKernel(iftAdjRel *A, int nbands, int nkernels);
iftMMKernel* iftV1likeMMKernel2D(int sizex,int sizey, int nbands, int norients, int nfreqs, float *freqs);
iftMMKernel *iftRandomZMMMKernel(iftAdjRel *A, int nbands, int nkernels);
iftMMKernel *iftOnesMMKernel(iftAdjRel *A, int nbands, int nkernels);
void         iftDestroyMMKernel(iftMMKernel **k_bank);

iftMKernel *iftGaborMKernel2D(int xsize, int ysize, float u0, float v0, float P, int x0, int y0, float a, float b, float theta);
iftMMKernel *iftGaborMMKernelBank2D(int xsize, int ysize, int n_kernels, float **gabor_parameters, int nbands);
//! swig(newobject, stable)
iftKernel   *iftGaussianKernel(float radius, float stdev);

//! swig(newobject, stable)
iftKernel   *iftGaussianKernel2D(float radius, float stdev);

//! swig(newobject, stable)
iftKernel   *iftSobelXKernel(void);

//! swig(newobject, stable)
iftKernel   *iftSobelYKernel(void);

//! swig(newobject, stable)
iftKernel   *iftSobelZKernel(void);

//! swig(newobject, stable)
iftKernel   *iftSobelXKernel2D(void);

//! swig(newobject, stable)
iftKernel   *iftSobelYKernel2D(void);
iftKernel   *iftDoGKernel(float stdev1, float stdev2);
iftKernel   *iftDoGKernel2D(float stdev1, float stdev2);

//! swig(newobject, stable)
iftKernel   *iftGabor2D(float gw,float gh,float gx0,float gy0,float wfreq,float worient,float wphase,iftAdjRel* A);
iftKernel   *iftUniformKernel(iftAdjRel *A); /* kernel with all weights equal to 1 */

iftHessianKernel *iftGaussianHessianKernels(float radius, float stdev);


#ifdef __cplusplus
}
#endif

#endif
