#ifndef _iftMSPS_H_
#define _iftMSPS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftMatrix.h"
#include "iftFitnessFunc.h"
#include "iftMemory.h"

struct ift_MSPS;

/*Function to generate the perturbation, given a matrix error in the ith dimension, jth scale*/
typedef void (*iftPerturbationFunc)(struct ift_MSPS *msps, int i, int j, float *perturbation);

typedef void (*iftMinMaxPerturbationFunc)(struct ift_MSPS *msps, float *cur, int i, int j, float *posPerturbation,
                                          float *negPerturbation);

/* Data structure for Multiscale Parameter Search */
typedef struct ift_MSPS {
  iftFitnessFunc iftFitness;   /* fitness function for optimization. */
  void   *problem;             /* the context of the problem */
  int     n;                   /* number of parameters for optimization */
  int     m;                   /* number of scales */
  float     *theta;            /* parameter vector to be optimized */
  iftMatrix *delta;            /* matrix of displacements for each
				  parameter and all scales */

  iftMatrix* sigma;		/*
				* Random deviation of the matrix error
				* for each parameter and scale
				*/

  float *max; /** Maximum value for each parameter, NULL if there are no boundaries */
  float *min; /** Minimum value for each parameter, NULL if there are no boundaries */

  int     niters;              /* maximum number of iterations */

  float stopcriteria;/*If the gain between iterations is smaller than this stop*/

  iftMinMaxPerturbationFunc iftMinMaxPerturbation;
  iftPerturbationFunc iftPerturbation;/* generator of the perturbation vector */
  int iterstar;/*iteration that found the optimum value*/
  bool verbose;
} iftMSPS;

/**
 * @brief Creates an MSPS structure to optimize a black box function
 * @param n Number of parameters to optimize
 * @param m Number of scales
 * @param iftFitnessFunc Black box function
 * @param Pointer to a data-structure with other parameters specific to the problem
 *
 */
iftMSPS *iftCreateMSPS(int n, int m, iftFitnessFunc f, void *prob);

void iftDestroyMSPS(iftMSPS **msps);

float iftMSPSMax(iftMSPS *msps);

float iftMSPSMin(iftMSPS *msps);

float iftMSPSMinTest(iftMSPS *msps);

void iftMSPSRandomMinMaxPerturbation(iftMSPS *msps, float *cur, int row, int col, float *posout, float *negout);

void iftMSPSDeterministicMinMaxPerturbation(iftMSPS *msps, float *cur, int row, int col, float *posout, float *negout);

void iftMSPSLinearRandomPerturbation(iftMSPS *, int row, int col, float *out);

void iftMSPSGaussianRandomPerturbation(iftMSPS *, int row, int col, float *out);

void iftMSPSDeterministicPerturbation(iftMSPS *, int row, int col, float *out);

void iftMSPSRandomDirectionPerturbation(iftMSPS *, int row, int col, float *out);

#ifdef __cplusplus
}
#endif


#endif
