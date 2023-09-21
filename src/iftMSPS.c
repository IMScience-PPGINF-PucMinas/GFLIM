#include "iftMSPS.h"

#include "ift/core/io/Stream.h"


/* ------------------------- Public ----------------------------------------*/

void iftMSPSLinearRandomPerturbation(iftMSPS* msps, int row, int col, float *out)
{
  iftMatrix *delta = msps->delta;
  iftMatrix* sigma = msps->sigma;

  int i=0;
  float mean, deviation;

  for (i = 0; i < delta->ncols; ++i) {
    out[i] = 0.0;
  }

  i=iftGetMatrixIndex(delta, col, row);
  mean = delta->val[i];
  deviation = sigma->val[i];

  float f = ((float)iftRandomInteger(0,1000))/1000.0;
  out[col] = mean-deviation+f*2*deviation;
}

void iftMSPSRandomDirectionPerturbation(iftMSPS* msps, int row, int col, float *out)
{
  int i=0, index;
  float cur, prev;

  iftMatrix *delta = msps->delta;

  for (i = 0; i < delta->ncols; ++i) {
    index = iftGetMatrixIndex(delta, col, row);
    cur = delta->val[index];
    if(row>0)
    {
      index = iftGetMatrixIndex(delta, col, row-1);
    }
    else
    {
      prev = 0.0f;
    }
    out[i] = prev + (iftRandomInteger(0, 1000)/1000.0)*(cur - prev);
    if(rand()%2==0)
    {
      out[i] = -out[i];
    }
  }
}


void iftMSPSGaussianRandomPerturbation(iftMSPS* msps, int row, int col, float *out)
{
  iftMatrix *delta = msps->delta;
  iftMatrix *sigma = msps->sigma;
  int i=0;
  for (i = 0; i < delta->ncols; ++i) {
    out[i] = 0.0;
  }
  i=iftGetMatrixIndex(delta, col, row);
  out[col] = delta->val[i] + sigma->val[i]*iftRandomNormalFloat(0.0, 1.0);
}

void iftMSPSDeterministicPerturbation(iftMSPS* msps, int row, int col, float *out)
{
  iftMatrix* delta = msps->delta;

  int i=0;
  for (i = 0; i < delta->ncols; ++i) {
    out[i] = 0.0;
  }
  i=iftGetMatrixIndex(delta, col, row);
  out[col] = delta->val[i];
}

void iftMSPSRandomMinMaxPerturbation(iftMSPS* msps, float* cur, int row, int col, float *posout, float* negout)
{
  iftMatrix* delta = msps->delta;

  int i=0;
  for (i = 0; i < delta->ncols; ++i) {
    posout[i] = 0.0;
    negout[i] = 0.0;
  }

  float pos = iftRandomInteger(1,1000)/1000.0;

  i=iftGetMatrixIndex(delta, col, row);
  posout[col] = (msps->max[row] - cur[col])*pos;
  negout[col] = (msps->min[row] - cur[col])*pos;
}

void iftMSPSDeterministicMinMaxPerturbation(iftMSPS* msps, float* cur, int row, int col, float *posout, float* negout)
{
  iftMatrix* delta = msps->delta;

  int i=0;
  for (i = 0; i < delta->ncols; ++i) {
    posout[i] = 0.0;
    negout[i] = 0.0;
  }

  float pos = (col+1.0)/(delta->nrows + 1.0);

  i=iftGetMatrixIndex(delta, col, row);
  posout[col] = (msps->max[row] - cur[col])*pos;
  negout[col] = (msps->min[row] - cur[col])*pos;
}

iftMSPS *iftCreateMSPS(int n, int m, iftFitnessFunc f, void *prob)
{

  iftMSPS *msps=(iftMSPS *)iftAlloc(1,sizeof(iftMSPS));
  msps->n          = n;
  msps->m          = m;
  msps->theta      = iftAllocFloatArray(n);

  msps->min = iftAllocFloatArray(n);
  msps->max = iftAllocFloatArray(n);

  msps->delta      = iftCreateMatrix(n,m);
  msps->sigma      = iftCreateMatrix(n,m);
  msps->iftFitness = f;
  msps->problem    = prob;
  msps->iftPerturbation = iftMSPSDeterministicPerturbation;
  msps->iftMinMaxPerturbation = iftMSPSRandomMinMaxPerturbation;
  msps->niters     = IFT_INFINITY_INT-1;
  msps->iterstar   = 0;
  msps->stopcriteria = 0.0f;
  msps->verbose = false;
  return(msps);
}

void iftDestroyMSPS(iftMSPS **msps)
{
  iftMSPS *aux=*msps;

  if (aux != NULL) {
    iftFree(aux->theta);
    iftDestroyMatrix(&aux->delta);
    iftDestroyMatrix(&aux->sigma);
    iftFree(aux);
    *msps = NULL;
  }
}

float iftMSPSMax(iftMSPS *msps)
{
  int i,j,n=msps->n,m=msps->m;
  float V_star,*theta_star=iftAllocFloatArray(n);
  float  V0,V,Vplus,Vminus;
  float *Vp = iftAllocFloatArray(n); // best value for each
  // parameter along the scales
  float *Dp = iftAllocFloatArray(n); // best displacement for each
  // parameter along the scales
  float *Ds = iftAllocFloatArray(n); // best displacement for each parameter
  // in a given scale
  float *Dplus  = iftAllocFloatArray(n); // current positive displacement
  // for a given parameter
  float *Dminus = iftAllocFloatArray(n); // current negative displacement
  // for a given parameter
  float *D_star = iftAllocFloatArray(n); // best current displacement for
  // a given parameter
  float acquisition;
  float *theta=msps->theta;
  float *theta_plus  = iftAllocFloatArray(n); // positive perturbation
  float *theta_minus = iftAllocFloatArray(n); // negative perturbation
  iftFitnessFunc F=msps->iftFitness;
  int     iter=1;

  // initialize MSPS optimization
  iftCopyFloatArray(theta_star,theta,n);
  V_star = F(msps->problem,theta_star);


  do { // initialize for the next iteration

    V0 = V_star;
    iftCopyFloatArray(theta,theta_star,n);
    iftInitFloatArray(Dp,0.0,n);
    iftInitFloatArray(Vp,V0,n);
    for (j=0; j < m; j++) { // for each scale
      iftInitFloatArray(Ds,0.0,n);
      for (i=0; i < n; i++) { // for each parameter
        // evaluate local perturbation

        msps->iftPerturbation(msps, j, i, Dplus);
        iftScaleFloatArray(Dplus, -1.0f, Dminus, n);

        //	iftInitFloatArrayPos(Dplus,i,error->val[iftGetMatrixIndex(error,i,j)],n);
        //	iftInitFloatArrayPos(Dminus,i,-error->val[iftGetMatrixIndex(error,i,j)],n);
        V = V0;
        iftInitFloatArray(D_star,0.0,n);
        iftSumFloatArrays(theta,Dplus,theta_plus,n);
        iftSumFloatArrays(theta,Dminus,theta_minus,n);
        Vplus  = F(msps->problem,theta_plus);
        Vminus = F(msps->problem,theta_minus);
        if (Vplus > V) {
          V = Vplus; iftCopyFloatArray(D_star,Dplus,n);
        }
        if (Vminus > V) {
          V = Vminus; iftCopyFloatArray(D_star,Dminus,n);
        }
        if (V > V_star) { // update the best theta so far
          msps->iterstar = iter;
          V_star = V; iftSumFloatArrays(theta,D_star,theta_star,n);
        }
        if (Vp[i] < V) { // update the best displacement for
          // this parameter along the scales
          Vp[i] = V; Dp[i] = D_star[i];
        }
        iftSumFloatArrays(Ds,D_star,Ds,n); // update the best displacement
        // for each parameter in this scale
      }
      // update the best theta based on
      // the best displacement vector
      // of the current scale
      iftSumFloatArrays(theta,Ds,theta_plus,n);
      V = F(msps->problem,theta_plus);
      // printf("V = %f\nV_star = %f\n\n", V, V_star);
      if (V > V_star) {
        msps->iterstar = iter;
        V_star = V; iftCopyFloatArray(theta_star,theta_plus,n);
      }
    }
    // update the best theta based on the displacement vector
    // with the best parameters along the scales
    iftSumFloatArrays(theta,Dp,theta_plus,n);
    V = F(msps->problem,theta_plus);
    if (V > V_star) {
      V_star = V; iftCopyFloatArray(theta_star,theta_plus,n);
    }
    iter++;
    acquisition = fabs(V_star - V0);
  } while ((acquisition > msps->stopcriteria)&&(iter < msps->niters));

  if (iter >= msps->niters)
    iftWarning("Reached maximum number of iterations","iftMSPSMax");

  if(msps->verbose) {
    printf("Final fitness value %f, parameters: ", V_star);
    for (i = 0; i < n; i++)
      printf("%f ", theta_star[i]);
    printf("\n");
  }

  iftCopyFloatArray(theta,theta_star,n);

  iftFree(theta_star);
  iftFree(Vp);
  iftFree(Dp);
  iftFree(Ds);
  iftFree(Dplus);
  iftFree(Dminus);
  iftFree(D_star);
  iftFree(theta_plus);
  iftFree(theta_minus);

  return(V0);
}

float iftMSPSMinTest(iftMSPS* msps)
{
  int i,j,n=msps->n,m=msps->m;
  float V_star,*theta_star=iftAllocFloatArray(n);
  float  V0,V,Vplus,Vminus;
  float *Vp = iftAllocFloatArray(n); // best value for each
  // parameter along the scales
  float *Dp = iftAllocFloatArray(n); // best displacement for each
  // parameter along the scales
  float *Ds = iftAllocFloatArray(n); // best displacement for each parameter
  // in a given scale
  float *Dplus  = iftAllocFloatArray(n); // current positive displacement
  // for a given parameter
  float *Dminus = iftAllocFloatArray(n); // current negative displacement
  // for a given parameter
  float *D_star = iftAllocFloatArray(n); // best current displacement for
  // a given parameter
  float acquisition;
  float *theta=msps->theta;
  float *theta_plus  = iftAllocFloatArray(n); // positive perturbation
  float *theta_minus = iftAllocFloatArray(n); // negative perturbation
  iftFitnessFunc F=msps->iftFitness;
  int iter=1;

  // initialize MSPS optimization
  iftCopyFloatArray(theta_star,theta,n);
  V_star = F(msps->problem,theta_star);


  do { // initialize for the next iteration

    V0 = V_star;
    iftCopyFloatArray(theta,theta_star,n);
    iftInitFloatArray(Dp,0.0,n);
    iftInitFloatArray(Vp,V0,n);
    for (j=0; j < m; j++) { // for each scale
      iftInitFloatArray(Ds,0.0,n);
      for (i=0; i < n; i++) { // for each parameter
        // evaluate local perturbation

        msps->iftMinMaxPerturbation(msps, theta, j, i, Dplus, Dminus);
//        iftScaleFloatArray(Dplus, -1.0, Dminus, n);
        // iftInitFloatArrayPos(Dplus,i,error->val[iftGetMatrixIndex(error,i,j)],n);
        // iftInitFloatArrayPos(Dminus,i,-error->val[iftGetMatrixIndex(error,i,j)],n);

        V = V0;
        iftInitFloatArray(D_star,0.0,n);
        iftSumFloatArrays(theta,Dplus,theta_plus,n);
        iftSumFloatArrays(theta,Dminus,theta_minus,n);

        Vplus  = F(msps->problem,theta_plus);
        Vminus = F(msps->problem,theta_minus);

        if (Vplus < V) {
          V = Vplus; iftCopyFloatArray(D_star,Dplus,n);
        }
        if (Vminus < V) {
          V = Vminus; iftCopyFloatArray(D_star,Dminus,n);
        }
        if (V < V_star) { // update the best theta so far
          msps->iterstar = iter;
          V_star = V; iftSumFloatArrays(theta,D_star,theta_star,n);
        }
        if (Vp[i] > V) { // update the best displacement for
          // this parameter along the scales
          Vp[i] = V; Dp[i] = D_star[i];
        }
        iftSumFloatArrays(Ds,D_star,Ds,n); // update the best displacement
        // for each parameter in this scale
      }
      // update the best theta based on
      // the best displacement vector
      // of the current scale
      iftSumFloatArrays(theta,Ds,theta_plus,n);
      V = F(msps->problem,theta_plus);
      if (V < V_star) {
        msps->iterstar = iter;
        V_star = V; iftCopyFloatArray(theta_star,theta_plus,n);
      }
    }
    // update the best theta based on the displacement vector
    // with the best parameters along the scales
    iftSumFloatArrays(theta,Dp,theta_plus,n);
    V = F(msps->problem,theta_plus);
    if (V < V_star) {
      V_star = V; iftCopyFloatArray(theta_star,theta_plus,n);
    }
    iter++;
    acquisition = fabs(V_star - V0);
  } while (( acquisition > msps->stopcriteria) &&(iter < msps->niters));

  if (iter >= msps->niters)
    iftWarning("Reached maximum number of iterations","iftMSPSMin");


  if(msps->verbose) {
    printf("Final fitness value %f, parameters: ", V_star);
    for (i = 0; i < n; i++)
      printf("%f ", theta_star[i]);
    printf("\n");
  }

  iftCopyFloatArray(theta,theta_star,n);

  iftFree(theta_star);
  iftFree(Vp);
  iftFree(Dp);
  iftFree(Ds);
  iftFree(Dplus);
  iftFree(Dminus);
  iftFree(D_star);
  iftFree(theta_plus);
  iftFree(theta_minus);

  return(V0);
}

float iftMSPSMin(iftMSPS *msps)
{
  int i,j,n=msps->n,m=msps->m;
  float V_star,*theta_star=iftAllocFloatArray(n);
  float  V0,V,Vplus,Vminus;
  float *Vp = iftAllocFloatArray(n); // best value for each
  // parameter along the scales
  float *Dp = iftAllocFloatArray(n); // best displacement for each
  // parameter along the scales
  float *Ds = iftAllocFloatArray(n); // best displacement for each parameter
  // in a given scale
  float *Dplus  = iftAllocFloatArray(n); // current positive displacement
  // for a given parameter
  float *Dminus = iftAllocFloatArray(n); // current negative displacement
  // for a given parameter
  float *D_star = iftAllocFloatArray(n); // best current displacement for
  // a given parameter
  float acquisition;
  float *theta=msps->theta;
  float *theta_plus  = iftAllocFloatArray(n); // positive perturbation
  float *theta_minus = iftAllocFloatArray(n); // negative perturbation
  iftFitnessFunc F=msps->iftFitness;
  int iter=1;

  // initialize MSPS optimization
  iftCopyFloatArray(theta_star,theta,n);
  V_star = F(msps->problem,theta_star);


  do { // initialize for the next iteration

    V0 = V_star;
    iftCopyFloatArray(theta,theta_star,n);
    iftInitFloatArray(Dp,0.0,n);
    iftInitFloatArray(Vp,V0,n);
    for (j=0; j < m; j++) { // for each scale
      iftInitFloatArray(Ds,0.0,n);
      for (i=0; i < n; i++) { // for each parameter
        // evaluate local perturbation

        msps->iftPerturbation(msps, j, i, Dplus);
        iftScaleFloatArray(Dplus, -1.0, Dminus, n);

        // iftInitFloatArrayPos(Dplus,i,error->val[iftGetMatrixIndex(error,i,j)],n);
        // iftInitFloatArrayPos(Dminus,i,-error->val[iftGetMatrixIndex(error,i,j)],n);

        V = V0;
        iftInitFloatArray(D_star,0.0,n);
        iftSumFloatArrays(theta,Dplus,theta_plus,n);
        iftSumFloatArrays(theta,Dminus,theta_minus,n);

        Vplus  = F(msps->problem,theta_plus);
        Vminus = F(msps->problem,theta_minus);

        if (Vplus < V) {
          V = Vplus; iftCopyFloatArray(D_star,Dplus,n);
        }
        if (Vminus < V) {
          V = Vminus; iftCopyFloatArray(D_star,Dminus,n);
        }
        if (V < V_star) { // update the best theta so far
          msps->iterstar = iter;
          V_star = V; iftSumFloatArrays(theta,D_star,theta_star,n);
        }
        if (Vp[i] > V) { // update the best displacement for
          // this parameter along the scales
          Vp[i] = V; Dp[i] = D_star[i];
        }
        iftSumFloatArrays(Ds,D_star,Ds,n); // update the best displacement
        // for each parameter in this scale
      }
      // update the best theta based on
      // the best displacement vector
      // of the current scale
      iftSumFloatArrays(theta,Ds,theta_plus,n);
      V = F(msps->problem,theta_plus);
      if (V < V_star) {
        msps->iterstar = iter;
        V_star = V; iftCopyFloatArray(theta_star,theta_plus,n);
      }
    }
    // update the best theta based on the displacement vector
    // with the best parameters along the scales
    iftSumFloatArrays(theta,Dp,theta_plus,n);
    V = F(msps->problem,theta_plus);
    if (V < V_star) {
      V_star = V; iftCopyFloatArray(theta_star,theta_plus,n);
    }
    iter++;
    acquisition = fabs(V_star - V0);
  } while (( acquisition > msps->stopcriteria) &&(iter < msps->niters));

  if (iter >= msps->niters)
    iftWarning("Reached maximum number of iterations","iftMSPSMin");


  if(msps->verbose) {
    printf("Final fitness value %f, parameters: ", V_star);
    for (i = 0; i < n; i++)
      printf("%f ", theta_star[i]);
    printf("\n");
  }

  iftCopyFloatArray(theta,theta_star,n);

  iftFree(theta_star);
  iftFree(Vp);
  iftFree(Dp);
  iftFree(Ds);
  iftFree(Dplus);
  iftFree(Dminus);
  iftFree(D_star);
  iftFree(theta_plus);
  iftFree(theta_minus);

  return(V0);
}

