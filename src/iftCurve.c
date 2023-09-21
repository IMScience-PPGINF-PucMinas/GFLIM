#include "iftCurve.h"

#include "ift/core/io/Stream.h"


/* ------------------------ 2D Curve ----------------------------- */

iftCurve2D *iftCreateCurve2D(int npts)
{
  iftCurve2D *curve = (iftCurve2D *) iftAlloc(1,sizeof(iftCurve2D));

  curve->npts = npts;
  curve->x    = iftAllocDoubleArray(npts);
  curve->y    = iftAllocDoubleArray(npts);

  return(curve);
}

void      iftDestroyCurve2D(iftCurve2D **curve)
{
  iftCurve2D *tmp = *curve;

  iftFree(tmp->x);
  iftFree(tmp->y);
  iftFree(tmp);
  *curve = NULL;

}

iftCurve2D *iftCopyCurve2D(iftCurve2D *curve) 
{
  iftCurve2D *ccurve = iftCreateCurve2D(curve->npts);
  
  for (int i=0; i < curve->npts; i++) {
    ccurve->x[i] = curve->x[i];
    ccurve->y[i] = curve->y[i];
  }
  
  return(ccurve);
}

void iftMeansOf2DPoints(iftCurve2D *curve, double *mean_of_x, double *mean_of_y)
{
  *mean_of_x = 0.0;
  *mean_of_y = 0.0;

  for (int i=0; i < curve->npts; i++) {
    *mean_of_x += curve->x[i];
    *mean_of_y += curve->y[i];
  }

  *mean_of_x /= curve->npts;
  *mean_of_y /= curve->npts;
}

void iftStdevsOf2DPoints(iftCurve2D *curve, double mean_of_x, double mean_of_y, double *stdev_of_x, double *stdev_of_y)
{
  *stdev_of_x = 0.0;
  *stdev_of_y = 0.0;

  for (int i=0; i < curve->npts; i++) {
    *stdev_of_x += (curve->x[i]-mean_of_x)*(curve->x[i]-mean_of_x);
    *stdev_of_y += (curve->y[i]-mean_of_y)*(curve->y[i]-mean_of_y);
  }

  *stdev_of_x = sqrt(*stdev_of_x / (curve->npts-1));
  *stdev_of_y = sqrt(*stdev_of_y / (curve->npts-1));
}

iftCurve2D *iftCentralizeCurve2D(iftCurve2D *curve, double mean_of_x,  double mean_of_y)
{
  iftCurve2D *ccurve = iftCopyCurve2D(curve);

  for (int i=0; i < ccurve->npts; i++) {
      ccurve->x[i] = (ccurve->x[i] - mean_of_x); 
      ccurve->y[i] = (ccurve->y[i] - mean_of_y);
  }

  return(ccurve);
}

double iftCorrCoefficient(iftCurve2D *curve, double mean_of_x, double mean_of_y) 
{
  iftCurve2D *ccurve = iftCentralizeCurve2D(curve,mean_of_x,mean_of_y);
  double sum_of_xy=0.0, sum_of_x2 = 0.0, sum_of_y2 = 0.0;
  
  for (int i=0; i < ccurve->npts; i++) {
    sum_of_xy += ccurve->x[i]*ccurve->y[i];
    sum_of_x2 += ccurve->x[i]*ccurve->x[i];
    sum_of_y2 += ccurve->y[i]*ccurve->y[i];
  }
  iftDestroyCurve2D(&ccurve); 

  return(sum_of_xy/sqrt(sum_of_x2*sum_of_y2));
}

void  iftLinearRegression2D(iftCurve2D *curve, float *slope, float *intercept, double *sqr_corr_coef)
{
  double mean_of_x,  mean_of_y, stdev_of_x, stdev_of_y, corr_coef;
  
  iftMeansOf2DPoints(curve,&mean_of_x,&mean_of_y);
  iftStdevsOf2DPoints(curve,mean_of_x,mean_of_y,&stdev_of_x,&stdev_of_y);
  corr_coef = iftCorrCoefficient(curve,mean_of_x,mean_of_y);

  if (fabs(stdev_of_x) > IFT_EPSILON)
    *slope = (float)(corr_coef * (stdev_of_y / stdev_of_x)); 
  else
    *slope = IFT_INFINITY_FLT;

  if (*slope != IFT_INFINITY_FLT)
    *intercept = (float)(mean_of_y - ((*slope) * mean_of_x));
  else
    *intercept = 0.0;

  *sqr_corr_coef = corr_coef*corr_coef;

}

void iftWriteCurve2D(iftCurve2D *curve, char *filename)
{
  FILE *fp = fopen(filename,"w");

  for (int i=0; i < curve->npts; i++) 
    fprintf(fp,"%lf %lf\n",curve->x[i],curve->y[i]);
  
  fclose(fp);
}

/* -------------------- 3D Curve ------------------------------- */

iftCurve   *iftCreateCurve(int npts)
{
  iftCurve *curve = (iftCurve *)iftAlloc(1,sizeof(iftCurve));

  curve->pt   = (iftPoint *)iftAlloc(npts,sizeof(iftPoint));
  curve->npts = npts;

  return(curve);
}

void        iftDestroyCurve(iftCurve **curve)
{
  iftFree((*curve)->pt);
  iftFree(*curve);
  (*curve)=NULL;
}

iftCurve   *iftReadCurve(char *filename)
{
  iftCurve *curve;
  int       npts;
  FILE     *fp = fopen(filename,"r");
  
  if (fscanf(fp,"%d",&npts)!=1) iftError(MSG_FILE_FORMAT_ERROR, "iftReadCurve");
  curve = iftCreateCurve(npts);
  for (int i=0; i < npts; i++) 
    if (fscanf(fp,"%f %f %f",&curve->pt[i].x,&curve->pt[i].y,&curve->pt[i].z)!=3)
        iftError(MSG_FILE_FORMAT_ERROR, "iftReadCurve");
    
  fclose(fp);

  return(curve);
}

void iftWriteCurve(iftCurve *curve, char *filename)
{
  FILE     *fp = fopen(filename,"w");
  
  fprintf(fp,"%d\n",curve->npts);
  for (int i=0; i < curve->npts; i++) 
    fprintf(fp,"%f %f %f\n",curve->pt[i].x,curve->pt[i].y,curve->pt[i].z);
    
  fclose(fp);
}

void iftPrintCurve(iftCurve *curve)
{
  printf("%d\n",curve->npts);
  for (int i=0; i < curve->npts; i++) 
    printf("%f %f %f\n",curve->pt[i].x,curve->pt[i].y,curve->pt[i].z);
}

