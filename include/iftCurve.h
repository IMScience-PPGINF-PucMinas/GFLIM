#ifndef IFT_Curve_H_
#define IFT_Curve_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <iftCommon.h>

typedef struct ift_curve2D {
  double *x, *y;
  int     npts;
} iftCurve2D;

typedef struct ift_curve {
  iftPoint *pt;
  int       npts;
} iftCurve;

  /* ----------------- 2D curve ---------------------------- */

iftCurve2D *iftCreateCurve2D(int npts);
void        iftDestroyCurve2D(iftCurve2D **curve);
iftCurve2D *iftCopyCurve2D(iftCurve2D *curve);
void        iftMeansOf2DPoints(iftCurve2D *curve, double *mean_of_x, double *mean_of_y);
void        iftStdevsOf2DPoints(iftCurve2D *curve, double mean_of_x, double mean_of_y, double *stdev_of_x, double *stdev_of_y);
iftCurve2D *iftCentralizeCurve2D(iftCurve2D *curve, double mean_of_x,  double mean_of_y);
double      iftCorrCoefficient2D(iftCurve2D *curve, double mean_of_x, double mean_of_y); 
void        iftLinearRegression2D(iftCurve2D *curve, float *slope, float *intercept, double *sqr_corr_coef);
void        iftWriteCurve2D(iftCurve2D *curve, char *filename);

  /* ----------------- 3D curve ---------------------------- */

iftCurve   *iftCreateCurve(int npts);
void        iftDestroyCurve(iftCurve **curve);
iftCurve   *iftReadCurve(char *filename);
void        iftWriteCurve(iftCurve *curve, char *filename);
void        iftPrintCurve(iftCurve *curve);


#ifdef __cplusplus
}
#endif

#endif
