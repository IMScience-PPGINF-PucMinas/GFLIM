/**
 * @file
 * @brief  Fast Fourier Transform and spectrum manipulation functions.
 * @author Alexandre Falcao
 * @date   Feb 5th 2022
 */

#ifndef IFT_SPECTRUM_H
#define IFT_SPECTRUM_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftImage.h"
#include "iftMImage.h"
  
void       iftFFT1D(int dir, long nn, double *x, double *y);
iftMImage *iftFFT2D(iftMImage *mimg);
iftMImage *iftFFT3D(iftMImage *mimg); 
iftMImage *iftInvFFT2D(iftMImage *spec);
iftMImage *iftInvFFT3D(iftMImage *spec);
// Adding my ViewMagnitude
iftImage  *iftViewMagnitude(iftMImage *spec);


iftImage  *iftViewLogMagnitude(iftMImage *spec);
iftImage  *iftViewPhase(iftMImage *spec);
iftMImage *iftMultSpectra(iftMImage *spec, iftMImage *filt);
iftMImage *iftFFTImageToMImage(iftImage *img);
iftImage  *iftFFTMImageToImage(iftMImage *mimg, iftImage *orig);

/**
 * @brief Creates a complex Gabor Filter 2D as a iftMIMage with same dimensions as
 * spec. The Gabor parameters are available at section 2.3 of https://inc.ucsd.edu/mplab/75/media//gabor.pdf
 * 
 * @date Mar 22, 2022
 *
 * @param u0 X axis spatial frequency
 * @param v0 Y axis spatial frequency
 * @param x0 X location of gaussian peak
 * @param y0 Y location of gaussian peak
 * @param a Scale X axis of gaussian
 * @param b Scale Y axis of gaussian

 * @return The new Gabor iftMImage.
 */
//! swig(newobject)

iftMImage *iftGaborFilter2D(iftMImage *spec, float u0, float v0, float P, float x0, float y0, float a, float b, float theta);
  /* n is a low number in [1,5] and Dh in (0,1) is a percentage of the
     semi-diagnonal of the image */  
  iftMImage *iftButterworthHighPass(iftMImage *spec, int n, float Dh);
  /* n is a low number in [1,5] and Dl in (0,1) is a percentage of the
     semi-diagnonal of the image */  
  iftMImage *iftButterworthLowPass(iftMImage *spec, int n, float Dl);
  /* n is a low number in [1,5] and Dl > Dh in (0,1) are percentages of the
     semi-diagnonal of the image */  
  iftMImage *iftButterworthBandPass(iftMImage *spec, int n, float Dl, float Dh);
  /* n is a low number in [1,5] and Dl < Dh in (0,1) are percentages of the
     semi-diagnonal of the image */  
  iftMImage *iftButterworthBandReject(iftMImage *spec, int n, float Dl, float Dh);

#ifdef __cplusplus
}
#endif


#endif


