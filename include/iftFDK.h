#ifndef FDK_H
#define FDK_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftFImage.h"

int iftFDKConeBeam(   
	double *proj,                     //projection data
	int numProj,                      //number of projections
	int nSignalXPts,                  //pixels detector
	int nSignalYPts,                  //lines detector
	double dSignalXInc,               //pixel spacing
	double dSignalYInc,               //line spacing
	double AOR,			  //column on detector that
	                                  // corresponds to the axis of rotation
	double midPlane,                  //row on detector that corresponds to the midplane
	double dFocalLength,              //source to actual AOR distance
	double dSourceDetectorLength,     //source to detector distance
	double theta[356],                //angle data
	double *vox,                      //slice data
	int numSlices,                    //number of slices,
	int pixels,			  //cols per slice
	int lines,                        //rows per slice
	double subSample,                 // sub sampling factor
	int extProj,                      //extend projection by extProj*nSignalXPts 
	int bWindow                       //use smoothing window?
	);


#ifdef __cplusplus
}
#endif

#endif
