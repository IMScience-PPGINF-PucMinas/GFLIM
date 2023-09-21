//
// Created by deangeli on 5/16/17.
//

#ifndef _DISTANCEFUNCTIONS_H_
#define _DISTANCEFUNCTIONS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftArgumentList.h"

typedef double (*iftDistanceFunction)(float* vector1, float* vector2, size_t numberElements,iftArgumentList* argumentList);

double iftComputeL1Norm(float* vector1,float* vector2, size_t  numberElements, iftArgumentList* argumentList);

double iftComputeNormalizedL1Norm(float* vector1,float* vector2, size_t  numberElements,iftArgumentList* argumentList);

double iftComputeL2Norm(float* vector1,float* vector2, size_t  numberElements,iftArgumentList* argumentList);

double iftComputeNormalizedL2Norm(float* vector1,float* vector2, size_t  numberElements,iftArgumentList* argumentList);



float iftEarthMoversDistance(float *vector1, float* vector2, size_t num);



//typedef double (*InverseDistanceFunction)(double distance);
//
//inline double computeInverseDistanceUsingExponential(double distance){
//    return exp(-distance);
//}
//
//inline double computeInverseDistanceUsingInverseX(double distance){
//    return 1.0/(distance + 0.00000001);
//}

#ifdef __cplusplus
}
#endif

#endif //LIBFL_DISTANCEFUNCTIONS_H
