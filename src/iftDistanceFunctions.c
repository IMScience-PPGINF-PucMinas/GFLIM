#include "iftDistanceFunctions.h"

#include "ift/core/io/Stream.h"


double iftComputeL1Norm(float* vector1,float* vector2, size_t  numberElements, iftArgumentList* argumentList){
    iftUnusedParameter(argumentList);
    float diff = 0;
    double total = 0;
    for (size_t i = 0; i < numberElements; ++i) {
        diff = vector1[i] - vector2[i];
        diff = (diff < 0) ? -diff : diff;
        total += diff;
    }
    return  total;
}

double iftComputeNormalizedL1Norm(float* vector1,float* vector2, size_t  numberElements,iftArgumentList* argumentList){
    iftUnusedParameter(argumentList);
    float diff = 0;
    double total = 0;
    for (size_t i = 0; i < numberElements; ++i) {
        diff = vector1[i] - vector2[i];
        diff = (diff < 0) ? -diff : diff;
        total += diff;
    }
    total /= (numberElements+0.000001);
    return  total;
}

double iftComputeL2Norm(float* vector1,float* vector2, size_t  numberElements,iftArgumentList* argumentList){
    iftUnusedParameter(argumentList);
    float diff = 0;
    double total = 0;
    for (size_t i = 0; i < numberElements; ++i) {
        diff = vector1[i] - vector2[i];
        diff = (diff*diff);
        total += diff;
    }
    total = sqrt(total);
    return  total;
}
//
double iftComputeNormalizedL2Norm(float* vector1,float* vector2, size_t  numberElements,iftArgumentList* argumentList){
    iftUnusedParameter(argumentList);
    float diff = 0;
    double total = 0;
    for (size_t i = 0; i < numberElements; ++i) {
        diff = vector1[i] - vector2[i];
        diff = (diff*diff);
        total += diff;
    }
    total /= (numberElements+0.000001);
    total = sqrt(total);
    return  total;
}
//
//typedef double (*InverseDistanceFunction)(double distance);
//
//inline double computeInverseDistanceUsingExponential(double distance){
//    return exp(-distance);
//}
//
//inline double computeInverseDistanceUsingInverseX(double distance){
//    return 1.0/(distance + 0.00000001);
//}


float iftEarthMoversDistance(float *vector1, float* vector2, size_t num) {

    float emd=0;
    float *arr = (float*) iftAllocFloatArray(num+1);

    arr[0]=0;
    for (int i=1; i<num+1;i++)
        arr[i] = vector1[i-1] + arr[i-1] - vector2[i-1];
    for (int i=0; i<num+1;i++)
        emd += arr[i];

    iftFree(arr);
    return(emd);
}

