//
// Created by peixinho on 6/23/15.
//

#include "iftFunctions.h"

#include "ift/core/io/Stream.h"
#include <math.h>

float iftSphereFunction(float *x, int n) {

    float sum = 0.0;
    int i;

    for (i = 0; i<n; ++i) {
        sum += x[i]*x[i];
    }

    return sum;
}

float iftRosembrockFunction(float *x, int n) {

    float sum = 0.0;
    int i;

    for (i = 0; i < n-1; ++i) {
        sum += 100.0*pow(x[i+1] - x[i]*x[i], 2.0) + pow(x[i]-1, 2.0);/** We could use multiplication instead of pow function, but speed is not a major problem here. */
    }

    return sum;
}

//TODO: check if its correct
float iftAckleyFunction(float *x, int n) {

    int i;
    float sum1 = 0.0, sum2 = 0.0;
    float term1, term2;
    float a = 20.0, b = 0.2, c = 2*M_PI;/** Common values for parameters a, b and c */

    for (i = 0; i < n; ++i) {
        sum1 += x[i]*x[i];
        sum2 += cos(c*x[i]);
    }

    term1 = -a*exp(-b*sqrt(sum1/n));
    term2 = -exp(sum2/n);

    return term1 + term2 + a + M_E;
}

float iftGrienwankFunction(float *x, int n) {

    int i;
    float sum = 0.0, prod = 1.0;

    for (i = 0; i < n; ++i) {
        sum += x[i]*x[i]/4000.0;
        prod *= cos(x[i]/sqrt(i+1));
    }

    return sum - prod + 1.0;
}

float iftRastringinFunction(float *x, int n) {
    float sum = 0.0;
    int i;

    for (i = 0; i < n; ++i) {
        sum += x[i]*x[i] - 10*cos(2*M_PI*x[i]);
    }

    return 10*n + sum;
}
