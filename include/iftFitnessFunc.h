/**
 * @file
 * @brief Fitness Function definition file.
 * @author Peixinho
 * @date 28 July 2015
 */
#ifndef IFT_IFTFITNESSFUNC_H
#define IFT_IFTFITNESSFUNC_H

/* @brief Fitness function used for optimization problems */
typedef float (*iftFitnessFunc)   (void *problem, float *theta);

/**
 * @brief Fitness function used for optimization problems using double precision. See also iftFitnessFunc
 */
typedef double (*iftDFitnessFunc)   (void *problem, double *theta);

#endif //IFT_IFTFITNESSFUNC_H
