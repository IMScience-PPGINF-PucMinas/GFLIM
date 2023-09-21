//
// Created by Samuel Martins on 10/12/18.
//

#ifndef IFT_RANGE_H
#define IFT_RANGE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/DblArray.h"
#include "ift/core/dtypes/IntArray.h"

/**
 * @brief Creates an integer array in the specified range.
 *
 * @note Inpired in range function from python.
 *
 * @param begin Starting value
 * @param end Ending value (inclusive)
 * @param inc Increment ammount
 * @return Range Integer Array
 *
 * @author Samuka
 * @date Jan 3, 2017
 */
iftIntArray *iftIntRange(int begin, int end, int inc);


/**
 * @brief Creates an integer array by repeating <n_repetitions> the value <x>.
 *
 * @note The size of the integer array will be <n_repetitions>
 *
 * @param x Value to be repeated.
 * @param n_repetitions Number of repetitions.
 * @return  Integer array with <n_repetitions> of the value <x>.
 *
 * @author Samuka
 * @date Aug 23, 2018
 */
iftIntArray *iftIntRepeat(int x, int n_repetitions);


/**
 * @author Peixinho
 * @date Nov, 2016
 * @brief Creates a double array in the specified range
 * @param begin Starting value
 * @param end Ending value (inclusive)
 * @param inc Increment ammount
 * @return
 */
iftDblArray *iftRange(double begin, double end, double inc);


/**
 * @author Peixinho
 * @date Nov, 2016
 * @brief Creates a double array in the specified geometric range
 * @param begin Starting value
 * @param end Ending value (not inclusive)
 * @param inc Geommetric increment ammount
 * @return
 */
iftDblArray *iftGeomRange(double begin, double end, double mul);


#ifdef __cplusplus
}
#endif

#endif //IFT_RANGE_H
