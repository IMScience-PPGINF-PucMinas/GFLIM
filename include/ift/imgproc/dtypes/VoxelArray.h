//
// Created by Cesar Castelo on Mar 21, 2019.
//

#ifndef IFT_VOXELARRAY_H
#define IFT_VOXELARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/io/Stream.h"

/**
 * @brief Array of Voxels.
 * @author Samuka Martins
 * @date Aug 23, 2018
 * @ingroup VoxelArray
 */
//! swig(destroyer = iftDestroyVoxelArray, extend = VoxelArrayExt.i)
typedef struct ift_voxel_arr {
    /** Number of Elements */
    long n;
    /** Array of Voxels */
    iftVoxel *val;
} iftVoxelArray;

/**
 * @brief Creates an array of Voxels.
 * 
 * @param n Number of elements.
 * @return Array of Voxels
 * 
 * @author Samuka Martins
 * @date Aug 23, 2018
 */
iftVoxelArray *iftCreateVoxelArray(long n);

/**
 * @brief Destroys a Voxel Array.
 * 
 * @param varr Voxel array
 * 
 * @author Samuka Martins
 * @date Aug 23, 2018
 */
void iftDestroyVoxelArray(iftVoxelArray **varr);

//! swig()
int iftVoxelArrayFurthestPair(const iftVoxelArray *a, const iftVoxelArray *b);


void iftInsertVoxel(iftVoxelArray *arr, int index, iftVoxel *u);


#ifdef __cplusplus
}
#endif

#endif //IFT_VOXELARRAY_H
