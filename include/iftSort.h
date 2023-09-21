#ifndef IFT_SORT_H_
#define IFT_SORT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"

/* Suppose you have a list of objects, each with a given value. These
   methods can sort the indices of the objects in that list in the
   increasing (decreasing) order of values. They will return in index,
   the sorted indices, and in value, the sorted values. Note that
   value[i] will no longer correspond to value[index[i]] after
   sorting. The objects can be accessed in order by index[i] (e.g.,
   object[index[i]]), but their values must remain the same (e.g.,
   object[index[i]].value). The sorted list of values (i.e.,
   value[i]=object[index[i]].value), on the other hand, can be used
   when we want to sort values only, with no underlying object
   concept. In such a case, index[i] has no meaning. See examples in
   demo/iftSortObjects.c. 

*/

void iftBucketSort(int *value, int *index, int nelems, uchar order);
void iftFHeapSort(float *value, int *index, int nelems, uchar order);

/**
 * @brief Sort an integer array.
 * @param value char array
 * @param index index on the original array
 * @param i0 starting position (should be [0,n))
 * @param i1 ending position (should be [0,n))
 * @param order IFT_INCREASING or IFT_DECREASING order.
 */
void iftQuickSort( int *value, int *index, int i0, int i1, uchar order ); 
void iftSQuickSort( char **value, int *index, int i0, int i1, uchar order, int size ); 
void iftFQuickSort( float *value, int *index, int i0, int i1, uchar order); 
void iftDQuickSort( double *value, int *index, int i0, int i1, uchar order); 

/**
 * @brief   Selects the value of the kth element in an unsorted array of values
 * @author  Alexandre Falcão (this is the quick-selection algorithm from Cormen)
 * @date    February, 10th 2017
 * @ingroup 
 *
 * @param value   Array of values 
 * @param nelems  Size of the array 
 * @param k       The desired value of k
 * @return The value of the kth element
 */

int iftSelectTheKthElem(int *value, int nelems, int k);

/**
 * @brief   Selects the index (e.g., index of a voxel) of the kth element in an unsorted array of values
 * @author  Alexandre Falcão (this is the quick-selection algorithm from Cormen)
 * @date    February, 10th 2017
 * @ingroup 
 *
 * @param value   Array of values 
 * @param index   Array of indices 
 * @param nelems  Size of the array 
 * @param k       The desired value of k
 * @return The value of the kth element
 */
  
  int iftSelectIndexOfTheKthElem(int *value, int *index, int nelems, int k);

  
#ifdef __cplusplus
}
#endif

#endif
