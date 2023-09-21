//
// Created by Samuel Martins on 2019-01-04.
//

#ifndef IFT_DHEAP_H
#define IFT_DHEAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


#define MINVALUE   0 /* define heap to remove node with minimum value */
#define MAXVALUE   1 /* define heap to remove node with maximum value */

#define iftDad(i) ((i - 1) / 2)
#define iftLeftSon(i) (2 * i + 1)
#define iftRightSon(i) (2 * i + 2)
#define iftSetRemovalPolicyDHeap(a,b) a->removal_policy = b

/**
 * @brief Struct that defines a Heap of double values
 * @author Samuel Martins (documentation)
 * @date Jan 4, 2019
 */
typedef struct ift_dheap {
    /** Pointer to an array o double values */
    double *value;
    /** Color of each Heap element */
    char  *color;
    /** Node of each Heap element */
    int   *node;
    /** Position of each Heap element */
    int   *pos;
    /** Index of the last Heap element */
    int    last;
    /** Number of elements from the heap. */
    int    n;
    /** Removal Policy */
    char removal_policy; /* 0 is MINVALUE and 1 is MAXVALUE */
} iftDHeap;

/**
 * @brief Allocates a Heap of <n> Double values.
 *
 * @warning The input array of values passed is only assigned to the corresponding pointer of the Heap. Be careful.
 *
 * @param n Number of elements.
 * @param value Array of double values whose reference will be stored by the heap.
 * @return Double Heap.
 *
 * @author Samuel Martins (documentation)
 * @date Jan 4, 2019
 */
iftDHeap *iftCreateDHeap(int n, double *value);
void      iftDestroyDHeap(iftDHeap **H);
char      iftFullDHeap(iftDHeap *H);
char      iftEmptyDHeap(iftDHeap *H);
char      iftInsertDHeap(iftDHeap *H, int pixel);
int       iftRemoveDHeap(iftDHeap *H);
void      iftRemoveDHeapElem(iftDHeap *H, int pixel);
void      iftGoUpDHeap(iftDHeap *H, int i);
void      iftGoDownDHeap(iftDHeap *H, int i);
void      iftResetDHeap(iftDHeap *H);

void  iftGoUpDHeap_(iftDHeap *H, int i, bool show);

#ifdef __cplusplus
}
#endif

#endif //IFT_DHEAP_H
