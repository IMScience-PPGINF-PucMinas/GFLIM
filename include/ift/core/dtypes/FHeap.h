//
// Created by Samuel Martins on 2019-01-04.
//

#ifndef IFT_FHEAP_H
#define IFT_FHEAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


#define IFT_MINVALUE   0 /* define heap to remove node with minimum value */
#define IFT_MAXVALUE   1 /* define heap to remove node with maximum value */

#define iftDad(i) ((i - 1) / 2)
#define iftLeftSon(i) (2 * i + 1)
#define iftRightSon(i) (2 * i + 2)
#define iftSetRemovalPolicyFHeap(a,b) a->removal_policy = b

typedef struct ift_fheap {
    float *value; /* values of each node*/
    char  *color; /* colors of each node (IFT_WHITE=never inserted, IFT_GRAY=inserted,
  IFT_BLACK=removed)*/
    int   *node;  /* array of nodes of the binary heap */
    int   *pos;   /* positions of each node in the array of the binary heap.*/
    int    last;  /* quantity of inserted elements in the binary heap minus one*/
    int    n;     /* size of the binary heap */
    char removal_policy; /* removal policy used in the binary heap. 0 is MINVALUE and 1 is MAXVALUE */
} iftFHeap;

/**
 * @brief Creates a binary heap
 * @param n Number of elements
 * @param value Values of elements
 * @return A binary heap
 */
iftFHeap *iftCreateFHeap(int n, float *value);

/**
 * @brief Destroys a binary heap
 * @param H The binary heap
 */
void   iftDestroyFHeap(iftFHeap **H);

/**
 * @brief Checks if a binary heap is full
 * @param H The binary heap
 * @return 1 if it's full, 0 otherwise
 */
char   iftFullFHeap(iftFHeap *H);

/**
 * @brief Checks if a binary heap is empty
 * @param H A binary heap
 * @return 1 if it's empty, 0 otherwise
 */
char   iftEmptyFHeap(iftFHeap *H);

/**
 * @brief Inserts an element into a binary heap
 * @param H A binary heap
 * @param pixel The position of the element to be inserted, the value of this element is found in H->val[pixel]
 * @return 1 if the element was inserted successfully, 0 otherwise
 */
char   iftInsertFHeap(iftFHeap *H, int pixel);

/**
 * @brief Returns and removes the top element in the binary heap
 * @param H A binary heap
 * @return The top element of the binary heap.
 */
int    iftRemoveFHeap(iftFHeap *H);

/**
 * @brief Removes the <br> pixel </br> element from the binary heap
 * @param H A binary heap
 * @param pixel The element to be removed
 */
void   iftRemoveFHeapElem(iftFHeap *H, int pixel);

/**
 * @brief Goes up an element in a binary heap
 * @param H A binary heap
 * @param i The element to be promoted
 */
void   iftGoUpFHeap(iftFHeap *H, int i);

/**
 * @brief Goes down an element in a binary heap
 * @param H A binary heap
 * @param i The element that goes down
 */
void   iftGoDownFHeap(iftFHeap *H, int i);

/**
 * @brief Resets a binary heap to its default values
 * @param H A binary heap
 */
void   iftResetFHeap(iftFHeap *H);


#ifdef __cplusplus
}
#endif

#endif //IFT_FHEAP_H
