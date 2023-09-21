#ifndef IFT_ADJSET_H_
#define IFT_ADJSET_H_

/**
 * @file
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftMatrix.h"

/**
 * Defines an adjacent node set of a node in a not complete graph
 */
typedef struct ift_adjset {
  int   node;   /* current node */
  float arcw;   /* arc weight respecting to the current node edge*/
  struct ift_adjset *next; /* next adjacent node */
} iftAdjSet;

/**
 * @brief Insert a node into the beggining of an adjacent node set
 * @param S A node adjacent set
 * @param node The node to be inserted
 * @param arcw The arc value
 */
void iftInsertAdjSet(iftAdjSet **S, int node, float arcw);
int  iftRemoveAdjSet(iftAdjSet **S, float *arcw);
void iftRemoveAdjSetNode(iftAdjSet **S, int node);
void iftDestroyAdjSet(iftAdjSet **S);

// Uses recursion to copy the set in order
iftAdjSet* iftAdjSetCopyOrdered(iftAdjSet *S);


/**
 * @brief Get the size of an adjacency set.
 * @author Samuel Martins
 * @date May 15, 2018
 */
int iftAdjSetSize(const iftAdjSet *S);

/**
 * @brief Convert an Adjacency Set into a Double Matrix [n, 2].
 * 
 * Each matrix row is a node from the AdjSet, containing 2 columns.
 * 
 * @author Samuel Martins
 * @date May 15, 2018
 */
iftMatrix *iftAdjSetToMatrix(const iftAdjSet *S);

/**
 * @brief Convert a Matrix [n, 2] into an Adjacency Set.
 * 
 * Each matrix row is a node from the AdjSet, containing 2 columns.
 * 
 * @author Samuel Martins
 * @date May 15, 2018
 */
iftAdjSet *iftMatrixToAdjSet(const iftMatrix *M);

int iftAdjSetHasElement(iftAdjSet *S, int elem);

void iftReverseAdjSet(iftAdjSet **S);

#ifdef __cplusplus
}
#endif

#endif

