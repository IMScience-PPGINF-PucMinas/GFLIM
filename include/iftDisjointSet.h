#ifndef IFT_DISJOINT_SET_H_
#define IFT_DISJOINT_SET_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"

/**
 * Disjoint-set data-structure
 * @author Adan Echemendia
 * @date Jun 10, 2016
 */
typedef struct{
    /** List of set representatives. Ex. parent[i] corresponds to the representative of the set containing item i */
    int *parent;
    /** List of set ranks. Ex. rank[i] corresponds to the rank of the set containing item i */
    int *rank;
    /** Number of sets */
    int n; 
} iftDisjointSet;

/**
 * @brief Creates a disjoint-set data-structure with n sets
 * @details [long description]
 * 
 * @param n Number of sets
 * @return Disjoint-set data-structure
 */
iftDisjointSet* iftCreateDisjointSet(int n);

/**
 * @brief Destroys a disjoint-set data-structure
 * @details [long description]
 * 
 * @param dset Pointer to pointer to a disjoint-set data-structure
 */
void iftDestroyDisjointSet(iftDisjointSet **dset);

/**
 * @brief Joins the sets containing elem1 and elem2 and returns the representative of this new set
 * @details [long description]
 * 
 * @param dset Disjoint-set data-structure
 * @param elem1 First item
 * @param elem2 Second item
 * @return Representative of the resulting set
 */
int iftDisjointSetUnion(iftDisjointSet* dset, int elem1, int elem2);


/**
 * @brief Finds the representative for the set containing elem
 * 
 * @param dset Disjoint-set data-structure
 * @param elem Searched item
 * 
 * @return Representative for the set containing elem
 */
int iftDisjointSetFind(iftDisjointSet* dset, int elem);

#ifdef __cplusplus
}
#endif

#endif

