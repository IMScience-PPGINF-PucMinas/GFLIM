//
// Created by Samuel Martins on 2019-01-04.
//

#ifndef IFT_SET_H
#define IFT_SET_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"

/**
 * Set of integer elements. Functions as a single linked list
 */
//! swig(extend = iftSetExt.i, destroyer = iftDestroySet)
typedef struct ift_set {
    int elem;
    struct ift_set *next;
} iftSet;


/**
 * @brief Print elements in the set
 * @param S A linked list
 */
void    iftPrintSet(iftSet *S);

/**
 * @brief Inserts an element in the begining of the linked list
 * @param S A linked list
 * @param elem The element to be inserted
 */
void    iftInsertSet(iftSet **S, int elem);

/**
 * @brief Removes and returns the first element of the linked list
 * @param S A linked list
 * @return The first element
 */

//! swig(newobject)
int     iftRemoveSet(iftSet **S);

void    iftRemoveSetElem(iftSet **S, int elem);
void    iftDestroySet(iftSet **S);
iftSet* iftSetUnion(iftSet *S1,iftSet *S2);

//! swig(newobject)
iftSet* iftSetConcat(iftSet *S1,iftSet *S2);
char    iftUnionSetElem(iftSet **S, int elem);
void    iftInvertSet(iftSet **S);
void iftReverseSet(iftSet **S);

/**
 * @brief Returns the size of an iftSET.
 */
int 	iftSetSize(const iftSet* S);
iftSet* iftSetCopy(iftSet* S);
int     iftSetHasElement(iftSet *S, int elem);

iftSet* iftSetCopyOrdered(const iftSet *S);

/**
 * @brief Converts an iftSet to a 1D-array.
 *
 * @author Samuel Martins
 * @date September 14, 2015
 *
 * @param S iftSet to be copied.
 * @return The resulting integer 1D-array.
 */

//! swig(newobject)
iftIntArray *iftSetToArray(iftSet *S);

/**
 * @brief Converts a 1D-array to a iftSet.
 *
 * @author Samuel Martins
 * @date September 14, 2015
 *
 * @param array The integer 1D-array to be copied.
 * @param n_elems Array size.
 * @return The resulting iftSet.
 */
iftSet *iftIntArrayToSet(iftIntArray *array);


#ifdef __cplusplus
}
#endif

#endif //IFT_SET_H
