/**
 * @file iftDList.h
 * @brief Definitions and functions of Doubly Linked List of doubles.
 * @author Samuel Martins
 * @date Sep 19, 2015
 *
 * @note Programs:
 * @ref iftTestDList.c (demo/Miscellaneous/iftTestDList.c)
 *
 */

#ifndef IFT_DLIST_H
#define IFT_DLIST_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/DblArray.h"


/**
 * @brief Node of the Doubly Linked List of Doubles.
 * @author Samuel Martins
 * @date Sep 19, 2015
 */
typedef struct _ift_dnode {
    double elem;
    struct _ift_dnode *previous;
    struct _ift_dnode *next;
} iftDNode;


/**
 * @brief Doubly Linked List of Doubles.
 * @author Samuel Martins
 * @date Sep 19, 2015
 */
typedef struct _ift_dlist {
    /** Number of Nodes of the Doubly Linked List of Doubles */
    int n;
    /** Pointer to the Head (first) Node */
    iftDNode *head;
    /** Pointer to the Tail (last) Node */
    iftDNode *tail;
} iftDList;



/**
 * @brief Creates an iftDList.
 * @author Samuel Martins
 * @date Sep 19, 2015
 */
iftDList *iftCreateDList();


/**
 * @brief Destroys an iftDList.
 * @author Samuel Martins
 * @date Sep 19, 2015
 *
 * @param DL Reference to the iftDList to be destroyed.
 */
void iftDestroyDList(iftDList **DL);


/**
 * @brief Insert a double into the head of the iftDList.
 *
 * @author Samuel Martins
 * @date Sep 19, 2015
 *
 * Insert the double <elem> into the head of the iftDList.
 * If the iftDList is not allocated, an error will be raised.
 *
 * @param DL The iftDList.
 * @param elem The double to be inserted into the head of the iftDList.
 */
void iftInsertDListIntoHead(iftDList *DL, double elem);


/**
 * @brief Insert a double into the tail of the iftDList.
 *
 * @author Samuel Martins
 * @date Sep 19, 2015
 *
 * Insert the double <elem> into the tail of the iftDList.
 * If the iftDList is not allocated, an error will be raised.
 *
 * @param DL The iftDList.
 * @param elem The double to be inserted into the tail of the iftDList.
 */
void iftInsertDListIntoTail(iftDList *DL, double elem);


/**
 * @brief Removes the Head element from the iftDList and removes its node.
 *
 * @author Samuel Martins
 * @date Sep 19, 2015
 *
 * Removes the Head double node of the iftDList and removes its node.
 * If there is no element, the NIL value is returned.
 * If the iftDList is not allocated, an error will be raised.
 *
 * @param DL The iftDList.
 * @return The Head double from the iftDList.
 */
double iftRemoveDListHead(iftDList *DL);


/**
 * @brief Removes the Tail element from the iftDList and removes its node.
 *
 * @author Samuel Martins
 * @date Sep 19, 2015
 *
 * Removes the Tail double node of the iftDList and removes its node.
 * If there is no element, the NIL value is returned.
 * If the iftDList is not allocated, an error will be raised.
 *
 * @param DL The iftDList.
 * @return The Tail double from the iftDList.
 */
double iftRemoveDListTail(iftDList *DL);


/**
 * @brief Prints an iftDList (from the Node to the Tail node).
 *
 * @author Samuel Martins
 * @date Sep 19, 2015
 *
 * Prints an iftDList from the Head to Tail node.
 *
 * @param DL The iftDList.
 */
void iftPrintDList(const iftDList *DL);


/**
 * @brief Prints an iftDList in Reveresed Order (from the Tail to the Head node).
 *
 * @author Samuel Martins
 * @date Sep 19, 2015
 *
 * Prints an iftDList from the Tail to Head node.
 *
 * @param DL The iftDList.
 */
void iftPrintReversedDList(const iftDList *DL);



/**
 * @brief Converts a DList to a Double Array.
 * @author Samuka Martins
 * @date July 12, 2019
 */
iftDblArray *iftDListToDblArray(const iftDList *DL);


/**
 * @brief Converts a Double Array to a (Double) List.
 * @author Samuka Martins
 * @date July 12, 2019.
 */
iftDList *iftDblArrayToDList(const iftDblArray *arr);


#ifdef __cplusplus
}
#endif

#endif //IFT_DLIST_H

