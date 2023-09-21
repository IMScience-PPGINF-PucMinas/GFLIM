/**
 * @file iftSList.h
 * @brief String Doubly Linked List definitions and functions.
 * @author Samuel Martins
 * @date Sep 16, 2015
 * @ingroup StringDoublyLinkedList
 *
 * @note Programs:
 * * @ref iftTestSList.c = It shows how to use the String Doubly Linked List.
 */

#ifndef IFT_SLIST_H
#define IFT_SLIST_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Node of the String Doubly Linked List.
 * @author Samuel Martins
 * @date Sep 16, 2015
 * @ingroup StringDoublyLinkedList
 */
typedef struct ift_snode {
    char *elem; 
    struct ift_snode *prev;
    struct ift_snode *next;
} iftSNode;


/**
 * @brief String Doubly Linked List. 
 * @author Samuel Martins
 * @date Sep 16, 2015
 * @ingroup StringDoublyLinkedList
 */
typedef struct ift_slist {
    /** Number of Nodes of the String Doubly Linked List */
    long n;
    /** Pointer to the Head (first) Node */
    iftSNode *head;
    /** Pointer to the Tail (last) Node */
    iftSNode *tail;
} iftSList;



/**
 * @brief Creates an iftSList.
 * @author Samuel Martins
 * @date Sep 16, 2015
 * @ingroup StringDoublyLinkedList
 */
iftSList *iftCreateSList();


/**
 * @brief Destroys an iftSList.
 * @author Samuel Martins
 * @date Sep 16, 2015
 * @ingroup StringDoublyLinkedList
 *
 * @param SL Reference to the iftSList to be destroyed.
 */
void iftDestroySList(iftSList **SL);


/**
 * @brief Checks ift a String List is empty.
 * @author Samuel Martins
 * @date Sep 16, 2015
 */
static inline bool iftIsSListEmpty(const iftSList *SL) {
    return (SL->n == 0);
}


/**
 * @brief Insert a string into the head of the iftSList.
 * @author Samuel Martins
 * @date Sep 18, 2015
 * @ingroup StringDoublyLinkedList
 *
 * Insert the string <elem> into the head of the iftSList. The string <elem> is ALWAYS COPIED.
 * If the iftSList is not allocated, an error will be raised.
 *
 * @param SL The iftSList.
 * @param elem The string to be inserted into the head of the iftSList.
 */
void iftInsertSListIntoHead(iftSList *SL, const char *elem);


/**
 * @brief Insert a string into the tail of the iftSList.
 * @author Samuel Martins
 * @date Sep 18, 2015
 * @ingroup StringDoublyLinkedList
 *
 * Insert the string <elem> into the tail of the iftSList. The string <elem> is ALWAYS COPIED.
 * If the iftSList is not allocated, an error will be raised.
 *
 * @param SL The iftSList.
 * @param elem The string to be inserted into the tail of the iftSList.
 */
void iftInsertSListIntoTail(iftSList *SL, const char *elem);


/**
 * @brief Removes the Head element from the iftSList and removes its node. If there is no element, a NULL pointer will be returned.
 * @author Samuel Martins
 * @date Sep 18, 2015
 * @ingroup StringDoublyLinkedList
 *
 * Removes the Head string node of the iftSList and removes its node.
 * If there is no element, a NULL pointer will be returned.
 * If the iftSList is not allocated, an error will be raised.
 *
 * @param SL The iftSList.
 * @return The Head string from the iftSList.
 */
char *iftRemoveSListHead(iftSList *SL);


/**
 * @brief Removes the Tail element from the iftSList and removes its node. If there is no element, a NULL pointer will be returned.
 * @author Samuel Martins
 * @date Sep 18, 2015
 * @ingroup StringDoublyLinkedList
 *
 * Removes the Tail string node of the iftSList and removes its node.
 * If there is no element, a NULL pointer will be returned.
 * If the iftSList is not allocated, an error will be raised.
 *
 * @param SL The iftSList.
 * @return The Tail string from the iftSList.
 */
char *iftRemoveSListTail(iftSList *SL);


/**
 * @brief Prints an iftSList (from the Node to the Tail node).
 * @author Samuel Martins
 * @date Sep 18, 2015
 * @ingroup StringDoublyLinkedList
 *
 * Prints an iftSList from the Head to Tail node.
 *
 * @param SL The iftSList.
 */
void iftPrintSList(iftSList *SL);


/**
 * @brief Prints an iftSList in Reveresed Order (from the Tail to the Head node).
 * @author Samuel Martins
 * @date Sep 18, 2015
 * @ingroup StringDoublyLinkedList
 *
 * Prints an iftSList from the Tail to Head node.
 *
 * @param SL The iftSList.
 */
void iftPrintReversedSList(iftSList *SL);


/**
 * @brief Copies a String List.
 * @author Samuka Martins
 * @date Oct 30, 2017
 */
iftSList *iftCopySList(const iftSList *SL);


/**
 * @brief Convert a String List into a String Array.
 * @author Samuka Martins
 * @date Jul 31, 2018
 */
iftStrArray *iftSListToStrArray(const iftSList *SL);


#ifdef __cplusplus
}
#endif

#endif //IFT_SLIST_H

