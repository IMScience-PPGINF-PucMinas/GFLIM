//
// Created by Samuel Martins on 2019-01-04.
//

/*
  This program is a collection of functions to create, destroy, and
  manipulate a priority queue.

  A priority queue Q consists of two data structures: a circular
  queue C and a table L that encodes all possible doubly-linked
  lists.

  Q requires that the maximum possible increment along the paths be a
  non-negative integer less than the number of buckets in C. An extra
  bucket is created to store infinity values (positive and negative)
  for the LIFO policy. The queue size increases dynamically whenever
  (maxvalue-minvalue) > (nbuckets-1).

  Q->C.first[i] gives the first element that is in bucket i.
  Q->C.last[i]  gives the last  element that is in bucket i.
  Q->C.nbuckets gives the number of buckets in C.
  Q->C.minvalue  gives the minimum value of a node in queue.
  Q->C.maxvalue  gives the maximum value of a node in queue.
  Q->C.tiebreak gives the FIFO or LIFO tie breaking policy
  Q->C.removal_policy gives the MINVALUE or MAXVALUE removal policy

  All possible doubly-linked lists are represented in L. Each bucket
  contains a doubly-linked list that is treated as a FIFO.

  Q->L.elem[i].next: the next element to i
  Q->L.elem[i].prev: the previous element to i
  Q->L.elem[i].color: the color of i (IFT_WHITE=never inserted, IFT_GRAY=inserted,
  IFT_BLACK=removed)
  Q->L.nelems: gives the total number of elements that can be
  inserted in Q (It is usually the number of pixels in a given image
  or the number of nodes in a graph)
  Q->L.value[i]: gives the value of element i in the graph.

  Insertions and updates are done in O(1).
  Removal may take O(K+1), where K+1 is the number of buckets.
 */

#ifndef IFT_GQUEUE_H
#define IFT_GQUEUE_H


#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"

/** define queue to remove node with minimum value */
#define MINVALUE   0
/** define queue to remove node with maximum value */
#define MAXVALUE   1
/** define queue to solve ambiguity by FIFO */
#define FIFOBREAK  0
/** define queue to solve ambiguity by LIFO */
#define LIFOBREAK  1
/** define maximum size of the queue*/
#define IFT_QSIZE      65536

#define iftSetTieBreak(a,b) a->C.tiebreak=b
#define iftSetRemovalPolicy(a,b) a->C.removal_policy=b


typedef struct ift_gqnode {
    int  next;  /* next node */
    int  prev;  /* prev node */
    char color; /* IFT_WHITE=0, IFT_GRAY=1, IFT_BLACK=2 (IFT_WHITE=never inserted, IFT_GRAY=inserted,
  IFT_BLACK=removed)*/
} iftGQNode;

struct ift_gqueue;

typedef struct ift_gdoublylinkedlists {
    iftGQNode *elem;  /* all possible doubly-linked lists of the circular queue */
    int nelems;  /* total number of elements */
    int *value;   /* DEPRECATED: the value of the nodes in the graph.*/

    void *value_data;
    int (*value_function)(struct ift_gqueue *Q, int p);
} iftGDoublyLinkedLists;

typedef struct ift_gcircularqueue {
    int  *first;   /* list of the first elements of each doubly-linked list */
    int  *last;    /* list of the last  elements of each doubly-linked list  */
    int  nbuckets; /* number of buckets in the circular queue */
    int  minvalue;  /* minimum value of a node in queue */
    int  maxvalue;  /* maximum value of a node in queue */
    char tiebreak; /* 1 is LIFO, 0 is FIFO (default) */
    char removal_policy; /* 0 is MINVALUE and 1 is MAXVALUE */
} iftGCircularQueue;

typedef struct ift_gqueue { /*  Priority queue by Dial implemented as
                           proposed by A. Falcao */
    iftGCircularQueue C;
    iftGDoublyLinkedLists L;
} iftGQueue;

/**
 * @brief Creates a priority queue with predefined cost function
 * @details [long description]
 *
 * @param nbuckets Number of buckets
 * @param nelems Total number of elements
 * @param value Auxiliary data used by the cost function
 * @param cost_function Cost function predefined by the user
 * @return A priority queue
 */
iftGQueue *iftCreateGQueueWithCostFunction(int nbuckets, int nelems, void *value, int (*cost_function)(iftGQueue *Q, int p));

/**
 * @brief Creates a priority queue
 * @details [long description]
 *
 * @param nbuckets Number of buckets
 * @param nelems Total number of elements
 * @param value List of element values
 * @return A priority queue
 */
iftGQueue *iftCreateGQueue(int nbuckets, int nelems, int *value);

/**
 * @brief Destroys a priority queue
 * @details [long description]
 *
 * @param Q Queue to be destroyed
 */
void   iftDestroyGQueue(iftGQueue **Q);

/**
 * @brief Resets a priority queue to its default values
 * @details [long description]
 *
 * @param Q Priority queue to be reseted
 */
void   iftResetGQueue(iftGQueue *Q);

/**
 * @brief Resets a priority queue to its default values for a voxel
 * list. It must be used to speed up the process when the IFT
 * algorithm is executed too many times.
 * @details [long description]
 *
 * @param Q Priority queue to be reseted
 */
void iftResetGQueueForVoxelList(iftGQueue *Q, int *voxel, int nelems);

/**
 * @brief Verifies if a priority queue is empty or not
 * @details [long description]
 *
 * @param Q Target priority queue
 * @return Returns 1 if the priority queue is empty and 0 if not.
 */
int    iftEmptyGQueue(iftGQueue *Q);

/**
 * @brief Inserts an element into the priority queue
 * @details [long description]
 *
 * @param Q The priority queue
 * @param elem Element to insert
 */
void   iftInsertGQueue(iftGQueue **Q, int elem);

/**
 * @brief Removes the next element from the priority queue
 * @param The priority queue
 * @return The element removed
 */
int    iftRemoveGQueue(iftGQueue *Q);

/**
 * @brief Removes the specified element from the priority queue
 * @param Q The priority queue
 * @param elem Element to be removed
 */
void   iftRemoveGQueueElem(iftGQueue *Q, int elem);

/**
 * @brief Grows a priority queue
 * @param Old priority queue
 * @param nbuckets Number of buckets
 * @return The new priority queue
 */
iftGQueue *iftGrowGQueue(iftGQueue **Q, int nbuckets);

#ifdef __cplusplus
}
#endif

#endif //IFT_GQUEUE_H
