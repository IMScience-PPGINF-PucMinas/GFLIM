//
// Created by Samuel Martins on 2019-01-04.
//

#ifndef IFT_INTQUEUE_H
#define IFT_INTQUEUE_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Circular Integer Queue (FIFO).
 * @author Samuel Martins
 * @date Jun 14, 2016
 * @ingroup DataStructures
 */
typedef struct ift_int_queue {
    /** Index of the first element of the queue */
    int first;
    /** Index of the last element of the queue */
    int last;
    /** Integer Queue: array of integer elements */
    int *val;
    /** Number of elements inserted into queue */
    int n;
    /** Maximum number of elements of the queue (array size) */
    int max_size;
} iftIntQueue;


/**
 * @brief Creates a Static Circular Integer Queue.
 * @author Samuel Martins
 * @date Jun 14, 2016
 * @ingroup DataStructures
 *
 * @param  max_size Maximum number of elements of the Integer Queue.
 * @return          The allocated FIFO.
 */
iftIntQueue *iftCreateIntQueue(int max_size);


/**
 * @brief Destroys a Circular Integer Queue.
 * @author Samuel Martins
 * @date Jun 14, 2016
 * @ingroup DataStructures
 */
void iftDestroyIntQueue(iftIntQueue **Q);


/**
 * @brief Checks if a Circular Integer Queue is empty.
 * @author Samuel Martins
 * @date Jun 14, 2016
 * @ingroup DataStructures
 */
bool iftIsIntQueueEmpty(const iftIntQueue *Q);


/**
 * @brief Checks if a Circular Integer Queue is full.
 * @author Samuel Martins
 * @date Jun 14, 2016
 * @ingroup DataStructures
 */
bool iftIsIntQueueFull(const iftIntQueue *Q);


/**
 * @brief Inserts an element into Integer Queue.
 * @author Samuel Martins
 * @date Jun 14, 2016
 * @ingroup DataStructures
 *
 * @param  Q    Target Integer Queue.
 * @param  elem Element to be inserted into Integer Queue.
 * @return      True if the insertion was succeed, false otherwise.
 */
bool iftInsertIntQueue(iftIntQueue *Q, int elem);


/**
 * @brief Removes an element from the Integer Queue.
 * @author Samuel Martins
 * @date Jun 14, 2016
 * @ingroup DataStructures
 *
 * @param  Q    Target Integer Queue.
 * @param  elem Return by Reference the element removed from the Queue.
 * @return      True if the remotion was succeed, false otherwise.
 */
bool iftRemoveIntQueue(iftIntQueue *Q, int *elem);


/**
 * @brief Prints an Integer Queue.
 * @author Samuel Martins
 * @date Jun 14, 2016
 * @ingroup DataStructures
 */
void iftPrintIntQueue(const iftIntQueue *Q);


#ifdef __cplusplus
}
#endif

#endif //IFT_INTQUEUE_H
