//
// Created by Samuel Martins on 2019-01-04.
//

#ifndef IFT_FIFO_H
#define IFT_FIFO_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"

/** Simple Static Queue/FIFO struct. */
typedef struct ift_fifo {
    /** Array of elements */
    int *FIFO;
    /** Maximum number of elements of the queue (array size) */
    int n;
    /** Index of the first element of the queue */
    int first;
    /** Index of the last element of the queue */
    int last;
    /** Array of flags, one for each element from the queue */
    char *color;
} iftFIFO;


/**
 * @brief Creates a Static Queue/FIFO.
 * @param  n Maximum number of elements of the Queue.
 * @return   The allocated FIFO.
 */
iftFIFO *iftCreateFIFO(int n);


/**
 * @brief Destroys a FIFO.
 */
void iftDestroyFIFO(iftFIFO **F);


/**
 * @brief Inserts an element into Queue.
 * @param  F    Target Queue/FIFO.
 * @param  elem Element to be inserted into Queue.
 * @return      True if the insertion was succeed.
 */
char iftInsertFIFO(iftFIFO *F, int elem);


/**
 * @brief Removes an element from Queue or IFT_NILL if the remotion was not succeeded.
 * @param  F Target Queue.
 * @return   Element removed or IFT_NILL if a problem occurred.
 */
int iftRemoveFIFO(iftFIFO *F);


/**
 * @brief Checks if a Queue/FIFO is FULL.
 */
bool iftFullFIFO(iftFIFO *F);


/**
 * @brief Checks if a Queue/FIFO is Empty.
 */
bool iftEmptyFIFO(iftFIFO *F);


/**
 * @brief Resets the Queue/FIFO releasing all position and setting their colors to IFT_WHITE.
 */
void iftResetFIFO(iftFIFO *F);


/**
 * @brief Returns the Color from the position <b>pos</b> from the Queue/FIFO
 * @param  F   Target Queue/FIFO.
 * @param  pos Target position.
 * @return     The color of the position <b>pos</b>.
 */
int iftColorFIFO(iftFIFO *F, int pos);


#ifdef __cplusplus
}
#endif


#endif //IFT_FIFO_H
