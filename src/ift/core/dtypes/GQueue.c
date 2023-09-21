#include "ift/core/dtypes/GQueue.h"

/*
  Copyright (C) <2003> <Alexandre Xavier Falc�o>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  please see full copyright in COPYING file.
  -------------------------------------------------------------------------
  written by A.X. Falc�o <afalcao@ic.unicamp.br>, May 13th 2007

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

#include "ift/core/dtypes/GQueue.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"

int default_value_function(iftGQueue *Q, int p) {
    int *value = (int*) Q->L.value_data;
    return value[p];
}

iftGQueue *iftCreateGQueue(int nbuckets, int nelems, int *value)
{
    iftGQueue *Q =  iftCreateGQueueWithCostFunction(nbuckets, nelems, value,
                                                    default_value_function);
    Q->L.value = value;

    return Q;
}

iftGQueue *iftCreateGQueueWithCostFunction(int nbuckets, int nelems, void *value, int (*cost_function)(iftGQueue *Q, int p)) {
    iftGQueue *Q=NULL;

    Q = (iftGQueue *) iftAlloc(1, sizeof(iftGQueue));

    if (Q != NULL)
    {
        Q->C.first = (int *)iftAlloc((nbuckets+1), sizeof(int));
        Q->C.last  = (int *)iftAlloc((nbuckets+1), sizeof(int));
        Q->C.nbuckets = nbuckets;
        if ( (Q->C.first != NULL) && (Q->C.last != NULL) )
        {
            Q->L.elem = (iftGQNode *)iftAlloc(nelems, sizeof(iftGQNode));
            Q->L.nelems = nelems;
            Q->L.value  = NULL;
            Q->L.value_data  = value;
            Q->L.value_function = cost_function;
            if (Q->L.elem != NULL)
            {
                iftResetGQueue(Q);
            }
            else
                iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateGQueue");
        }
        else
            iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateGQueue");
    }
    else
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateGQueue");

    /* default */

    iftSetTieBreak(Q,FIFOBREAK);
    iftSetRemovalPolicy(Q,MINVALUE);

    return(Q);
}

void iftResetGQueue(iftGQueue *Q)
{
    int i;

    Q->C.minvalue = IFT_INFINITY_INT;
    Q->C.maxvalue = IFT_INFINITY_INT_NEG;
    /* No need for that, since the programmer might have changed them  */
    //    iftSetTieBreak(Q,FIFOBREAK);
    //    iftSetRemovalPolicy(Q,MINVALUE);
    for (i=0; i < Q->C.nbuckets+1; i++)
        Q->C.first[i]=Q->C.last[i]= IFT_NIL;

    for (i=0; i < Q->L.nelems; i++)
    {
        Q->L.elem[i].next =  Q->L.elem[i].prev = IFT_NIL;
        Q->L.elem[i].color = IFT_WHITE;
    }

}

void iftResetGQueueForVoxelList(iftGQueue *Q, int *voxel, int nelems)
{
    int i;

    Q->C.minvalue = IFT_INFINITY_INT;
    Q->C.maxvalue = IFT_INFINITY_INT_NEG;
    for (i=0; i < Q->C.nbuckets+1; i++)
        Q->C.first[i]=Q->C.last[i]= IFT_NIL;

    for (i=0; i < nelems; i++)
    {
        Q->L.elem[voxel[i]].next  =  Q->L.elem[voxel[i]].prev = IFT_NIL;
        Q->L.elem[voxel[i]].color =  IFT_WHITE;
    }

}

void iftDestroyGQueue(iftGQueue **Q)
{
    iftGQueue *aux=*Q;

    if (aux != NULL)
    {
        if (aux->C.first != NULL) iftFree(aux->C.first);
        if (aux->C.last  != NULL) iftFree(aux->C.last);
        if (aux->L.elem  != NULL) iftFree(aux->L.elem);
        iftFree(aux);
        *Q = NULL;
    }
}

iftGQueue *iftGrowGQueue(iftGQueue **Q, int nbuckets)
{
    iftGQueue *Q1=iftCreateGQueueWithCostFunction(nbuckets,(*Q)->L.nelems,
                                                  (*Q)->L.value_data,
                                                  (*Q)->L.value_function);
    int i,bucket;

    Q1->C.minvalue  = (*Q)->C.minvalue;
    Q1->C.maxvalue  = (*Q)->C.maxvalue;
    Q1->C.tiebreak = (*Q)->C.tiebreak;
    Q1->C.removal_policy = (*Q)->C.removal_policy;
    for (i=0; i<(*Q)->C.nbuckets; i++)
        if ((*Q)->C.first[i] != IFT_NIL)
        {
            bucket = iftSafeMod((*Q)->L.value_function((*Q), (*Q)->C.first[i]),Q1->C.nbuckets);
            Q1->C.first[bucket] = (*Q)->C.first[i];
            Q1->C.last[bucket]  = (*Q)->C.last[i];
        }
    if ((*Q)->C.first[(*Q)->C.nbuckets] != IFT_NIL)
    {
        bucket = Q1->C.nbuckets;
        Q1->C.first[bucket] = (*Q)->C.first[(*Q)->C.nbuckets];
        Q1->C.last[bucket]  = (*Q)->C.last[(*Q)->C.nbuckets];
    }

    for (i=0; i < (*Q)->L.nelems; i++)
        Q1->L.elem[i]  = (*Q)->L.elem[i];

    iftDestroyGQueue(Q);
    return(Q1);
}


void iftInsertGQueue(iftGQueue **Q, int elem)
{
    int bucket,minvalue=(*Q)->C.minvalue,maxvalue=(*Q)->C.maxvalue;

    if (((*Q)->L.value_function((*Q), elem) == IFT_INFINITY_INT) || ((*Q)->L.value_function((*Q), elem) == IFT_INFINITY_INT_NEG))
        bucket=(*Q)->C.nbuckets;
    else
    {
        if ((*Q)->L.value_function((*Q), elem) < minvalue)
            minvalue = (*Q)->L.value_function((*Q), elem);
        if ((*Q)->L.value_function((*Q), elem) > maxvalue)
            maxvalue = (*Q)->L.value_function((*Q), elem);
        if ((maxvalue-minvalue) > ((*Q)->C.nbuckets-1))
        {
            (*Q) = iftGrowGQueue(Q,2*(maxvalue-minvalue)+1);
            fprintf(stdout,"Warning: Doubling queue size\n");
        }
        if ((*Q)->C.removal_policy==MINVALUE)
        {
            bucket=iftSafeMod((*Q)->L.value_function((*Q), elem),(*Q)->C.nbuckets);
        }
        else
        {
            bucket=(*Q)->C.nbuckets-1-(iftSafeMod((*Q)->L.value_function((*Q), elem),(*Q)->C.nbuckets));
        }
        (*Q)->C.minvalue = minvalue;
        (*Q)->C.maxvalue = maxvalue;
    }
    if ((*Q)->C.first[bucket] == IFT_NIL)
    {
        (*Q)->C.first[bucket]   = elem;
        (*Q)->L.elem[elem].prev = IFT_NIL;
    }
    else
    {
        (*Q)->L.elem[(*Q)->C.last[bucket]].next = elem;
        (*Q)->L.elem[elem].prev = (*Q)->C.last[bucket];
    }

    (*Q)->C.last[bucket]     = elem;
    (*Q)->L.elem[elem].next  = IFT_NIL;
    (*Q)->L.elem[elem].color = IFT_GRAY;
}

int iftRemoveGQueue(iftGQueue *Q)
{
    int elem= IFT_NIL, next, prev;
    int last, current;

    if (Q->C.removal_policy==MINVALUE)
        current=iftSafeMod(Q->C.minvalue,Q->C.nbuckets);
    else
        current=Q->C.nbuckets-1-iftSafeMod(Q->C.maxvalue,Q->C.nbuckets);

    /** moves to next element **/

    if (Q->C.first[current] == IFT_NIL)
    {
        last = current;

        current = iftSafeMod(current + 1, Q->C.nbuckets);

        while ((Q->C.first[current] == IFT_NIL) && (current != last))
        {
            current = iftSafeMod(current + 1, Q->C.nbuckets);
        }

        if (Q->C.first[current] != IFT_NIL)
        {
            if (Q->C.removal_policy==MINVALUE)
                Q->C.minvalue = Q->L.value_function(Q, Q->C.first[current]);
            else
                Q->C.maxvalue = Q->L.value_function(Q, Q->C.first[current]);
        }
        else
        {
            if (Q->C.first[Q->C.nbuckets] != IFT_NIL)
            {
                current = Q->C.nbuckets;
                if (Q->C.removal_policy==MINVALUE)
                    Q->C.minvalue = Q->L.value_function(Q, Q->C.first[current]);
                else
                    Q->C.maxvalue = Q->L.value_function(Q, Q->C.first[current]);
            }
            else
            {
                iftError("iftGQueue is empty\n", "iftRemoveGQueue");
            }
        }
    }

    if (Q->C.tiebreak == LIFOBREAK)
    {
        elem = Q->C.last[current];
        prev = Q->L.elem[elem].prev;
        if (prev == IFT_NIL)           /* there was a single element in the list */
        {
            Q->C.last[current] = Q->C.first[current]  = IFT_NIL;
        }
        else
        {
            Q->C.last[current]   = prev;
            Q->L.elem[prev].next = IFT_NIL;
        }
    }
    else   /* Assume FIFO policy for breaking ties */
    {
        elem = Q->C.first[current];
        next = Q->L.elem[elem].next;
        if (next == IFT_NIL)           /* there was a single element in the list */
        {
            Q->C.first[current] = Q->C.last[current]  = IFT_NIL;
        }
        else
        {
            Q->C.first[current] = next;
            Q->L.elem[next].prev = IFT_NIL;
        }
    }

    Q->L.elem[elem].color = IFT_BLACK;

    return elem;
}

void iftRemoveGQueueElem(iftGQueue *Q, int elem)
{
    int prev,next,bucket;

    if ((Q->L.value_function(Q, elem) == IFT_INFINITY_INT) || (Q->L.value_function(Q, elem) == IFT_INFINITY_INT_NEG))
        bucket = Q->C.nbuckets;
    else
    {
        if (Q->C.removal_policy == MINVALUE)
            bucket = iftSafeMod(Q->L.value_function(Q, elem),Q->C.nbuckets);
        else
            bucket = Q->C.nbuckets-1-iftSafeMod(Q->L.value_function(Q, elem),Q->C.nbuckets);
    }

    prev = Q->L.elem[elem].prev;
    next = Q->L.elem[elem].next;

    /* if elem is the first element */
    if (Q->C.first[bucket] == elem)
    {
        Q->C.first[bucket] = next;
        if (next == IFT_NIL) /* elem is also the last one */
            Q->C.last[bucket] = IFT_NIL;
        else
            Q->L.elem[next].prev = IFT_NIL;
    }
    else    /* elem is in the middle or it is the last */
    {
        Q->L.elem[prev].next = next;
        if (next == IFT_NIL) /* if it is the last */
            Q->C.last[bucket] = prev;
        else
            Q->L.elem[next].prev = prev;
    }

    Q->L.elem[elem].color = IFT_WHITE;

}

int iftEmptyGQueue(iftGQueue *Q)
{
    int last,current;

    if (Q->C.removal_policy == MINVALUE)
        current=iftSafeMod(Q->C.minvalue,Q->C.nbuckets);
    else
        current=Q->C.nbuckets - 1 - (iftSafeMod(Q->C.maxvalue,Q->C.nbuckets));

    if (Q->C.first[current] != IFT_NIL)
        return 0;

    last = current;

    current = iftSafeMod(current + 1, Q->C.nbuckets);

    while ((Q->C.first[current] == IFT_NIL) && (current != last))
    {
        current = iftSafeMod(current + 1, Q->C.nbuckets);
    }

    if (Q->C.first[current] == IFT_NIL)
    {
        if (Q->C.first[Q->C.nbuckets] == IFT_NIL)
        {
            //Changed by Falcao and Nikolas
            // iftResetGQueue(Q);
            return(1);
        }
    }

    return (0);
}

