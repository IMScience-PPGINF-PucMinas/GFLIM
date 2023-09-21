#include "ift/core/dtypes/IntQueue.h"

#include "ift/core/io/Stream.h"
#include "iftMemory.h"


iftIntQueue *iftCreateIntQueue(int max_size) {
    iftIntQueue *Q = (iftIntQueue*) iftAlloc(1, sizeof(iftIntQueue));
    Q->max_size    = max_size;
    
    Q->val   = iftAllocIntArray(max_size);
    Q->first = Q->last = 0;
    Q->n     = 0;
    
    return Q;
}


void iftDestroyIntQueue(iftIntQueue **Q) {
    if (Q != NULL) {
        iftIntQueue *Qaux = *Q;
        
        if (Qaux != NULL) {
            iftFree(Qaux->val);
            iftFree(Qaux);
            *Q = NULL;
        }
    }
}


bool iftIsIntQueueEmpty(const iftIntQueue *Q) {
    return (Q->n == 0);
}


bool iftIsIntQueueFull(const iftIntQueue *Q) {
    return (Q->n == Q->max_size);
}


bool iftInsertIntQueue(iftIntQueue *Q, int elem) {
    if (iftIsIntQueueFull(Q)) {
        iftWarning("Queue is Full: size: %d\n", "iftInsertIntQueue", Q->n);
        return false;
    }
    else {
        Q->val[Q->last] = elem;
        Q->last         = (Q->last+1) % Q->max_size;
        Q->n++;
    }
    
    return true;
}


bool iftRemoveIntQueue(iftIntQueue *Q, int *elem) {
    if (iftIsIntQueueEmpty(Q)) {
        iftWarning("Queue is Empty", "iftRemoveIntQueue");
        return false;
    }
    else {
        *elem    = Q->val[Q->first];
        Q->first = (Q->first+1) % Q->max_size;
        Q->n--;
    }
    
    return true;
}


void iftPrintIntQueue(const iftIntQueue *Q) {
    if (!iftIsIntQueueEmpty(Q)) {
        int i = Q->first;
        
        do {
            printf("[%d] = %d\n", i, Q->val[i]);
            i = (i+1) % Q->max_size;
        } while (i != Q->last);
        puts("");
    }
}

