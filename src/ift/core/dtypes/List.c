#include "ift/core/dtypes/List.h"

#include "ift/core/io/Stream.h"


iftList *iftCreateList() {
    iftList *L = (iftList *) iftAlloc(1, sizeof(iftList));
    
    // just to really force
    L->n    = 0;
    L->head = NULL;
    L->tail = NULL;
    
    return L;
}


void iftDestroyList(iftList **L) {
    iftList *L_aux = *L;
    
    if (L_aux != NULL) {
        iftNode *snode = L_aux->head;
        iftNode *tnode = NULL;
        
        while (snode != NULL) {
            tnode = snode;
            snode = snode->next;
            iftFree(tnode);
        }
        iftFree(L_aux);
        *L = NULL;
    }
}


bool iftIsEmptyList(const iftList *L) {
    return (L->n == 0);
}


void iftInsertListIntoHead(iftList *L, int elem) {
    if (L == NULL)
        iftError("The Integer Doubly Linked List L is NULL. Allocated it firstly", "iftInsertListIntoHead");
    
    iftNode *node  = (iftNode*) iftAlloc(1, sizeof(iftNode));
    node->elem     = elem;
    node->previous = NULL;
    node->next     = NULL;
    
    
    // The String Linked List is empty
    if (L->head == NULL) {
        L->head = node;
        L->tail = node;
        L->n    = 1;
    }
    else {
        node->next        = L->head;
        L->head->previous = node;
        L->head           = node;
        L->n++;
    }
}


void iftInsertListIntoTail(iftList *L, int elem) {
    if (L == NULL)
        iftError("The Integer Linked List L is NULL. Allocated it firstly", "iftInsertListIntoTail");
    
    iftNode *node  = (iftNode*) iftAlloc(1, sizeof(iftNode));
    node->elem     = elem;
    node->previous = NULL;
    node->next     = NULL;
    
    
    // The String Linked List is empty
    if (L->head == NULL) {
        L->head = node;
        L->tail = node;
        L->n    = 1;
    }
    else {
        L->tail->next  = node;
        node->previous = L->tail;
        L->tail        = node;
        L->n++;
    }
}


int iftRemoveListHead(iftList *L) {
    if (L == NULL)
        iftError("The Doubly Linked List L is NULL. Allocated it firstly", "iftRemoveListHead");
    
    int elem      = IFT_NIL;
    iftNode *node = NULL;
    
    // if there are elements
    if (L->head != NULL) {
        node    = L->head;
        L->head = L->head->next;
        
        // checks if the list is empty now
        if (L->head == NULL)
            L->tail = NULL;
        else
            L->head->previous = NULL;
        
        L->n--;
        elem = node->elem;
        iftFree(node);
    }
    
    return elem;
}


int iftRemoveListTail(iftList *L) {
    if (L == NULL)
        iftError("The Integer Doubly Linked List L is NULL. Allocated it firstly", "iftRemoveListTail");
    
    int elem      = IFT_NIL;
    iftNode *node = NULL;
    
    // if there are elements
    if (L->head != NULL) {
        node    = L->tail;
        L->tail = L->tail->previous;
        
        // checks if the list is empty now
        if (L->tail == NULL)
            L->head = NULL;
        else
            L->tail->next = NULL;
        
        L->n--;
        elem = node->elem;
        iftFree(node);
    }
    
    return elem;
}


void iftPrintList(const iftList *L) {
    if (L != NULL) {
        iftNode *node = NULL;
        printf("Number of nodes: %d\n", L->n);
        
        node = L->head;
        
        while (node != NULL) {
            printf("%d ", node->elem);
            node = node->next;
        }
        puts("");
    }
}


void iftPrintReversedList(const iftList *L) {
    if (L != NULL) {
        iftNode *node = NULL;
        printf("Number of nodes: %d\n", L->n);
        
        node = L->tail;
        
        while (node != NULL) {
            printf("%d\n", node->elem);
            node = node->previous;
        }
        puts("");
    }
}



iftIntArray *iftListToIntArray(const iftList *L) {
    iftIntArray *arr = iftCreateIntArray(L->n);
    
    iftNode *node = L->head;
    
    for (int i = 0; i < L->n; i++) {
        arr->val[i] = node->elem;
        node = node->next;
    }
    
    return arr;
}


iftList *iftIntArrayToList(const iftIntArray *iarr) {
    iftList *L = iftCreateList();
    
    for (long i = 0; i < iarr->n; i++)
        iftInsertListIntoTail(L, iarr->val[i]);
    
    return L;
}

