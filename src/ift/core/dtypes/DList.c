#include "ift/core/dtypes/DList.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"


iftDList *iftCreateDList() {
    iftDList *DL = (iftDList *) iftAlloc(1, sizeof(iftDList));

    // just to really force
    DL->n    = 0;
    DL->head = NULL;
    DL->tail = NULL;

    return DL;
}


void iftDestroyDList(iftDList **DL) {
    iftDList *DL_aux = *DL;

    if (DL_aux != NULL) {
        iftDNode *snode = DL_aux->head;
        iftDNode *tnode = NULL;

        while (snode != NULL) {
            tnode = snode;
            snode = snode->next;
            iftFree(tnode);
        }
        iftFree(DL_aux);
        *DL = NULL;
    }
}


void iftInsertDListIntoHead(iftDList *DL, double elem) {
    if (DL == NULL)
        iftError("The Doubly Linked List DL is NULL. Allocated it firstly", "iftInsertDListIntoHead");

    iftDNode *node  = (iftDNode*) iftAlloc(1, sizeof(iftDNode));
    node->elem     = elem;
    node->previous = NULL;
    node->next     = NULL;


    // The String Linked List is empty
    if (DL->head == NULL) {
        DL->head = node;
        DL->tail = node;
        DL->n    = 1;
    }
    else {
        node->next         = DL->head;
        DL->head->previous = node;
        DL->head           = node;
        DL->n++;
    }
}


void iftInsertDListIntoTail(iftDList *DL, double elem) {
    if (DL == NULL)
        iftError("The Linked List DL is NULL. Allocated it firstly", "iftInsertDListIntoTail");

    iftDNode *node  = (iftDNode*) iftAlloc(1, sizeof(iftDNode));
    node->elem     = elem;
    node->previous = NULL;
    node->next     = NULL;


    // The String Linked List is empty
    if (DL->head == NULL) {
        DL->head = node;
        DL->tail = node;
        DL->n    = 1;
    }
    else {
        DL->tail->next = node;
        node->previous = DL->tail;
        DL->tail       = node;
        DL->n++;
    }
}


double iftRemoveDListHead(iftDList *DL) {
    if (DL == NULL)
        iftError("The Doubly Linked List DL is NULL. Allocated it firstly", "iftRemoveDListHead");

    double  elem  = IFT_NIL;
    iftDNode *node = NULL;

    // if there are elements
    if (DL->head != NULL) {
        node     = DL->head;
        DL->head = DL->head->next;

        // checks if the list is empty now
        if (DL->head == NULL)
            DL->tail = NULL;
        else
            DL->head->previous = NULL;

        DL->n--;
        elem = node->elem;
        iftFree(node);
    }

    return elem;
}


double iftRemoveDListTail(iftDList *DL) {
    if (DL == NULL)
        iftError("The Doubly Linked List DL is NULL. Allocated it firstly", "iftRemoveListTail");

    double elem   = IFT_NIL;
    iftDNode *node = NULL;

    // if there are elements
    if (DL->head != NULL) {
        node     = DL->tail;
        DL->tail = DL->tail->previous;

        // checks if the list is empty now
        if (DL->tail == NULL)
            DL->head = NULL;
        else
            DL->tail->next = NULL;

        DL->n--;
        elem = node->elem;
        iftFree(node);
    }

    return elem;
}


void iftPrintDList(const struct _ift_dlist *DL) {
    if (DL != NULL) {
        iftDNode *node = NULL;
        printf("Number of nodes: %d\n", DL->n);

        node = DL->head;

        while (node != NULL) {
            printf("%lf\n", node->elem);
            node = node->next;
        }
        puts("");
    }
}


void iftPrintReversedDList(const struct _ift_dlist *DL) {
    if (DL != NULL) {
        iftDNode *node = NULL;
        printf("Number of nodes: %d\n", DL->n);

        node = DL->tail;

        while (node != NULL) {
            printf("%lf\n", node->elem);
            node = node->previous;
        }
        puts("");
    }
}


iftDblArray *iftDListToDblArray(const iftDList *DL) {
    iftDblArray *arr = iftCreateDblArray(DL->n);
    
    iftDNode *node = DL->head;
    
    for (int i = 0; i < DL->n; i++) {
        arr->val[i] = node->elem;
        node = node->next;
    }
    
    return arr;
}


iftDList *iftDblArrayToDList(const iftDblArray *arr) {
    iftDList *DL = iftCreateDList();
    
    for (long i = 0; i < arr->n; i++)
        iftInsertDListIntoTail(DL, arr->val[i]);
    
    return DL;
}


