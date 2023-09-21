#include "ift/core/dtypes/SList.h"

#include "ift/core/dtypes/StrArray.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftMemory.h"
 
iftSList *iftCreateSList() {
    iftSList *SL = (iftSList *) iftAlloc(1, sizeof(iftSList));

    // just to really force
    SL->n    = 0;
    SL->head = NULL;
    SL->tail = NULL;

    return SL;
}



void iftDestroySList(iftSList **SL) {
    if (SL != NULL) {
        iftSList *SL_aux = *SL;

        if (SL_aux != NULL) {
            iftSNode *snode = SL_aux->head;
            iftSNode *tnode = NULL;

            while (snode != NULL) {
                tnode = snode;
                snode = snode->next;

                if (tnode->elem != NULL)
                    iftFree(tnode->elem);
                iftFree(tnode);
            }
            iftFree(SL_aux);
            *SL = NULL;
        }
    }
}



void iftInsertSListIntoHead(iftSList *SL, const char *elem) {
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftInsertSListIntoHead");

    iftSNode *snode = (iftSNode*) iftAlloc(1, sizeof(iftSNode));
    snode->elem     = iftAllocCharArray(512);
    snode->prev = NULL;
    snode->next     = NULL;
    strcpy(snode->elem, elem);


    // The String Linked List is empty
    if (SL->head == NULL) {
        SL->head = snode;
        SL->tail = snode;
        SL->n    = 1;
    }
    else {
        snode->next        = SL->head;
        SL->head->prev = snode;
        SL->head           = snode;
        SL->n++;
    }
}


void iftInsertSListIntoTail(iftSList *SL, const char *elem) {
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftInsertSListIntoTail");
    if (elem == NULL)
        iftError("The Element to be Inserted is NULL", "iftInsertSListIntoTail");

    iftSNode *snode = (iftSNode*) iftAlloc(1, sizeof(iftSNode));
    snode->prev     = NULL;
    snode->next     = NULL;
    snode->elem     = iftAllocCharArray(strlen(elem) + 1);
    strcpy(snode->elem, elem);


    // The String Linked List is empty
    if (SL->head == NULL) {
        SL->head = snode;
        SL->tail = snode;
        SL->n    = 1;
    }
    else {
        SL->tail->next = snode;
        snode->prev    = SL->tail;
        SL->tail       = snode;
        SL->n++;
    }
}


char *iftRemoveSListHead(iftSList *SL) {
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftRemoveSListHead");
    
    char *elem      = NULL;
    iftSNode *snode = NULL;

    // if there are elements
    if (SL->head != NULL) {
        snode    = SL->head;
        SL->head = SL->head->next;

        // checks if the list is empty now
        if (SL->head == NULL)
            SL->tail = NULL;
        else
            SL->head->prev = NULL;

        SL->n--;
        elem = snode->elem;
        snode->elem = NULL;
        iftFree(snode); // it does not deallocates snode->free
    }

    return elem;
}


char *iftRemoveSListTail(iftSList *SL) {
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftRemoveSListTail");

    char *elem      = NULL;
    iftSNode *snode = NULL;

    // if there are elements
    if (SL->head != NULL) {
        snode    = SL->tail;
        SL->tail = SL->tail->prev;

        // checks if the list is empty now
        if (SL->tail == NULL)
            SL->head = NULL;
        else
            SL->tail->next = NULL;

        SL->n--;
        elem = snode->elem;
        snode->elem = NULL;
        iftFree(snode); // it does not deallocates snode->free
    }

    return elem;
}


void iftPrintSList(iftSList *SL) {
    if (SL != NULL) {
        iftSNode *snode = NULL;
        printf("Number of nodes: %lu\n", SL->n);

        snode = SL->head;

        puts("---");
        while (snode != NULL) {
            printf("%s\n", snode->elem);
            snode = snode->next;
        }
        puts("---");
    }
}


void iftPrintReversedSList(iftSList *SL) {
    if (SL != NULL) {
        iftSNode *snode = NULL;
        printf("Number of nodes: %lu\n", SL->n);

        snode = SL->tail;

        while (snode != NULL) {
            printf("%s\n", snode->elem);
            snode = snode->prev;
        }

        puts("");
    }
}


iftSList *iftCopySList(const iftSList *SL) {
    iftSList *copy = iftCreateSList();

    for (iftSNode *snode = SL->head; snode != NULL; snode = snode->next)
        iftInsertSListIntoTail(copy, snode->elem);

    return copy;
}


iftStrArray *iftSListToStrArray(const iftSList *SL) {
    iftStrArray *sarr = iftCreateStrArray(SL->n);

    iftSNode *node = SL->head;
    int i = 0;
    while (node != NULL) {
        strcpy(sarr->val[i++], node->elem);
        node = node->next;
    }


    return sarr;
}


















