//
// Created by deangeli on 5/1/17.
//

#ifndef _IFTGENETICLINKEDLIST_H_
#define _IFTGENETICLINKEDLIST_H_

#include "iftCommon.h"
#include "iftDataTransference.h"

typedef bool (*iftComparatorFunction)(void *,void *);


typedef struct _iftGenericLinkedListNode {
    void *data;
    struct _iftGenericLinkedListNode *next;
    struct _iftGenericLinkedListNode *previous;
} iftGenericLinkedListNode;

typedef struct _iftGenericlinkedList{
    size_t length;
    size_t dataSize;
    iftGenericLinkedListNode *head;
    iftGenericLinkedListNode *tail;
    iftGenericLinkedListNode *currentNode;
    int currentNodeIndex;
    bool isCircular;
    FreeFunction freeFunction;
    iftComparatorFunction comparatorFunction;
} iftGenericLinkedList;

typedef struct _iftGenericlistIterator {
    iftGenericLinkedListNode *node;
} iftGenericListIterator;

static inline void  iftSetCursorPosition(iftGenericLinkedList* list, int index,iftGenericLinkedListNode* node){
    list->currentNodeIndex = index;
    list->currentNode = node;
}



//LinkedList* createLinkedList(size_t dataSize);
iftGenericLinkedList* iftcreateLinkedList(size_t dataSize);
static inline void iftAppendElementInList(iftGenericLinkedList *list, void *element){
    iftGenericLinkedListNode *node = (iftGenericLinkedListNode*)calloc(1,sizeof(iftGenericLinkedListNode));
    node->data = calloc(1,list->dataSize);
    node->next = NULL;

    TRANSFER_DATA_COPY(node->data, element, list->dataSize);
    //memcpy(node->data, element, list->dataSize);

    if(list->length == 0) {
        list->head = list->tail = node;
    } else {
        list->tail->next = node;
        node->previous = list->tail;
        list->tail = node;
    }
    if(list->isCircular){
        list->tail->next = list->head;
        list->head->previous = list->tail;
    }

    list->length++;
}
#define LIST_PUSH_BACK(type,list,element) \
    LinkedListNode *node = (LinkedListNode*)calloc(1,sizeof(LinkedListNode)); \
    node->data = calloc(1,list->dataSize); \
    node->next = NULL; \
    ((type*)node->data)[0] = element; \
    if(list->length == 0) {\
        list->head = list->tail = node;\
    }\
    else {\
        list->tail->next = node;\
        node->previous = list->tail;\
        list->tail = node;\
    }\
    if(list->isCircular){\
        list->tail->next = list->head; \
        list->head->previous = list->tail;\
    }\
    list->length++;

#define LIST_GET_ELEMENT_AS(type,list,index) (*((type*)getElementInList(list,index)))

static inline void iftPushBackElementInList(iftGenericLinkedList *list, void *element){
    iftAppendElementInList(list,element);
}
static inline void iftPrependElementInList(iftGenericLinkedList *list, void *element){
    iftGenericLinkedListNode *node = (iftGenericLinkedListNode*)malloc(sizeof(iftGenericLinkedListNode));
    node->data = malloc(list->dataSize);
    TRANSFER_DATA_COPY(node->data, element, list->dataSize);

    if(list->length == 0) {
        list->head = node;
        list->tail = list->head;
    }else{
        node->next = list->head;
        list->head->previous = node;
        list->head = node;
    }
    if(list->isCircular){
        list->tail->next = list->head;
        list->head->previous = list->tail;
    }
    list->length++;
}

static inline void iftPushFrontElementInList(iftGenericLinkedList *list, void *element){
    iftPrependElementInList(list,element);
}
void iftInsertElementInListAt(iftGenericLinkedList *list, void *element, size_t index);
void iftRemoveGenericListHead(iftGenericLinkedList *list);
void iftRemoveGenericListTail(iftGenericLinkedList *list);
void iftRemoveElementInGenericListAt(iftGenericLinkedList *list, size_t index);
void iftRemoveElementInGenericListByReference(iftGenericLinkedList *list, void *element);
void iftRemoveElementInGenericListGivenValue(iftGenericLinkedList *list, void *element);
void iftRemoveElementsInGenericListGivenValue(iftGenericLinkedList *list, void *element);
void iftResetGenericListIterator(iftGenericLinkedList *list);
void* iftGetNextElementInGenericList(iftGenericLinkedList *list);
void* iftGetPreviousElementInGenericList(iftGenericLinkedList *list);
static inline void* iftGetElementIngenericList(iftGenericLinkedList *list, size_t index){
    if(index >= list->length){
        if(list->isCircular){
            index = index % list->length;
        }else{
            printf("[getElement] invalid position %lu. The list length is %lu (indexing start from 0)\n", index,list->length);
            return NULL;
        }
    }

    if (index == 0){
        return list->head->data;
    }
    if(index == list->length-1){
        return list->tail->data;
    }

    int distance2Head = index;
    int distance2Tail = list->length -index;
    int distance2Current = index - list->currentNodeIndex;
    int currentDirection = 0; //foward
    if(distance2Current <= 0){
        currentDirection = 1;//backward
        distance2Current = -distance2Current;
    }

    if(distance2Head <= distance2Tail && distance2Head <= distance2Current){//head 2 element
        iftGenericLinkedListNode *currentNode = list->head;
        for (size_t i = 0; i < list->length; ++i) {
            if(i == index){
                iftSetCursorPosition(list,i,currentNode);
                return currentNode->data;
            }else{
                currentNode = currentNode->next;
            }

        }
    }else if(distance2Tail <= distance2Current) {//tail 2 element
        iftGenericLinkedListNode *currentNode = list->tail;
        for (int i = list->length-1; i >= 0; --i) {
            if(i == (int)index){
                iftSetCursorPosition(list,i,currentNode);
                return currentNode->data;
            }else{
                currentNode = currentNode->previous;
            }
        }
    }else{//current 2 element
        if(currentDirection){//element is back
            iftGenericLinkedListNode *currentNode = list->currentNode;
            for (int i = list->currentNodeIndex; i >= 0; --i) {
                if(i == (int)index){
                    iftSetCursorPosition(list,i,currentNode);
                    return currentNode->data;
                }else{
                    currentNode = currentNode->previous;
                }
            }
        }else{//element is front
            iftGenericLinkedListNode *currentNode = list->currentNode;
            for (size_t i = list->currentNodeIndex; i < list->length; ++i) {
                if(i == index){
                    iftSetCursorPosition(list,i,currentNode);
                    return currentNode->data;
                }else{
                    currentNode = currentNode->next;
                }
            }
        }
    }

    //unlikely hit this bit
    return NULL;
}



iftGenericLinkedListNode* iftGetLinkedListNode(iftGenericLinkedList *list, size_t index);
void iftDestroyLinkedList(iftGenericLinkedList **list);
void iftDestroyNodeList(iftGenericLinkedList* list,iftGenericLinkedListNode **node);
void iftClearLinkedList(iftGenericLinkedList* list);

//ITERATOR
iftGenericListIterator* iftGetListIteratorBegin(iftGenericLinkedList *list);
iftGenericListIterator* iftGetListIteratorEnd(iftGenericLinkedList *list);

#define LISTITERATOR_GET_NEXT_AS(type,iterator) (*((type*)getNextValueInListIterator(iterator)))
#define LISTITERATOR_GET_PREVIOUS_AS(type,iterator) (*((type*)getPreviousValueInListIterator(iterator)))
#define LISTITERATOR_GET_CURRENT_AS(type,iterator) (*((type*)getValueInListIterator(iterator)))

static inline void* iftGetValueInListIterator(iftGenericListIterator* iterator){
    return iterator->node->data;
}

static inline void* iftGetNextValueInListIterator(iftGenericListIterator* iterator){
    void* data = iterator->node->data;
    iterator->node = iterator->node->next;
    return data;
}

static inline void* iftGetPreviousValueInListIterator(iftGenericListIterator* iterator){
    void* data = iterator->node->data;
    iterator->node = iterator->node->previous;
    return data;
}

static inline void* iftHasNextInListIterator(iftGenericListIterator* iterator){
    return iterator->node;
}



size_t iftGetIteratorIndexInList(iftGenericLinkedList* list, iftGenericListIterator* iterator);


#endif //_LIST_H_
