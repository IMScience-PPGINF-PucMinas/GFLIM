//
// Created by deangeli on 4/30/17.
//

#include "iftGenericLinkedList.h"

#include "ift/core/io/Stream.h"


iftGenericLinkedList* iftcreateLinkedList(size_t dataSize){
    if(dataSize == 0){
        printf("[createList] data size is 0\n");
        return NULL;
    }
    iftGenericLinkedList* list = (iftGenericLinkedList*)calloc(1,sizeof(iftGenericLinkedList));
    list->dataSize = dataSize;
    list->head = list->tail = NULL;
    list->freeFunction = NULL;
    list->length = 0;
    list->currentNode = NULL;
    list->currentNodeIndex = -1;
    list->isCircular = false;
    list->comparatorFunction = NULL;
    return list;
}

void iftRemoveGenericListHead(iftGenericLinkedList *list){
    if(list->head == NULL){
        return;
    }
    iftGenericLinkedListNode *currentNode = NULL;
    currentNode = list->head;
    list->head = currentNode->next;
    iftDestroyNodeList(list,&currentNode);
    iftSetCursorPosition(list,0,list->head);
    if(list->isCircular){
        list->tail->next = list->head;
        list->head->previous = list->tail;
    }
}

void iftRemoveGenericListTail(iftGenericLinkedList *list){
    if(list->head == NULL){
        return;
    }

    iftGenericLinkedListNode *previousNode = list->tail->previous;
    iftGenericLinkedListNode *currentTail = list->tail;
    iftDestroyNodeList(list,&currentTail);
    list->tail = previousNode;
    iftSetCursorPosition(list,list->length,list->tail);
    //setCursorPosition(list,-1,NULL);
    if(list->isCircular){
        list->tail->next = list->head;
        list->head->previous = list->tail;
    }

    return;
}



void iftInsertElementInListAt(iftGenericLinkedList *list, void *element, size_t index){
    if(index > list->length){
        printf("[insertElementInListAt] invalid position %lu. The list length is %lu (indexing start from 0)\n", index,list->length);
        return;
    }

    if (index == 0){
        iftPrependElementInList(list, element);
        return;
    }
    if(index == list->length){
        iftAppendElementInList(list, element);
        return;
    }

    iftGenericLinkedListNode *node = (iftGenericLinkedListNode*)malloc(sizeof(iftGenericLinkedListNode));
    node->data = malloc(list->dataSize);
    //memcpy(node->data, element, list->dataSize);
    TRANSFER_DATA_COPY(node->data, element, list->dataSize);

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
                node->previous = currentNode->previous;
                node->next = currentNode;
                currentNode->previous->next = node;
                currentNode->previous = node;
                list->length++;
                iftSetCursorPosition(list,i,node);
                return;
            }else{
                currentNode = currentNode->next;
            }

        }
    }else if(distance2Tail <= distance2Current){//tail 2 element
        iftGenericLinkedListNode *currentNode = list->tail;
        for (int i = list->length-1; i >= 0; --i) {
            if(i == (int)index){
                node->previous = currentNode->previous;
                node->next = currentNode;
                currentNode->previous->next = node;
                currentNode->previous = node;
                list->length++;
                iftSetCursorPosition(list,i,node);
                return;
            }else{
                currentNode = currentNode->previous;
            }
        }
    }else{//current 2 element
        if(currentDirection){//element is back
            iftGenericLinkedListNode *currentNode = list->currentNode;
            for (int i = list->currentNodeIndex; i >= 0; --i) {
                if(i == (int)index){
                    node->previous = currentNode->previous;
                    node->next = currentNode;
                    currentNode->previous->next = node;
                    currentNode->previous = node;
                    list->length++;
                    iftSetCursorPosition(list,i,node);
                    return;
                }else{
                    currentNode = currentNode->previous;
                }
            }
        }else{//element is front
            iftGenericLinkedListNode *currentNode = list->currentNode;
            for (size_t i = list->currentNodeIndex; i < list->length; ++i) {
                if(i == index){
                    node->previous = currentNode->previous;
                    node->next = currentNode;
                    currentNode->previous->next = node;
                    currentNode->previous = node;
                    list->length++;
                    iftSetCursorPosition(list,i,node);
                    return;
                }else{
                    currentNode = currentNode->next;
                }
            }
        }

    }
}

void iftRemoveElementInGenericListAt(iftGenericLinkedList *list, size_t index){
    if(index >= list->length){
        printf("[insertElementInListAt] invalid position %lu. The list length is %lu (indexing start from 0)\n", index,list->length);
        return;
    }

    if (index == 0){
        iftRemoveGenericListHead(list);
        return;
    }
    if(index == list->length-1){
        iftRemoveGenericListTail(list);
        return;
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
                currentNode->previous->next = currentNode->next;
                currentNode->next->previous = currentNode->previous;
                iftSetCursorPosition(list,i,currentNode->next);
                iftDestroyNodeList(list,&currentNode);//already decreasing list size
                return;
            }else{
                currentNode = currentNode->next;
            }

        }
    }else if(distance2Tail <= distance2Current) {//tail 2 element
        iftGenericLinkedListNode *currentNode = list->tail;
        for (int i = list->length-1; i >= 0; --i) {
            if(i == (int)index){
                currentNode->previous->next = currentNode->next;
                currentNode->next->previous = currentNode->previous;
                //setCursorPosition(list,-1,NULL);
                iftDestroyNodeList(list,&currentNode);//already decreasing list size
                return;
            }else{
                currentNode = currentNode->previous;
            }
        }
    }else{//current 2 element
        if(currentDirection){//element is back
            iftGenericLinkedListNode *currentNode = list->currentNode;
            for (int i = list->currentNodeIndex; i >= 0; --i) {
                if(i == (int)index){
                    currentNode->previous->next = currentNode->next;
                    currentNode->next->previous = currentNode->previous;
                    iftSetCursorPosition(list,i,currentNode->next);
                    iftDestroyNodeList(list,&currentNode);//already decreasing list size
                    return;
                }else{
                    currentNode = currentNode->previous;
                }
            }
        }else{//element is front
            iftGenericLinkedListNode *currentNode = list->currentNode;
            for (size_t i = list->currentNodeIndex; i < list->length; ++i) {
                if(i == index){
                    currentNode->previous->next = currentNode->next;
                    currentNode->next->previous = currentNode->previous;
                    iftSetCursorPosition(list,i,currentNode->next);
                    iftDestroyNodeList(list,&currentNode);//already decreasing list size
                    return;
                }else{
                    currentNode = currentNode->next;
                }
            }
        }
    }

}

void iftRemoveElementInGenericListByReference(iftGenericLinkedList *list, void *element){
    if(list->head->data == element){
        iftRemoveGenericListHead(list);
        return;
    }
    if(list->tail->data == element){
        iftRemoveGenericListTail(list);
        return;
    }

    iftGenericLinkedListNode *currentNode = list->head;
    for (size_t i = 0; i < list->length; ++i) {
        if(currentNode->data == element){
            currentNode->previous->next = currentNode->next;
            currentNode->next->previous = currentNode->previous;
            iftSetCursorPosition(list,i,currentNode->next);
            iftDestroyNodeList(list,&currentNode);//already decreasing list size
            return;
        }else{
            currentNode = currentNode->next;
        }
    }
    printf("[removeElement] could not find the pointer in the list\n");
}

void iftRemoveElementInGenericListGivenValue(iftGenericLinkedList *list, void *element){
    if(!list->comparatorFunction){
        printf("[removeElementInListGivenValue] compare function not defined\n");
    }

    if(list->head->data == element){
        iftRemoveGenericListHead(list);
        return;
    }
    if(list->tail->data == element){
        iftRemoveGenericListTail(list);
        return;
    }

    iftGenericLinkedListNode *currentNode = list->head;
    for (size_t i = 0; i < list->length; ++i) {
        if( list->comparatorFunction(currentNode->data,element)){
            currentNode->previous->next = currentNode->next;
            currentNode->next->previous = currentNode->previous;
            iftSetCursorPosition(list,i,currentNode->next);
            iftDestroyNodeList(list,&currentNode);//already decreasing list size
            return;
        }else{
            currentNode = currentNode->next;
        }
    }
}

void iftRemoveElementsInGenericListGivenValue(iftGenericLinkedList *list, void *element){
    if(!list->comparatorFunction){
        printf("[removeElementInListGivenValue] compare function not defined\n");
    }

    if(list->head->data == element){
        iftRemoveGenericListHead(list);
        return;
    }
    if(list->tail->data == element){
        iftRemoveGenericListTail(list);
        return;
    }

    iftGenericLinkedListNode *currentNode = list->head;
    for (size_t i = 0; i < list->length; ++i) {
        if( list->comparatorFunction(currentNode->data,element)){
            currentNode->previous->next = currentNode->next;
            currentNode->next->previous = currentNode->previous;
            iftSetCursorPosition(list,i,currentNode->next);
            iftDestroyNodeList(list,&currentNode);//already decreasing list size
        }else{
            currentNode = currentNode->next;
        }
    }
}

iftGenericLinkedListNode* iftGetLinkedListNode(iftGenericLinkedList *list, size_t  index){
    if(index >= list->length){
        if(list->isCircular){
            index = index % list->length;
        }else{
            printf("[getElement] invalid position %lu. The list length is %lu (indexing start from 0)\n", index,list->length);
            return NULL;
        }
    }

    if (index == 0){
        return list->head;
    }
    if(index == list->length-1){
        return list->tail;
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
                return currentNode;
            }else{
                currentNode = currentNode->next;
            }

        }
    }else if(distance2Tail <= distance2Current) {//tail 2 element
        iftGenericLinkedListNode *currentNode = list->tail;
        for (int i = list->length-1; i >= 0; --i) {
            if(i == (int)index){
                iftSetCursorPosition(list,i,currentNode);
                return currentNode;
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
                    return currentNode;
                }else{
                    currentNode = currentNode->previous;
                }
            }
        }else{//element is front
            iftGenericLinkedListNode *currentNode = list->currentNode;
            for (size_t i = list->currentNodeIndex; i < list->length; ++i) {
                if(i == index){
                    iftSetCursorPosition(list,i,currentNode);
                    return currentNode;
                }else{
                    currentNode = currentNode->next;
                }
            }
        }
    }

    //unlikely hit this bit
    return NULL;
}

void iftDestroyNodeList(iftGenericLinkedList* list,iftGenericLinkedListNode **node){
    iftGenericLinkedListNode *aux = *node;
    if(aux == NULL){
        return;
    }
    if(list->freeFunction) {
        list->freeFunction(aux->data);
        aux->data = NULL;
    }
    free(aux);
    aux = NULL;
    list->length--;
}

void iftClearLinkedList(iftGenericLinkedList* list){
    iftGenericLinkedListNode *currentNode = NULL;
    if(list->isCircular){
        list->tail->next = NULL;
        list->head->previous = NULL;
    }
    while(list->head != NULL) {
        currentNode = list->head;
        list->head = currentNode->next;
        iftDestroyNodeList(list,&currentNode);
    }
    iftSetCursorPosition(list,-1,NULL);
}

void iftDestroyLinkedList(iftGenericLinkedList **list)
{

    iftGenericLinkedList* aux = *list;
    iftGenericLinkedListNode *currentNode = NULL;
    if(aux == NULL){
        return;
    }

    aux->tail->next = NULL;
    aux->head->previous = NULL;

    while(aux->head != NULL) {
        currentNode = aux->head;
        aux->head = currentNode->next;
        iftDestroyNodeList(aux,&currentNode);
    }
    free(aux);
    aux = NULL;
}

iftGenericListIterator* iftGetListIteratorBegin(iftGenericLinkedList *list){
    iftGenericListIterator* iterator = (iftGenericListIterator*)calloc(1,sizeof(iftGenericListIterator));
    iterator->node = list->head;
    return iterator;
}

iftGenericListIterator* iftGetListIteratorEnd(iftGenericLinkedList *list){
    iftGenericListIterator* iterator = (iftGenericListIterator*)calloc(1,sizeof(iftGenericListIterator));
    iterator->node = list->tail;
    return iterator;
}

size_t iftGetIteratorIndexInList(iftGenericLinkedList* list, iftGenericListIterator* iterator){
    iftGenericLinkedListNode* node = list->head;
    for (size_t i = 0; i < list->length; ++i) {
        if(iterator->node == node){
            return i;
        }else{
            node = node->next;
        }
    }
    return -1;
}
