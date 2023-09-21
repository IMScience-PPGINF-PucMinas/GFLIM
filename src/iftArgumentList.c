#include "iftArgumentList.h"

#include "ift/core/io/Stream.h"


iftArgumentList* iftCreateArgumentList(){
    iftArgumentList* list = (iftArgumentList*)calloc(1,sizeof(iftArgumentList));
    list->head = list->tail = NULL;
    list->length = 0;
    list->currentNode = NULL;
    list->currentNodeIndex = -1;
    return list;
}

iftArgumentListNode* iftCreateArgumentNode(void* element, size_t dataSize){

    iftArgumentListNode* node = (iftArgumentListNode*)calloc(1,sizeof(iftArgumentListNode));
    node->data = calloc(1,dataSize);
    node->dataSize = dataSize;
    node->freeFunction = NULL;
    if(dataSize < 100){
        TRANSFER_DATA_COPY(node->data, element, dataSize);
    }else{
        memcpy(node->data, element, dataSize);
    }
    return node;
}

iftArgumentListNode* iftCreateRawArgumentNode(size_t dataSize){
    iftArgumentListNode* node = (iftArgumentListNode*)calloc(1,sizeof(iftArgumentListNode));
    node->data = calloc(1,dataSize);
    node->dataSize = dataSize;
    node->freeFunction = NULL;
    node->next = NULL;
    node->previous = NULL;
    return node;
}

void iftDestroyArgumentList(iftArgumentList **list){
    iftArgumentList* aux = *list;
    iftArgumentListNode *currentNode = NULL;
    if(*list == NULL || list == NULL){
        return;
    }

    aux->tail->next = NULL;
    aux->head->previous = NULL;

    while(aux->head != NULL) {
        currentNode = aux->head;
        aux->head = currentNode->next;
        iftDestroyNodeArgumentList(&currentNode);
        aux->length--;
    }
    free(aux);
    aux = NULL;
}

void iftDestroyNodeArgumentList(iftArgumentListNode **node){
    iftArgumentListNode *aux = *node;
    if(aux == NULL){
        return;
    }
    if(aux->freeFunction) {
        aux->freeFunction(aux->data);
    }else{
        free(aux->data);
    }
    aux->data = NULL;
    free(*node);
    aux = NULL;
}

void iftClearArgumentList(iftArgumentList* list){
    if(list == NULL){
        return;
    }
    iftArgumentListNode *currentNode = NULL;
    list->tail->next = NULL;
    list->head->previous = NULL;
    while(list->head != NULL) {
        currentNode = list->head;
        list->head = currentNode->next;
        iftDestroyNodeArgumentList(&currentNode);
        list->length--;
    }
}

