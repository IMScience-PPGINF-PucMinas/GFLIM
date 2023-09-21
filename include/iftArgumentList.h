//
// Created by deangeli on 5/1/17.
//

#ifndef _IFTARGUMENTLIST_H_
#define _IFTARGUMENTLIST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftDataTransference.h"

typedef struct _iftArgumentListNode {
    void *data;
    size_t dataSize;
    struct _iftArgumentListNode *next;
    struct _iftArgumentListNode *previous;
    char typeName[50];
    FreeFunction freeFunction;
} iftArgumentListNode;

typedef struct _iftArgumentList{
    size_t length;
    iftArgumentListNode *head;
    iftArgumentListNode *tail;
    iftArgumentListNode *currentNode;
    int currentNodeIndex;
} iftArgumentList;


iftArgumentList* iftCreateArgumentList();
iftArgumentListNode* iftCreateArgumentNode(void* element, size_t dataSize);
iftArgumentListNode* iftCreateRawArgumentNode(size_t dataSize);

static inline void  iftSetCursorPositionInArgumentList(iftArgumentList* list, int index,iftArgumentListNode* node){
    list->currentNodeIndex = index;
    list->currentNode = node;
}

static inline void iftAppendElementInArgumentList(iftArgumentList *argumentList, void *element, size_t dataSize){
    iftArgumentListNode* node = iftCreateArgumentNode(element,dataSize);

    if(argumentList->length == 0) {
        argumentList->head = argumentList->tail = node;
    } else {
        argumentList->tail->next = node;
        node->previous = argumentList->tail;
        argumentList->tail = node;
    }

    argumentList->length++;
}


static inline void iftPushBackElementInArgumentList(iftArgumentList *argumentList, void *element, size_t dataSize){
    iftAppendElementInArgumentList(argumentList,element,dataSize);
}

static inline void iftPushBackNodeInArgumentList(iftArgumentList *argumentList, iftArgumentListNode* node){
    if(argumentList->length == 0) {
        argumentList->head = argumentList->tail = node;
    } else {
        argumentList->tail->next = node;
        node->previous = argumentList->tail;
        argumentList->tail = node;
    }

    argumentList->length++;
}
static inline void iftPrependElementInArgumentList(iftArgumentList *argumentList, void *element, size_t dataSize){
    iftArgumentListNode* node = iftCreateArgumentNode(element,dataSize);

    if(argumentList->length == 0) {
        argumentList->head = node;
        argumentList->tail = argumentList->head;
    }else{
        node->next = argumentList->head;
        argumentList->head->previous = node;
        argumentList->head = node;
    }

    argumentList->length++;
}

static inline void iftPushFrontElementInArgumentList(iftArgumentList *argumentList, void *element, size_t dataSize){
    iftPrependElementInArgumentList(argumentList,element,dataSize);
}

#define ARGLIST_PUSH_BACK_AS(type,list,element) {\
    ArgumentListNode* node = createRawArgumentNode(sizeof(type)); \
    *( (type*)node->data) = element; \
    pushBackElementInArgumentList(list,node); \
    }

#define ARGLIST_PUSH_BACK_NAMED_NODE_AS(type,list,element,name) {\
    ArgumentListNode* node = createRawArgumentNode(sizeof(type)); \
    size_t stringSize = 0; \
    while(name[stringSize] != '\0'){\
        stringSize++;\
    }\
    for (size_t i = 0; i < stringSize; ++i) {\
        node->typeName[i] = name[i];\
    }\
    node->typeName[stringSize] = '\0';\
    *( (type*)node->data) = element; \
    pushBackElementInArgumentList(list,node); \
}


static inline void* iftGetElementInArgumentList(iftArgumentList *list, size_t index){
    if(index >= list->length){
        iftError("invalid position %lu. The list length is %lu (indexing start from 0)\n","iftGetElementInArgumentList",index,list->length);
        return NULL;
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
        iftArgumentListNode *currentNode = list->head;
        for (size_t i = 0; i < list->length; ++i) {
            if(i == index){
                iftSetCursorPositionInArgumentList(list,i,currentNode);
                return currentNode->data;
            }else{
                currentNode = currentNode->next;
            }

        }
    }else if(distance2Tail <= distance2Current) {//tail 2 element
        iftArgumentListNode *currentNode = list->tail;
        for (int i = list->length-1; i >= 0; --i) {
            if(i == (int)index){
                iftSetCursorPositionInArgumentList(list,i,currentNode);
                return currentNode->data;
            }else{
                currentNode = currentNode->previous;
            }
        }
    }else{//current 2 element
        if(currentDirection){//element is back
            iftArgumentListNode *currentNode = list->currentNode;
            for (int i = list->currentNodeIndex; i >= 0; --i) {
                if(i == (int)index){
                    iftSetCursorPositionInArgumentList(list,i,currentNode);
                    return currentNode->data;
                }else{
                    currentNode = currentNode->previous;
                }
            }
        }else{//element is front
            iftArgumentListNode *currentNode = list->currentNode;
            for (size_t i = list->currentNodeIndex; i < list->length; ++i) {
                if(i == index){
                    iftSetCursorPositionInArgumentList(list,i,currentNode);
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

static inline void* iftGetElementInArgumentListByName(iftArgumentList *list, const char* name){
    iftArgumentListNode *currentNode = list->head;
    for (size_t i = 0; i < list->length; ++i) {
        if( strcmp(name,currentNode->typeName) == 0){
            printf("%s\n",currentNode->typeName);
            break;
        }
        currentNode = currentNode->next;
    }
    if(currentNode == NULL){
        iftWarning("Could not find name","iftGetElementInArgumentListByName");
    }
    return currentNode->data;
}

#define ARGLIST_GET_ELEMENT_AS(type,list,index) (*((type*)iftGetElementInArgumentList(list,index)))

#define ARGLIST_GET_NAMED_ELEMENT_AS(type,list,name) \
    (*((type*)iftGetElementInArgumentListByName(list,name)))

void iftDestroyArgumentList(iftArgumentList **list);
void iftDestroyNodeArgumentList(iftArgumentListNode **node);
void iftClearArgumentList(iftArgumentList* list);

#ifdef __cplusplus
}
#endif


#endif //_ARGUMENTLIST_H_
