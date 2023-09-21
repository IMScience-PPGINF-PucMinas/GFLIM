//
// Created by deangeli on 5/1/17.
//

#ifndef _IFTGENERICVECTOR_H_
#define _IFTGENERICVECTOR_H_

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "iftCommon.h"


typedef struct _iftGenericVector {
    size_t size;
    size_t capacity;
    size_t elementSize;
  float growthFactor;
  FreeFunction freeFunction;

    void* data;
} iftGenericVector;



typedef struct _iftGenericVectorIterator {
    void* pointer;
    size_t elementSize;
} iftGenericVectorIterator;


/*Creation and Destruction*/
iftGenericVector* iftCreateGenericVector(size_t capacity, size_t elementSize);
iftGenericVector* iftCreateGenericNullVector(size_t capacity, size_t elementSize);
void iftDestroyGenericVector(iftGenericVector **vector);
void iftDestroyVectorVoidPointer(void *vector);
iftGenericVector* iftCopyVector(iftGenericVector *vector);
iftGenericVector* iftCopyVector(iftGenericVector *vector);
/*******************************************/


/*Capacity and bounds*/
#define IS_INDEX_OUT_BOUNDS(vector,i) (i >=vector->size)
#define IS_INDEX_OUT_BOUNDS_INSERT(vector,i) (i > vector->size)
#define SHOULD_GROW(vector) (vector->size >= vector->capacity)



//private
//used to grow or Shrink the vector
void iftSetVectorCapacity(iftGenericVector* vector, size_t newCapacity);
void iftResizeVector(iftGenericVector* vector, size_t newSize);
void iftShrinkToFit(iftGenericVector* vector);
/*******************************************/

/*Element access*/
#define iftVectorAt(type, vectorPointer, index) ((type*)vectorPointer->data)[index]
#define iftVectorBack(type, vectorPointer) ((type*)vectorPointer->data)[vectorPointer->size-1]

/*******************************************/

/*Modifiers*/
#define iftGenericVectorPushBack(type, vector, element) \
        if (SHOULD_GROW(vector)){iftSetVectorCapacity(vector,vector->size*vector->growthFactor);} \
        ((type*)vector->data)[vector->size] = element; \
        vector->size++;

#define iftGenericVectorInsertAt(type, vector, element) \
        if (SHOULD_GROW(vector)){iftSetVectorCapacity(vector,vector->size*vector->growthFactor);} \
        ((type*)vector->data)[vector->size] = element; \
        vector->size++;

void iftAssignElementInVectorAt(iftGenericVector* vector,void* element, size_t index);
void iftAssignElementInVector(iftGenericVector* vector,void* element, size_t indexBegin,size_t indexEnd);
inline void iftPushBackElementInVector(iftGenericVector* vector,void* element){
    if( SHOULD_GROW(vector) ){
        iftSetVectorCapacity(vector,vector->size*vector->growthFactor);
    }
    unsigned char* offset = (unsigned char*)vector->data + (vector->elementSize*vector->size);//last position
    TRANSFER_DATA_COPY(offset,element,vector->elementSize);
    //transferData(offset,data,vector->elementSize);
    //memcpy(offset, element, vector->elementSize);
    vector->size++;
}



void iftPopBackElementInVector(iftGenericVector* vector);
void iftInsertElementInVectorAt(iftGenericVector* vector,void* element, size_t index);
void iftRemoveElementInVectorAt(iftGenericVector* vector, size_t index);
void iftRemoveElementsInVector(iftGenericVector* vector, size_t indexBegin,size_t indexEnd);
void iftSwapVectors(iftGenericVector *vector1,iftGenericVector *vector2);
void iftClearGenericVector(iftGenericVector *vector);

//private
void iftShiftVectorToRightAt(iftGenericVector* vector,size_t index);
void iftShiftVectorToRightAt2(iftGenericVector* vector,size_t index);
void iftShiftVectorToLeftAt(iftGenericVector* vector,size_t index);
void iftShiftVectorToLeft(iftGenericVector* vector,size_t indexBegin,size_t indexEnd);

/*******************************************/

/*Iterator*/
iftGenericVectorIterator* iftGetVectorIteratorBegin(iftGenericVector* vector);
iftGenericVectorIterator* iftGetVectorIteratorEnd(iftGenericVector* vector);
inline void* iftGetValueInVectorIterator(iftGenericVectorIterator* iterator){
    return iterator->pointer;
}
inline void iftIncrementVectorIterator(iftGenericVectorIterator* iterator){
    iterator->pointer = (char*)iterator->pointer + iterator->elementSize;
}
inline void iftDecrementVectorIterator(iftGenericVectorIterator* iterator){
    iterator->pointer = (char*)iterator->pointer - iterator->elementSize;
}

inline void* iftGetNextValueInVectorIterator(iftGenericVectorIterator* iterator){
    void* data = iterator->pointer;
    iftIncrementVectorIterator(iterator);
    return data;
}

inline void* iftGetPreviousValueInVectorIterator(iftGenericVectorIterator* iterator){
    void* data = iterator->pointer;
    iftDecrementVectorIterator(iterator);
    return data;
}

size_t iftGetIteratorIndexInVector(iftGenericVector* vector, iftGenericVectorIterator* iterator);

#define VECTOR_GET_ITERATOR_NEXT_AS(type, iterator) *((type*)getNextValueInVectorIterator((iterator))
#define VECTOR_GET_ITERATOR_PREVIOUS_AS(type, iterator) *((type*)getPreviousValueInVectorIterator((iterator))
#define VECTOR_GET_ITERATOR_VALUE_AS(type, iterator) *((type*)getValueInVectorIterator((iterator))

#define iftGenericVectorPrint(type, symbol, vector) \
    for(size_t currentIndexVector = 0; currentIndexVector < vector->size; currentIndexVector++){\
        printf(symbol,  VECTOR_GET_ELEMENT_AS(type,vector,currentIndexVector) ); \
    }\
    printf("\n");

/*******************************************/



#endif //BITBUCKET_VECTOR_H
