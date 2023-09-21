#include "ift/core/dtypes/FSet.h"

#include "ift/core/io/Stream.h"
#include "iftCommon.h"


void iftInsertFSet(iftFSet **S, float elem)
{
    iftFSet *p=NULL;
    
    p = (iftFSet *) iftAlloc(1,sizeof(iftFSet));
    if (p == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertFSet");
    if (*S == NULL){
        p->elem  = elem;
        p->next  = NULL;
    }else{
        p->elem  = elem;
        p->next  = *S;
    }
    *S = p;
}

float iftRemoveFSet(iftFSet **S)
{
    iftFSet *p;
    float elem = IFT_NIL;
    
    if (*S != NULL){
        p    =  *S;
        elem = p->elem;
        *S   = p->next;
        iftFree(p);
    }
    
    return(elem);
}


void iftRemoveFSetElem(iftFSet **S, float elem){
    iftFSet *tmp = NULL, *remove;
    
    tmp = *S;
    if ( tmp->elem == elem ) {
        *S = tmp->next;
        iftFree( tmp );
    }
    else {
        while( tmp->next->elem != elem ) tmp = tmp->next;
        remove = tmp->next;
        tmp->next = remove->next;
        iftFree( remove );
    }
}


void iftDestroyFSet(iftFSet **S)
{
    iftFSet *p;
    while(*S != NULL){
        p = *S;
        *S = p->next;
        iftFree(p);
    }
}


