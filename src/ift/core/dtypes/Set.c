#include "ift/core/dtypes/Set.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Stream.h"
#include "iftCommon.h"


void iftInsertSet(iftSet **S, int elem)
{
    iftSet *p=NULL;
    
    p = (iftSet *) iftAlloc(1,sizeof(iftSet));
    if (p == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertSet");
    if (*S == NULL){
        p->elem  = elem;
        p->next  = NULL;
    }else{
        p->elem  = elem;
        p->next  = *S;
    }
    *S = p;
}


int iftRemoveSet(iftSet **S)
{
    iftSet *p;
    int elem = IFT_NIL;
    
    if (*S != NULL){
        p    =  *S;
        elem = p->elem;
        *S   = p->next;
        iftFree(p);
    }
    
    return(elem);
}


void iftRemoveSetElem(iftSet **S, int elem){
    if (S == NULL || *S == NULL)
        return;
    
    iftSet *tmp = *S;
    
    if (tmp->elem == elem) {
        *S = tmp->next;
        iftFree(tmp);
    } else {
        while (tmp->next != NULL && tmp->next->elem != elem)
            tmp = tmp->next;
        if (tmp->next == NULL)
            return;
        iftSet *remove = tmp->next;
        tmp->next = remove->next;
        iftFree(remove);
    }
}


void iftDestroySet(iftSet **S)
{
    iftSet *p;
    while(*S != NULL){
        p = *S;
        *S = p->next;
        iftFree(p);
    }
    *S = NULL;
}


// Could be made efficient by sorting both sets and merging (Theta(n + mlg(m))).
// This implementation is Theta(mn)
iftSet* 	iftSetUnion(iftSet *S1,iftSet *S2){
    iftSet *S = 0;
    
    iftSet *s = S1;
    while(s){
        iftInsertSet(&S,s->elem);
        s = s->next;
    }
    
    s = S2;
    while(s){
        iftUnionSetElem(&S,s->elem);
        s = s->next;
    }
    
    return S;
}


iftSet* iftSetCopy(iftSet* S){
    return iftSetUnion(S,0);
}


// Uses recursion to copy the set in order
iftSet* iftSetCopyOrdered(const iftSet *S)
{
    iftSet *newS = NULL, *part_cpy = NULL;
    
    // For the empty set, return NULL
    if(S == NULL)
        return NULL;
    
    // Recrusively copy the n-1 elements of the set
    // excluding the current one
    part_cpy = iftSetCopyOrdered(S->next);
    
    // Copy the current element to the new set in the
    // appropriate position, given that newS
    // for now is NULL and S->elem will be inserted
    // at the beginning of the new set. Since
    // iftInsertSet always inserts at the beginning
    // we are ensuring that each element will be copied
    // backwards
    iftInsertSet(&newS, S->elem);
    
    // Attached the copied set after the newly added
    // element of newS and return it
    newS->next = part_cpy;
    
    return newS;
}


//Doesn't check for duplicates (faster than iftSetUnion)
iftSet* 	iftSetConcat(iftSet *S1,iftSet *S2){
    iftSet *S = 0;
    
    iftSet *s = S1;
    while(s){
        iftInsertSet(&S,s->elem);
        s = s->next;
    }
    
    s = S2;
    while(s){
        iftInsertSet(&S,s->elem);
        s = s->next;
    }
    
    return S;
}


char iftUnionSetElem(iftSet **S, int elem)
{
    iftSet *aux=*S;
    
    while (aux != NULL) {
        if (aux->elem == elem)
            return(0);
        aux = aux->next;
    }
    iftInsertSet(S,elem);
    return(1);
}


void iftInvertSet(iftSet **S)
{
    iftSet *set=NULL;
    
    while (*S != NULL)
        iftInsertSet(&set,iftRemoveSet(S));
    *S = set;
}


void iftReverseSet(iftSet **S)
{
    iftSet *prev = NULL, *next = NULL, *cur = *S;
    while (cur)
    {
        next = cur->next;
        cur->next = prev;
        prev = cur;
        cur = next;
    }
    *S = prev;
}



int iftSetHasElement(iftSet *S, int elem){
    iftSet *s = S;
    while(s){
        if(s->elem == elem)
            return 1;
        
        s = s->next;
    }
    
    return 0;
}


int iftSetSize(const iftSet* S){
    const iftSet *s = S;
    
    int i = 0;
    while (s != NULL){
        i++;
        s = s->next;
    }
    
    return i;
    
}


iftIntArray *iftSetToArray(iftSet *S) {
    int n_elems = iftSetSize(S);
    iftIntArray *array = iftCreateIntArray(n_elems);
    
    iftSet *sp = S;
    int i      = 0;
    while (sp != NULL) {
        array->val[i++] = sp->elem;
        sp = sp->next;
    }
    
    return array;
}


iftSet *iftIntArrayToSet(iftIntArray *array) {
    iftSet *S = NULL;
    
    // it is copied from back to front to keep the original order of the elements
    // because the Set insertion is always in the head
    for (long i = array->n - 1; i >= 0; i--) {
        iftInsertSet(&S, array->val[i]);
    }
    
    return S;
}


void iftPrintSet(iftSet *S)
{
    while (S != NULL) {
        printf("elem %d \n",S->elem);
        S = S->next;
    }
    
}


