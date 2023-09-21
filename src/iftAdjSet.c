#include "iftAdjSet.h"

#include "ift/core/io/Stream.h"

/**
@file
@brief A description is missing here
*/
void iftInsertAdjSet(iftAdjSet **S, int node, float arcw)
{
  iftAdjSet *p=NULL;

  p = (iftAdjSet *) iftAlloc(1,sizeof(iftAdjSet));
  if (p == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertAdjSet");
  if (*S == NULL){
    p->node  = node;
    p->arcw  = arcw;
    p->next  = NULL;
  }else{
    p->node  = node;
    p->arcw  = arcw;
    p->next  = *S;
  }
  *S = p;
}

int iftRemoveAdjSet(iftAdjSet **S, float *arcw)
{
  iftAdjSet *p;
  int node = IFT_NIL;

  if (*S != NULL){
    p    =  *S;
    node = p->node;
    (*arcw) = p->arcw;
    *S   = p->next;
    iftFree(p);
  }

  return(node);
}


void iftRemoveAdjSetNode(iftAdjSet **S, int node){
  iftAdjSet *tmp = NULL, *remove;

  tmp = *S;
  if ( tmp->node == node ) {
    *S = tmp->next;
    iftFree( tmp );
  }
  else {
    while( tmp->next->node != node ) tmp = tmp->next;
    remove = tmp->next;
    tmp->next = remove->next;
    iftFree( remove );
  }
}


void iftDestroyAdjSet(iftAdjSet **S)
{
  iftAdjSet *p;
  while(*S != NULL){
    p = *S;
    *S = p->next;
    iftFree(p);
  }
}


// Uses recursion to copy the set in order
iftAdjSet* iftAdjSetCopyOrdered(iftAdjSet *S)
{
	iftAdjSet *newS = NULL, *part_cpy = NULL;

	// For the empty set, return NULL
	if(S == NULL)
		return NULL;

	// Recrusively copy the n-1 elements of the set
	// excluding the current one
	part_cpy = iftAdjSetCopyOrdered(S->next);

	// Copy the current element to the new set in the
	// appropriate position, given that newS
	// for now is NULL and S->elem will be inserted
	// at the beginning of the new set. Since
	// iftInsertSet always inserts at the beginning
	// we are ensuring that each element will be copied
	// backwards
	iftInsertAdjSet(&newS, S->node, S->arcw);

	// Attached the copied set after the newly added
	// element of newS and return it
	newS->next = part_cpy;

	return newS;
}



int iftAdjSetSize(const iftAdjSet *S) {
    int n = 0;
    for (const iftAdjSet *T = S; T != NULL; T = T->next, n++){}

    return n;
}


iftMatrix *iftAdjSetToMatrix(const iftAdjSet *S) {
    int nrows = iftAdjSetSize(S);

    iftMatrix *M = iftCreateMatrix(2, nrows);

    const iftAdjSet *T = S;

    for (int r = 0; r < nrows; r++) {
        iftMatrixElem(M, 0, r) = T->node;
        iftMatrixElem(M, 1, r) = T->arcw;
        T = T->next;
    }

    return M;
}


iftAdjSet *iftMatrixToAdjSet(const iftMatrix *M) {
    iftAdjSet *S = NULL;

    // copying it backwards to keep the order in the Set
    for (int r = M->nrows-1; r >= 0; r--) {
        int node = (int) iftMatrixElem(M, 0, r);
        float arcw = iftMatrixElem(M, 1, r);

        iftInsertAdjSet(&S, node, arcw);
    }

    return S;
}

int iftAdjSetHasElement(iftAdjSet *S, int elem) {
    iftAdjSet *s = S;
    while (s)
    {
        if (s->node == elem) return 1;
        s = s->next;
    }
    return 0;
}

void iftReverseAdjSet(iftAdjSet **S)
{
    iftAdjSet *prev = NULL, *next = NULL, *cur = *S;
    while (cur)
    {
        next = cur->next;
        cur->next = prev;
        prev = cur;
        cur = next;
    }
    *S = prev;
}



