#include "ift/core/dtypes/LabeledSet.h"

#include "ift/core/io/Stream.h"


int  iftNumberOfLabels(const iftLabeledSet *S)
{
  int c=IFT_NIL; /* number of objects */
  int incr = 0;
  
  while (S != NULL) {
    int label = S->label;
    if (label == 0)
      incr = 1;
    if (label > c)
      c = label;
    S = S->next;
  }

  
  return(c+incr);
}

int  iftNumberOfMarkers(const iftLabeledSet *S)
{
	int c=IFT_NIL; /* number of objects */

	while (S != NULL) {
		int marker = S->marker;
		if (marker > c)
			c = marker;
		S = S->next;
	}

	return(c);
}


void iftInsertLabeledSet(iftLabeledSet **S, int elem, int label)
{
  iftLabeledSet *p=NULL;

  p = (iftLabeledSet *) iftAlloc(1,sizeof(iftLabeledSet));
  if (p == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertLabeledSet");
  if (*S == NULL){
	p->elem     = elem;
	p->label    = label;
	p->marker   = IFT_NIL;
	p->next     = NULL;    
	p->handicap = 0;
  }else{
	p->elem     = elem;
	p->label    = label;
	p->marker   = IFT_NIL;
	p->next     = *S;
	p->handicap = 0;
  }
  *S = p;
}
void iftInsertLabeledSetMarkerAndHandicap(iftLabeledSet **S, int elem, int label, int marker, int handicap)
{
  iftLabeledSet *p=NULL;

  p = (iftLabeledSet *) iftAlloc(1,sizeof(iftLabeledSet));
  if (p == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertLabeledSet");
  if (*S == NULL)
  {
    p->elem     = elem;
    p->label    = label;
    p->marker   = marker;
    p->handicap = handicap;
    p->next     = NULL;
  }
  else{
    p->elem     = elem;
    p->label    = label;
    p->next     = *S;
    p->marker   = marker;
    p->handicap = handicap;
  }
  *S = p;
}

int iftRemoveLabeledSet(iftLabeledSet **S, int *label)
{
  iftLabeledSet *p;
  int elem = IFT_NIL;

  if (*S != NULL){
    p    =  *S;
    elem = p->elem;
    (*label) = p->label;
    *S   = p->next;
    iftFree(p);
  }

  return(elem);
}


void iftRemoveLabeledSetElem(iftLabeledSet **S, int elem){
  iftLabeledSet *tmp = NULL, *remove;

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


void iftDestroyLabeledSet(iftLabeledSet **S)
{
  iftLabeledSet *p;
  while(*S != NULL){
    p = *S;
    *S = p->next;
    iftFree(p);
  }
}


void iftInsertSetIntoLabeledSet(iftSet **S, int label, iftLabeledSet **T) {
    iftSet *node = *S;
    while (node != NULL) {
        int p = iftRemoveSet(&node);
        iftInsertLabeledSet(T, p, label);
    }
}


iftLabeledSet *iftTranslateLabeledSet(const iftLabeledSet *S, const iftImage *img, iftVector disp_vec) {
    iftLabeledSet *T = NULL;

    const iftLabeledSet *node = S;
    while (node != NULL) {
        iftVoxel v       = iftGetVoxelCoord(img, node->elem);
        iftVoxel trans_v = iftVectorSum(v, disp_vec);
        
        if (iftValidVoxel(img, trans_v)) {
            int p = iftGetVoxelIndex(img, trans_v);
            iftInsertLabeledSet(&T, p, node->label);
        }
        node = node->next;
    }

    return T;
}


iftSet *iftTranslateSet(const iftSet *S, const iftImage *img, iftVector disp_vec) {
    iftSet *T = NULL;

    const iftSet *node = S;
    while (node != NULL) {
        iftVoxel v       = iftGetVoxelCoord(img, node->elem);
        iftVoxel trans_v = iftVectorSum(v, disp_vec);

        if (iftValidVoxel(img, trans_v)) {
            int p = iftGetVoxelIndex(img, trans_v);
            iftInsertSet(&T, p);
        }
        node = node->next;
    }

    return T;
}


void iftConcatLabeledSet(iftLabeledSet **S1, iftLabeledSet **S2){
	if(*S2 == NULL)
		return;

	iftLabeledSet *i = *S2;
	while(i != NULL){
		iftInsertLabeledSetMarkerAndHandicap(S1,i->elem,i->label,i->marker, i->handicap);
		i = i->next;
	}
}

/* Warning: S2 must be a subset of S1! */

void iftRemoveSubsetLabeledSet(iftLabeledSet **S1,iftLabeledSet **S2){
	if(*S2 == NULL)
		return;

	iftLabeledSet *i = *S2;
	while(i != NULL){
		iftRemoveLabeledSetElem(S1,i->elem);
		i = i->next;
	}
}

iftLabeledSet* iftCopyLabeledSet(iftLabeledSet *s){
	iftLabeledSet *lset = NULL;

	while(s != NULL){
		iftInsertLabeledSetMarkerAndHandicap(&lset,s->elem,s->label,s->marker,s->handicap);
		s = s->next;
	}

	return lset;
}

iftLabeledSet* iftCopyOrderedLabeledSet(iftLabeledSet *s){
	iftLabeledSet *lset = NULL;
	if(s == NULL)
		return lset;

	iftLabeledSet *i = s;
	int nelem = 0;
	while(i != NULL){
		nelem++;
		i = i->next;
	}

	int *elems = iftAllocIntArray(nelem);
	int *labels = iftAllocIntArray(nelem);
	int *markers = iftAllocIntArray(nelem);
	int *handicap = iftAllocIntArray(nelem);
	i = s;
	int index = 0;
	while(i != NULL){
		elems[index] = i->elem;
		labels[index] = i->label;
		markers[index] = i->marker;
		handicap[index] = i->handicap;
		index++;
		i = i->next;
	}

	for(index = nelem - 1; index >= 0; index--){
		iftInsertLabeledSetMarkerAndHandicap(&lset,elems[index],labels[index],markers[index], handicap[index]);
	}
	iftFree(markers);
	iftFree(handicap);
	iftFree(elems);
	iftFree(labels);

	return lset;
}

iftLabeledSet* iftReverseLabeledSet(iftLabeledSet *s){
	iftLabeledSet *lset = NULL;

	while(s != NULL){
		iftInsertLabeledSetMarkerAndHandicap(&lset,s->elem,s->label,s->marker,s->handicap);
		s = s->next;
	}

	return lset;
}

int iftLabeledSetSize(const iftLabeledSet *s)
{
	const iftLabeledSet *aux = s;
	int counter = 0;
	while (aux != NULL)
	{
		counter++;
		aux = aux->next;
	}
	return counter;
}


iftSet* iftLabeledSetToSet(iftLabeledSet *S, int lb){
	iftSet *Snew = NULL;
	iftLabeledSet *s = S;

	while(s != NULL){
		if(lb == s->label)
			iftInsertSet(&Snew, s->elem);
		s = s->next;
	}

	return Snew;
}

iftSet* iftLabeledSetElemsToSet(iftLabeledSet *S) {
	iftSet *Snew = NULL;
	iftLabeledSet *s = S;

	while(s != NULL){
		iftInsertSet(&Snew, s->elem);
		s = s->next;
	}

	return Snew;
}

iftLabeledSet* iftCopyLabels(iftLabeledSet *S, int lb){
	iftLabeledSet *Snew = NULL;
	iftLabeledSet *s = S;

	while(s != NULL){
		if(lb == s->label)
			iftInsertLabeledSet(&Snew, s->elem, lb);
		s = s->next;
	}

	return Snew;
}

iftSet* iftLabeledSetMarkersToSet(iftLabeledSet *S, int marker){
	iftSet *Snew = NULL;
	iftLabeledSet *s = S;

	while(s != NULL){
		if(marker == s->marker)
			iftInsertSet(&Snew, s->elem);
		s = s->next;
	}

	return Snew;
}

iftLabeledSet* iftCopyLabeledSetMarkers(iftLabeledSet *S, int marker) {
	iftLabeledSet *Snew = NULL;
	iftLabeledSet *s = S;

	while(s != NULL){
		if(marker == s->marker)
			iftInsertLabeledSetMarkerAndHandicap(&Snew, s->elem, s->label, s->marker, s->handicap);
		s = s->next;
	}

	return Snew;
}

int iftLabeledSetHasElement(iftLabeledSet *S, int elem) {
  iftLabeledSet *s = S;
  while(s){
    if(s->elem == elem)
      return 1;
    
    s = s->next;
  }

  return 0;
}

void iftPrintLabeledSet(iftLabeledSet *S)
{
  while (S != NULL) {
    printf("elem %d label %d\n",S->elem,S->label);
    S = S->next;
  }
}


void iftInsertLabeledSetIntoImage(const iftLabeledSet *S, iftImage *img) {
    while (S != NULL) {
        img->val[S->elem] = S->label;
        S = S->next;
    }
}

char iftUnionLabeledSetElem(iftLabeledSet **S, int elem, int label)
{
    iftLabeledSet *aux=*S;
    
    while (aux != NULL) {
        if (aux->elem == elem)
            return(0);
        aux = aux->next;
    }
    iftInsertLabeledSet(S,elem, label);
    return(1);
}

char iftUnionLabeledSetElemMarkerAndHandicap(iftLabeledSet **S, int elem, int label, int marker, int handicap)
{
    iftLabeledSet *aux=*S;
    
    while (aux != NULL) {
        if (aux->elem == elem)
            return(0);
        aux = aux->next;
    }
    iftInsertLabeledSetMarkerAndHandicap(S,elem, label, marker, handicap);
    return(1);
}









