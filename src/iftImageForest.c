#include "iftImageForest.h"

#include "ift/core/dtypes/FIFO.h"
#include "ift/core/io/Stream.h"


iftImageForest  *iftCreateImageForest(iftImage *img, iftAdjRel *A)
{
  iftImageForest *fst=(iftImageForest *)iftAlloc(1,sizeof(iftImageForest));

  if (fst != NULL) {
    fst->pathval        = iftCreateImage(img->xsize,img->ysize,img->zsize);
    fst->label          = iftCreateImage(img->xsize,img->ysize,img->zsize);
    fst->marker         = iftCreateImage(img->xsize,img->ysize,img->zsize);
    fst->root           = iftCreateImage(img->xsize,img->ysize,img->zsize);
    fst->pred           = iftCreateImage(img->xsize,img->ysize,img->zsize);
    fst->processed      = iftCreateBMap(img->n);
    fst->Q              = iftCreateGQueue(iftMax(IFT_QSIZE, iftMaximumValue(img) + 1), img->n, fst->pathval->val);
    fst->img            = img;
    fst->A              = A;
  } else {
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateImageForest");
  }

  // Default
  iftResetImageForest(fst);

  iftCopyVoxelSize(img, fst->pathval);
  iftCopyVoxelSize(img, fst->label);
  iftCopyVoxelSize(img, fst->marker);
  iftCopyVoxelSize(img, fst->root);
  iftCopyVoxelSize(img, fst->pred);

  return(fst);
}


void iftResetImageForest(iftImageForest *fst)
{
  int p;

  iftResetGQueue(fst->Q);

  if (fst->Q->C.removal_policy == MINVALUE){
    for (p=0; p < fst->img->n; p++) {
      fst->pathval->val[p] = IFT_INFINITY_INT;
      fst->label->val[p]   = 0;
      fst->marker->val[p]  = 0;
      fst->pred->val[p]    = IFT_NIL;
      fst->root->val[p]    = p;
    }
  } else { // MAXVALUE
    for (p=0; p < fst->img->n; p++) {
      fst->pathval->val[p] = IFT_INFINITY_INT_NEG;
      fst->label->val[p]   = 0;
      fst->marker->val[p]  = 0;
      fst->pred->val[p]    = IFT_NIL;
      fst->root->val[p]    = p;
    }
  }
}

void iftDestroyImageForest(iftImageForest **fst) {
    iftImageForest *aux = *fst;

    if (aux != NULL) {
        if (aux->processed != NULL)
            iftDestroyBMap(&aux->processed);
        if (aux->pathval != NULL)
            iftDestroyImage(&aux->pathval);
        if (aux->label != NULL)
            iftDestroyImage(&aux->label);
        if (aux->marker != NULL)
            iftDestroyImage(&aux->marker);
        if (aux->pred != NULL)
            iftDestroyImage(&aux->pred);
        if (aux->root != NULL)
            iftDestroyImage(&aux->root);
        if (aux->Q != NULL)
            iftDestroyGQueue(&aux->Q);

        iftFree(aux);    
        *fst = NULL;
    }
}


iftSet *iftCompRemoval(iftImageForest *fst, iftLabeledSet *seed)
{
  int        i, p, q, n=fst->img->n, V0;
  iftVoxel   u,v;
  iftAdjRel *A=fst->A;
  iftLabeledSet *S=NULL;
  iftSet   *Frontier=NULL;
  iftBMap  *inFrontier = iftCreateBMap(n);
  iftImage *pathval=fst->pathval,*pred=fst->pred;
  iftImage *img=fst->img, *label=fst->label;
  iftFIFO  *F1=iftCreateFIFO(n),*F2=iftCreateFIFO(n);

  if (fst->Q->C.removal_policy == MINVALUE)
    V0 = IFT_INFINITY_INT;
  else // MAXVALUE
    V0 = IFT_INFINITY_INT_NEG;

  /* Remove all connected components with at least one removal seed
     and insert the removal seed in another set to find the frontier
     voxels afterwards. */

  S = seed;
  while (S != NULL) {
    if ((S->label == IFT_NIL) && (pathval->val[S->elem] != V0)) { // Removal marker must be in regions not removed yet
      iftInsertFIFO(F2,S->elem);
      // Remove its largest connected component
      pathval->val[S->elem] = V0;  pred->val[S->elem] = IFT_NIL; // Remove element
      iftInsertFIFO(F1,S->elem);
      while (!iftEmptyFIFO(F1)) {
	p = iftRemoveFIFO(F1);
	u = iftGetVoxelCoord(img,p);
	for (i=1; i < A->n; i++) {
	  v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(img,v)){
	    q   = iftGetVoxelIndex(img,v);
	    if (pathval->val[q] != V0){
	      if (label->val[q]==label->val[p]){
		pathval->val[q] = V0;  pred->val[q] = IFT_NIL; // Remove element
		iftInsertFIFO(F1,q);
	      }
	    }
	  }
	}
      }
      // Reset FIFO without reinitializing the color array
      F1->first = F1->last = 0;
    }
    S = S->next;
  }
  iftDestroyFIFO(&F1);


  /* Find the frontier voxels of non-removed regions */

  while (!iftEmptyFIFO(F2)) {
    p = iftRemoveFIFO(F2);
    u = iftGetVoxelCoord(img,p);
    for (i=1; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){
	q   = iftGetVoxelIndex(img,v);
	if (pathval->val[q] == V0){
	  if (F2->color[q] == IFT_WHITE){
	    iftInsertFIFO(F2,q);
	  }
	}else{ // q is a frontier voxel
	  if (iftBMapValue(inFrontier,q)==0){
	    iftInsertSet(&Frontier,q);
	    iftBMapSet1(inFrontier,q);
	  }
	}
      }
    }
  }

  iftDestroyFIFO(&F2);
  iftDestroyBMap(&inFrontier);

  return(Frontier);
}

void iftCompRemovalWithoutFrontier(iftImageForest *fst, iftLabeledSet *seed)
{
  int      i, r, p, q, n=fst->img->n, V0;
  iftVoxel u,v;
  iftAdjRel *A=fst->A;
  iftLabeledSet *S;
  iftImage *pathval=fst->pathval,*root=fst->root,*pred=fst->pred;
  iftImage *img=fst->img, *label=fst->label;
  iftFIFO  *F=iftCreateFIFO(n);

  if (fst->Q->C.removal_policy == MINVALUE)
    V0 = IFT_INFINITY_INT;
  else // MAXVALUE
    V0 = IFT_INFINITY_INT_NEG;

  S = seed;
  while (S != NULL) {

    if (label->val[S->elem] > 0) { // object region
      if ((label->val[S->elem]!=S->label)&&
	  (S->label >= 0) &&
	  (pathval->val[S->elem] != V0)){ // remove only trees of
	                                  // distinct labels in regions
                                          // already conquered
	r = root->val[S->elem]; // get a removal root r
	pathval->val[r] = V0;  pred->val[r] = IFT_NIL;
	iftInsertFIFO(F,r);
	while (!iftEmptyFIFO(F)) {
	  p = iftRemoveFIFO(F);
	  u = iftGetVoxelCoord(img,p);
	  for (i=1; i < A->n; i++) {
	    v = iftGetAdjacentVoxel(A,u,i);
	    if (iftValidVoxel(img,v)){	      
	      q   = iftGetVoxelIndex(img,v);
	      if (pathval->val[q]!=V0){
		if (label->val[q]==label->val[p]){ 
		  pathval->val[q] = V0;  pred->val[q] = IFT_NIL;
		  iftInsertFIFO(F,q);
		}
	      }
	    }
	  }
	}
	F->first = F->last = 0;
      }
    }else{ // then label->val[S->elem] == 0) => background region
      if ((S->label > 0)&&
	  (pathval->val[S->elem] != V0)){ // remove only the marked tree	
	r = root->val[S->elem]; // get a removal root r
	pathval->val[r] = V0;  pred->val[r] = IFT_NIL; 
	iftInsertFIFO(F,r);
	while (!iftEmptyFIFO(F)) {
	  p = iftRemoveFIFO(F);
	  u = iftGetVoxelCoord(img,p);
	  for (i=1; i < A->n; i++) {
	    v = iftGetAdjacentVoxel(A,u,i);
	    if (iftValidVoxel(img,v)){
	      q   = iftGetVoxelIndex(img,v);
	      if (pathval->val[q]!=V0){
		if (pred->val[q]==p){ // visit tree
		  iftInsertFIFO(F,q);
		  pathval->val[q] = V0;  pred->val[q] = IFT_NIL;
		}
	      }
	    }
	  }
	}
	F->first = F->last = 0;
      }
    }
    S = S->next;
  }

  iftDestroyFIFO(&F);
}


iftImage *iftSwitchTreeLabels(iftImageForest *fst, iftLabeledSet *seed, iftImage *tmp_label)
{
  int        i, r, p, q, n=fst->img->n;
  iftVoxel   u,v;
  iftAdjRel *A=fst->A;
  iftLabeledSet *S;
  iftImage *root=fst->root;
  iftImage *img=fst->img, *label=fst->label;
  iftFIFO  *F=iftCreateFIFO(n);
  iftBMap  *SelComp=iftCreateBMap(n);

  if (seed==NULL){
      iftError("Seeds must be specified", "iftSwitchTreeLabels");
  }

  /* Copy labels of the selected components */

  if (tmp_label == NULL) 
    tmp_label = iftCreateImage(label->xsize, label->ysize, label->zsize);


  /* Switch label of each marked tree to either 0 or to the
     corresponding seed label in tmp_label */

  S = seed;
  while (S != NULL) {
    r = root->val[S->elem]; 
    if (iftBMapValue(SelComp,r)==0){
      iftInsertFIFO(F,r);
      iftBMapSet1(SelComp,r);
      while (!iftEmptyFIFO(F)) {
	p = iftRemoveFIFO(F);	
	u = iftGetVoxelCoord(img,p);
	if (tmp_label->val[p]==0){
	  tmp_label->val[p]=S->label;
	}else{
	  tmp_label->val[p]=0;
	}
	for (i=1; i < A->n; i++) {
	  v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(img,v)){
	    q   = iftGetVoxelIndex(img,v);
	    if ((label->val[q]==label->val[r])&&
		(iftBMapValue(SelComp,q)==0)){

	      iftBMapSet1(SelComp,q);
	      iftInsertFIFO(F,q);
	    }
	  }
	}
      }
      iftResetFIFO(F);
    }
    S = S->next;
  }
  
  iftDestroyFIFO(&F);
  iftDestroyBMap(&SelComp);

  return(tmp_label);
}

void iftRelabelBorderTreesAsBkg(iftImageForest *fst)
{
  int        p, n=fst->img->n;
  iftSet    *S=iftImageBorderSet(fst->img);
  iftImage  *label=fst->label;
  iftBMap   *RemLabel=iftCreateBMap(iftMaximumValue(label)+1);

  // Relabel the connected components that touch the image border as
  // background regions

  while (S != NULL) {
    p=iftRemoveSet(&S);
    if (iftBMapValue(RemLabel,label->val[p])==0){
      iftBMapSet1(RemLabel,label->val[p]);
    }
  }

  for (p=0; p < n; p++) 
    if (iftBMapValue(RemLabel,label->val[p])==1){
      label->val[p]=0;
    }

  iftDestroyBMap(&RemLabel);

}


iftSet *iftTreeRemoval(iftImageForest *fst, iftSet *trees_for_removal)
{
    int        i, p, q, r, n = fst->img->n, V0;
    iftVoxel   u, v;
    iftAdjRel *A = fst->A;
    iftSet   *Frontier = NULL;
    iftBMap  *inFrontier = iftCreateBMap(n);
    iftImage *pathval = fst->pathval, *pred = fst->pred, *root = fst->root;
    iftImage *img = fst->img;
    iftSet   *T1 = NULL, *T2 = NULL;

    if (fst->Q->C.removal_policy == MINVALUE)
        V0 = IFT_INFINITY_INT;
    else // MAXVALUE
        V0 = IFT_INFINITY_INT_NEG;

    /* Remove all marked trees and find the frontier voxels
       afterwards. */

    while (trees_for_removal != NULL)
      {
	p = trees_for_removal->elem;
    
	if (pathval->val[root->val[p]] != V0) //tree not marked yet
	  {
	    r = root->val[p];
	    pathval->val[r] = V0; // mark removed node
	    pred->val[r] = IFT_NIL;
	    iftInsertSet(&T1, r);
	    while (T1 != NULL)
	      {
		p = iftRemoveSet(&T1);
		iftInsertSet(&T2, p);
		
		u = iftGetVoxelCoord(img, p);
		for (i = 1; i < A->n; i++)
		  {
		    v = iftGetAdjacentVoxel(A, u, i);
		    if (iftValidVoxel(img, v))
		      {
			q   = iftGetVoxelIndex(img, v);
			if (pathval->val[q] != V0)
			  {
			    if (pred->val[q] == p)
			      {
				iftInsertSet(&T1, q);
				
				pathval->val[q] = V0; // mark removed node
				pred->val[q]    = IFT_NIL;
			      }
			  }
		      }
		  }
	      }
	  }
	trees_for_removal = trees_for_removal->next;
      }

  /* Find the frontier voxels of non-removed trees */

  while (T2 != NULL)
  {
    p = iftRemoveSet(&T2);
    u = iftGetVoxelCoord(img, p);
    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);
      if (iftValidVoxel(img, v))
      {
        q   = iftGetVoxelIndex(img, v);
        if (pathval->val[q] != V0)
        {
          if (iftBMapValue(inFrontier, q) == 0)
          {
            iftInsertSet(&Frontier, q);
            iftBMapSet1(inFrontier, q);
          }
        }
      }
    }
  }
  iftDestroyBMap(&inFrontier);

  return (Frontier);
}
iftSet *iftMarkerRemoval(iftImageForest *fst, iftSet *removal_markers)
{
    int        i, p, q, r, n = fst->img->n, V0;
    iftVoxel   u, v;
    iftAdjRel *A = fst->A;
    iftSet   *Frontier = NULL;
    iftBMap  *inFrontier = iftCreateBMap(n);
    iftImage *pathval = fst->pathval, *pred = fst->pred, *root = fst->root, *markers = fst->marker;
    iftImage *img = fst->img;
    iftSet   *T1 = NULL, *T2 = NULL;

    if (fst->Q->C.removal_policy == MINVALUE)
        V0 = IFT_INFINITY_INT;
    else // MAXVALUE
        V0 = IFT_INFINITY_INT_NEG;

    /* Remove all marked trees and find the frontier voxels
       afterwards. */

  while (removal_markers != NULL)
  {
    p = removal_markers->elem;

    if (pathval->val[root->val[p]] != V0) //tree not marked yet
    {
      r = root->val[p];
      pathval->val[r] = V0; // mark removed node
      pred->val[r] = IFT_NIL;
      iftInsertSet(&T1, r);
      while (T1 != NULL)
      {
        p = iftRemoveSet(&T1);
        iftInsertSet(&T2, p);

        u = iftGetVoxelCoord(img, p);
        for (i = 1; i < A->n; i++)
        {
          v = iftGetAdjacentVoxel(A, u, i);
          if (iftValidVoxel(img, v))
          {
            q   = iftGetVoxelIndex(img, v);
            if (pathval->val[q] != V0)
            {
              if (markers->val[q] == markers->val[r])
              {
                iftInsertSet(&T1, q);

                pathval->val[q] = V0; // mark removed node
                pred->val[q]    = IFT_NIL;
              }
            }
          }
        }
      }
    }
    removal_markers = removal_markers->next;
  }

  /* Find the frontier voxels of non-removed trees */

  while (T2 != NULL)
  {
    p = iftRemoveSet(&T2);
    u = iftGetVoxelCoord(img, p);
    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);
      if (iftValidVoxel(img, v))
      {
        q   = iftGetVoxelIndex(img, v);
        if (pathval->val[q] != V0)
        {
          if (iftBMapValue(inFrontier, q) == 0)
          {
            iftInsertSet(&Frontier, q);
            iftBMapSet1(inFrontier, q);
          }
        }
      }
    }
  }
  iftDestroyBMap(&inFrontier);

  return (Frontier);
}


void iftPropagateLabelMarkerAndRootToSubtree(iftImageForest *fst, iftAdjRel *A, iftImage *new_label, iftImage *new_marker, int r)
{
  iftImage *pred = fst->pred, *label = fst->label, *root = fst->root, *marker = fst->marker;
  int       p, q, i;
  iftVoxel  u, v;
  iftSet   *T = NULL;

  iftInsertSet(&T, r);

  while (T != NULL)
  {
    p              = iftRemoveSet(&T);
    u              = iftGetVoxelCoord(pred, p);
    new_label->val[p] = label->val[p]  = new_label->val[r];
    root->val[p]   = r;
    new_marker->val[p] = marker->val[p] = new_marker->val[r];

    for (i = 1; i < A->n; i++)
    {
      v   = iftGetAdjacentVoxel(A, u, i);
      if (iftValidVoxel(pred, v))
      {
        q = iftGetVoxelIndex(pred, v);
        if (pred->val[q] == p)
        {
          iftInsertSet(&T, q);
        }
      }
    }
  }
}

int iftPathLength(iftImage *pred, int p)
{
  if (pred->val[p] == IFT_NIL){
    return(1);
  } else {
    return(iftPathLength(pred, pred->val[p])+1);
  }
}

iftDataSet *iftImageForestToDataSet(iftImageForest* fst, iftAdjRel *A)
{
    int n_roots = 0;
//#pragma omp parallel for
    for (int i = 0; i < fst->pred->n; i++) {
        if (fst->pred->val[i] == IFT_NIL) {
            n_roots++;
        }
    }

    iftDataSet *Z = iftCreateDataSet(n_roots, A->n);

    int idx = 0;
    for (int p = 0; p < fst->pred->n; p++)
    {
        iftVoxel u = iftGetVoxelCoord(fst->pred, p);
        if (fst->pred->val[p] == IFT_NIL)
        {
            Z->sample[idx].id = fst->label->val[p];
            for (int j = 0; j < A->n; j++)
            {
                iftVoxel v = iftGetAdjacentVoxel(A, u, j);
                if (iftValidVoxel(fst->pred, v))
                {
                    int q = iftGetVoxelIndex(fst->pred, v);
                    Z->sample[idx].feat[j] = fst->img->val[q];
                }
            }
            idx++;
        }
    }

    return Z;
}
