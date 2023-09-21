#include "iftRepresentation.h"

#include "ift/core/dtypes/FIFO.h"
#include "ift/core/dtypes/LIFO.h"
#include "ift/imgproc/basic/Histogram.h"



/*----------------------------------------- Private methods -------------------------------------------------------------------*/

typedef struct ift_glen_node {
  int     voxel;  /* surface voxel */
  iftSet *adj;    /* its list to adjacent nodes */
  iftSet *closest_voxel; /* its closest voxels in the EDT */
} iftGLenNode;


typedef struct ift_GLen {
  iftGLenNode *node;        /* list of surface voxels */
  int          nnodes;      /* number of surface voxels */
  int         *index;       /* index to the closest surface voxel */
  iftAdjRel *A;             /* adjacency relation to compute geodesic lengths */
  float     *arcw;          /* arc weights for geodesic length computation */
  float     *cost;          /* optimum path cost of each node */
  float     *geolen;        /* geodesic length at each node */
  iftFHeap  *Q;             /* priority queue to compute geodesic lengths on a given surface */
  iftSet    *visited_node; /* list of nodes that have been inserted in Q */
  iftImage  *label;         /* label image used to create it */
  iftImage  *root;          /* root  image used to create it */
} iftGLen;



iftGLen *iftCreateGLenDataStruct(iftImage *root, iftImage *label)
{
  iftGLen *gl=(iftGLen *) iftAlloc(1,sizeof(iftGLen));
  iftVoxel u, v;
  int p, q, i, s, t; 
  
  gl->A        = iftSpheric(sqrtf(3.0));
  gl->arcw     = iftAllocFloatArray(gl->A->n);
  gl->root     = root;
  gl->label    = label;

  gl->arcw[0] = 0.0;  
  for (i=1; i < gl->A->n; i++) {
    float w0 = (gl->A->dx[i]*gl->A->dx[i])+
      (gl->A->dy[i]*gl->A->dy[i])+(gl->A->dz[i]*gl->A->dz[i]);
    gl->arcw[i]  = iftSmoothEuclideanDistance(w0);
  }

  gl->index  = iftAllocIntArray(root->n);
  gl->nnodes = 0;
  for (p=0; p < gl->root->n; p++){    
    if (gl->label->val[p]!=0){
      if (gl->root->val[p]==p){
	gl->nnodes++;
      }
    }
  }
  gl->cost    = iftAllocFloatArray(gl->nnodes);
  gl->geolen  = iftAllocFloatArray(gl->nnodes);
  gl->Q       = iftCreateFHeap(root->n,gl->cost);
  gl->visited_node = NULL;

  gl->node = (iftGLenNode *)iftAlloc(gl->nnodes,sizeof(iftGLenNode));
  s = 0;
  for (p=0; p < gl->root->n; p++)
    {
      if (gl->label->val[p]!=0){
	gl->index[p]= IFT_NIL;
	if (gl->root->val[p]==p){
	  gl->cost[s]= IFT_INFINITY_FLT;
	  gl->geolen[s]=0;
	  gl->node[s].voxel = p;
	  gl->index[p] = s;
	  gl->node[s].adj = NULL;
	  gl->node[s].closest_voxel = NULL;
	  s++;
	}
      }
    }
  
  for (p=0; p < gl->root->n; p++){
    if (gl->label->val[p]!=0){ 
      if (gl->root->val[p]==p){
	u = iftGetVoxelCoord(gl->root,p);
	s = gl->index[p];
	for (i=1; i < gl->A->n; i++) {
	  v = iftGetAdjacentVoxel(gl->A,u,i);
	  if (iftValidVoxel(gl->root,v)){
	    q = iftGetVoxelIndex(gl->root,v);
	    if ((gl->label->val[q]==gl->label->val[p])&&
		(gl->root->val[q]==q)){
	      t = gl->index[q];
	      iftInsertSet(&gl->node[s].adj,t);
	    }
	  }
	}
      } else {
	q = gl->root->val[p];
	s = gl->index[q];
	iftInsertSet(&gl->node[s].closest_voxel,p);
      }
    }
  }

  return(gl);
}

void iftDestroyGLenDataStruct(iftGLen **gl)
{
  iftGLen *aux = *gl;

  if (aux != NULL) {
    iftDestroyAdjRel(&aux->A);
    iftFree(aux->arcw);
    iftDestroyFHeap(&aux->Q);
    iftFree(aux->cost);
    iftFree(aux->geolen);
    iftFree(aux->index);
    for (int s=0; s < aux->nnodes; s++) {
      iftDestroySet(&aux->node[s].adj);
      iftDestroySet(&aux->node[s].closest_voxel);
    }
    iftFree(aux->node);
    if (aux->visited_node != NULL)
      iftDestroySet(&aux->visited_node);
    iftFree(aux);
    *gl = NULL;
  }
}

void iftResetGLenDataStruct(iftGLen *gl)
{

  while(gl->visited_node != NULL) { /* reset visited voxels and queue */
    int s           = iftRemoveSet(&gl->visited_node);
    gl->cost[s]     = IFT_INFINITY_FLT;
    gl->geolen[s]   = 0;
    gl->Q->color[s] = IFT_WHITE;
    gl->Q->pos[s]   = IFT_NIL;
    gl->Q->node[s]  = IFT_NIL;
    gl->Q->last     = IFT_NIL;
  }

}

float iftGeodesicLength(iftGLen *gl, int end)
{
  int         p,q,s,t;
  float       tmp;
  iftSet     *adj;
  iftVoxel    u, v;

  /* Geodesic Image Foresting Transform up to the end point */

  while(!iftEmptyFHeap(gl->Q)) {

    s = iftRemoveFHeap(gl->Q);
    p = gl->node[s].voxel;
    u = iftGetVoxelCoord(gl->root,p);

    for (adj=gl->node[s].adj; adj != NULL; adj = adj->next){

      t = adj->elem;

      /* must be a boundary voxel and p should be able to improve
	 its geodesic value */
      
      if (gl->Q->color[t] != IFT_BLACK){
	
	q = gl->node[t].voxel;
	v = iftGetVoxelCoord(gl->root,q);

	tmp = gl->cost[s] + iftSmoothEuclideanDistance((float) iftSquaredVoxelDistance(u, v));

	if (tmp < gl->cost[t]){

	  gl->cost[t]   = tmp;

	  if(gl->Q->color[t] == IFT_WHITE){
	    iftInsertFHeap(gl->Q, t);
	    iftInsertSet(&gl->visited_node,t);
	  }else{
	    iftGoUpFHeap(gl->Q, gl->Q->pos[t]);
	  }
	}
      }
    }
      
    if (p==end) {
      return(gl->cost[s]);
    }
  }

  return(IFT_INFINITY_FLT);
}

float iftGeodesicLength_Astar(iftGLen *gl, int end)
{
  int         p,q,s,t;
  float       tmp1, tmp2;
  iftSet     *adj;
  iftVoxel    u, v, e = iftGetVoxelCoord(gl->root,end);

  /* Geodesic Image Foresting Transform up to the end point */

  while(!iftEmptyFHeap(gl->Q)) {

    s = iftRemoveFHeap(gl->Q);
    p = gl->node[s].voxel;
    u = iftGetVoxelCoord(gl->root,p);
    gl->cost[s] = gl->geolen[s];

    for (adj=gl->node[s].adj; adj != NULL; adj = adj->next){

      t = adj->elem;

      /* must be a boundary voxel and p should be able to improve
	 its geodesic value */
      
      if (gl->Q->color[t] != IFT_BLACK){
	
	q = gl->node[t].voxel;
	v = iftGetVoxelCoord(gl->root,q);

	tmp1 = gl->cost[s] + iftSmoothEuclideanDistance((float) iftSquaredVoxelDistance(u, v));
	tmp2 = tmp1 + iftVoxelDistance(u,e);

	if (tmp2 < gl->cost[t]){

	  gl->cost[t]   = tmp2;
	  gl->geolen[t] = tmp1;

	  if(gl->Q->color[t] == IFT_WHITE){
	    iftInsertFHeap(gl->Q, t);
	    iftInsertSet(&gl->visited_node,t);
	  }else{
	    iftGoUpFHeap(gl->Q, gl->Q->pos[t]);
	  }
	}
      }
    }
      
    if (p==end) {
      
      iftSet *S = NULL;
      while(!iftEmptyFHeap(gl->Q)) {
	t = iftRemoveFHeap(gl->Q);
	gl->cost[t] = gl->geolen[t];	
	iftInsertSet(&S,t);
      }
      while (S != NULL) {
	t = iftRemoveSet(&S);
	iftInsertFHeap(gl->Q, t);
      }	

      return(gl->cost[s]);
    }
  }

  return(IFT_INFINITY_FLT);
}


float iftMaxGeodesicMeasure(iftGLen *gl, int *R, int nroots)
{
  int         i;
  float       auxgl,maxgl=0, d;
  iftVoxel    u, v;

  for (i=1; i < nroots; i++) {
    if (gl->cost[gl->index[R[i]]] == IFT_INFINITY_FLT){
      u = iftGetVoxelCoord(gl->root,R[0]);
      v = iftGetVoxelCoord(gl->root,R[i]);
      d = iftSquaredVoxelDistance(u,v);
      if (d > 3)
	//auxgl = iftGeodesicLength(gl, R[i]);	
	auxgl = iftGeodesicLength_Astar(gl, R[i]);	
      else
	auxgl = iftSmoothEuclideanDistance(d);
    }else{    
      auxgl = gl->cost[gl->index[R[i]]];
    }

    if (auxgl > maxgl){
      maxgl = auxgl;
    }
  }

  return(maxgl);
}

/* For each object, select as seed a point with maximum importance
   value. Compute an optimum-path forest, for the connectivity
   function defined by the last value of the path with origin in that
   seed set, and then compute the number of paths that pass to each
   point. This number is the size of the subtrees rooted at each
   point. Such importance measure is by itself a ms-skeleton. Another
   option is to compute the number of paths that pass to each point
   with terminus at a leaf node. This measure should be the geodesic
   surface area of Rennier. */

void iftAddMSCurveSkeletons(iftFImage *skel, iftImage *root, iftImage *label, iftAdjRel *A)
{
  iftFImage *pathval = iftCreateFImage(label->xsize,label->ysize,label->zsize);
  iftImage  *pred    = iftCreateImage(label->xsize,label->ysize,label->zsize);
  iftFHeap  *Q = iftCreateFHeap(label->n,pathval->val);

  iftSetRemovalPolicyFHeap(Q,IFT_MAXVALUE);

  /* Initialize trivial-path value map */

  float maxval_skel = iftFMaximumValue(skel);

  for (int p=0; p < label->n; p++) {
    if ((label->val[p] > 0)&&(skel->val[p]!=IFT_NIL)) { /* avoid SKIZ
							   points */ 
      pathval->val[p] = skel->val[p]-maxval_skel-1;
      pred->val[p]    = IFT_NIL;
      iftInsertFHeap(Q,p);
    }
  }

 /* Image Foresting Transform from the root in the curve skeleton to
    the object's boundary */

  while(!iftEmptyFHeap(Q)) {

    int p=iftRemoveFHeap(Q);
    iftVoxel u = iftGetVoxelCoord(label,p);
    if (pred->val[p]==IFT_NIL) {
      pathval->val[p] = skel->val[p];
    }

    for (int i=1; i < A->n; i++){
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(label,v)){
        int q = iftGetVoxelIndex(label,v);
	if ((label->val[q]>0)&&(skel->val[q]!=IFT_NIL)){ /* avoid SKIZ points */
	  if ((pathval->val[q] < skel->val[q]) && (Q->color[q] != IFT_BLACK)){
            pathval->val[q]  = skel->val[q];
	    pred->val[q]     = p;
	    iftGoUpFHeap(Q, Q->pos[q]);
	  }
	}
      }
    }
  }

  iftDestroyFHeap(&Q);
  iftDestroyFImage(&pathval);

  iftImage   *leaf     = iftLeafVoxels(pred,A);
  iftFImage  *ms_cskel = iftCreateFImage(label->xsize,label->ysize,label->zsize);
  int        nleafs=0;

  /* Compute the geodesic area for each curve skeleton voxel */

  for (int p=0; p < label->n; p++) {
    if ((leaf->val[p]>0)&&(root->val[p]==p)){
      nleafs++;
      int q = p;
      while (pred->val[q]!=IFT_NIL){
	ms_cskel->val[pred->val[q]]++;
	q = pred->val[q];
      }
    }
  }
  iftDestroyImage(&pred);
  printf("nleafs %d\n",nleafs);

  float maxval_ms_cskel = iftFMaximumValue(ms_cskel);

  /* Alex, this part of combining both skeletons needs more
     attention */

  for (int p=0; p < label->n; p++) {
    if (skel->val[p]>0){ /* avoid the SKIZ */
      if (ms_cskel->val[p] > 0.005*maxval_ms_cskel){
	ms_cskel->val[p] = (ms_cskel->val[p]/maxval_ms_cskel); 
	skel->val[p]     = ms_cskel->val[p]+maxval_skel;
      }
      else
	ms_cskel->val[p] = 0;
    }
  }
  /* iftFWriteVTKImage(ms_cskel,"multiscale-curve-skeletons.vtk"); */
  /* iftImage *cskel = iftFImageToImage(ms_cskel,4095); */
  /* iftWriteImage(cskel,"multiscale-curve-skeletons.scn"); */

  iftDestroyImage(&leaf);
  iftDestroyFImage(&ms_cskel);
}



int iftIdentifyRootsForSkelAndSkizComp(iftImage *root, iftImage *label, iftAdjRel *B, int p, int *R, char *skiz)
{
  int      i, j=0, k, q;
  iftVoxel u, v;

  R[0]=root->val[p]; j=1;
  u = iftGetVoxelCoord(root,p);
  *skiz = 0;

  for (i=1; i < B->n; i++) { // visit 6 neigbors
    v = iftGetAdjacentVoxel(B,u,i);
    if (iftValidVoxel(label,v)){
      q = iftGetVoxelIndex(label,v);
      if (label->val[p]<label->val[q]){ // skiz
	*skiz=1;
	break;
      }else{
	if (label->val[p]==label->val[q]){ // skeleton
	  if (root->val[p]!=root->val[q]){

	    // Compute union between root and root list

	    for (k=0; k < j; k++)
	      if (R[k]==root->val[q])
		break;
	    if (k==j){
	      R[j]=root->val[q]; j++;
	    }
	  }
	}
      }
    }
  }

  return(j); /* return number of identified roots */
}

float iftMaxAngleMeasure(int *R, int nroots, iftImage *root, int p)
{
  int         i;
  float       auxAng,maxAng=0,iprod;
  iftVector   v1,v2;
  float       m1,m2;
  iftVoxel    u1,u2;

  u1 = iftGetVoxelCoord(root,p);
  u2 = iftGetVoxelCoord(root,R[0]);

  v1.x = u2.x - u1.x;
  v1.y = u2.y - u1.y;
  v1.z = u2.z - u1.z;
  v1.t = v2.t = 0;
  
  m1   = iftVectorMagnitude(v1);

  for (i=1; i < nroots; i++) {

      u2     = iftGetVoxelCoord(root,R[i]);
      v2.x   = u2.x - u1.x;
      v2.y   = u2.y - u1.y;
      v2.z   = u2.z - u1.z;
      m2     = iftVectorMagnitude(v2);
      iprod  = iftVectorInnerProduct(v1, v2);
      auxAng = acos(iprod/(m1*m2)) * 180.0 / IFT_PI;

      if (auxAng > maxAng){
	maxAng = auxAng;
      }

  }


  return(maxAng);
}



/*------------------------------------------ Public Methods ---------------------------------------------------------------*/
int iftIntializeDistTransCostInterior(int val) {
    return (val > 0) ? IFT_INFINITY_INT : 0;
}

int iftIntializeDistTransCostExterior(int val) {
    return (val == 0) ? IFT_INFINITY_INT : 0;
}

int iftIntializeDistTransCostBoth(int val) {
    return IFT_INFINITY_INT;
}


iftIntializeDistTransCost iftGetIntializeDistTransCost(iftSide side) {
    switch (side) {
        case IFT_INTERIOR:
            return iftIntializeDistTransCostInterior;
        case IFT_EXTERIOR:
            return iftIntializeDistTransCostExterior;
        case IFT_BOTH:
        default:
            return iftIntializeDistTransCostBoth;
    }
}


iftFImage *iftGeodesicDistTrans(const iftSet *S, const iftImage *mask, const iftAdjRel *A)
{
  iftFImage *dist=NULL;
  iftFHeap *Q=NULL;
  int i,p,q;
  float tmp;
  iftVoxel u,v;
  const iftSet *Saux;

  dist  = iftCreateFImage(mask->xsize,mask->ysize,mask->zsize);
  Q     = iftCreateFHeap(mask->n,dist->val);


  for (p=0; p < mask->n; p++) {
    if (mask->val[p] != 0)
      dist->val[p]= IFT_INFINITY_FLT;
  }

  Saux = S;
  while (Saux != NULL) {
    p = Saux->elem;
    if (mask->val[p] == 0)
        iftError("Seed set must be inside mask", "iftGeodesicDistTrans");
    dist->val[p]=0;
    iftInsertFHeap(Q,p);
    Saux = Saux->next;
  }

  // Image Foresting Transform

  while(!iftEmptyFHeap(Q)) {

    p=iftRemoveFHeap(Q);
    //Gets the voxel.
    u = iftGetVoxelCoord(mask,p);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(mask,v)){
        q = iftGetVoxelIndex(mask,v);
        if ((dist->val[q] > dist->val[p]) && (mask->val[q] != 0)){
          tmp = dist->val[p] + iftSmoothEuclideanDistance((float) iftSquaredVoxelDistance(u, v));
          if (tmp < dist->val[q]){
            dist->val[q]  = tmp;
            if(Q->color[q] == IFT_WHITE){
              iftInsertFHeap(Q, q);
            }else{
              iftGoUpFHeap(Q, Q->pos[q]);
            }
          }
        }
      }
    }
  }

  iftDestroyFHeap(&Q);

  dist->dx = mask->dx;
  dist->dy = mask->dy;
  dist->dz = mask->dz;
  return(dist);
}

iftImage  *iftBorderDistTrans(const iftImage *label, iftAdjRel *A)
{
  iftImage *dist=NULL,*root=NULL;
  iftGQueue *Q=NULL;
  int i,p,q,tmp;
  iftVoxel u,v,r;
  iftSet *S=NULL;
  iftAdjRel *B;

  // Initialization

  if (iftIs3DImage(label)) {
    B = iftSpheric(1.0);
  }
  else { 
    B = iftCircular(1.0);
  }

  dist  = iftCreateImage(label->xsize,label->ysize,label->zsize);
  root  = iftCreateImage(label->xsize,label->ysize,label->zsize);
  Q     = iftCreateGQueue(IFT_QSIZE,label->n,dist->val);

  for (p=0; p < label->n; p++) {
      dist->val[p]= IFT_INFINITY_INT;
  }

  S     = iftObjectBorderSet(label, B); // seeds
  iftDestroyAdjRel(&B);
  while (S != NULL) {
    p = iftRemoveSet(&S);
    dist->val[p]=0;
    root->val[p]=p;
    iftInsertGQueue(&Q,p);
  }
  
  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {

    p=iftRemoveGQueue(Q);
    //Gets the voxel and its root.
    u = iftGetVoxelCoord(label,p);
    r = iftGetVoxelCoord(label,root->val[p]);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(label,v)){
	q = iftGetVoxelIndex(label,v);
	if (dist->val[q] > dist->val[p]){
	    tmp = iftSquaredVoxelDistance(v,r);
	    if (tmp < dist->val[q]){
	      if (dist->val[q] != IFT_INFINITY_INT)
		iftRemoveGQueueElem(Q, q);
	      dist->val[q]  = tmp;
	      root->val[q]  = root->val[p];
	      iftInsertGQueue(&Q, q);
	    }
	}
      }
    }
  }
  
  
  iftDestroyGQueue(&Q);
  iftDestroyImage(&root);

  iftCopyVoxelSize(label,dist);
  return(dist);
}

iftImage  *iftSetDistTrans(iftSet *S, int xsize, int ysize, int zsize)
{
  iftImage *dist=NULL,*root=NULL;
  iftGQueue *Q=NULL;
  int i,p,q,tmp;
  iftVoxel u,v,r;
  iftAdjRel *A = NULL;
  
  // Initialization

  if (zsize > 1){
    A = iftSpheric(1.74);
  }else{
    if (zsize == 0){
      iftWarning("zsize >= 1 is mandatory","iftSetDistTrans");
      zsize = 1;
    }
    A = iftCircular(1.5);
  }

  dist  = iftCreateImage(xsize,ysize,zsize);
  root  = iftCreateImage(xsize,ysize,zsize);
  Q     = iftCreateGQueue(IFT_QSIZE,dist->n,dist->val);

  for (p=0; p < dist->n; p++) {
      dist->val[p]= IFT_INFINITY_INT;
  }

  iftSet *Saux = S;
  while (Saux != NULL) {
    p = Saux->elem;
    iftVoxel u = iftGetVoxelCoord(dist,p);
    if (iftValidVoxel(dist,u)){
      dist->val[p]=0;
      root->val[p]=p;
      iftInsertGQueue(&Q,p);
    } else {
      iftError("Image domain must be bigger","iftSetDistTrans");
    }
    Saux = Saux->next;
  }
  
  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {

    p=iftRemoveGQueue(Q);
    //Gets the voxel and its root.
    u = iftGetVoxelCoord(dist,p);
    r = iftGetVoxelCoord(dist,root->val[p]);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(dist,v)){
	q = iftGetVoxelIndex(dist,v);
	if (dist->val[q] > dist->val[p]){
	    tmp = iftSquaredVoxelDistance(v,r);
	    if (tmp < dist->val[q]){
	      if (dist->val[q] != IFT_INFINITY_INT)
		iftRemoveGQueueElem(Q, q);
	      dist->val[q]  = tmp;
	      root->val[q]  = root->val[p];
	      iftInsertGQueue(&Q, q);
	    }
	}
      }
    }
  }

  
  iftDestroyAdjRel(&A);
  iftDestroyGQueue(&Q);
  iftDestroyImage(&root);

  return(dist);
}

/* This function will output infinity values outside the shell. */

iftImage  *iftShellDistTrans(iftImage *label, iftAdjRel *A, char side, float max_dist)
{
  iftImage *dist=NULL,*root=NULL;
  iftGQueue *Q=NULL;
  int i,p,q,tmp;
  iftVoxel u,v,r;
  iftSet *S=NULL;
  iftAdjRel *B=NULL;

  max_dist *= max_dist;

  // Initialization

  if (label->zsize==1) // 2D
    B = iftCircular(1.0);
  else // 3D
    B = iftSpheric(1.0);

  S = iftObjectBorderSet(label,B);
  iftDestroyAdjRel(&B);

  dist  = iftCreateImage(label->xsize,label->ysize,label->zsize);
  root  = iftCreateImage(label->xsize,label->ysize,label->zsize);
  Q  = iftCreateGQueue(IFT_QSIZE,label->n,dist->val);


  switch (side) {
  case IFT_INTERIOR:
    for (p=0; p < label->n; p++) {
      if (label->val[p] > 0)
	dist->val[p]= IFT_INFINITY_INT;
    }
    break;
  case IFT_EXTERIOR:
    for (p=0; p < label->n; p++) {
      if (label->val[p] == 0)
	dist->val[p]= IFT_INFINITY_INT;
    }
    break;
  case IFT_BOTH:
  default:
    for (p=0; p < label->n; p++) {
      dist->val[p]= IFT_INFINITY_INT;
    }
    break;
  }

  while (S != NULL) {
    p = iftRemoveSet(&S);
    dist->val[p]=0;
    root->val[p]=p;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    if (dist->val[p] <= max_dist){

      u = iftGetVoxelCoord(label,p);
      r = iftGetVoxelCoord(label,root->val[p]);

      for (i=1; i < A->n; i++){
	v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(label,v)){
	  q = iftGetVoxelIndex(label,v);
	  if (dist->val[q] > dist->val[p]){
	    tmp = (v.x-r.x)*(v.x-r.x) + (v.y-r.y)*(v.y-r.y) + (v.z-r.z)*(v.z-r.z);
	    if (tmp < dist->val[q]){
	      if (dist->val[q] != IFT_INFINITY_INT)
		iftRemoveGQueueElem(Q, q);
	      dist->val[q]  = tmp;
	      root->val[q]  = root->val[p];
	      iftInsertGQueue(&Q, q);
	    }
	  }
	}
      }
    }
  }

  iftDestroyGQueue(&Q);
  iftDestroyImage(&root);

  iftCopyVoxelSize(label,dist);
  return(dist);
}


void iftDistTransRootMap(iftImage *bin, iftAdjRel *A, char side, iftImage **dist, iftImage **root)
{
  iftGQueue *Q=NULL;
  int i,p,q,tmp;
  iftVoxel u,v,r;
  iftSet *S=NULL;
  iftAdjRel *B=NULL;

  // Initialization

  if (bin->zsize==1) // 2D
    B = iftCircular(1.0);
  else // 3D
    B = iftSpheric(1.0);

  S = iftObjectBorderSet(bin,B);
  iftDestroyAdjRel(&B);


  *dist  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  *root  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  Q  = iftCreateGQueue(IFT_QSIZE,bin->n,(*dist)->val);


  switch (side) {
  case IFT_INTERIOR:
    for (p=0; p < bin->n; p++) {
      if (bin->val[p] > 0)
  (*dist)->val[p]= IFT_INFINITY_INT;
    }
    break;
  case IFT_EXTERIOR:
    for (p=0; p < bin->n; p++) {
      if (bin->val[p] == 0)
  (*dist)->val[p]= IFT_INFINITY_INT;
    }
    break;
  case IFT_BOTH:
  default:
    for (p=0; p < bin->n; p++) {
      (*dist)->val[p]= IFT_INFINITY_INT;
    }
    break;
  }

  while (S != NULL) {
    p = iftRemoveSet(&S);
    (*dist)->val[p]=0;
    (*root)->val[p]=p;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    u = iftGetVoxelCoord(bin,p);
    r = iftGetVoxelCoord(bin,(*root)->val[p]);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(bin,v)){
  q = iftGetVoxelIndex(bin,v);
  if ((*dist)->val[q] > (*dist)->val[p]){
    tmp = (v.x-r.x)*(v.x-r.x) + (v.y-r.y)*(v.y-r.y) + (v.z-r.z)*(v.z-r.z);
    if (tmp < (*dist)->val[q]){
      if ((*dist)->val[q] != IFT_INFINITY_INT)
        iftRemoveGQueueElem(Q, q);
      (*dist)->val[q]  = tmp;
      (*root)->val[q]  = (*root)->val[p];
      iftInsertGQueue(&Q, q);
    }
  }
      }
    }
  }
  iftDestroyGQueue(&Q);

}


iftFImage  *iftSignedDistTrans(iftImage *bin, iftAdjRel *A)
{
  iftImage  *dist;
  iftFImage *sdist=NULL;
  int        p;

  dist  = iftEuclDistTrans(bin, A, IFT_BOTH, NULL, NULL, NULL);
  sdist = iftCreateFImage(dist->xsize,dist->ysize,dist->zsize);

  for (p=0; p < dist->n; p++)
    if (bin->val[p]==0){ /* exterior */
      sdist->val[p] = -sqrtf(dist->val[p]);
    }else{ /* interior */
      sdist->val[p] = +sqrtf(dist->val[p]);
    }

  iftDestroyImage(&dist);

  return(sdist);
}

iftFImage  *iftShellSignedDistTrans(iftImage *bin, iftAdjRel *A, float max_dist)
{
  iftImage  *dist;
  iftFImage *sdist=NULL;
  int        p;

  dist  = iftShellDistTrans(bin, A, IFT_BOTH, max_dist);
  sdist = iftCreateFImage(dist->xsize,dist->ysize,dist->zsize);

  for (p=0; p < dist->n; p++)
    if (bin->val[p]==0){ /* exterior */
      if (dist->val[p] != IFT_INFINITY_INT)
	sdist->val[p] = -sqrtf(dist->val[p]);
      else
	sdist->val[p] = 0;
    }else{ /* interior */
      if (dist->val[p] != IFT_INFINITY_INT)
	sdist->val[p] = +sqrtf(dist->val[p]);
      else
	sdist->val[p] = 0;
    }

  iftDestroyImage(&dist);

  iftFCopyVoxelSizeFromImage(bin, sdist);

  return(sdist);
}



/* In order to implement the curve skeleton, it is important to paint
   on the surface a line (thickness=3) with the optimum paths between
   each pair of roots of a given voxel and its 6-neighbors (maybe
   26?). Then if this line is closed, it divides the surface into two
   parts and you should find the geodesic area of this division by
   flood filling each part. The iftMaxGeodesicMeasure(gl,R,nroots)
   below contains optimizations to avoid computation of all pairs. You
   need to evaluate this idea in another function. One option is to
   use the function below to compute a surface skeleton first and then
   choose voxels above some threshold (they will be on the surface
   skeleton) to compute the curve skeleton.
*/

iftFImage *iftMSSkel(iftImage *bin) {
    iftFImage *skel;
    char      skiz = 0;
    int       s, p;
    int       nroots, R[7]; /* number of extended roots and extended root set */
    iftImage  *root, *label, *dist, *pred;
    iftAdjRel *B   = iftSpheric(1.0), *A = iftSpheric(sqrtf(3.0));
    iftGLen   *gl;
    iftSet    *closest_voxel;

    dist = iftEuclDistTrans(bin, A, IFT_INTERIOR, &root, &label, &pred);
    iftDestroyImage(&dist);

    iftDestroyImage(&pred);
    gl   = iftCreateGLenDataStruct(root, label);

    skel = iftCreateFImage(root->xsize, root->ysize, root->zsize);

    for (s = 0; s < gl->nnodes; s++) {
      closest_voxel = gl->node[s].closest_voxel; /* list of closest
						    voxels to the
						    surface voxel s */

      /* initialize data structure for a new geodesic computation */
      gl->cost[s] = 0.0;
      iftInsertSet(&gl->visited_node, s);
      iftInsertFHeap(gl->Q, s);

      while (closest_voxel != NULL) {

	p = closest_voxel->elem;

	nroots = iftIdentifyRootsForSkelAndSkizComp(root, label, B, p, R, &skiz);	  
	/* compute the maximum distance between roots */
	
	if (!skiz) {
	  
	  skel->val[p] = iftMaxGeodesicMeasure(gl, R, nroots);
	  
	} else {
	  skel->val[p] = IFT_NIL;
	}
       
	closest_voxel = closest_voxel->next;
      }

      iftResetGLenDataStruct(gl);
    }
    

    iftDestroyGLenDataStruct(&gl);


    /* Compute the multiscale curve skeletons */

    iftAddMSCurveSkeletons(skel,root,label,A);
    
    iftDestroyImage(&root);
    iftDestroyImage(&label);
    iftDestroyAdjRel(&A);
    iftDestroyAdjRel(&B);

    float skel_max_val = iftFMaximumValue(skel);

    for (p = 0; p < skel->n; p++) {
        if (skel->val[p] == IFT_NIL) {
            skel->val[p] = skel_max_val + 1;
        }
    }

    printf("%f\n", skel_max_val + 1);

    iftFCopyVoxelSizeFromImage(bin, skel);


    return (skel);
}


iftFImage *iftMSSkel2DDistMap(iftImage *label_img, iftAdjRel *A, iftSide side, iftImage **dist_out,
                              iftImage **relabel_img_out) {
    iftImage *root = NULL, *dist = NULL, *relabel_img = NULL, *pred = NULL;

    dist = iftEuclDistTrans(label_img, A, side, &root, &relabel_img, &pred);
    iftDestroyImage(&pred);

    int n = iftMaximumValue(relabel_img) + 1; // number of relabeled objects + bg
    iftFloatArray *perim = iftCreateFloatArray(n);

    iftFImage *glen = iftLabelContourPixelByGeoLen(label_img); // this is equivalent to Lpx from the Falcao's lectures
    iftFImage *skel = iftCreateFImage(root->xsize, root->ysize, root->zsize);

    // compute the perimeter for each label, which corresponds to the value of the last labeled contour pixel
    for (int p = 0; p < relabel_img->n; p++)
        if (relabel_img->val[p] != 0) {
            if (perim->val[relabel_img->val[p]] < glen->val[p])
                perim->val[relabel_img->val[p]] = glen->val[p];
        }

    iftAdjRel *B = iftCircular(1.0);

    iftVoxel u, v;
    u.z = v.z = 0;
    for (u.y = 0; u.y < relabel_img->ysize; u.y++)
        for (u.x = 0; u.x < relabel_img->xsize; u.x++) {
            int p = iftGetVoxelIndex(relabel_img, u);

            // if it is object pixel but not in the contour
            if ((relabel_img->val[p] != 0) && (root->val[p] != p)) {
                // Identify distinct roots for skeletonization and skiz
                float maxgl = 0.0;
                for (int i = 1; i < B->n; i++) {
                    v.x = u.x + B->dx[i];
                    v.y = u.y + B->dy[i];
                
                    if (iftValidVoxel(relabel_img, v)) {
                        int q = iftGetVoxelIndex(relabel_img, v);
                        
                        if (relabel_img->val[p] < relabel_img->val[q]) { // skiz
                            maxgl = IFT_NIL;
                            break;
                        } else {
                            if (relabel_img->val[p] == relabel_img->val[q]) { // skeleton
                                if (root->val[p] != root->val[q]) {
				                    float gl = (glen->val[root->val[q]] - glen->val[root->val[p]]);
                                    
                                    if (perim->val[relabel_img->val[p]] - gl < gl)
                                        gl = perim->val[relabel_img->val[p]] - gl;
                                    
                                    if (gl > maxgl) {
                                        maxgl = gl;
                                    }
                                }
                            }
                        }
                    }
                }
                skel->val[p] = maxgl;
            }
        }

    // compute the maximum geodesic length for each object (including the bg, i.e., external skeleton)
    float skel_max_val = 0.0f;
    for (int p = 0; p < relabel_img->n; p++) {
        if(skel_max_val < skel->val[p]) skel_max_val = skel->val[p];
    }

    for (int p = 0; p < skel->n; p++) {
        if (skel->val[p] == IFT_NIL)
            skel->val[p] = 100.0;
        else
            skel->val[p] = 100.0 * (skel->val[p] / skel_max_val);
    }

    iftDestroyFloatArray(&perim);
    iftDestroyFImage(&glen);
    iftDestroyImage(&root);
    iftDestroyAdjRel(&B);

    if (dist_out == NULL)
        iftDestroyImage(&dist);
    else *dist_out = dist;
    
    if (relabel_img_out == NULL)
        iftDestroyImage(&relabel_img);
    else *relabel_img_out = relabel_img;

    iftFCopyVoxelSizeFromImage(label_img, skel);

    return skel;
}

iftFImage *iftMSSkel2D(iftImage *label_img, iftAdjRel *A, iftSide side, iftImage **dist_out,
                       iftImage **relabel_img_out) {
    return iftMSSkel2DDistMap(label_img, A, side, dist_out, relabel_img_out);
}


iftImage *iftIntMSSkel2D(iftImage *bin, iftAdjRel *A, iftSide side) {
    iftImage  *skel;
    iftVoxel  u, v;
    int       p, q, i;
    int       maxgl, gl, *perim;
    iftImage  *root, *label, *cont, *dist, *pred;
    iftAdjRel *B = iftCircular(1.0);

    dist = iftEuclDistTrans(bin, A, side, &root, &label, &pred);
    iftDestroyImage(&dist);
    iftDestroyImage(&pred);

    perim = iftAllocIntArray(iftMaximumValue(label) + 1);
    cont  = iftLabelContPixel(bin);
    skel = iftCreateImage(root->xsize, root->ysize, root->zsize);

    for (p = 0; p < label->n; p++)
        if (label->val[p] != 0) {
            if (perim[label->val[p]] < cont->val[p])
                perim[label->val[p]] = cont->val[p];
        }

    u.z = v.z = 0;
    for (u.y = 0; u.y < label->ysize; u.y++)
        for (u.x = 0; u.x < label->xsize; u.x++) {

            p = iftGetVoxelIndex(label, u);

            if ((label->val[p] != 0) && (root->val[p] != p)) {

                // Identify distinct roots for skeletonization and skiz

                maxgl  = 0;
                for (i = 1; i < B->n; i++) {
                    v.x = u.x + B->dx[i];
                    v.y = u.y + B->dy[i];
                    if (iftValidVoxel(label, v)) {
                        q = iftGetVoxelIndex(label, v);
                        if (label->val[p] < label->val[q]) { // skiz
                            maxgl = IFT_NIL;
                            break;
                        } else {
                            if (label->val[p] == label->val[q]) { // skeleton
                                if (root->val[p] != root->val[q]) {
                                    gl     = (cont->val[root->val[q]] - cont->val[root->val[p]]);
                                    if (perim[label->val[p]] - gl < gl)
                                        gl = perim[label->val[p]] - gl;
                                    if (gl > maxgl) {
                                        maxgl = gl;
                                    }
                                }
                            }
                        }
                    }
                }
                skel->val[p] = maxgl;
            }
        }

    iftFree(perim);
    iftDestroyImage(&cont);
    iftDestroyImage(&root);
    iftDestroyImage(&label);
    iftDestroyAdjRel(&B);

    iftCopyVoxelSize(bin, skel);

    return (skel);
}


void  iftMultiLabelDistTransFromBorders(iftImage *label, iftAdjRel *A, char side,
											iftImage **dist, iftImage **root)
{
	iftSet *S=NULL;
	iftAdjRel *B=NULL;
	iftImage *_dist = NULL, *_root = NULL;

	if (label->zsize==1) // 2D
		B = iftCircular(1.0);
	else // 3D
		B = iftSpheric(1.0);

	S = iftObjectBorderSet(label, B);

	iftMultiLabelDistTransFromSet(S, label, A, side, &_dist, &_root);

	iftDestroySet(&S);
	iftDestroyAdjRel(&B);

	if(dist != NULL)
		*dist = _dist;
	else
		iftDestroyImage(&_dist);

	if(root != NULL)
		*root = _root;
	else
		iftDestroyImage(&_root);
}


void iftMultiLabelDistTransFromSet(iftSet *S, iftImage *label, iftAdjRel *A, char side,
				   iftImage **dist, iftImage **root){
  iftMultiLabelShellDistTransFromSet(S, label, A, side, IFT_INFINITY_DBL, dist, root);
}

void iftMultiLabelShellDistTransFromSet(iftSet *S, iftImage *label, iftAdjRel *A, char side, double max_dist,
					iftImage **dist, iftImage **root){
  iftGQueue *Q=NULL;
  int i,p,q,tmp;
  iftVoxel u,v,r;

  max_dist *= max_dist;

  // Initialization
  *dist  = iftCreateImage(label->xsize,label->ysize,label->zsize);
  *root  = iftCreateImage(label->xsize,label->ysize,label->zsize);
  Q  = iftCreateGQueue(IFT_QSIZE,label->n,(*dist)->val);

  switch (side) {
  case IFT_INTERIOR:
    for (p=0; p < label->n; p++) {
      if (label->val[p] > 0)
	(*dist)->val[p]= IFT_INFINITY_INT;
    }
    break;
  case IFT_EXTERIOR:
    for (p=0; p < label->n; p++) {
      if (label->val[p] == 0)
	(*dist)->val[p]= IFT_INFINITY_INT;
    }
    break;
  case IFT_BOTH:
  default:
    for (p=0; p < label->n; p++) {
      (*dist)->val[p]= IFT_INFINITY_INT;
    }
    break;
  }

  while (S != NULL) {
    p = S->elem;
    (*dist)->val[p]=0;
    (*root)->val[p]=p;
    iftInsertGQueue(&Q,p);
    S = S->next;
  }

  // Image Foresting Transform

  
  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);
    
    if ((*dist)->val[p] <= max_dist){
      
      u = iftGetVoxelCoord(label,p);
      r = iftGetVoxelCoord(label,(*root)->val[p]);
      
      for (i=1; i < A->n; i++){
	v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(label,v)){
	  q = iftGetVoxelIndex(label,v);
	  if ((*dist)->val[q] > (*dist)->val[p]){
	    tmp = (v.x-r.x)*(v.x-r.x) + (v.y-r.y)*(v.y-r.y) + (v.z-r.z)*(v.z-r.z);
	    if (tmp < (*dist)->val[q]){
	      if ((*dist)->val[q] != IFT_INFINITY_INT)
		iftRemoveGQueueElem(Q, q);
	      (*dist)->val[q]  = tmp;
	      (*root)->val[q]  = (*root)->val[p];
	      iftInsertGQueue(&Q, q);
	    }
	  }
	}
      }
    }
  }


  iftDestroyGQueue(&Q);
  
  iftCopyVoxelSize(label, *dist);
  iftCopyVoxelSize(label, *root);
}



iftImage *iftEuclDistTrans(const iftImage *label_img, const iftAdjRel *Ain, iftSide side, iftImage **root_out,
                           iftImage **edt_label_out, iftImage **pred_out) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = (iftIs3DImage(label_img)) ? iftSpheric(1.0) : iftCircular(1.0);
    else A = iftCopyAdjacency(Ain);

    // Initialization
    iftImage *dist      = iftCreateImageFromImage(label_img);
    iftImage *root      = iftCreateImageFromImage(label_img);
    iftImage *pred      = iftCreateImageFromImage(label_img);
    iftImage *edt_label = iftFindAndLabelObjectBorders(label_img, A);

    iftGQueue *Q = iftCreateGQueue(IFT_QSIZE, label_img->n, dist->val);
    iftIntializeDistTransCost initializationCostFunction = iftGetIntializeDistTransCost(side);

    for (int p = 0; p < label_img->n; p++) {
        dist->val[p] = initializationCostFunction(label_img->val[p]);

        if (edt_label->val[p] != 0) {
            dist->val[p] = 0;
            root->val[p] = p;
            pred->val[p] = IFT_NIL;
            iftInsertGQueue(&Q, p);
        }
    }


    // Image Foresting Transform
    while(!iftEmptyGQueue(Q)) {
        int p = iftRemoveGQueue(Q);

        iftVoxel u = iftGetVoxelCoord(label_img, p);
        iftVoxel r = iftGetVoxelCoord(label_img, root->val[p]);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            
            if (iftValidVoxel(label_img, v)) {
                int q = iftGetVoxelIndex(label_img, v);
    
                if (dist->val[q] > dist->val[p]) {
                    float tmp = (v.x-r.x)*(v.x-r.x) + (v.y-r.y)*(v.y-r.y) + (v.z-r.z)*(v.z-r.z);
                    
                    if (tmp < dist->val[q]) {
                        if (dist->val[q] != IFT_INFINITY_INT)
                            iftRemoveGQueueElem(Q, q);
                        dist->val[q]      = tmp;
                        root->val[q]      = root->val[p];
                        edt_label->val[q] = edt_label->val[p];
                        pred->val[q]      = p;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyGQueue(&Q);


    if (root_out == NULL)
        iftDestroyImage(&root);
    else *root_out = root;
    
    if (edt_label_out == NULL)
        iftDestroyImage(&edt_label);
    else *edt_label_out = edt_label;
    
    if (pred_out == NULL)
        iftDestroyImage(&pred);
    else *pred_out = pred;

    iftDestroyAdjRel(&A);

    return dist;
}


iftImage *iftGeodesicDistTransFromSet(const iftImage *mask, const iftSet *set, iftImage **root_out)
{
    iftAdjRel *A = iftIs3DImage(mask) ? iftSpheric(1.0) : iftCircular(1.0);

    // Initialization
    iftImage *dist      = iftCreateImageFromImage(mask);
    iftImage *root      = iftCreateImageFromImage(mask);

    iftGQueue *Q = iftCreateGQueue(IFT_QSIZE, mask->n, dist->val);

    for (int p = 0; p < mask->n; p++)
    {
        root->val[p] = p;
        if (mask->val[p])
            dist->val[p] = IFT_INFINITY_INT;
        else
            dist->val[p] = 0;
    }

    for (const iftSet *s = set; s; s = s->next)
    {
        int p = s->elem;
        if (mask->val[p]) {
            dist->val[p] = 0;
            iftInsertGQueue(&Q, p);
        }
    }


    // Image Foresting Transform
    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);

        iftVoxel u = iftGetVoxelCoord(mask, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(mask, v)) {
                int q = iftGetVoxelIndex(mask, v);

                if (mask->val[q] && dist->val[q] > (dist->val[p] + 1))
                {
                    if (Q->L.elem[q].color == IFT_GRAY)
                        iftRemoveGQueueElem(Q, q);
                    dist->val[q] = dist->val[p] + 1;
                    root->val[q] = root->val[p];
                    iftInsertGQueue(&Q, q);
                }
            }
        }
    }

    iftDestroyGQueue(&Q);

    if (root_out == NULL)
        iftDestroyImage(&root);
    else *root_out = root;

    iftDestroyAdjRel(&A);

    return dist;
}


iftFImage * iftMultiLabelSignedDistTransFromBorders(iftImage *label, iftAdjRel *A, char side, iftImage **root) {
    iftFImage *sdist = NULL;
    iftImage *dist = NULL;
    int p;

    iftMultiLabelDistTransFromBorders(label, A, side, &dist, root);
    sdist = iftCreateFImage(dist->xsize, dist->ysize, dist->zsize);

    for (p = 0; p < dist->n; p++)
        if (label->val[p] == 0) { /* exterior */
            sdist->val[p] = -sqrtf(dist->val[p]);
        } else { /* interior */
            sdist->val[p] = sqrtf(dist->val[p]);
        }

    iftDestroyImage(&dist);

    return sdist;
}

iftImage *iftRootPropagation(iftImage *bin, iftAdjRel *A, char side, float max_dist)
{
  iftImage *dist=NULL,*root=NULL;
  iftGQueue *Q=NULL;
  int i,p,q,tmp;
  iftVoxel u,v,r;
  iftSet *S=NULL;
  iftAdjRel *B=NULL;

  // Initialization

  dist     = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  root     = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  Q        = iftCreateGQueue(IFT_QSIZE,bin->n,dist->val);

  switch (side) {
  case IFT_INTERIOR:
    for (p=0; p < bin->n; p++) {
      if (bin->val[p] > 0)
	dist->val[p]= IFT_INFINITY_INT;
    }
    break;
  case IFT_EXTERIOR:
    for (p=0; p < bin->n; p++) {
      if (bin->val[p] == 0)
	dist->val[p]= IFT_INFINITY_INT;
    }
    break;
  case IFT_BOTH:
  default:
    for (p=0; p < bin->n; p++) {
      dist->val[p]= IFT_INFINITY_INT;
    }
    break;
  }

  // Seed computation

  if (bin->zsize == 1)
    B = iftCircular(1.0);
  else
    B = iftSpheric(1.0);

  S = iftObjectBorderSet(bin,B);
  iftDestroyAdjRel(&B);

  while (S != NULL) {
    p = iftRemoveSet(&S);
    dist->val[p]=0;
    root->val[p]=p;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    if (dist->val[p] > max_dist)
      break;

    u = iftGetVoxelCoord(bin,p);
    r = iftGetVoxelCoord(bin,root->val[p]);

    for (i=1; i < A->n; i++){
      v  = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(bin,v)){
	q = iftGetVoxelIndex(bin,v);
	if (dist->val[q] > dist->val[p]){
	  tmp = iftSquaredVoxelDistance(v,r);
	  if (tmp < dist->val[q]){
	    if (dist->val[q] != IFT_INFINITY_INT)
	      iftRemoveGQueueElem(Q, q);
	    dist->val[q]  = tmp;
	    root->val[q]  = root->val[p];
	    iftInsertGQueue(&Q, q);
	  }
	}
      }
    }
  }
  iftDestroyImage(&dist);
  iftDestroyGQueue(&Q);

  iftCopyVoxelSize(bin, root);

  return(root);
}


iftFImage *iftLabelContourPixelByGeoLen(iftImage *bin)
{
  int p,q,i;
  iftVoxel u,v;  
  
  if (iftIs3DImage(bin))
    iftError("Input image must be a 2D image", "iftLabelContourPixelByGeoLen");

  iftSet *B     = iftObjectBorderSet(bin, NULL);
  iftBMap *bndr = iftCreateBMap(bin->n);
  
  iftSet *Baux = B;
  while (Baux != NULL) {
    int p = Baux->elem;
    iftBMapSet1(bndr,p);
    Baux = Baux->next;
  }
    
  iftAdjRel *A = iftClockCircular(sqrtf(2.0));
  float *arcw = iftAllocFloatArray(A->n);
  
  arcw[0] = 0.0;  
  
  for (int i = 1; i < A->n; i++) {
    float w0 = iftPowerOfTwo(A->dx[i]) + iftPowerOfTwo(A->dy[i]) + iftPowerOfTwo(A->dz[i]);
    arcw[i]  = iftSmoothEuclideanDistance(w0);
  }
  
  iftAdjRel *L = iftLeftSide(A, 1.0);
  iftAdjRel *R = iftRightSide(A, 1.0);
  
  iftFImage *glen = iftCreateFImage(bin->xsize, bin->ysize, bin->zsize);
  iftImage *pred = iftCreateImageFromImage(bin);
  iftSetImage(pred, IFT_NIL);
  iftLIFO *LIFO = iftCreateLIFO(bin->n);
  
  Baux = B;
  while (Baux != NULL) {
    int r = Baux->elem; 
    int nvalidarcs=0;
    if (iftValidStartingPixel(bndr,pred, bin, A, L, R, r, &nvalidarcs)){
      iftInsertLIFO(LIFO,r);
      pred->val[r] = r;
      glen->val[r] = 0.0;
      
      while(!iftEmptyLIFO(LIFO)){
	
	p   = iftRemoveLIFO(LIFO);
	u   = iftGetVoxelCoord(bin,p);
	
	for (i=1; i < A->n; i++){
	    
	  v = iftGetAdjacentVoxel(A,u,i);
	  
	  if (iftValidVoxel(bin,v)){
	    q = iftGetVoxelIndex(bin,v);
	      
	    if ((q==r)&&(pred->val[p]!=r)){ /* A pixel p is the last
					      visited in the contour
					      when its neighbor is r,
					      but the precessor of p
					      is not r */	      
	      if (pred->val[pred->val[p]]==r){ /* if it is a silly
	      					  loop
	      					  r=(x,y),p=(x+1,y),
	      					  and pred[p]=(x,y+1),
	      					  then erase the glen
	      					  computation */

		while(pred->val[p]!=p){
		  glen->val[p] = 0.0;
		  p    = pred->val[p];
		}
		glen->val[p] = 0.0;
	      } 	 	  
	      iftResetLIFO(LIFO);
	      break;
	    }
	    
	    if (iftValidArc(bndr,pred,bin,A,L,R,u,i,q)){
	      pred->val[q] = p;
	      glen->val[q] = glen->val[p] + arcw[i];
	      iftInsertLIFO(LIFO,q);
	    }
	  }
	}
      }
    }
    Baux = Baux->next;
  }

  iftDestroyAdjRel(&A);
  iftDestroyAdjRel(&L);
  iftDestroyAdjRel(&R);
  iftDestroySet(&B);
  iftDestroyBMap(&bndr);
  iftDestroyImage(&pred);
  iftDestroyLIFO(&LIFO);
  iftFree(arcw);
  
  iftFCopyVoxelSizeFromImage(bin, glen);
  
  return(glen);
}

iftImage *iftLiftImage(iftImage *img)
{
  iftImage *bin = iftCreateImage(img->xsize, img->ysize, 255 + 1);
  int p,q;
  iftVoxel u;

  if (img->zsize != 1)
      iftError("It requires a 2D image", "iftLiftImage");

  for (p=0; p < img->n; p++) {
    u   = iftGetVoxelCoord(img,p);
    for (u.z = img->val[p]; u.z >= 0; u.z--){
      q   = iftGetVoxelIndex(bin,u);
      bin->val[q] = 1;
    }
  }

  return(bin);
}

iftImage *iftDropImage(iftImage *bin)
{
  iftImage *img = iftCreateImage(bin->xsize, bin->ysize, 1);
  int p,q;
  iftVoxel u,v;

  if ((bin->zsize == 1)||(iftMaximumValue(bin) != 1))
      iftError("It requires a 3D binary image", "iftDropImage");

  for (u.y=0; u.y < bin->ysize; u.y++)
    for (u.x=0; u.x < bin->xsize; u.x++){
      for (u.z=bin->zsize-1; u.z >= 0; u.z--){
	p = iftGetVoxelIndex(bin,u);
	if (bin->val[p]==1)
	  break;
      }
      v.x=u.x; v.y=u.y; v.z=0;
      q = iftGetVoxelIndex(img,v);
      img->val[q] = u.z;
    }

  return(img);
}

iftFImage *iftIntegralImage(iftImage *img)
{
  iftFImage *integ=NULL;
  iftVoxel u;
  int p;

  integ = iftCreateFImage(img->xsize,img->ysize,img->zsize);
  iftFCopyVoxelSizeFromImage(img, integ);

  if (img->zsize==1) {
    iftVoxel v[3];
    int q;
    u.z = v[0].z = v[1].z = v[2].z = 0;
    for (u.y=0; u.y < img->ysize; u.y++)
      for (u.x=0; u.x < img->xsize; u.x++){
	p = iftFGetVoxelIndex(integ,u);
	v[0].x = u.x - 1; v[0].y = u.y - 1;
	v[1].x = u.x;     v[1].y = u.y - 1;
	v[2].x = u.x - 1; v[2].y = u.y;
	integ->val[p] = img->val[p];
	if (iftFValidVoxel(integ,v[0])){
	  q = iftFGetVoxelIndex(integ,v[0]);
	  integ->val[p] -= integ->val[q];
	}
	if (iftFValidVoxel(integ,v[1])){
	  q = iftFGetVoxelIndex(integ,v[1]);
	  integ->val[p] += integ->val[q];
	}
	if (iftFValidVoxel(integ,v[2])){
	  q = iftFGetVoxelIndex(integ,v[2]);
	  integ->val[p] += integ->val[q];
	}
      }
  } else { // 3D
    iftVoxel v[7];
    int q;
    for (u.z=0; u.z < img->zsize; u.z++)
      for (u.y=0; u.y < img->ysize; u.y++)
	for (u.x=0; u.x < img->xsize; u.x++){
	  p = iftFGetVoxelIndex(integ,u);
	  v[0].x = u.x - 1; v[0].y = u.y - 1; v[0].z = u.z - 1;
	  v[1].x = u.x;     v[1].y = u.y - 1; v[1].z = u.z - 1;
	  v[2].x = u.x - 1; v[2].y = u.y;     v[2].z = u.z - 1;
	  v[3].x = u.x;     v[3].y = u.y;     v[3].z = u.z - 1;
	  v[4].x = u.x - 1; v[4].y = u.y - 1; v[4].z = u.z;
	  v[5].x = u.x;     v[5].y = u.y - 1; v[5].z = u.z;
	  v[6].x = u.x - 1; v[6].y = u.y;     v[6].z = u.z;

	  integ->val[p] = img->val[p];

	  if (iftFValidVoxel(integ,v[0])){
	    q = iftFGetVoxelIndex(integ,v[0]);
	    integ->val[p] += integ->val[q];
	  }
	  if (iftFValidVoxel(integ,v[1])){
	    q = iftFGetVoxelIndex(integ,v[1]);
	    integ->val[p] -= integ->val[q];
	  }
	  if (iftFValidVoxel(integ,v[2])){
	    q = iftFGetVoxelIndex(integ,v[2]);
	    integ->val[p] -= integ->val[q];
	  }
	  if (iftFValidVoxel(integ,v[3])){
	    q = iftFGetVoxelIndex(integ,v[3]);
	    integ->val[p] += integ->val[q];
	  }
	  if (iftFValidVoxel(integ,v[4])){
	    q = iftFGetVoxelIndex(integ,v[4]);
	    integ->val[p] -= integ->val[q];
	  }
	  if (iftFValidVoxel(integ,v[5])){
	    q = iftFGetVoxelIndex(integ,v[5]);
	    integ->val[p] += integ->val[q];
	  }
	  if (iftFValidVoxel(integ,v[6])){
	    q = iftFGetVoxelIndex(integ,v[6]);
	    integ->val[p] += integ->val[q];
	  }
	}
  }

  return(integ);
}

/* npts is the number of points that define the region.

   In 3D: npts = 8
   v[0] is the vertex (0,0,0) of the region in image coordinates,
   v[1] is the vertex (1,0,0) of the region in image coordinates,
   v[2] is the vertex (0,1,0) of the region in image coordinates,
   v[3] is the vertex (1,1,0) of the region in image coordinates,
   v[4] is the vertex (0,0,1) of the region in image coordinates,
   v[5] is the vertex (1,0,1) of the region in image coordinates,
   v[6] is the vertex (0,1,1) of the region in image coordinates,
   v[7] is the vertex (1,1,1) of the region in image coordinates,

   In 2D: npts = 4
   v[0] is the vertex (0,0,0) of the region in image coordinates,
   v[1] is the vertex (1,0,0) of the region in image coordinates,
   v[2] is the vertex (0,1,0) of the region in image coordinates,
   v[3] is the vertex (1,1,0) of the region in image coordinates,
*/

float iftGetIntegralValueInRegion(const iftFImage *integ, iftVoxel *v, int npts)
{
  int i, p[npts];

  if (integ->zsize==1){
    if (npts != 4)
        iftError("2D Image requires 4 points", "iftGetIntegralValueInRegion");

    for (i=0; i < 4; i++)
      if (iftFValidVoxel(integ,v[i]))
	p[i] = iftFGetVoxelIndex(integ,v[i]);
      else
	return(0);
    
    return(integ->val[p[3]]-integ->val[p[1]]-integ->val[p[2]]+integ->val[p[0]]);

  }else{

    if (npts != 8)
        iftError("3D Image requires 8 points", "iftGetIntegralValueInRegion");

    for (i=0; i < 8; i++)
      if (iftFValidVoxel(integ,v[i]))
	p[i] = iftFGetVoxelIndex(integ,v[i]);
      else
	return(0);
    
    return(-integ->val[p[0]] + integ->val[p[1]] + integ->val[p[2]] - integ->val[p[3]] + integ->val[p[4]] - integ->val[p[5]] - integ->val[p[6]] + integ->val[p[7]]);
  }

}

iftImage *iftBorderSize(iftImage *bin)
{
    int p;
    iftImage *border,*size;
    iftHist *hist;
    iftAdjRel *A;

    if (iftIs3DImage(bin))
        A=iftSpheric(1.0);
    else
        A=iftCircular(1.0);

    border = iftFindAndLabelObjectBorders(bin,A);
    size   = iftCreateImage(border->xsize,border->ysize,border->zsize);

    int max_val = iftMaximumValue(size);
    int nbins = max_val + 1;
    hist   = iftCalcGrayImageHist(border, NULL, nbins, max_val, 0);
    for (p=0; p < border->n; p++)
        if (border->val[p] != 0)
            size->val[p] = hist->val[border->val[p]];

    iftDestroyHist(&hist);
    iftDestroyImage(&border);

    return(size);
}

iftImage   *iftMarkGeometricCenters(iftImage *bin)
{
  int p,i,Lmax;
  iftImage *label,*geom;
  int *N;
  iftPoint *P;
  iftAdjRel *A;
  iftVoxel u;

  if (bin->zsize==1){
    A     = iftCircular(sqrtf(2.0));
  }else{
    A     = iftSpheric(sqrtf(3.0));
  }

  label = iftFastLabelComp(bin,A);
  geom  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  Lmax  = iftMaximumValue(label);
  P     = (iftPoint *)iftAlloc(Lmax+1, sizeof(iftPoint));
  if (P == NULL)
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftGeometricCenter");

  N     = iftAllocIntArray(Lmax+1);

  for (p=0; p < label->n; p++){
    P[label->val[p]].x += (float)(iftGetXCoord(label,p));
    P[label->val[p]].y += (float)(iftGetYCoord(label,p));
    P[label->val[p]].z += (float)(iftGetZCoord(label,p));
    N[label->val[p]]++;
  }

  for (i=1; i <= Lmax; i++){
    P[i].x /= N[i];
    P[i].y /= N[i];
    P[i].z /= N[i];
  }
  for (i=1; i <= Lmax; i++) {
    u.x = (int) P[i].x;
    u.y = (int) P[i].y;
    u.z = (int) P[i].z;
    p   = iftGetVoxelIndex(label,u);
    geom->val[p] = 255;
  }
  iftDestroyImage(&label);
  iftDestroyAdjRel(&A);
  iftFree(P);
  iftFree(N);

  return(geom);
}

iftImage *iftComponentSizes(iftImage *bin, iftAdjRel *A)
{
  iftImage  *label;
  iftAdjRel *B;
  int       *size;

  if (iftIs3DImage(bin))
    B=iftSpheric(sqrtf(3.0));
  else
    B=iftCircular(sqrtf(2.0));

  label = iftFastLabelComp(bin,B);
  iftDestroyAdjRel(&B);

  size = iftAllocIntArray(iftMaximumValue(label)+1);

  for (int p=0; p < label->n; p++)
    size[label->val[p]]++;

  for (int p=0; p < label->n; p++)
    label->val[p] = size[label->val[p]];

  iftFree(size);

  return(label);
}

void iftObjectAreaFromPixel(const iftImage* comp, int index, iftImage *area_img)
{
    iftFIFO *F = iftCreateFIFO(comp->n);
    iftAdjRel *A = iftCircular(1.0);

    int area = 0;
    iftInsertFIFO(F, index);
    while (!iftEmptyFIFO(F))
    {
        int p = iftRemoveFIFO(F);
        area++;
        iftVoxel u = iftGetVoxelCoord(comp, p);

        for(int i = 0; i <A->n; i++){
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if(iftValidVoxel(comp, v)){
                int q = iftGetVoxelIndex(comp, v);
                if(comp->val[q] == comp->val[p] && F->color[q] == IFT_WHITE)
                    iftInsertFIFO(F, q);
            }
        }
    }

    iftResetFIFO(F);
    iftInsertFIFO(F, index);

    while (!iftEmptyFIFO(F)) {
        int p = iftRemoveFIFO(F);
        area_img->val[p] = area;
        iftVoxel u = iftGetVoxelCoord(comp, p);

        for(int i = 0; i <A->n; i++){
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if(iftValidVoxel(comp, v)){
                int q = iftGetVoxelIndex(comp, v);
                if(comp->val[q] == comp->val[p] && F->color[q] == IFT_WHITE)
                    iftInsertFIFO(F, q);
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftDestroyFIFO(&F);
}


iftImage *iftSurfaceSkeleton(iftFImage *skel, float thres, int number_of_components)
{
  iftImage  *surf_skel, *aux;
  int        p;
  iftAdjRel *A=iftSpheric(sqrtf(3.0));

  aux = iftCreateImage(skel->xsize,skel->ysize,skel->zsize);
  for (p=0; p < skel->n; p++)
    if (skel->val[p]>=thres)
      aux->val[p]=1;
  surf_skel    = iftSelectKLargestComp(aux,A,number_of_components);

  iftDestroyImage(&aux);
  iftDestroyAdjRel(&A);

  iftFCopyVoxelSizeToImage(skel, surf_skel);

  return(surf_skel);
}

iftImage *iftMedialAxisTrans2D(iftImage *bin, float scale_thres, iftSide side)
{
  iftAdjRel *A      = iftCircular(sqrtf(2.0)); 
  iftImage  *dist   = NULL;
  iftFImage *msskel = iftMSSkel2DDistMap(bin, A, side, &dist, NULL);
  iftImage  *skel   = NULL, *mat = NULL;

  if ((scale_thres <= 0.0)&&(scale_thres > 100.0))
      iftError("Scale threshold must be within (0,100]", "iftMedialAxisTrans2D");

  skel = iftFThreshold(msskel,scale_thres,100.0,1);

  mat = iftMult(skel,dist);

  iftDestroyAdjRel(&A);
  iftDestroyImage(&dist);
  iftDestroyImage(&skel);
  iftDestroyFImage(&msskel);

  return(mat);
}

/* We need to substitute this naive code by the reverse EDT, which
   computes the EDT from the medial axis and stops propagation when
   the voxel removed from queue is at the object's boundary (i.e., at
   the distance to the medial axis equal to the value of the medial
   axis stored in its root voxel) */


iftImage *iftShapeReconstruction(iftImage *medial_axis, int value)
{
  iftImage  *bin = iftCreateImage(medial_axis->xsize,medial_axis->ysize,medial_axis->zsize);

  if (iftIs3DImage(medial_axis)){
#pragma omp parallel for shared(medial_axis,bin)
    for (int p=0; p < medial_axis->n; p++)
      if (medial_axis->val[p] > 0){
	bin->val[p]=1;
	iftVoxel u   = iftGetVoxelCoord(medial_axis,p);
	iftAdjRel *A = iftSpheric(sqrtf(medial_axis->val[p]));
	for (int i=0; i < A->n; i++) {
	  iftVoxel v = iftGetAdjacentVoxel(A,u,i); // By definition,
						      // all voxels v are valid
	  int q = iftGetVoxelIndex(medial_axis,v);
	  bin->val[q]=value;
	}
	iftDestroyAdjRel(&A);
      }
  }else{ // 2D image
#pragma omp parallel for shared(medial_axis,bin)
    for (int p=0; p < medial_axis->n; p++)
      if (medial_axis->val[p] > 0){
	bin->val[p]=1;
	iftVoxel u   = iftGetVoxelCoord(medial_axis,p);
	iftAdjRel *A = iftCircular(sqrtf(medial_axis->val[p]));
	for (int i=0; i < A->n; i++) {
	  iftVoxel v = iftGetAdjacentVoxel(A,u,i); // By definition,
						      // all voxels v are valid
	  int q = iftGetVoxelIndex(medial_axis,v);
	  bin->val[q]=value;
	}
	iftDestroyAdjRel(&A);
      }
  }

  return(bin);
}


/*
iftImage *iftShapeReconstruction(iftImage *medial_axis, int value)
{
    iftImage *bin = iftCreateImage(medial_axis->xsize,medial_axis->ysize,medial_axis->zsize);
    iftImage *root = iftCreateImage(medial_axis->xsize, medial_axis->ysize, medial_axis->zsize);
    iftImage *pathval = iftCreateImage(medial_axis->xsize, medial_axis->ysize, medial_axis->zsize);

    const int int2float = 1000;
    iftGQueue *Q = iftCreateGQueue((int) (sqrt(iftMaximumValue(medial_axis)) * int2float) + 1,
                                    medial_axis->n, pathval->val);

    for (int i = 0; i < medial_axis->n; i++) {
        pathval->val[i] = IFT_INFINITY_INT;
        if (medial_axis->val[i] > 0) {
            bin->val[i] = value;
            root->val[i] = i;
            pathval->val[i] = 0;
            iftInsertGQueue(&Q, i);
        }
    }

    iftAdjRel *A = iftCircular(sqrtf(2.0f));

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        int r = root->val[p];
        iftVoxel v = iftGetVoxelCoord(pathval, p);
        iftVoxel t = iftGetVoxelCoord(pathval, r);

        int max_dist = (int) (sqrt(medial_axis->val[r]) * int2float);

        for (int i = 1; i < A->n; i++) {
            iftVoxel u = iftGetAdjacentVoxel(A, v, i);
            if (!iftValidVoxel(pathval, u))
                continue;
            int q = iftGetVoxelIndex(pathval, u);
            if (Q->L.elem[q].color != IFT_BLACK) {
                int dist = (int) (iftVoxelDistance(u, t) * int2float);
                if (dist < max_dist && dist < pathval->val[q])
                {
                    if (Q->L.elem[q].color == IFT_GRAY)
                        iftRemoveGQueueElem(Q, q);
                    bin->val[q] = value;
                    root->val[q] = r;
                    pathval->val[q] = dist;
                    iftInsertGQueue(&Q, q);
                }
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftDestroyGQueue(&Q);
    iftDestroyImage(&pathval);
    iftDestroyImage(&root);

    return bin;
}
*/

iftSet *iftTerminalPointSet2D(iftImage *skel)
{
  iftAdjRel *A=iftCircular(sqrtf(2.0));
  iftSet    *S=NULL;

  for (int p=0; p < skel->n; p++){
    if (skel->val[p]!=0){
      iftVoxel u = iftGetVoxelCoord(skel,p);
      int n=0;
      for (int i=1; i < A->n; i++) {
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(skel,v)){
	  int q = iftGetVoxelIndex(skel,v);
	  if (skel->val[q]!=0){
	    n++;
	  }
	} else {n=3; break; /* invalid terminal point of the
			       background skeleton */ }
      }
      if (n==1)
	iftInsertSet(&S,p);
    }
  }
  
  iftDestroyAdjRel(&A);
  
  return(S);
}

iftImage *iftTerminalPoints2D(iftImage *skel)
{
  int        i, p, q, n;		
  iftVoxel   u, v;
  iftAdjRel *A=iftCircular(sqrtf(2.0));
  iftImage  *termpt = iftCreateImage(skel->xsize,skel->ysize,skel->zsize);
  
  if (iftIs3DImage(skel))
      iftError("It works for 2D images only", "iftTerminalPoints2D");

  for (p=0; p < skel->n; p++) {
    if (skel->val[p]){
      u = iftGetVoxelCoord(skel,p);
      n = 0;
      for (i=1; i < A->n; i++) {
	v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(skel,v)){
	  q = iftGetVoxelIndex(skel,v);
	  if (skel->val[q]){
	    n++;
	  }
	}else {n=3; break; /* invalid terminal point of the
			       background skeleton */ }      	
      }
      if (n==1) 
	termpt->val[p]=1;      
    }
  }
  iftDestroyAdjRel(&A);
    
  return(termpt);
}

iftImage *iftBranchPoints2D(iftImage *skel)
{
  iftImage  *bpts = iftCreateImageFromImage(skel);
  iftAdjRel *A    = iftCircular(1.5);
  
  for (int p=0; p < skel->n; p++) {
    if (skel->val[p]){
      iftVoxel u    = iftGetVoxelCoord(skel,p);
      int   npts    = 0; /* number of neighbors in the skeleton */
      for (int i=1; i < A->n; i++) { 
	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(skel,v)){
	  int q = iftGetVoxelIndex(skel,v);
	  if (skel->val[q]){
	      npts++;
	  }
	}
      }
      if (npts > 2) {
	bpts->val[p] = 255;
      }
    }
  }
  iftDestroyAdjRel(&A);
  
  return(bpts);
}


iftImage *iftNumberOfDescendantsInTarget(iftImage *pred, iftImage *target_bin)
{
  iftImage *ndesc = iftCreateImage(pred->xsize,pred->ysize,pred->zsize);
  int       p,q;

  iftVerifyImageDomains(pred,target_bin,"iftNumberOfDescendantsInTarget");

  for (p=0; p < target_bin->n; p++) {
    if (target_bin->val[p]>0){
      q = pred->val[p]; 
      while (q != IFT_NIL){
	ndesc->val[q]++;
	q = pred->val[q];
      }
    }
  }

  iftCopyVoxelSize(pred,ndesc);

  return(ndesc);
}

iftImage *iftNumberOfAscendants(iftImage *pred)
{
  iftImage *nasc = iftCreateImage(pred->xsize,pred->ysize,pred->zsize);
  int       p,q,n;

  for (p=0; p < pred->n; p++) {
    q = pred->val[p]; n=0;
    while (q != IFT_NIL){
	n++;
	q = pred->val[q];
      }
    nasc->val[p]=n;
  }

  iftCopyVoxelSize(pred,nasc);

  return(nasc);
}

iftImage *iftNumberOfDescendants(iftImage *pred)
{
  iftImage *ndesc = iftCreateImage(pred->xsize,pred->ysize,pred->zsize);

  for (int p=0; p < pred->n; p++) {    
    int q = pred->val[p];
    while (q != IFT_NIL){
      ndesc->val[q]++;
      q = pred->val[q];
    }
  }

  iftCopyVoxelSize(pred,ndesc);

  return(ndesc);
}

float iftWeightedIntegralValue(const iftFImage *int_img, iftVoxel center, int width, int height, float weight)
{
    iftVoxel pt[4];

    pt[0].x = center.x - width/2;
    pt[0].y = center.y - height/2;
    pt[0].z = 0;

    pt[1].x = center.x + width/2;
    pt[1].y = center.y - height/2;
    pt[1].z = 0;

    pt[2].x = center.x - width/2;
    pt[2].y = center.y + height/2;
    pt[2].z = 0;

    pt[3].x = center.x + width/2;
    pt[3].y = center.y + height/2;
    pt[3].z = 0;

    float int_value = weight * iftGetIntegralValueInRegion(int_img, pt, 4);

    return int_value;
}


iftMImage *iftExtractRegionOfInterest(iftMImage *morig, iftMImage *mlabel, float radius)
{
    iftMImage *mroi = iftCreateMImage(2 * radius + 1,2 * radius + 1,1,morig->m);
    iftAdjRel *A=iftCircular(sqrtf(2.0)); /* 8-neighbors */

    for (int b=0; b < mlabel->m; b++) {
        iftImage *label  =  iftMImageToImage(mlabel, 255, b); /* get label image b */
        int       maxval =  iftMaximumValue(label); /* select is brightest objects */
        for (int p=0; p < label->n; p++){
            if (label->val[p]!=maxval)
                label->val[p]=0;
        }
        iftImage *mask = iftSelectLargestComp(label,A);
        iftVoxel c = iftGeometricCenterVoxel(mask);
        iftDestroyImage(&label);
        iftDestroyImage(&mask);

        iftVoxel u; u.z = 0;
        int q = 0;
        for (u.y = c.y - radius; u.y <= c.y+radius; u.y++)
        {
            for (u.x = c.x - radius; u.x <= c.x+radius; u.x++)
            {
                if (iftMValidVoxel(morig,u))
                {
                    float dist = sqrtf((u.y-c.y)*(u.y-c.y)+(u.x-c.x)*(u.x-c.x));

                    if (dist <= radius) {
                        int p = iftMGetVoxelIndex(morig,u);
                        mroi->val[q][b] = morig->val[p][b];
                    } else {
                        mroi->val[q][b] = 0;
                    }
                    q++;
                }
            }
        }
    }

    iftDestroyAdjRel(&A);

    return(mroi);
}


double iftPerpendicularDist(const iftVoxel *cur, const iftVoxel *first, const iftVoxel *last)
{
    if (first->x == last->x)
        return abs(cur->x - first->x);

    double slope = (last->y - first->y) / (last->x - first->x);
    double intercept = first->y - (slope * first->x);
    return fabs(slope * cur->x - cur->y + intercept) / sqrt(slope * slope + 1);
}


double iftPointToLineSqrDist(const iftVoxel *cur, const iftVoxel *first, const iftVoxel *last)
{
    double len = iftSquaredVoxelDistance((*first), (*last)); // First, we need the length of the line segment.
    if (len == 0) // if it's 0, the line is actually just a point.
        return iftSquaredVoxelDistance((*first), (*cur));

    // t is very important. t is a number that essentially compares the individual coordinates
    // distances between the point and each point on the line.

    double t = ((cur->x - first->x) * (last->x - first->x) +
            (cur->y - first->y) * (last->y - first->y)) / len;

    if (t < 0) { 	// if t is less than 0, the point is behind i, and closest to i.
        return iftSquaredVoxelDistance((*first), (*cur));
    } else if (t > 1) { // if greater than 1, it's closest to j.
        return iftSquaredVoxelDistance((*last), (*cur));
    }
    // this figure represents the point on the line that p is closest to.
    double tmp_x = first->x + t * (last->x - first->x);
    double tmp_y = first->y + t * (last->y - first->y);
    return (cur->x - tmp_x) * (cur->x - tmp_x) + (cur->y - tmp_y) * (cur->y - tmp_y);
}


// Ramer Douglas Peucker
iftVoxel *iftRDP(const iftVoxel *array, long len, double epsilon, long *out_len)
{
    const iftVoxel *first = &array[0], *last = &array[len - 1];
    if (len <= 2)
    {
        *out_len = 2;
        iftVoxel *res = iftAlloc(2, sizeof *res);
        res[0].x = first->x;
        res[0].y = first->y;
        res[0].z = first->z;
        res[1].x = last->x;
        res[1].y = last->y;
        res[1].z = last->z;
        return res;
    }

    int idx = -1;
    double dist = 0;
    for (int i = 1; i < len - 1; i++)
    {
//        double cur_dist = iftPerpendicularDist(&array[i], first, last);
        double cur_dist = iftPointToLineSqrDist(&array[i], first, last);
        if (cur_dist > dist) {
            dist = cur_dist;
            idx = i;
        }
    }

    if (dist > epsilon)
    {
        long left_len = 0, right_len = 0;
        iftVoxel *left = iftRDP(array, idx + 1, epsilon, &left_len);
        iftVoxel *right = iftRDP(array + idx, (len - idx), epsilon, &right_len);
        *out_len = left_len + right_len - 1;
        iftVoxel *res = iftAlloc(*out_len, sizeof *res);
        int i = 0;
        for (; i < left_len; i++) {
            res[i].x = left[i].x;
            res[i].y = left[i].y;
            res[i].z = left[i].z;
        }
        // left end and right start are the same
        for (int j = 1; j < right_len; j++, i++) {
            res[i].x = right[j].x;
            res[i].y = right[j].y;
            res[i].z = right[j].z;
        }
        iftFree(left);
        iftFree(right);
        return res;
    } else {
        *out_len = 2;
        iftVoxel *res = iftAlloc(2, sizeof *res);
        res[0].x = first->x;
        res[0].y = first->y;
        res[0].z = first->z;
        res[1].x = last->x;
        res[1].y = last->y;
        res[1].z = last->z;
        return res;
    }
}


iftVoxelArray *iftContourToArray(const iftImage *contour)
{
    // this is necessary to keep the array on clockwise direction
    iftAdjRel *A = iftClockCircular(sqrtf(2.0f));

    iftImage *pred = iftCreateImage(contour->xsize, contour->ysize, contour->zsize);
    iftSetImage(pred, IFT_NIL);

    iftLIFO *Q = iftCreateLIFO(contour->n);

    int src = 1;
    for (int i = 0; i < contour->n; i++) {
        if (contour->val[i]) {
            src = i;
            break;
        }
    }

    iftVoxelArray *array = iftCreateVoxelArray(contour->n);
    int size = 0;
    iftInsertLIFO(Q, src);
    while (!iftEmptyLIFO(Q))
    {
        int p = iftRemoveLIFO(Q);
        iftVoxel u = iftGetVoxelCoord(contour, p);

        iftCopyVoxel(&u, &array->val[size]);
        size += 1;

        int j = 1; // because of the selection of the starting point (no predecessor)
        if (pred->val[p] != IFT_NIL) {
            iftVoxel t = iftGetVoxelCoord(pred, pred->val[p]);
            j = ift6NeighborsClockFromLeft(&t, &u);
        }

        for (int i = 0; i < A->n; i++, j = (j + 1) % A->n)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, j);
            if (!iftValidVoxel(contour, v))
                continue;
            int q = iftGetVoxelIndex(contour, v);
            if (contour->val[q] && Q->color[q] != IFT_BLACK) {
                pred->val[q] = p;
                if (Q->color[q] == IFT_WHITE) {
                    iftInsertLIFO(Q, q);
                    break;
                }
            }
        }
    }

    iftDestroyLIFO(&Q);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&pred);

    array->val = iftRealloc(array->val, size *(sizeof *array));
    array->n = size;
    return array;
}


void iftSplitByMaxDist(const iftVoxelArray *array, iftVoxelArray **out1, iftVoxelArray **out2)
{
    int max_i = -1, max_j = -1;
    double max_dist = 0;

    for (int i = 0; i < array->n - 1; i++) {
        for (int j = i + 1; j < array->n; j++) {
            double dist = iftSquaredVoxelDistance(array->val[i], array->val[j]);
            if (dist > max_dist) {
                max_i = i;
                max_j = j;
                max_dist = dist;
            }
        }
    }

    *out1 = iftCreateVoxelArray(max_j - max_i + 1);
    *out2 = iftCreateVoxelArray(max_i + array->n - max_j + 1);

    for (int i = 0, j = max_i; i < (*out1)->n; i++, j++)
        iftCopyVoxel(&array->val[j], &(*out1)->val[i]);

    int i = 0;
    for (int j = max_j; j < array->n; i++, j++)
        iftCopyVoxel(&array->val[j], &(*out2)->val[i]);

    for (int j = 0; j < max_i + 1; i++, j++)
        iftCopyVoxel(&array->val[j], &(*out2)->val[i]);
}


iftVoxelArray *iftApproxContour(const iftImage *contour, double epsilon)
{
    iftVoxelArray *array = iftContourToArray(contour);

    iftVoxelArray *approx_pts = iftApproxVoxelArray(array, epsilon);

    iftDestroyVoxelArray(&array);

    return approx_pts;
}


iftVoxelArray *iftApproxVoxelArray(const iftVoxelArray *array, double epsilon)
{
    iftVoxelArray *v1 = NULL, *v2 = NULL;
    iftSplitByMaxDist(array, &v1, &v2);

    iftVoxelArray *approx_v1 = iftAlloc(1, sizeof *approx_v1);
    approx_v1->n = 0;
    approx_v1->val = iftRDP(v1->val, v1->n, epsilon, &approx_v1->n);

    iftVoxelArray *approx_v2 = iftAlloc(1, sizeof *approx_v2);
    approx_v2->n = 0;
    approx_v2->val = iftRDP(v2->val, v2->n, epsilon, &approx_v2->n);

    iftDestroyVoxelArray(&v1);
    iftDestroyVoxelArray(&v2);

    iftVoxelArray *approx_pts = iftCreateVoxelArray(approx_v2->n + approx_v1->n - 2);
    int i = 0;
    for (; i < approx_v1->n; i++)
        iftCopyVoxel(&approx_v1->val[i], &approx_pts->val[i]);
    for (int j = 1; j < approx_v2->n - 1; i++, j++)
        iftCopyVoxel(&approx_v2->val[j], &approx_pts->val[i]);

    iftDestroyVoxelArray(&approx_v1);
    iftDestroyVoxelArray(&approx_v2);

    return approx_pts;
}


iftSet *iftVoxelArrayToSet(const iftImage *image, const iftVoxelArray *array)
{
    iftSet *set = NULL;
    for (int i = array->n - 1; i >= 0; i--)
    {
        int p = iftGetVoxelIndex(image, array->val[i]);
        iftInsertSet(&set, p);
    }
    return set;
}


iftSet *iftNearestInContour(const iftImage *contour, const iftSet *set, const iftImage *mask)
{
    // roots marks closest
    iftImage *roots = iftCreateImageFromImage(contour);
    iftImage *dist = iftGeodesicDistTransFromSet(mask, set, &roots);
    iftImage *index = iftCreateImageFromImage(contour);

    // resetting closest
    for (const iftSet *s = set; s; s = s->next) {
        dist->val[s->elem] = IFT_INFINITY_INT;
        index->val[s->elem] = s->elem;
    }

    for (int p = 0; p < contour->n; p++) {
        if (contour->val[p])
        {
            int r = roots->val[p];
            if (dist->val[p] < dist->val[r]) {
                index->val[r] = p;
                dist->val[r] = dist->val[p];
            }
        }
    }

    iftSet *closest = NULL;
    for (const iftSet *s = set; s; s = s->next)
    {
        iftInsertSet(&closest, index->val[s->elem]);
    }
    iftReverseSet(&closest);

    iftDestroyImage(&roots);
    iftDestroyImage(&dist);
    iftDestroyImage(&index);

    return closest;
}
