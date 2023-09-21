#include "iftSeeds.h"

#include "ift/core/dtypes/FIFO.h"
#include "ift/core/dtypes/GQueue.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/LIFO.h"
#include "ift/core/dtypes/ULongArray.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/Numerical.h"


/*************************** PRIVATE *****************************/
/**
 * @brief Find the best adjacency relation for the grid sampling, in which all its extreme displacements
 * (boundaries) have size/norm of at least <radius>.
 *
 * From an initial spheric adjacency relation of radius <radius>, it will increase its radius until
 * that all extreme displacements (boundaries) have size/norm of at least <radius>.
 * 
 * @param  radius    Radius in which all extreme displacements (boundaries) of the resulting adjacency
 *                   have size/norm of at least <radius>.
 * @param  is_3D_adj If true, a 3D spheric adjacency is applied. Otherwise, a 2D circular adjacency is used.
 * @return           Resulting adjacency relation.
 * 
 * @author Samuel Martins
 * @date Apr 27, 2018
 */
iftAdjRel *_iftFindBestAdjRelForGridSampling(float radius, bool is_3D_adj) {
    iftAdjRel *A = NULL;

    // Adjacency to compute the boundaries of the best adjacency relation for the grid sampling
    // We considered 8-neighborhood and 26-neighborhood to avoid possible leaks during sampling
    iftAdjRel *B = (is_3D_adj) ? iftSpheric(sqrtf(3.0)) : iftCircular(sqrtf(2.0));

    int it = 0;
    float epsilon = 0.1;
    float min_dist = 0.0;

    do {
        iftDestroyAdjRel(&A);

        A = (is_3D_adj) ? iftSpheric(radius + (it * epsilon)) : iftCircular(radius + (it * epsilon));
        it++;

        iftAdjRel *Abound = iftAdjacencyBoundaries(A, B);

        min_dist = IFT_INFINITY_FLT;
        for (int i = 0; i < Abound->n; i++) {
            float dist = sqrtf(Abound->dx[i]*Abound->dx[i] + Abound->dy[i]*Abound->dy[i] + Abound->dz[i]*Abound->dz[i]);
            if (dist < min_dist)
                min_dist = dist;
        }
        iftDestroyAdjRel(&Abound);
    } while (min_dist < radius);

    iftDestroyAdjRel(&B);    

    return A;
}


/**
 * @brief Performs the Differential IFT (DIFT) whose resulting cost map (pathvals) will be used for grid sampling.
 * 
 * It performs DIFT in the input image <img> with the labeled seeds <seeds> and adjacency relation <A>.
 * The cost function is the sum of the absolute difference of voxel values along the path.
 * 
 * It requires that the priority queue <Q>, cost map <pathvals>, and labeled image are previously allocated.
 * If the cost map already has precomputed path values (resulting from other iterations, for example),
 * it considers such costs for the current delineations;
 * 
 * The max priority queue <max_pathvals_queue> stores the voxels and their pathvalues as a max-heap.
 * It is used for rapidly selection of the voxels with highest costs in the map.
 * If it is NULL, nothing is done. 
 * 
 * @param img Input image used for segmentation.
 * @param seeds Seeds for segmentation.
 * @param A Adjacency Relation considered for the delineation.
 * @param Q Pre-alocated Priority Queue.
 * @param pathvals (Pre-computed) Cost Map.
 * @param label_img Label Image where the delineation is assigned.
 * @param max_pathvals_queue (Optional) Max Priority Queue to store the voxels according to the highest path values;
 * 
 * @author Samuel Martins
 * @date Sep 21, 2018
 */
void _iftPerformDIFTForGridSampling(const iftImage *img, const iftLabeledSet *seeds, const iftAdjRel *A,
                                    iftGQueue *Q, iftImage *pathvals, iftImage *label_img, iftGQueue *max_pathvals_queue) {
    if (Q == NULL)
        iftError("Priority Queue (Q) must be != NULL", "iftPerformIFT");
    if (pathvals == NULL)
        iftError("Cost Map (pathvals) must be != NULL", "iftPerformIFT");
    if (Q->L.value != pathvals->val)
        iftError("Priority Queue (Q) does not point to Cost Map (pathvals)", "iftPerformIFT");
    if (label_img == NULL)
        iftError("Label Image must be != NULL", "iftPerformIFT");

    // initialization
    const iftLabeledSet *node = seeds;
    while (node != NULL) {
        int p = node->elem;
        label_img->val[p] = node->label;
        pathvals->val[p] = 0;
        iftInsertGQueue(&Q, p);
        node = node->next;
    }


    while (!iftEmptyGQueue(Q)) {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftGetVoxelCoord(img, p);

        if (max_pathvals_queue != NULL) {
            if (max_pathvals_queue->L.elem[p].color == IFT_GRAY)
                iftRemoveGQueueElem(max_pathvals_queue, p);
            iftInsertGQueue(&max_pathvals_queue, p);
        }

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(img, v)) {
                int q = iftGetVoxelIndex(img, v);
                long tmp = pathvals->val[p] + abs(img->val[q] - img->val[p]);

                if (tmp < pathvals->val[q]) {
                    if (Q->L.elem[q].color == IFT_GRAY)
                        iftRemoveGQueueElem(Q, q);

                    pathvals->val[q] = tmp;
                    iftInsertGQueue(&Q, q);
                    label_img->val[q] = label_img->val[p];
                }
            }
        }
    }
}



/*************************** PUBLIC *****************************/

iftLabeledSet *     iftReadSeeds(const iftImage *img, const char *_filename, ...)
{
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];
    va_start(args, _filename);
    vsprintf(filename, _filename, args);
    va_end(args);

    
    FILE *fp=fopen(filename,"r");
    int i,l,mk,nseeds,handicap;
    int w,h,d;
    iftVoxel v;
    iftLabeledSet *seed=NULL;

    if (!iftIs3DImage(img)) {
        if (fscanf(fp,"%d %d %d",&nseeds, &w, &h)!=3)
            iftError("Reading error", "iftReadSeeds");

        if ((w != img->xsize)||(h != img->ysize)||(img->zsize != 1))
            iftError("Invalid 2D Seed File", "iftReadSeeds");

        v.z = 0;

        for (i=0; i < nseeds; i++){
            if (fscanf(fp,"%d %d %d %d",&v.x,&v.y,&mk,&l)!=4)
                iftError("Reading error", "iftReadSeeds");
            if (iftValidVoxel(img,v)){
                iftInsertLabeledSetMarkerAndHandicap(&seed, iftGetVoxelIndex(img,v), l, mk, 0);
            }
        }
    } else {
        if (fscanf(fp,"%d %d %d %d",&nseeds, &w, &h, &d)!=4) {
            iftError("Reading error", "iftReadSeeds");
        }

        for (i=0; i < nseeds; i++){
            if (fscanf(fp,"%d %d %d %d %d %d",&v.x,&v.y,&v.z,&mk,&l, &handicap)!=6) {
                iftError("Reading error", "iftReadSeeds");
            }
            if (iftValidVoxel(img,v)){
                iftInsertLabeledSetMarkerAndHandicap(&seed, iftGetVoxelIndex(img,v), l, mk, handicap);
            }
        }
    }

    fclose(fp);
    return(seed);
}

iftLabeledSet *iftMReadSeeds(iftMImage *img, char *filename)
{
    FILE *fp=fopen(filename,"r");
    int i,l,mk,nseeds,w,h,d,handicap;
    iftVoxel v;
    iftLabeledSet *seed=NULL;

    if (!iftIs3DMImage(img)) {
        if (fscanf(fp,"%d %d %d",&nseeds, &w, &h)!=3)
            iftError("Reading error", "MReadSeeds");

        if ((w != img->xsize)||(h != img->ysize)||(img->zsize != 1))
            iftError("Invalid 2D Seed File", "MReadSeeds");

        v.z = 0;


        for (i=0; i < nseeds; i++){
            if (fscanf(fp,"%d %d %d %d",&v.x,&v.y,&mk,&l)!=4)
                iftError("Reading error", "MReadSeeds");
            if (iftMValidVoxel(img,v)){
                iftInsertLabeledSetMarkerAndHandicap(&seed, iftMGetVoxelIndex(img,v), l, mk, 0);
            }
        }
    } else {
        if (fscanf(fp,"%d %d %d %d",&nseeds, &w, &h, &d)!=4) {
            iftError("Reading error", "MReadSeeds");
        }

        for (i=0; i < nseeds; i++){
            if (fscanf(fp,"%d %d %d %d %d %d",&v.x,&v.y,&v.z,&mk,&l, &handicap)!=6)
                iftError("Reading error", "MReadSeeds");
            if (iftMValidVoxel(img,v)){
                iftInsertLabeledSetMarkerAndHandicap(&seed, iftMGetVoxelIndex(img,v), l, mk, handicap);
            }
        }
    }

    fclose(fp);
    return(seed);
}



iftLabeledSet *iftReadSeedsGraph(char *filename)
{
    FILE *fp=fopen(filename,"r");
    int i,l,mk,id,nseeds;
    iftLabeledSet *seed=NULL;

    if (fscanf(fp,"%d\n",&nseeds)!=1)
        iftError("Reading error", "iftReadSeedsGraph");

    for (i=0; i < nseeds; i++){
        if (fscanf(fp,"%d %d %d\n",&id,&mk,&l)!=3)
            iftError("Reading error", "MReadSeeds");
        iftInsertLabeledSetMarkerAndHandicap(&seed, id, l, mk, 0);
    }

    fclose(fp);
    return(seed);
}

void iftWriteSeedsGraph(iftLabeledSet* seed, const char *_filename, ...)
{
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];
    va_start(args, _filename);
    vsprintf(filename, _filename, args);
    va_end(args);

    FILE *file = fopen(filename,"w");
    if(file == NULL)
        iftError("Invalid destination file", "iftWriteSeedsGraph");

    iftLabeledSet *s = seed;
    int nseeds = 0;

    while(s != NULL){
        nseeds++;
        s = s->next;
    }
    
    fprintf(file,"%d\n", nseeds);
    s = seed;
    while(s != NULL){
      fprintf(file, "%d %d %d\n", s->elem, s->marker, s->label);
      s = s->next;
    }

    fclose(file);
}



void iftMWriteSeeds(iftLabeledSet* seed, iftMImage* img, char *filename)
{
  FILE *file = fopen(filename,"w");
  if(file == NULL)
    iftError("Invalid destination file", "iftMWriteSeeds");
  
  iftLabeledSet *s = seed;
  int nseeds = 0;
  
  while(s != NULL){
    nseeds++;
    s = s->next;
  }
  
  if (!iftIs3DMImage(img)) {
    fprintf(file,"%d %d %d\n", nseeds, img->xsize, img->ysize);
    s = seed;
    while(s != NULL){
      iftVoxel voxel = iftMGetVoxelCoord(img,s->elem);
      fprintf(file, "%d %d %d %d\n", voxel.x, voxel.y, s->marker, s->label);
      s = s->next;
    }
  } else {
    fprintf(file,"%d %d %d %d\n", nseeds, img->xsize, img->ysize, img->zsize);
    s = seed;
    while(s != NULL){
      iftVoxel voxel = iftMGetVoxelCoord(img,s->elem);
      fprintf(file, "%d %d %d %d %d %d\n", voxel.x, voxel.y, voxel.z, s->marker, s->label, s->handicap);
      s = s->next;
    }
  }
  
  fclose(file);
}

void iftWriteSeeds(iftLabeledSet* seed, const iftImage* img, const char *_filename, ...)
{
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];
    va_start(args, _filename);
    vsprintf(filename, _filename, args);
    va_end(args);

    FILE *file = fopen(filename,"w");
    if(file == NULL)
        iftError("Invalid destination file", "iftWriteSeeds");

    iftLabeledSet *s = seed;
    int nseeds = 0;

    while(s != NULL){
        nseeds++;
        s = s->next;
    }

    if (!iftIs3DImage(img)) {
        fprintf(file,"%d %d %d\n", nseeds, img->xsize, img->ysize);
        s = seed;
        while(s != NULL){
            iftVoxel voxel = iftGetVoxelCoord(img,s->elem);
            fprintf(file, "%d %d %d %d\n", voxel.x, voxel.y, s->marker, s->label);
            s = s->next;
        }
    } else {
        fprintf(file,"%d %d %d %d\n", nseeds, img->xsize, img->ysize, img->zsize);
        s = seed;
        while(s != NULL){
            iftVoxel voxel = iftGetVoxelCoord(img,s->elem);
            fprintf(file, "%d %d %d %d %d %d\n", voxel.x, voxel.y, voxel.z, s->marker, s->label, s->handicap);
            s = s->next;
        }
    }

    fclose(file);
}

iftSet * iftExtractRemovalMarkers(iftLabeledSet **s)
{
  iftSet *removal_markers = NULL;
  iftSet *set_aux;
  iftLabeledSet *aux = *s;

  while (aux != NULL)
  {
    if (aux->label == IFT_NIL)
      iftInsertSet(&removal_markers, aux->elem);
    aux = aux->next;
  }
  set_aux = removal_markers;
  while (set_aux != NULL)
  {
    iftRemoveLabeledSetElem(s, set_aux->elem);
    set_aux = set_aux->next;
  }
  return removal_markers;

}



iftSet *iftImageBorderSet(const iftImage *img)
{
  iftSet *border=NULL;
  iftVoxel u;
  int p;

  if (img->zsize > 1){

    u.z=0;
    for (u.y=0; u.y < img->ysize; u.y++)
      for (u.x=0; u.x < img->xsize; u.x++){
	p = iftGetVoxelIndex(img,u);
	iftInsertSet(&border,p);
      }

    u.z=img->zsize-1;
    for (u.y=0; u.y < img->ysize; u.y++)
      for (u.x=0; u.x < img->xsize; u.x++){
	p = iftGetVoxelIndex(img,u);
	iftInsertSet(&border,p);
      }

  u.x=0;
  for (u.z=1; u.z < img->zsize-1; u.z++)
    for (u.y=0; u.y < img->ysize; u.y++){
      p = iftGetVoxelIndex(img,u);
      iftInsertSet(&border,p);
    }
  u.x=img->xsize-1;
  for (u.z=1; u.z < img->zsize-1; u.z++)
    for (u.y=0; u.y < img->ysize; u.y++){
      p = iftGetVoxelIndex(img,u);
      iftInsertSet(&border,p);
    }

  u.y=0;
  for (u.z=1; u.z < img->zsize-1; u.z++)
    for (u.x=1; u.x < img->xsize-1; u.x++){
      p = iftGetVoxelIndex(img,u);
      iftInsertSet(&border,p);
    }
  u.y=img->ysize-1;
  for (u.z=1; u.z < img->zsize-1; u.z++)
    for (u.x=1; u.x < img->xsize-1; u.x++){
      p = iftGetVoxelIndex(img,u);
      iftInsertSet(&border,p);
    }

  }else{ // 2D

    u.z=0;

    u.x=0;
    for (u.y=0; u.y < img->ysize; u.y++){
      p = iftGetVoxelIndex(img,u);
      iftInsertSet(&border,p);
    }
    u.x=img->xsize-1;
    for (u.y=0; u.y < img->ysize; u.y++){
      p = iftGetVoxelIndex(img,u);
      iftInsertSet(&border,p);
    }

    u.y=0;
    for (u.x=1; u.x < img->xsize-1; u.x++){
      p = iftGetVoxelIndex(img,u);
      iftInsertSet(&border,p);
    }

    u.y=img->ysize-1;
    for (u.x=1; u.x < img->xsize-1; u.x++){
      p = iftGetVoxelIndex(img,u);
      iftInsertSet(&border,p);
    }
  }

  return(border);
}

iftLabeledSet *iftImageBorderLabeledSet(iftImage *img)
{
  iftLabeledSet *border=NULL;
  iftVoxel u;
  int p;

  if (img->zsize > 1){

    u.z=0;
    for (u.y=0; u.y < img->ysize; u.y++)
      for (u.x=0; u.x < img->xsize; u.x++){
	p = iftGetVoxelIndex(img,u);
	iftInsertLabeledSet(&border,p,0);
      }

    u.z=img->zsize-1;
    for (u.y=0; u.y < img->ysize; u.y++)
      for (u.x=0; u.x < img->xsize; u.x++){
	p = iftGetVoxelIndex(img,u);
	iftInsertLabeledSet(&border,p,0);
      }

  u.x=0;
  for (u.z=1; u.z < img->zsize-1; u.z++)
    for (u.y=0; u.y < img->ysize; u.y++){
      p = iftGetVoxelIndex(img,u);
      iftInsertLabeledSet(&border,p,0);
    }
  u.x=img->xsize-1;
  for (u.z=1; u.z < img->zsize-1; u.z++)
    for (u.y=0; u.y < img->ysize; u.y++){
      p = iftGetVoxelIndex(img,u);
      iftInsertLabeledSet(&border,p,0);
    }

  u.y=0;
  for (u.z=1; u.z < img->zsize-1; u.z++)
    for (u.x=1; u.x < img->xsize-1; u.x++){
      p = iftGetVoxelIndex(img,u);
      iftInsertLabeledSet(&border,p,0);
    }
  u.y=img->ysize-1;
  for (u.z=1; u.z < img->zsize-1; u.z++)
    for (u.x=1; u.x < img->xsize-1; u.x++){
      p = iftGetVoxelIndex(img,u);
      iftInsertLabeledSet(&border,p,0);
    }

  }else{ // 2D

    u.z=0;

    u.x=0;
    for (u.y=0; u.y < img->ysize; u.y++){
      p = iftGetVoxelIndex(img,u);
      iftInsertLabeledSet(&border,p,0);
    }
    u.x=img->xsize-1;
    for (u.y=0; u.y < img->ysize; u.y++){
      p = iftGetVoxelIndex(img,u);
      iftInsertLabeledSet(&border,p,0);
    }

    u.y=0;
    for (u.x=1; u.x < img->xsize-1; u.x++){
      p = iftGetVoxelIndex(img,u);
      iftInsertLabeledSet(&border,p,0);
    }

    u.y=img->ysize-1;
    for (u.x=1; u.x < img->xsize-1; u.x++){
      p = iftGetVoxelIndex(img,u);
      iftInsertLabeledSet(&border,p,0);
    }
  }

  return(border);
}



iftSet *iftObjectBorderSet(const iftImage *label_img, iftAdjRel *Ain) {
    iftAdjRel *A = Ain;
    if (A == NULL) {
        if (iftIs3DImage(label_img))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }

    iftSet *borders = NULL;

    for (int p = 0; p < label_img->n; p++) {
        if (label_img->val[p] != 0) {
            iftVoxel u = iftGetVoxelCoord(label_img, p);
            
            for (int i = 1; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                if (iftValidVoxel(label_img, v)) {
                    int q = iftGetVoxelIndex(label_img, v);
                    
                    if (label_img->val[q] != label_img->val[p]) {
                        iftInsertSet(&borders, p);
                        break;
                    }
                }
                else {
                    // p is an object pixel which is on the image's border
                    iftInsertSet(&borders, p);
                    break;
                }
            }
        }
    }


    if (Ain == NULL)
        iftDestroyAdjRel(&A);

    return borders;
}

/* It assumes that the object does not touch the image border */

iftSet *iftBackgroundBorderSet(const iftImage *label_img, iftAdjRel *Ain) {
    iftAdjRel *A = Ain;
    if (A == NULL) {
        if (iftIs3DImage(label_img))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }

    iftSet *borders = NULL;

    for (int p = 0; p < label_img->n; p++) {
        if (label_img->val[p] == 0) {
            iftVoxel u = iftGetVoxelCoord(label_img, p);
            
            for (int i = 1; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                if (iftValidVoxel(label_img, v)) {
                    int q = iftGetVoxelIndex(label_img, v);
                    
                    if (label_img->val[q] != label_img->val[p]) {
                        iftInsertSet(&borders, p);
                        break;
                    }
                }
            }
        }
    }


    if (Ain == NULL)
        iftDestroyAdjRel(&A);

    return borders;
}


iftBMap *iftObjectBorderBitMap(const iftImage *label, const iftAdjRel *A, int *n_border_spels) {
    iftBMap *border;
    iftVoxel u, v;
    int p, q;

    // A bitmap with the number of elements (bits) of the label image.
    // It will contain 1 if the object is on the border and 0 otherwise.
    border = iftCreateBMap(label->n);

    (*n_border_spels) = 0;
    for (u.z = 0; u.z < label->zsize; u.z++) {
        for (u.y = 0; u.y < label->ysize; u.y++) {
            for (u.x = 0; u.x < label->xsize; u.x++) {
                p = iftGetVoxelIndex(label,u);

                if (label->val[p] != 0) {
                    for (int i = 1; i < A->n; i++) {
                        v.x = u.x + A->dx[i];
                        v.y = u.y + A->dy[i];
                        v.z = u.z + A->dz[i];

                        if (iftValidVoxel(label, v)) {
                            q = iftGetVoxelIndex(label, v);
                            
                            if (label->val[q] != label->val[p]) {
                                // Sets this voxel as a border voxel and increments their quantitiy.
                                iftBMapSet1(border, p);
                                (*n_border_spels) = (*n_border_spels) + 1;
                                break;
                            }
                        }
                        else { 
                            // Sets this voxel as a border voxel and increments their quantitiy.
                            iftBMapSet1(border, p);
                            (*n_border_spels) = (*n_border_spels) + 1;
                            break;
                        }
                    }
                }
            }
        }
    }

    return border;
}


void iftMaskImageToSet(const iftImage *bin_img, iftSet **S) {
    if (bin_img == NULL)
        iftError("Image is NULL", "iftMaskImageToSet");
    if (S == NULL)
        iftError("Set is NULL", "iftMaskImageToSet");

    for (int p = 0; p < bin_img->n; p++) {
        if (bin_img->val[p] != 0) {
            iftInsertSet(S, p);
        }
    }
}

iftImage *iftMaskImageFromSet(iftSet *S, int xsize, int ysize, int zsize)
{
  iftImage *mask = iftCreateImage(xsize, ysize, zsize);
  iftSet *seed = S;
  while (seed != NULL) {
    int p = seed->elem;
    mask->val[p] = 1;
    seed = seed->next;
  }

  return(mask);
}


void iftSetToImage(const iftSet *S, iftImage *img, int obj) {
    for (const iftSet *node = S; node != NULL; node = node->next)
        img->val[node->elem] = obj;
}

void iftListToImage(const iftList *L, iftImage *img, int obj) {
    for (const iftNode *node = L->head; node != NULL; node = node->next)
        img->val[node->elem] = obj;
}


iftList *iftMaskToList(const iftImage *mask) {
    iftList *L = iftCreateList();

    for (int p = 0; p < mask->n; p++)
        if (mask->val[p] != 0)
            iftInsertListIntoTail(L, p);

    return L;
}


void iftLabeledSetToImage(const iftLabeledSet *S, iftImage *img, bool increment_label) {
    for (const iftLabeledSet *node = S; node != NULL; node = node->next)
        img->val[node->elem] = node->label + increment_label;
}

void iftIntArrayToImage(const iftIntArray *iarr, iftImage *img, int obj) {
    #pragma omp parallel for
    for (long i = 0; i < iarr->n; i++)
        img->val[iarr->val[i]] = obj;
}


iftIntArray *iftMaskToIntArray(const iftImage *mask) {
    iftList *L = iftMaskToList(mask);

    iftIntArray *iarr = iftCreateIntArray(L->n);

    for (int i = 0; L->n != 0; i++) {
        int p = iftRemoveListHead(L);
        iarr->val[i] = p;
    }

    iftDestroyList(&L);

    return iarr;
}



iftLabeledSet *iftMultiObjectBorderLabeledSet(iftImage *img, iftAdjRel *A)
{
  iftLabeledSet *border=NULL;
  int i,p,q;
  iftVoxel u,v;

  for (u.z=0; u.z < img->zsize; u.z++)
    for (u.y=0; u.y < img->ysize; u.y++)
      for (u.x=0; u.x < img->xsize; u.x++)
      {
        p = iftGetVoxelIndex(img,u);
        if (img->val[p]!=0)
        {
          for (i=1; i < A->n; i++)
          {
            v.x = u.x + A->dx[i];
            v.y = u.y + A->dy[i];
            v.z = u.z + A->dz[i];
            if (iftValidVoxel(img,v))
            {
              q = iftGetVoxelIndex(img,v);
              if (img->val[q]!=img->val[p])
              {
                iftInsertLabeledSet(&border, p, img->val[p]);
                break;
              }
            }
            else
            {
              iftInsertLabeledSet(&border, p, img->val[p]);
              break;
            }
          }
        }
      }
  return(border);
}


iftImage *iftObjectBorders(const iftImage *label_img, const iftAdjRel *Ain, bool keep_border_labels,
                           bool include_image_frame) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = (iftIs3DImage(label_img)) ? iftSpheric(1.0) : iftCircular(1.0);
    else A = iftCopyAdjacency(Ain);

    iftImage *border_img = iftCreateImageFromImage(label_img);

    for (int p = 0; p < label_img->n; p++) {
        if (label_img->val[p] != 0) {
            iftVoxel u = iftGetVoxelCoord(label_img, p);

            for (int i = 1; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                if (iftValidVoxel(label_img, v)) {
                    if (iftImgVoxelVal(label_img, u) != iftImgVoxelVal(label_img, v)) {
                        iftImgVoxelVal(border_img, u) = iftImgVoxelVal(label_img, u);
                        break;
                    }
                }
                // voxel u belongs to the image's frame
                else if (include_image_frame) {
                    iftImgVoxelVal(border_img, u) = iftValidVoxel(label_img, u);
                    break;
                }
            }
        }
    }

    if (!keep_border_labels) {
        iftImage *aux = border_img;
        border_img = iftBinarize(border_img);
        iftDestroyImage(&aux);
    }

    iftDestroyAdjRel(&A);

    return border_img;
}


iftLabeledSet *iftRegionsToLabeledSet(const iftImage *label, bool decrement_label)
{
    iftLabeledSet *S = NULL;
    for (int p = 0; p < label->n; p++) {
        if (label->val[p] > 0) {
            iftInsertLabeledSet(&S, p, label->val[p] - decrement_label);
        }
    }

    return S;
}


iftSet *iftObjectToSet(const iftImage *label_img) {
    iftSet *S = NULL;

    for (int p = 0; p < label_img->n; p++) {
        if (label_img->val[p] > 0) {
            iftInsertSet(&S, p);
        }
    }

    return S;
}

iftLabeledSet *iftObjectToLabeledSet(const iftImage *label_img) {
    iftLabeledSet *S = NULL;

    for (int p = 0; p < label_img->n; p++) {
        if (label_img->val[p] > 0) {
            iftInsertLabeledSet(&S, p,label_img->val[p]);
        }
    }

    return S;
}


iftList *iftObjectsToList(const iftImage *label_img) {
    iftList *L = iftCreateList();

    for (int p = 0; p < label_img->n; p++) {
        if (label_img->val[p] > 0) {
            iftInsertListIntoTail(L, label_img->val[p]);
        }
    }

    return L;   
}


iftLabeledSet *iftLabelObjBorderSet(iftImage *bin, iftAdjRel *A)
{
  iftLabeledSet *seed=NULL;
  iftImage *label;
  int p;

  label = iftFindAndLabelObjectBorders(bin,A);
  for (p=0; p < label->n; p++) {
    if (label->val[p]!=0){
      iftInsertLabeledSet(&seed,p,label->val[p]);
    }
  }
  iftDestroyImage(&label);

  return(seed);
}


iftImage *iftFindAndLabelObjectBorders(const iftImage *label_img, const iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = (iftIs3DImage(label_img)) ? iftSpheric(sqrtf(3)) : iftCircular(sqrtf(2));
    else A = iftCopyAdjacency(Ain);

    iftImage *border_img = iftObjectBorders(label_img, NULL, true, true);
    iftImage *border_label_img = iftFastLabelComp(border_img, A);

    iftDestroyImage(&border_img);
    
    iftDestroyAdjRel(&A);

    return border_label_img;
}


iftImage  *iftEasyLabelComp(iftImage *bin, iftAdjRel *A)
{
  iftImage *label=NULL;
  int i,j,p,q,l=1;
  iftVoxel u,v;
  iftFIFO *F=NULL;

  label  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  iftCopyVoxelSize(bin,label);
  F      = iftCreateFIFO(bin->n);

  for (j=0; j < bin->n; j++){
    if ((bin->val[j]!=0)&&(label->val[j]==0)){
      label->val[j]=l;
      iftInsertFIFO(F,j);
      while(!iftEmptyFIFO(F)){
	p = iftRemoveFIFO(F);
	u = iftGetVoxelCoord(bin,p);
	for (i=1; i < A->n; i++){
	  v = iftGetAdjacentVoxel(A,u,i);
	  if (iftValidVoxel(bin,v)){
	    q = iftGetVoxelIndex(bin,v);
	    if ((bin->val[q]!=0)&&(label->val[q] == 0)){
	      label->val[q] = label->val[p];
	      iftInsertFIFO(F,q);
	    }
	  }
	}
      }
      l++;
    }
  }

  iftDestroyFIFO(&F);
  return(label);
}

iftImage *iftFastLabelComp(const iftImage *bin, const iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(bin))
            A = iftSpheric(1.4);
        else A = iftCircular(1.5);
    }
    else A = iftCopyAdjacency(Ain);

  iftImage *label=NULL;
  int i,p,q,l=1, *cost;
  iftVoxel u,v;
  iftGQueue *Q;

  label  = iftCreateImageFromImage(bin);
  cost   = iftAllocIntArray(bin->n);
  Q      = iftCreateGQueue(2,bin->n,cost);

  for (p=0; p < bin->n; p++){
    if ((bin->val[p]!=0)&&(label->val[p]==0)){
      cost[p] = 1;
      iftInsertGQueue(&Q,p);
    }
  }

  while(!iftEmptyGQueue(Q)){
    p = iftRemoveGQueue(Q);
    if (cost[p]==1){ /* p is the first in a component */
      cost[p]=0;
      label->val[p] = l; l++;
    }
    u = iftGetVoxelCoord(bin,p);
    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(bin,v)){
	q = iftGetVoxelIndex(bin,v);
	if ((bin->val[p] == bin->val[q])&&(label->val[q] == 0)){
	      label->val[q] = label->val[p];
	      iftRemoveGQueueElem(Q,q);
	      cost[q]=0;
	      iftInsertGQueue(&Q,q);
	}
      }
    }
  }

  iftDestroyAdjRel(&A);

  iftDestroyGQueue(&Q);
  iftFree(cost);
  return(label);
}

int iftRootVoxel(iftImage *pred, int p)
{
  if (pred->val[p]==p)
    return(p);
  else
    return(pred->val[p] = iftRootVoxel(pred,pred->val[p])); /* path
							       compression */
}

/*
  Atribui o mesmo rótulo a sementes conectadas 
  no espaço da imagem conforme a adjacência A.
  Os rótulos são valores inteiros sequenciais a partir de 1.
  */
iftImage  *iftLabelComp(iftImage *bin, iftAdjRel *A)
{
  iftImage *label=NULL, *pred=NULL;
  int l = 1;

  pred  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);

  for (int p=0;  p < pred->n; p++)
    pred->val[p]=p;

  for (int p=0;  p < bin->n; p++) {
    if (bin->val[p]!=0) { // posições com sementes
      int rp = iftRootVoxel(pred,p);
      iftVoxel u = iftGetVoxelCoord(bin,p);

      for (int i=1; i < A->n; i++){
      	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      	if (iftValidVoxel(bin,v)){
      	  int q = iftGetVoxelIndex(bin,v);
      	  if (bin->val[q]!=0){ // posições vizinhas com sementes
      	    int rq = iftRootVoxel(pred,q);
	          if (rq != rp){
      	      if (rp < rq)
		            pred->val[rq] = rp;
      	      else
            		pred->val[rp] = rq;
      	    }
	        }
	      }
      }
    }
  }

  label = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);

  for (int p=0;  p < pred->n; p++) {
    if (bin->val[p]!=0){
      if (pred->val[p]==p){
      	label->val[p]=l; l++;
      }else{
      	if (pred->val[pred->val[p]]!=pred->val[p])
      	  pred->val[p]=iftRootVoxel(pred,p);
      }
    }
  }

#pragma omp parallel for shared(pred,bin,label)
  for (int p=0;  p < pred->n; p++)
    if ((bin->val[p]!=0)&&(pred->val[p]!=p)){
      label->val[p]=label->val[pred->val[p]];
    }

  iftCopyVoxelSize(bin,label);
  iftDestroyImage(&pred);

  return(label);
}


iftLabeledSet *iftLabelCompSet(iftImage *bin, iftAdjRel *A)
{
  iftImage *label=iftFastLabelComp(bin,A);
  iftLabeledSet *seed=NULL;
  int p;

  for (p=0; p < label->n; p++)
    if (label->val[p] > 0)
      iftInsertLabeledSet(&seed,p,label->val[p]);

  iftDestroyImage(&label);

  return(seed);
}


iftImage *iftSelectLargestComp(const iftImage *bin, const iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = (iftIs3DImage(bin)) ? iftSpheric(1.74) : iftCircular(1.45);
    else A = iftCopyAdjacency(Ain);

    iftImage *label_img = iftFastLabelComp(bin, A);

    int ncomps = iftMaximumValue(label_img);
    int *size  = iftAllocIntArray(ncomps+1);

    for (int p = 0; p < label_img->n; p++)
        if (bin->val[p])
            size[label_img->val[p]]++;

    int maxsize = 0;
    int imax = IFT_NIL;
    for (int i = 1; i <= ncomps; i++)
        if (size[i] > maxsize) {
            maxsize = size[i];
            imax = i;
        }

    for (int p = 0; p < label_img->n; p++)
        label_img->val[p] = ((label_img->val[p] == imax) ? bin->val[p] : 0);

    iftFree(size);
    if (Ain != NULL)
        iftDestroyAdjRel(&A);

    return label_img;
}

iftImage *iftSelectKLargestComp(iftImage *bin, iftAdjRel *A, int K)
{
  iftImage *label, *nbin;
  int       ncomps,p,i,*index,*size;

  if (K <= 0)
      iftError("Invalid number of components", "iftSelectKLargestComp");

  label  = iftFastLabelComp(bin,A);
  ncomps = iftMaximumValue(label);
  size   = iftAllocIntArray(ncomps+1);
  index  = iftAllocIntArray(ncomps+1);

  for (i=0; i <= ncomps; i++) {
    index[i]=i;
    size[i] =0.0;
  }

  for (p=0; p < label->n; p++)
    if (bin->val[p]){
      size[label->val[p]]++;
    }

  iftBucketSort(size, index, ncomps+1, IFT_DECREASING);

  iftFree(size);

  nbin = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);

  for (p=0; p < label->n; p++) {
    for (i=0; i < K; i++) {
      if (label->val[p]==index[i])
	nbin->val[p]=bin->val[p];
    }
  }

  iftCopyVoxelSize(bin,nbin);
  iftDestroyImage(&label);
  iftFree(index);

  return(nbin);
}

iftImage *iftSelectKSmallestComp(iftImage *bin, iftAdjRel *A, int K)
{
  iftImage *label, *nbin;
  int       ncomps,p,i,*index,*size;

  if (K <= 0)
      iftError("Invalid number of components", "iftSelectKSmallestComp");

  label  = iftFastLabelComp(bin,A);
  ncomps = iftMaximumValue(label);
  size   = iftAllocIntArray(ncomps+1);
  index  = iftAllocIntArray(ncomps+1);

  for (i=0; i <= ncomps; i++) {
    index[i]=i;
    size[i] =0.0;
  }

  for (p=0; p < label->n; p++)
    if (bin->val[p]){
      size[label->val[p]]++;
    }

  iftBucketSort(size, index, ncomps+1, IFT_INCREASING);

  iftFree(size);

  nbin = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);

  for (p=0; p < label->n; p++) {
    for (i=0; i < K; i++) {
      if (label->val[p]==index[i])
	nbin->val[p]=bin->val[p];
    }
  }

  iftCopyVoxelSize(bin,nbin);
  iftDestroyImage(&label);
  iftFree(index);

  return(nbin);
}

iftImage *iftSelectSmallestComp(iftImage *bin, iftAdjRel *A)
{
  iftImage *label;
  int *size,ncomps,p,i,imin,minsize;

  label  = iftFastLabelComp(bin,A);
  ncomps = iftMaximumValue(label);
  size   = iftAllocIntArray(ncomps+1);
  for (p=0; p < label->n; p++)
    if (bin->val[p])
      size[label->val[p]]++;

  minsize= IFT_INFINITY_INT; imin= IFT_NIL;
  for (i=1; i <= ncomps; i++)
    if (size[i] < minsize){
      minsize=size[i]; imin=i;
    }

  for (p=0; p < label->n; p++) {
    if (label->val[p]==imin)
      label->val[p]=bin->val[p];
    else
      label->val[p]=0;
  }

  iftFree(size);
  iftCopyVoxelSize(bin,label);
  return(label);
}


iftImage *iftComponentArea(iftImage *bin, iftAdjRel *A)
{
  iftImage *label;
  int      *size, ncomps, p;

  label  = iftFastLabelComp(bin,A);
  ncomps = iftMaximumValue(label);
  size   = iftAllocIntArray(ncomps+1);
  for (p=0; p < label->n; p++)
    if (label->val[p] != 0)
      size[label->val[p]]++;

  for (p=0; p < label->n; p++)
    label->val[p] = size[label->val[p]];

  iftFree(size);
  iftCopyVoxelSize(bin,label);

  return(label);
}

iftImage *iftSelectCompAboveArea(iftImage *bin, iftAdjRel *A, int thres)
{
  iftImage *area  = iftComponentArea(bin, A);

  for (int p=0; p < area->n; p++)
    if (area->val[p] >= thres)
      area->val[p] = bin->val[p];
    else
      area->val[p] = 0;

  iftCopyVoxelSize(bin,area);
  return(area);
}

iftImage *iftSelectCompBelowArea(iftImage *bin, iftAdjRel *A, int thres)
{
  iftImage *area  = iftComponentArea(bin, A);

  for (int p=0; p < area->n; p++)
    if (area->val[p] <= thres)
      area->val[p] = bin->val[p];
    else
      area->val[p] = 0;

  iftCopyVoxelSize(bin,area);
  return(area);
}

iftImage *iftSelectCompInAreaInterval(iftImage *bin, iftAdjRel *A, int thres_min, int thres_max)
{
  iftImage *area  = iftComponentArea(bin, A);
  for (int p=0; p < area->n; p++){
    if ((area->val[p] >= thres_min)&&
	(area->val[p] <= thres_max) )
      area->val[p] = bin->val[p];
    else
      area->val[p] = 0;
  }
  
  
  iftCopyVoxelSize(bin,area);
  return(area);
}

iftImage *iftRegionArea(iftImage *label)
{
  int      *size, ncomps, p;
  iftImage *area;

  area   = iftCreateImage(label->xsize,label->ysize,label->zsize);
  ncomps = iftMaximumValue(label);
  size   = iftAllocIntArray(ncomps+1);
  for (p=0; p < label->n; p++)
    if (label->val[p]!=0)
      size[label->val[p]]++;

  for (p=0; p < label->n; p++)
    area->val[p] = size[label->val[p]];

  iftFree(size);
  iftCopyVoxelSize(label,area);

  return(area);
}

iftImage *iftSelectLargestRegion(const iftImage *label_img) {
    iftImage *bin = iftCreateImageFromImage(label_img);
    int ncomps    = iftMaximumValue(label_img);
    int *size     = iftAllocIntArray(ncomps+1);
    
    for (int p = 0; p < label_img->n; p++)
        size[label_img->val[p]]++;

    int maxsize = 0;
    int imax    = IFT_NIL;
    for (int i = 1; i <= ncomps; i++)
        if (size[i] > maxsize) {
            maxsize = size[i];
            imax = i;
        }

    for (int p = 0; p < label_img->n; p++) {
        if (label_img->val[p] == imax)
            bin->val[p] = imax;
        else bin->val[p] = 0;
    }
  iftFree(size);

    return bin;
}

iftImage *iftSelectSmallestRegion(iftImage *label)
{
  iftImage *bin;
  int *size,ncomps,p,i,imin,minsize;

  bin    = iftCreateImage(label->xsize,label->ysize,label->zsize);
  ncomps = iftMaximumValue(label);
  size   = iftAllocIntArray(ncomps+1);
  for (p=0; p < label->n; p++)
    size[label->val[p]]++;

  minsize= IFT_INFINITY_INT; imin= IFT_NIL;
  for (i=1; i <= ncomps; i++)
    if (size[i] < minsize){
      minsize=size[i]; imin=i;
    }

  for (p=0; p < label->n; p++) {
    if (label->val[p]==imin)
      bin->val[p]=1;
    else
      bin->val[p]=0;
  }

  iftFree(size);
  iftCopyVoxelSize(label,bin);
  return(bin);
}

iftImage *iftSelectRegionsAboveAreaAndPropagateTheirLabels(iftImage *label, int thres) {
  int       *area   = NULL;
  iftImage  *nlabel = NULL;
  int       i, j;
  iftSet    *S      = NULL;
  iftAdjRel *A;

  if (iftIs3DImage(label))
    A = iftSpheric(sqrtf(3.0));
  else
    A = iftCircular(sqrtf(2.0));

  /* relabel components: it accounts for disconnected labels, which
     should be considered multiple components */

  nlabel = iftRelabelRegions(label, A);

  /* compute area of each region */
  int nlabel_max_val = iftMaximumValue(nlabel);
  area = iftAllocIntArray(nlabel_max_val + 1);
  for (int p = 0; p < nlabel->n; p++)
    area[nlabel->val[p]]++;

  /* eliminate regions below area and store the new labels in the area
     vector */

  area[0] = 0;
  for (i = 1, j = 1; i <= nlabel_max_val; i++) {
    if (area[i] < thres) {
      area[i] = 0;
    } else {
      area[i] = j; /* store new label */
      j++;
    }
  }
  
  /* set/get new labels and assign them to the new label image. Insert
     in the set the labeled voxels */

  for (int p = 0; p < nlabel->n; p++) {
    nlabel->val[p] = area[nlabel->val[p]];
    if (nlabel->val[p] != 0)
      iftInsertSet(&S, p);
  }
  
  iftCopyVoxelSize(label, nlabel);
  iftFree(area);

  /* fill regions with label zero from the labeled voxels in the
     set */

  while (S != NULL) {
    int      p = iftRemoveSet(&S);
    iftVoxel u = iftGetVoxelCoord(nlabel, p);
    for (i = 1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, i);
      if (iftValidVoxel(nlabel, v)) {
        int q = iftGetVoxelIndex(nlabel, v);
        if (nlabel->val[q] == 0) {
          nlabel->val[q] = nlabel->val[p];
          iftInsertSet(&S, q);
        }
      }
    }
  }

  iftDestroyAdjRel(&A);
  iftDestroySet(&S);

  return (nlabel);
}

iftImage *iftSelectRegionsAboveArea(iftImage *label, int thres) {
  int      label_max_value = iftMaximumValue(label);
  int      *area           = iftAllocIntArray(label_max_value + 1);
  iftImage *nlabel         = iftCopyImage(label);
  int      i, j;

  /* compute area of each region */

  for (int p = 0; p < label->n; p++)
    area[label->val[p]]++;

  /* eliminate regions below area and store the new labels in the area vector */

  area[0] = 0;
  for (i = 1, j = 1; i <= label_max_value; i++) {
    if (area[i] < thres)
      area[i] = 0;
    else {
      area[i] = j; /* store new label */
      j++;
    }
  }
  
  /* set/get new labels and assign them to the new label image */

  for (int p = 0; p < label->n; p++) {
    nlabel->val[p] = area[label->val[p]];
  }
  iftCopyVoxelSize(label, nlabel);

  iftFree(area);

  return (nlabel);
}

iftImage *iftSelectRegionsBelowArea(iftImage *label, int thres) {
  int      label_max_value = iftMaximumValue(label);
  int      *area           = iftAllocIntArray(label_max_value + 1);
  iftImage *nlabel         = iftCopyImage(label);

  /* compute area of each region */

  for (int p = 0; p < label->n; p++)
    area[label->val[p]]++;

  /* eliminate regions below area and store the new labels in the area vector */

  area[0] = 0;
  for (int i = 1, j = 1; i <= label_max_value; i++) {
    if (area[i] >= thres)
      area[i] = 0;
    else {
      area[i] = j; /* store new label */
      j++;
    }
  }

  /* set/get new labels and assign them to the new label image */

  for (int p = 0; p < label->n; p++) {
    nlabel->val[p] = area[label->val[p]];
  }
  iftCopyVoxelSize(label, nlabel);
  iftFree(area);

  return (nlabel);
}
/**
 * @brief Selects the regions with area size within the given interval.
 *
 * The regions are relabeled automatically in the <output_label> image.
 *
 * @param in_label The input label.
 * @param out_label The output label, which should be pre-allocated.
 * @param min_thres The minimum area size.
 * @param max_thres The maximum area size.
 *
 * @note Private function
 *
 */
void iftSelectRegionsInAreaIntervalAux(iftImage *in_label, iftImage *out_label, int min_thres, int max_thres) {
  iftULongArray *area  = NULL;
  int maxval;

  iftVerifyImageDomains(in_label, out_label, "iftSelectRegionsInAreaIntervalAux");

  /* compute area of each region */

  area = iftCreateULongArray(iftMaximumValue(in_label) + 1);
  for (ulong p=0; p < in_label->n; p++)
    area->val[in_label->val[p]]++;

  /* eliminate regions below area and store the new labels in the area vector */

  area->val[0]=0;
  maxval = iftMaximumValue(in_label);

  for (int i=1, j=1; i <= maxval; i++) {
    if ((area->val[i]>max_thres)||(area->val[i]<min_thres))
      area->val[i]=0;
    else{
      area->val[i]=j; /* store new label */
      j++;
    }
  }

  /* set/get new labels and assign them to the new label image */

  for (int p=0; p < in_label->n; p++){
    out_label->val[p]=area->val[in_label->val[p]];
  }
  iftCopyVoxelSize(in_label, out_label);
  iftDestroyULongArray(&area);
}

iftImage *iftSelectRegionsInAreaInterval(iftImage *label, int min_thres, int max_thres)
{
  iftImage *out_label =iftCopyImage(label);

  iftSelectRegionsInAreaIntervalAux(label, out_label, min_thres, max_thres);

  iftCopyVoxelSize(label, out_label);

  return(out_label);
}

void iftSelectRegionsInAreaIntervalInplace(iftImage *label, int min_thres, int max_thres)
{
  iftSelectRegionsInAreaIntervalAux(label, label, min_thres, max_thres);
}

iftImage *iftSelectKLargestRegions(iftImage *label, int K)
{
  iftImage *nlabel;
  int       ncomps,p,i,*index,*size;

  if (K <= 0)
      iftError("Invalid number of components", "iftSelectKLargestRegions");

  ncomps = iftMaximumValue(label);
  size   = iftAllocIntArray(ncomps+1);
  index  = iftAllocIntArray(ncomps+1);


  for (i=0; i <= ncomps; i++) {
    index[i]=i;
    size[i] =0.0;
  }

  /* Do not consider label 0 --- background */
  for (p=0; p < label->n; p++)
    size[label->val[p]]++;
  size[0]=0;

  iftBucketSort(size, index, ncomps+1, IFT_DECREASING);
  
  iftFree(size);

  nlabel = iftCreateImage(label->xsize,label->ysize,label->zsize);

  for (p=0; p < label->n; p++) {
    for (i=0; i < K; i++) {
      if (label->val[p]==index[i])
	nlabel->val[p]=i+1;
    }
  }

  iftCopyVoxelSize(label,nlabel);
  iftFree(index);

  return(nlabel);
}

iftImage *iftSelectKLargestRegionsAndPropagateTheirLabels(iftImage *label, iftAdjRel *A, int K)
{
  iftImage *nlabel[2];
  int       ncomps,p,i,*index,*size;
  iftSet *S = NULL;

  /* relabel components: it accounts for disconnected labels, which
     should be considered multiple components */

  nlabel[0] = iftRelabelRegions(label,A);

  ncomps = iftMaximumValue(nlabel[0]);
  size   = iftAllocIntArray(ncomps+1);
  index  = iftAllocIntArray(ncomps+1);

  if (K <= 0)
      iftError("Invalid number of components", "iftSelectKLargestRegions");

  if (K > ncomps)
    K = ncomps;

  for (i=0; i <= ncomps; i++) {
    index[i]=i;
    size[i] =0.0;
  }

  /* Do not consider label 0 --- background */
  for (p=0; p < nlabel[0]->n; p++)
    size[nlabel[0]->val[p]]++;
  size[0]=0;

  iftBucketSort(size, index, ncomps+1, IFT_DECREASING);
  iftFree(size);

  nlabel[1] = iftCreateImage(nlabel[0]->xsize,nlabel[0]->ysize,nlabel[0]->zsize);
  for (p=0; p < nlabel[0]->n; p++) {
    for (i=0; i < K; i++) {
      if (nlabel[0]->val[p]==index[i]){
	nlabel[1]->val[p]=i+1;
	iftInsertSet(&S,p);
      }
    }
  }

  iftCopyVoxelSize(label,nlabel[1]);
  iftFree(index);

  /* fill regions with label zero from the labeled voxels in the
     set */

  while(S != NULL) {
    int p      = iftRemoveSet(&S);
    iftVoxel u = iftGetVoxelCoord(nlabel[1],p);
    for (i=1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(nlabel[1],v)){
	int q = iftGetVoxelIndex(nlabel[1],v);
	if (nlabel[1]->val[q]==0){
	  nlabel[1]->val[q]=nlabel[1]->val[p];
	  iftInsertSet(&S,q);
	}
      }
    }
  }
	
  iftDestroySet(&S);
  iftDestroyImage(&nlabel[0]);

  return(nlabel[1]);
}

iftImage *iftSelectKSmallestRegions(iftImage *label, int K)
{
  iftImage *nlabel;
  int       ncomps,p,i,*index,*size;

  if (K <= 0)
      iftError("Invalid number of components", "iftSelectKSmallestRegions");

  ncomps = iftMaximumValue(label);
  size   = iftAllocIntArray(ncomps+1);
  index  = iftAllocIntArray(ncomps+1);

  for (i=0; i <= ncomps; i++) {
    index[i]=i;
    size[i] =0.0;
  }

  /* Do not consider label 0 --- background */
  for (p=0; p < label->n; p++)
    size[label->val[p]]++;
  size[0]=label->n;

  iftBucketSort(size, index, ncomps+1, IFT_INCREASING);
  
  iftFree(size);

  nlabel = iftCreateImage(label->xsize,label->ysize,label->zsize);

  for (p=0; p < label->n; p++) {
    for (i=0; i < K; i++) {
      if (label->val[p]==index[i])
	nlabel->val[p]=i+1;
    }
  }

  iftCopyVoxelSize(label,nlabel);
  iftFree(index);

  return(nlabel);
}

char iftValidArc(iftBMap *bndr, iftImage *pred, iftImage *bin, iftAdjRel *A, iftAdjRel *L, iftAdjRel *R, iftVoxel u, int i, int q)
{
  iftVoxel w; w.z = 0;
  int left, right;

  w.x = u.x + L->dx[i];
  w.y = u.y + L->dy[i];

  if (iftValidVoxel(bin,w))
    left = iftGetVoxelIndex(bin,w);
  else
    left = -1;

  w.x = u.x + R->dx[i];
  w.y = u.y + R->dy[i];
  
  if (iftValidVoxel(bin,w))
    right = iftGetVoxelIndex(bin,w);
  else
    right = -1;

  
  if (iftBMapValue(bndr,q) && (pred->val[q] == IFT_NIL) && (left != -1)){
    /* adjacent pixel on the boundary with a valid left pixel */
    if (bin->val[left]!=0){
      if (right != -1){
	/* the right pixel is valid and in the exterior */	
	if (bin->val[right]==0){
	  return(1);
	}else{
	  return(0);
	}
      } else {
	/* the object touches the image's border */	
	return(1);
      }
    } else {
      return(0);
    }
  }
  
  return(0);
}

char iftValidStartingPixel(iftBMap *bndr, iftImage *pred, iftImage *bin, iftAdjRel *A, iftAdjRel *L, iftAdjRel *R, int p, int *nvalidarcs)
{

  if ((*nvalidarcs) == 2)
    return(1);
    
  if (iftBMapValue(bndr,p) && (pred->val[p] == IFT_NIL)){
    iftVoxel u = iftGetVoxelCoord(bin,p);
    for (int i=1;(i < A->n); i++) {
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(bin,v)){
	int q        = iftGetVoxelIndex(bin,v);
	if (iftValidArc(bndr,pred,bin,A,L,R,u,i,q)){
	  (*nvalidarcs) += 1;
	  if (iftValidStartingPixel(bndr, pred, bin, A, L, R, q, nvalidarcs)){
	    return(1);
	  }
	}
      }
    }
  }

  return(0);
}

iftImage *iftLabelContPixel(iftImage *bin)
{
  iftImage *pred=NULL,*label=NULL;
  int p=0,q,r,i,l;
  iftLIFO *LIFO=NULL;
  iftAdjRel *A,*L,*R,*A4;
  iftVoxel   u,v;
  iftSet    *B=NULL,*Baux=NULL;
  iftBMap   *bndr;
  
  if (bin->zsize != 1)
      iftError("Input image must be a 2D image", "iftLabelContPixel");

  u.z=v.z=0;

  A4     = iftCircular(1.0);
  B      = iftObjectBorderSet(bin, A4);
  bndr   = iftCreateBMap(bin->n);
  Baux   = B;
  while (Baux != NULL) {
    p = Baux->elem;
    iftBMapSet1(bndr,p);
    Baux = Baux->next;
  }
  
  A      = iftClockCircular(sqrtf(2.0));
  L      = iftLeftSide(A,1.0);
  R      = iftRightSide(A,1.0);
  
  label  = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  pred   = iftCreateImage(bin->xsize,bin->ysize,bin->zsize);
  iftSetImage(pred, IFT_NIL);

  LIFO   = iftCreateLIFO(bin->n);
  
  Baux = B;
  while (Baux != NULL) {
    r = Baux->elem;
    int nvalidarcs=0;
    if (iftValidStartingPixel(bndr,pred, bin, A, L, R, r, &nvalidarcs)){
      iftInsertLIFO(LIFO,r);
      pred->val[r] = r;
      
      while(!iftEmptyLIFO(LIFO)){

	p   = iftRemoveLIFO(LIFO);
	u   = iftGetVoxelCoord(bin,p);
	
	for (i=1; i < A->n; i++){
	    
	  v = iftGetAdjacentVoxel(A,u,i);

	  if (iftValidVoxel(bin,v)){
	    q = iftGetVoxelIndex(bin,v);

	    if ((q==r)&&(pred->val[p]!=r)){ /*A pixel p is the last
					      visited in the contour
					      when its neighbor is r,
					      but the precessor of p
					      is not r */ 

	      if (pred->val[pred->val[p]]!=r){ /* if it is not
	      	      				  a silly loop
	      	      				  r=(x,y),p=(x+1,y),
	      	      				  and pred[p]=(x,y+1),
	      	      				  then label the
	      	      				  contour */
		
		l = 1;
		while(pred->val[p]!=p){ /* it will finish the labeling at r */
		  label->val[p] = l;
		  p    = pred->val[p];
		  l++;
		}
		label->val[p] = l;
	      }
	      iftResetLIFO(LIFO);
	      break;
	    }
	    
	    if (iftValidArc(bndr,pred,bin,A,L,R,u,i,q)){
	      pred->val[q] = p;
	      iftInsertLIFO(LIFO,q);
	    }
	  }
	}
      }
    }
    Baux = Baux->next;
  }

  iftDestroySet(&B);
  iftDestroyBMap(&bndr);
  iftDestroyAdjRel(&A);
  iftDestroyAdjRel(&A4);
  iftDestroyAdjRel(&L);
  iftDestroyAdjRel(&R);
  iftDestroyImage(&pred);
  iftDestroyLIFO(&LIFO);

  iftCopyVoxelSize(bin,label);

  return(label);
}




iftSet *iftEndPoints(iftImage *skel, iftAdjRel *A)
{
  int p,q,i,counter;
  iftVoxel u,v;
  iftSet *S=NULL;

  for (p=0; p < skel->n; p++) {
    if (skel->val[p]!=0){
      u.x = iftGetXCoord(skel,p);
      u.y = iftGetYCoord(skel,p);
      u.z = iftGetZCoord(skel,p);
      counter=0;
      for (i=1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	v.z = u.z + A->dz[i];
	if (iftValidVoxel(skel,v)){
	  q = iftGetVoxelIndex(skel,v);
	  if (skel->val[q]!=0){
	    counter++;
	  }
	}
      }
      if (counter==1){ // terminal point
	iftInsertSet(&S,p);
      }
    }
  }
  return(S);
}

iftSet *iftSkeletonPoints(iftImage *skel)
{
  iftSet *S=NULL;
  int p;

  for (p=0; p < skel->n; p++)
    if (skel->val[p]!=0)
      iftInsertSet(&S,p);

  return(S);
}

iftSet *iftFindPathOnSkeleton(iftImage *skel, iftAdjRel *A, int src, int dst)
{
  int i,p,q,*pred;
  iftVoxel u,v;
  iftSet  *path=NULL;
  iftFIFO *F=NULL;

  if ((skel->val[src]==0)||(skel->val[dst]==0))
      iftError("End points must be skeleton points", "iftFindPathOnSkeleton");

  F      = iftCreateFIFO(skel->n);
  pred   = iftAllocIntArray(skel->n);
  iftInsertFIFO(F,src);
  pred[src]= IFT_NIL;

  while(!iftEmptyFIFO(F)){
    p = iftRemoveFIFO(F);

    if (p == dst){
      while(p != IFT_NIL){
	u.x = iftGetXCoord(skel,p);
	u.y = iftGetYCoord(skel,p);
	u.z = iftGetZCoord(skel,p);
	iftInsertSet(&path,p);
	p = pred[p];
      }
      iftDestroyFIFO(&F);
      iftFree(pred);
      return(path);
    }

    u.x = iftGetXCoord(skel,p);
    u.y = iftGetYCoord(skel,p);
    u.z = iftGetZCoord(skel,p);

    for (i=1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      v.z = u.z + A->dz[i];
      if (iftValidVoxel(skel,v)){
	q = iftGetVoxelIndex(skel,v);
	if ((skel->val[q]!=0)&&(F->color[q] == IFT_WHITE)){
	  pred[q]=p;
	  iftInsertFIFO(F,q);
	}
      }
    }

  }

  iftDestroyFIFO(&F);
  iftFree(pred);
    iftError("Path was not found", "iftFindPathOnSkeleton");

  return(NULL);
}

iftLabeledSet *iftFuzzyModelToLabeledSet(iftImage *model)
{
  iftLabeledSet *S=NULL;
  int p;

  for (p=0; p < model->n; p++){
    if (model->val[p]==255) // interior
      iftInsertLabeledSet(&S,p,1);
    else
      if (model->val[p]==0) // exterior
	iftInsertLabeledSet(&S,p,0);
  }
  return(S);
}

void iftWriteSeeds2D(char* filename, iftLabeledSet* seed, iftImage* image){
	FILE *file = fopen(filename,"w");
	if(file == NULL)
        iftError("Invalid destination file", "iftWriteSeeds2D");

	iftLabeledSet *s = seed;
	int nseeds = 0;
	while(s != NULL){
		nseeds++;
		s = s->next;
	}

	fprintf(file,"%d %d %d\n", nseeds, image->xsize, image->ysize);

	s = seed;
	while(s != NULL){
		iftVoxel voxel = iftGetVoxelCoord(image,s->elem);

		fprintf(file, "%d %d 0 %d\n", voxel.x, voxel.y, s->label);
		s = s->next;
	}

	fclose(file);

}

iftImage* iftSeedImageFromLabeledSet(iftLabeledSet* labeled_set, iftImage *image){
	iftImage* seed_image = iftCreateImage(image->xsize,image->ysize,image->zsize);
	iftSetImage(seed_image, -1);

	iftWriteSeedsOnImage(seed_image,labeled_set);

	return seed_image;
}


void iftWriteSeedsOnImage(iftImage* image, iftLabeledSet* seed){
	iftLabeledSet* i = seed;

	while(i != NULL){
		image->val[i->elem] = i->label;
		i = i->next;
	}

}

iftLabeledSet *iftLabeledSetFromSeedImage(iftImage* seed_image, bool decrement){
	iftLabeledSet *seed = NULL;

	int p, label;
	if (decrement){
	  for(p = 0; p < seed_image->n; p++)
	    {
	      label = seed_image->val[p] - 1;
	      if(label >= 0)
		iftInsertLabeledSet(&seed, p, label);
	    }
	} else {
	  for(p = 0; p < seed_image->n; p++)
	    {
	      label = seed_image->val[p];
	      if(label > 0)
		iftInsertLabeledSet(&seed, p, label);
	    }
	}
	return seed;
}

iftLabeledSet *iftLabeledSetFromSeedImageMarkersAndHandicap(iftImage* seed_image, iftImage *marker, iftImage *handicap){
  /* seed_image needs to be initialized with -1, since 0 is a possible label */
  iftLabeledSet *seed = NULL;

  int p, label;
  for(p = 0; p < seed_image->n; p++)
  {
    label = seed_image->val[p];
    if(label != -1)
      iftInsertLabeledSetMarkerAndHandicap(&seed, p, label, marker->val[p], handicap->val[p]);
  }

  return seed;
}


/* Gt image needs to be: 0, 1, 2, ..., N */
iftLabeledSet *iftBorderMarkersForPixelSegmentation(iftImage *grad_image, iftImage *gt_image, float border_distance)
{
  iftAdjRel *A1,*A2;
  iftSet    *B=NULL;
  iftImage  *dist=NULL,*root=NULL,*bin=NULL;
  iftGQueue *Q=NULL;
  int        i,p,q,tmp,*grad = NULL;
  iftVoxel   u,v,r;
  float      dist_thres=border_distance*border_distance;
  iftLabeledSet *S=NULL;

  /* Initialization */

  if (iftIs3DImage(grad_image)){
    A1 = iftSpheric(1.0);
    A2 = iftSpheric(sqrtf(3.0));
  } else {
    A1 = iftCircular(1.0);
    A2 = iftCircular(1.5);
  }

  B = iftObjectBorderSet(gt_image, A1);

  dist  = iftCreateImage(gt_image->xsize,gt_image->ysize,gt_image->zsize);
  iftSetImage(dist, IFT_INFINITY_INT);
  root  = iftCreateImage(gt_image->xsize,gt_image->ysize,gt_image->zsize);
  bin   = iftCreateImage(gt_image->xsize,gt_image->ysize,gt_image->zsize);
  Q     = iftCreateGQueue(IFT_QSIZE,gt_image->n,dist->val);

  while (B != NULL) {
    p = iftRemoveSet(&B);
    dist->val[p]=0;
    root->val[p]=p;
    iftInsertGQueue(&Q,p);
  }

  /* Image Foresting Transform: Dilation of the border upto dist_thres */

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    if (dist->val[p] > dist_thres)
      break;
    else{
      bin->val[p]=1;
    }

    u = iftGetVoxelCoord(bin,p);
    r = iftGetVoxelCoord(bin,root->val[p]);

    for (i=1; i < A2->n; i++){
      v = iftGetAdjacentVoxel(A2,u,i);
      if (iftValidVoxel(bin,v)){
	q = iftGetVoxelIndex(bin,v);
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

  iftDestroyGQueue(&Q);
  iftDestroyImage(&dist);

  /* Compute markers in the decreasing order of gradient values (due to
     the LIFO nature of iftLabeledSet) */

  // This array holds the seed's root gradient value for sorting
  grad  = iftAllocIntArray(grad_image->n);
  Q     = iftCreateGQueue(IFT_QSIZE,gt_image->n, grad);
  iftSetRemovalPolicy(Q,MAXVALUE);

  B = iftObjectBorderSet(bin, A1);
  iftDestroyImage(&bin);

  iftDestroyAdjRel(&A2);

  if (iftIs3DImage(grad_image)){
    A2 = iftSpheric(1.0);
  }else
  {
    A2 = iftCircular(1.0);
  }

  while (B != NULL) {
    p = iftRemoveSet(&B);
    // Since the ground truth's border might be a little bit displaced, we compute the median gradient value around
    // the root of the voxel as the voxel's weight, aiming to minimize it along the boundary at each iteration
    // of the robot
    grad[p]=iftMedianValueInAdjacency(grad_image, root->val[p], A2);
    iftInsertGQueue(&Q,p);
  }

  iftDestroyImage(&root);

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);
    //XXX
    //if (gt_image->val[p]==0)
      //iftInsertLabeledSet(&S,p,1);
    //else
      //iftInsertLabeledSet(&S,p,2);
    iftInsertLabeledSet(&S,p,gt_image->val[p]);
  }

  iftDestroyGQueue(&Q);
  iftFree(grad);

  iftDestroyAdjRel(&A2);
  iftDestroyAdjRel(&A1);

  return(S);
}

iftImage *iftRelabelSegmentationErrorComponents(iftImage *gt_image, iftImage *label, iftAdjRel *adj_relabeling) {
    int p;

    iftImage* error_image = iftCreateImage(gt_image->xsize, gt_image->ysize, gt_image->zsize);
    iftImage *relabelled = NULL;

    // Computing error components. If the label image is not passed, then we
    // return an error component image for the entire gt
    if(label == NULL){
        for(p = 0; p < gt_image->n; p++){
            error_image->val[p] = gt_image->val[p]+1;
        }
    }else{
        for(p = 0; p < gt_image->n; p++){
            if (gt_image->val[p] != label->val[p])
            error_image->val[p] = gt_image->val[p]+1;
        }
    }

    relabelled = iftRelabelRegions(error_image, adj_relabeling);

    iftDestroyImage(&error_image);

    iftCopyVoxelSize(gt_image, relabelled);

    return relabelled;
}

/* Gt image needs to be: 0, 1, 2, ..., N */
iftLabeledSet *iftGeodesicMarkersForSegmentation(iftImage *gt_image, iftImage *label){
    iftAdjRel *adj_relabeling = NULL;
    iftImage *error_image;
    iftLabeledSet *geo_centers = NULL, *c = NULL;

    if (iftIs3DImage(gt_image))
        adj_relabeling = iftSpheric(1.0);
    else
        adj_relabeling = iftCircular(1.0);

    /////////////////////////////////////////////////////////////////////////////////////
//    error_image = iftCreateImage(gt_image->xsize, gt_image->ysize, gt_image->zsize);
//    int p;
//    if(label == NULL){
//        for(p = 0; p < gt_image->n; p++){
//            error_image->val[p] = gt_image->val[p] + 1;
//        }
//    }else{
//        for(p = 0; p < gt_image->n; p++){
//            if (gt_image->val[p] != label->val[p])
//                error_image->val[p] = gt_image->val[p] + 1;
//        }
//    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    error_image = iftRelabelSegmentationErrorComponents(gt_image, label, adj_relabeling);

    iftWriteImageByExt(error_image, "tmp/relabelled.pgm");
    geo_centers = iftGeodesicCenters(error_image);

    //////////////////////////////////////////////////////////////////
//    for(c = geo_centers; c != NULL; c = c->next){
//        printf("p %d, n %d\n", c->elem, gt_image->n);
//    }
    //////////////////////////////////////////////////////////////////
    c = geo_centers;

    //Redefines the labels for the geometric centers, since the image was rellabeled
    while(c != NULL){
        c->label = gt_image->val[c->elem];
        c = c->next;
    }

    iftDestroyImage(&error_image);
    iftDestroyAdjRel(&adj_relabeling);
    return geo_centers;
}


/* Gt image needs to be: 0, 1, 2, ..., N
    Not tested!! */
iftLabeledSet* iftBorderMarkersForSuperpixelSegmentation(iftImage* label_image,iftImage* gt_image, iftDataSet* dataset){

  iftVerifyImageDomains(label_image,gt_image,"iftBorderMarkersForSuperpixelSegmentation");

	iftAdjRel *adj_rel;
	if(iftIs3DImage(label_image))
		adj_rel = iftSpheric(1.0);
	else
		adj_rel = iftCircular(1.0);

	iftRegionGraph* region_graph = iftRegionGraphFromLabelImage(label_image,dataset, adj_rel);

	iftDestroyAdjRel(&adj_rel);

	//Counts how many pixels are from each class in a region
	int *n_obj_pixels = iftAllocIntArray(region_graph->nnodes);
	int *n_bkg_pixels = iftAllocIntArray(region_graph->nnodes);

	int p;
	for(p = 0; p < gt_image->n; p++){
		int region = label_image->val[p] - 1;

		if(gt_image->val[p] != 0)
			n_obj_pixels[region]++;
		else
			n_bkg_pixels[region]++;
	}

	int *region_gt = iftAllocIntArray(region_graph->nnodes);

	int r;
	for(r = 0; r < region_graph->nnodes; r++){
		if(n_obj_pixels[r] > n_bkg_pixels[r])
			region_gt[r] = 1;
		else
			region_gt[r] = 0;
	}

	iftFree(n_obj_pixels);
	iftFree(n_bkg_pixels);

	//Computes the number of borders between object and background
	int n_borders = 0;
	for(r = 0; r < region_graph->nnodes; r++){
	    iftSet *adj = region_graph->node[r].adjacent;
		while(adj != NULL){
			int v = adj->elem;
			if(region_gt[r] > region_gt[v]){ //From object to background
				n_borders++;
			}
			adj = adj->next;
		}
	}

	//Computes the cost between object and background
	int src_border[n_borders];
	int dst_border[n_borders];
	int border_index[n_borders];
	float border_distance[n_borders];

	int i;
	for(i = 0; i < n_borders; i++)
		border_index[i] = i;

	i = 0;
	for(r = 0; r < region_graph->nnodes; r++){
		iftSet *adj = region_graph->node[r].adjacent;
		while(adj != NULL){
			int v = adj->elem;
			if(region_gt[r] > region_gt[v]){ //From object to background
				float arcw = dataset->iftArcWeight(dataset->sample[r].feat,dataset->sample[v].feat, dataset->alpha,dataset->nfeats);

				src_border[i] = r;
				dst_border[i] = v;
				border_distance[i] = arcw;
				i++;
			}
			adj = adj->next;
		}
	}

	iftFree(region_gt);

	//Sorts the borders by decreasing distance, because iftLabeledSet is a LIFO structure
	iftFQuickSort(border_distance, border_index, 0,n_borders -1, IFT_DECREASING);

	int *centers = iftAllocIntArray(region_graph->nnodes);
	for(i = 0; i < region_graph->nnodes; i++)
		centers[i] = -1;

	iftDestroyRegionGraph(&region_graph);

  //iftImage * label_image_plus_one = iftCopyImage(label_image);
  //iftLabeledSet *geo_centers = iftGeodesicCenters(label_image_plus_one);
  //iftDestroyImage(&label_image_plus_one);
	iftLabeledSet *geo_centers = iftGeodesicCenters(label_image);
	iftLabeledSet *c = geo_centers;
	while(c != NULL){
		centers[c->label - 1] = c->elem;
		c = c->next;
	}
	iftDestroyLabeledSet(&geo_centers);

	iftLabeledSet* markers = NULL;
	//Select the markers
	for(i = 0; i < n_borders; i++){
		int border = border_index[i];

		int r_obj = src_border[border];
		int r_bkg = dst_border[border];

    //XXX
		iftInsertLabeledSet(&markers,centers[r_obj],1);
		iftInsertLabeledSet(&markers,centers[r_bkg],0);

	}

	iftFree(centers);

	return markers;
}

iftLabeledSet* iftGetSeeds(iftLabeledSet* S, int nelem, int label){
	iftLabeledSet* seeds = NULL;

	while(S != NULL && nelem > 0){
		if(S->label == label){
			iftInsertLabeledSet(&seeds,S->elem,S->label);
			nelem--;
		}

		S = S->next;
	}

	return seeds;
}

/* Gt image needs to be: 0, 1, 2, ..., N */
iftLabeledSet* iftGetMisclassifiedSeeds(iftLabeledSet* S, int nelem, int label, iftImage* gt_image, iftImage* cl_image){
	iftLabeledSet* seeds = NULL;

	while(S != NULL && nelem > 0){
		if(S->label == label){
			int p = S->elem;
			//If it is misclassified
      //XXX
			//if( (cl_image == NULL) || ((gt_image->val[p] == 0) && (cl_image->val[p] == 2) )
			//		|| ( (gt_image->val[p] != 0) && (cl_image->val[p] == 1) )){
      if( (cl_image == NULL) || (gt_image->val[p] != cl_image->val[p]) ){
				iftInsertLabeledSet(&seeds,S->elem,S->label);
				nelem--;
			}
		}

		S = S->next;
	}

	return seeds;
}

int iftCheckNewSeeds(int *nelem, int length)
{
    for (int i=0; i<length; i++)
        if (nelem[i] > 0)
            return 1;
    return 0;
}

// As of now, the iftBMap is redundant, but kept for future works.
int iftMarkersFromMisclassifiedSeeds(iftImage *seed_image, iftLabeledSet *all_seeds, iftBMap *used_seeds, int nseeds,
                                     int number_of_labels, iftImage *gt, iftImage *label, int dist_border,
                                     int max_marker_radius, int min_marker_radius) {
    //int nelem_1 = nseeds/2;
    //int nelem_2 = nelem_1;
    int i=0;
    int nseeds_per_object = (int) nseeds / number_of_labels;
    int *nelem = iftAllocIntArray(number_of_labels);
    iftLabeledSet *S = all_seeds;
    int total_seeds = 0;
    iftAdjRel *distance_border = NULL;

    for(i=0; i < number_of_labels; i++)
        nelem[i] = nseeds_per_object;

    if(iftIs3DImage(gt))
        distance_border = iftSpheric((float)(dist_border + max_marker_radius));
    else
        distance_border = iftCircular((float)(dist_border + max_marker_radius));

    while((S != NULL) && iftCheckNewSeeds(nelem, number_of_labels)) {
        int p = S->elem;

        int misclassified = ( (label == NULL) || (gt->val[p] != label->val[p]));

        int needs_seed = ( nelem[S->label] > 0 );  //If this object still need seed

        if ( (seed_image->val[p] < 0) && misclassified && needs_seed) {
            int closest_distance = IFT_INFINITY_INT;

            iftVoxel u = iftGetVoxelCoord(gt,p);

            // Calculating the minimum distance to the object's border in the ground truth image to determine the marker's
            // radius
            for(i = 1; i < distance_border->n; i++){
                iftVoxel v = iftGetAdjacentVoxel(distance_border,u,i);
                if(iftValidVoxel(gt,v)){
                    int q = iftGetVoxelIndex(gt,v);

                    if((gt->val[p] != gt->val[q]) && (closest_distance > iftVoxelDistance(u,v))){
                        closest_distance = floorf(iftVoxelDistance(u,v));
                    }
                }
            }

            // Computing the maximum possible marker radius
            closest_distance = iftMin(closest_distance - dist_border, max_marker_radius);

            // If the marker's radius is less than the acceptable minimum size, then we set the seed image with the
            // ground truth's label
            if(closest_distance >= min_marker_radius) {
                iftAdjRel *marker_size;
                if(iftIs3DImage(gt))
                    marker_size = iftSpheric(closest_distance);
                else
                    marker_size = iftCircular(closest_distance);

                for(i = 0; i < marker_size->n; i++){
                    iftVoxel v = iftGetAdjacentVoxel(marker_size,u,i);

                    if(iftValidVoxel(gt,v)){
                        int q = iftGetVoxelIndex(gt,v);

                        seed_image->val[q] = gt->val[q];
                    }
                }

                iftDestroyAdjRel(&marker_size);

                nelem[gt->val[p]]--;
                total_seeds++;

                iftBMapSet1(used_seeds,p);
            }
        }

        S = S->next;
    }

    iftDestroyAdjRel(&distance_border);
    iftFree(nelem);
    return total_seeds;

}

/* Seeds in the input image coordinate system are adjusted to the
   output image coordinate system (due to crop along convolutional
   neural network) */

iftLabeledSet *iftMAdjustSeedCoordinates(iftLabeledSet *Sin, iftMImage *input, iftMImage *output)
{
  int dx, dy, dz;
  iftLabeledSet *Sout=NULL, *S;

  dx = (output->xsize - input->xsize)/2;
  dy = (output->ysize - input->ysize)/2;
  dz = (output->zsize - input->zsize)/2;

  S = Sin;
  while(S != NULL) {
    int p = S->elem;
    iftVoxel u = iftMGetVoxelCoord(input,p);
    u.x += dx; u.y += dy; u.z += dz;
    if (iftMValidVoxel(output,u)){
      p = iftMGetVoxelIndex(output,u);
      iftInsertLabeledSet(&Sout,p,S->label);
    }
    S = S->next;
  }

  return(Sout);
}

/* Seeds in the input image coordinate system are adjusted to the
   output image coordinate system (due to crop along convolutional
   neural network) */

iftLabeledSet *iftAdjustSeedCoordinates(iftLabeledSet *Sin, iftImage *orig, iftMImage *output)
{
  int dx, dy, dz;
  iftLabeledSet *Sout=NULL, *S;

  dx = (output->xsize - orig->xsize)/2;
  dy = (output->ysize - orig->ysize)/2;
  dz = (output->zsize - orig->zsize)/2;

  S = Sin;
  while(S != NULL) {
    int p = S->elem;
    iftVoxel u = iftGetVoxelCoord(orig,p);
    u.x += dx; u.y += dy; u.z += dz;
    if (iftMValidVoxel(output,u)){
      p = iftMGetVoxelIndex(output,u);
      iftInsertLabeledSet(&Sout,p,S->label);
    }
    S = S->next;
  }

  return(Sout);
}

iftImage *iftRegionalMaxima(iftImage *img)
{
  iftImage   *pred    = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage   *pathval = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage   *regmax  = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftGQueue  *Q=NULL;
  int         i,p,q,tmp;
  iftVoxel    u,v;
  iftAdjRel  *A;

  if (iftIs3DImage(img))
    A = iftSpheric(1.0);
  else
    A = iftCircular(1.0);

  Q     = iftCreateGQueue(iftMaximumValue(img)+2,img->n,pathval->val);
  iftSetRemovalPolicy(Q,MAXVALUE);

  for (p=0; p < img->n; p++) {
    img->val[p]++;
    pred->val[p]    = IFT_NIL;
    pathval->val[p] = img->val[p]-1;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(img,p);
    
    if (pred->val[p] == IFT_NIL) {
      regmax->val[p]  = 255;
      pathval->val[p]++;// To output one seed per maximum

    } 
    img->val[p]--;

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){	
	q = iftGetVoxelIndex(img,v);
	if (pathval->val[q] < pathval->val[p]){
	  tmp = iftMin(pathval->val[p], img->val[q]);

	  if (tmp > pathval->val[q]){
	    iftRemoveGQueueElem(Q,q);
	    pathval->val[q] = tmp;
	    pred->val[q]    = p; 
	    iftInsertGQueue(&Q,q);
	  }
	}
      }
    }
  }
  
  iftDestroyImage(&pred);
  iftDestroyImage(&pathval);
  iftDestroyGQueue(&Q);
  iftDestroyAdjRel(&A);

  return(regmax);
}

iftImage *iftRegionalMaximaInRegion(iftImage *img, iftImage *mask)
{
  iftImage   *pred    = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage   *pathval = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage   *regmax  = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftGQueue  *Q=NULL;
  int         i,p,q,tmp;
  iftVoxel    u,v;
  iftAdjRel  *A;

  if (iftIs3DImage(img))
    A = iftSpheric(1.0);
  else
    A = iftCircular(1.0);

  Q     = iftCreateGQueue(iftMaximumValue(img)+2,img->n,pathval->val);
  iftSetRemovalPolicy(Q,MAXVALUE);

  for (p=0; p < img->n; p++) {

    img->val[p]++;

    if (mask->val[p]){
      pred->val[p]    = IFT_NIL;
      pathval->val[p] = img->val[p]-1;
      iftInsertGQueue(&Q,p);
    } 
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(img,p);

    img->val[p]--;

    if (pred->val[p] == IFT_NIL) {
      regmax->val[p]  = 255;
      pathval->val[p]++;// To output one seed per maximum
    } 

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){	
	q = iftGetVoxelIndex(img,v);
	if ((mask->val[q])&&(pathval->val[q] < pathval->val[p])){
	  tmp = iftMin(pathval->val[p], img->val[q]);

	  if (tmp > pathval->val[q]){
	    iftRemoveGQueueElem(Q,q);
	    pathval->val[q] = tmp;
	    pred->val[q]    = p; 
	    iftInsertGQueue(&Q,q);
	  }
	}
      }
    }
  }
  
  iftDestroyImage(&pred);
  iftDestroyImage(&pathval);
  iftDestroyGQueue(&Q);
  iftDestroyAdjRel(&A);

  return(regmax);
}

iftImage *iftRegionalMinima(iftImage *img)
{
  iftImage   *pred    = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage   *pathval = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage   *regmin  = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftGQueue  *Q=NULL;
  int         i,p,q,tmp;
  iftVoxel    u,v;
  iftAdjRel  *A;

  if (iftIs3DImage(img))
    A = iftSpheric(1.0);
  else
    A = iftCircular(1.0);

  Q     = iftCreateGQueue(iftMaximumValue(img)+2,img->n,pathval->val);

  for (p=0; p < img->n; p++) {
    pred->val[p]    = IFT_NIL;
    pathval->val[p] = img->val[p]+1;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(img,p);
    
    if (pred->val[p] == IFT_NIL) {
      regmin->val[p]  = 255;
      pathval->val[p]--;// To output one seed per minimum
    } 

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){	
	q = iftGetVoxelIndex(img,v);
	if (pathval->val[q] > pathval->val[p]){
	  tmp = iftMax(pathval->val[p], img->val[q]);

	  if (tmp < pathval->val[q]){
	    iftRemoveGQueueElem(Q,q);
	    pathval->val[q] = tmp;
	    pred->val[q]    = p; 
	    iftInsertGQueue(&Q,q);
	  }
	}
      }
    }
  }
  
  iftDestroyImage(&pred);
  iftDestroyImage(&pathval);
  iftDestroyGQueue(&Q);
  iftDestroyAdjRel(&A);

  return(regmin);
}

iftImage *iftRegionalMinimaInRegion(iftImage *img, iftImage *mask)
{
  iftImage   *pred    = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage   *pathval = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage   *regmin  = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftGQueue  *Q=NULL;
  int         i,p,q,tmp;
  iftVoxel    u,v;
  iftAdjRel  *A;

  if (iftIs3DImage(img))
    A = iftSpheric(1.0);
  else
    A = iftCircular(1.0);

  Q     = iftCreateGQueue(iftMaximumValue(img)+2,img->n,pathval->val);

  for (p=0; p < img->n; p++) {
    if (mask->val[p]){
      pred->val[p]    = IFT_NIL;
      pathval->val[p] = img->val[p]+1;
      iftInsertGQueue(&Q,p);
    } 
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(img,p);
    
    if (pred->val[p] == IFT_NIL) {
      regmin->val[p]  = 255;
      pathval->val[p]--;// To output one seed per minimum
    } 

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){	
	q = iftGetVoxelIndex(img,v);
	if ((mask->val[q])&&(pathval->val[q] > pathval->val[p])){
	  tmp = iftMax(pathval->val[p], img->val[q]);

	  if (tmp < pathval->val[q]){
	    iftRemoveGQueueElem(Q,q);
	    pathval->val[q] = tmp;
	    pred->val[q]    = p; 
	    iftInsertGQueue(&Q,q);
	  }
	}
      }
    }
  }
  
  iftDestroyImage(&pred);
  iftDestroyImage(&pathval);
  iftDestroyGQueue(&Q);
  iftDestroyAdjRel(&A);

  return(regmin);
}



iftSet *iftBinaryMaskToSet(iftImage *mask) {
  int p;
  iftSet *S = NULL;

  for(p = 0; p < mask->n; p++) {
    if(mask->val[p] != 0)
      iftInsertSet(&S, p);
  }

  return S;
}


iftImage *iftRootVoxels(iftImage *pred)
{
 iftImage *bin = iftCreateImage(pred->xsize,pred->ysize,pred->zsize);
 int      p;

 for (p=0; p < pred->n; p++) {
   if (pred->val[p] == IFT_NIL){
     bin->val[p]=255;
   }
 }

 return(bin);
}

iftImage *iftLeafVoxels(iftImage *pred, iftAdjRel *A)
{
  iftImage *bin = iftCreateImage(pred->xsize,pred->ysize,pred->zsize);
  int       p,q,i;
  iftVoxel  u, v;
  char      is_leaf;

  for (p=0; p < pred->n; p++) {
    if (pred->val[p]>0){
      u       = iftGetVoxelCoord(pred,p);
      is_leaf = 1;
      for (i=1; i < A->n; i++) {
    	v = iftGetAdjacentVoxel(A,u,i);
    	if (iftValidVoxel(pred,v)){
    	  q = iftGetVoxelIndex(pred,v);
    	  if (pred->val[q]==p){
    	    is_leaf=0;
    	    break;
    	  }
    	}
      }
      if (is_leaf)
    	bin->val[p]=255;
    }
  }

  iftCopyVoxelSize(pred,bin);

  return(bin);
}

iftImage *iftLeafVoxelsOnMask(iftImage *pred, iftAdjRel *A, iftImage *mask)
{
    iftImage *bin = iftCreateImage(pred->xsize,pred->ysize,pred->zsize);
    int       p,q,i;
    iftVoxel  u, v;
    char      is_leaf;

    for (p=0; p < pred->n; p++) {
        if ((mask->val[p] > 0) && (pred->val[p]>0)){
            u       = iftGetVoxelCoord(pred,p);
            is_leaf = 1;
            for (i=1; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A,u,i);
                if (iftValidVoxel(pred,v)){
                    q = iftGetVoxelIndex(pred,v);
                    if (pred->val[q]==p){
                        is_leaf=0;
                        break;
                    }
                }
            }
            if (is_leaf)
                bin->val[p]=255;
        }
    }

    iftCopyVoxelSize(pred,bin);

    return(bin);
}

iftImage *iftHBasins(iftImage *img, int H)
{
  iftImage *marker,*suprec,*hbasins;

  if (H <= 0)
      iftError("Invalid depth for basins", "iftHBasins");

  marker   = iftAddValue(img,H); 
  suprec   = iftSuperiorRec(img,marker,NULL);
  hbasins  = iftSub(suprec,img);
  
  iftDestroyImage(&suprec);
  iftDestroyImage(&marker);

  return(hbasins);
}

iftImage *iftHDomes(iftImage *img, int H)
{
  iftImage *hdomes,*aux;

  if (H <= 0)
      iftError("Invalid height for domes", "iftHDomes");

  aux    = iftComplement(img);
  hdomes = iftHBasins(aux,H);
  iftDestroyImage(&aux);

  return(hdomes);
}

iftIntArray *iftGridSamplingOnMask(const iftImage *bin_mask, float radius,
                                   int initial_obj_voxel_idx, long n_samples) {
    int first_obj_voxel = initial_obj_voxel_idx;

    if (initial_obj_voxel_idx >= 0) {
      if (bin_mask->val[initial_obj_voxel_idx] == 0) {
        iftError("Initial Voxel Index %d is not an object voxel",
                 "iftGridSamplingOnMask", initial_obj_voxel_idx);
      }
    }
    else {
      // finds the first object voxel from the binary mask
      int p = 0;
      for (p = 0; p < bin_mask->n && bin_mask->val[p] == 0; p++) {}
      first_obj_voxel = p;
    }
                                     

    if (iftAlmostZero(radius)) {
      iftIntArray *grid_chosen = iftCreateIntArray(1);
      grid_chosen->val[0] = first_obj_voxel;
      return grid_chosen;
    }

    iftAdjRel *A = _iftFindBestAdjRelForGridSampling(radius, iftIs3DImage(bin_mask));
    iftFloatArray *dist = iftCreateFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dist->val[i] = sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i]);
    
    iftImage *prior = iftCreateImageFromImage(bin_mask);
    iftImage *label_img = iftCreateImageFromImage(bin_mask);
    iftGQueue *Q = iftCreateGQueue(IFT_QSIZE, prior->n, prior->val);
    iftSetRemovalPolicy(Q, MAXVALUE);
    
    prior->val[first_obj_voxel] = 1;
    iftInsertGQueue(&Q, first_obj_voxel);

    iftList *grid = iftCreateList();

    int label = 1;

    while (!iftEmptyGQueue(Q)) {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftGetVoxelCoord(bin_mask, p);
        iftInsertListIntoTail(grid, p);

        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(bin_mask, v) && (iftImgVoxelVal(bin_mask, v) != 0) && (iftImgVoxelVal(label_img, v) == 0)) {
                int q = iftGetVoxelIndex(bin_mask, v);

                // q is inside the sphere
                if (dist->val[i] < radius) {
                    label_img->val[q] = label;
                    
                    // voxel was border of another ball
                    if (Q->L.elem[q].color == IFT_GRAY) {
                        iftRemoveGQueueElem(Q, q);
                        prior->val[q] = 0;
                    }
                }
                // q is on the border of the ball with center p
                else {
                    if (Q->L.elem[q].color == IFT_WHITE) {
                        prior->val[q]++;
                        iftInsertGQueue(&Q, q);
                    }
                    // q is on the border of the intersection of two balls 
                    else if (Q->L.elem[q].color == IFT_GRAY) {
                        iftRemoveGQueueElem(Q, q);
                        prior->val[q]++;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
        label++;
    }

    iftIntArray *grid_all = iftListToIntArray(grid);
    iftDestroyList(&grid);    

    iftIntArray *grid_chosen = NULL;
    if (n_samples <= 0)
        grid_chosen = grid_all;
    else if (n_samples >= grid_all->n) {
        printf("Warning: Number of required samples %ld is >= total number of sampling voxels %ld\n" \
               "All sampling points will be considered\n", n_samples, grid_all->n);
        grid_chosen = grid_all;
    }
    else {
        grid_chosen = iftCreateIntArray(n_samples);
        iftShuffleIntArray(grid_all->val, grid_all->n);

        #pragma omp parallel for
        for (long i = 0; i < n_samples; i++)
            grid_chosen->val[i] = grid_all->val[i];
        iftDestroyIntArray(&grid_all);
    }


    // cleaning up
    iftDestroyAdjRel(&A);
    iftDestroyImage(&prior);
    iftDestroyImage(&label_img);
    iftDestroyGQueue(&Q);
    iftDestroyFloatArray(&dist);

    return grid_chosen;
}

float iftEstimateGridOnMaskSamplingRadius
(const iftImage *binMask, int initialObjVoxelIdx, int nSamples)
{
  if (nSamples == 1) { return 0.0; }

  bool is3D = iftIs3DImage(binMask);

  // Compute radius if each seed actually covered its entire radius 
  int totalArea = 0;

  for (int p = 0; p < binMask->n; ++p)
    if (binMask->val[p] != 0)
      totalArea += 1;
  double baseR;
  if (is3D)
    baseR = pow(((double)totalArea * 3.0) / (4.0 * IFT_PI * (double)nSamples), 1.0/3.0);
  else
    baseR = sqrt((double)totalArea/(IFT_PI * (double)nSamples));

  // Optimization method based on binary search
  double lowerBound = 0.0;
  double upperBound = IFT_INFINITY_FLT;
  int bestError = IFT_INFINITY_INT;
  int maxError = iftMax(nSamples/20, 5);
  double bestR = -1.0;
  // Arbitrary initial estimate
  double rad = baseR * (IFT_PI / 2.0);
  while (bestError > maxError && upperBound - lowerBound > IFT_EPSILON) {
    iftIntArray *seeds = iftGridSamplingOnMask(binMask, rad, initialObjVoxelIdx, 0);

    long error = labs(nSamples - seeds->n);
    if (error < bestError) {
      bestR = rad;
      bestError = error;
    }

    if (seeds->n < nSamples && rad < upperBound)
      upperBound = rad;
    if (seeds->n > nSamples && rad > lowerBound)
      lowerBound = rad;

    if (upperBound != IFT_INFINITY_FLT)
      rad = (lowerBound + upperBound) / 2.0;
    else
      rad = rad * 2;

    iftDestroyIntArray(&seeds);
  }
  assert(bestR != -1.0);

  return rad;
}

iftSet *iftMultiLabelGridSamplingOnMaskByArea
(const iftImage *label, iftImage *mask, int nSeeds )
{
  // Variables
  int totalArea, numObj, max_num_seeds;
  iftImage *newLabels;
  iftAdjRel *A;
  iftSet *seeds;

  // Init
  seeds = NULL;

  // Assign
  if( iftIs3DImage(label) ) A = iftSpheric(1.74);
  else A = iftCircular(1.45);
  totalArea = 0;

  newLabels = iftFastLabelComp(label, A);
  numObj = iftMaximumValue(newLabels);

  // Compute total area (except background)
  #pragma omp parallel for reduction(+:totalArea)
  for(int p = 0; p < newLabels->n; p++ ) {
    if(newLabels->val[p] > 0) totalArea++;
  }
  
  max_num_seeds = iftMin(totalArea, nSeeds);
  
  // For each label
  #pragma omp parallel for
  for( int i = 1; i <= numObj; i++ ) {
    // Variable
    int objArea, amount_seeds;
    float objPerc;
    iftImage* objMask;

    // Assign
    objMask = iftThreshold(newLabels, i, i, 1);
    objArea = 0;

    // Compute respective object area
    for(int p = 0; p < objMask->n; p++ ) {
      if(objMask->val[p] > 0) objArea++;
    }

    // Compute the amount of seeds (proportionally to the object area)
    objPerc = objArea / (float)totalArea;
    amount_seeds = iftRound(max_num_seeds * objPerc);

    // If it is sufficiently big
    if( amount_seeds > 0 ) {
      // Variables
      float radius;
      iftIntArray *sampled;
      
      // Assign
      radius = iftEstimateGridOnMaskSamplingRadius(objMask, -1, amount_seeds);
      sampled = iftGridSamplingOnMask(objMask, radius, -1, 0);

      // Add seeds to the set
      for( int j = 0; j < sampled->n; j++ ) {
        if( mask->val[sampled->val[j]] != 0 ) {
          #pragma omp critical
          iftInsertSet(&seeds, sampled->val[j]);
        } 
      }

      // Free
      iftDestroyIntArray(&sampled);
    }

    // Free
    iftDestroyImage(&objMask);
  }
  
  // Free
  iftDestroyImage(&newLabels);
  iftDestroyAdjRel(&A);
  
  return seeds;
}

iftSet *iftMultiLabelCentroidSamplingOnMaskByArea
(const iftImage *label, iftImage *mask, float thresh)
{
  // Variables
  int totalArea, numObj;
  iftImage *newLabels;
  iftAdjRel *A;
  iftSet *seeds;

  // Init
  seeds = NULL;

  // Assign
  if( iftIs3DImage(label) ) A = iftSpheric(1.74);
  else A = iftCircular(1.45);

  newLabels = iftFastLabelComp(label, A);
  numObj = iftMaximumValue(newLabels);

  // Compute total area (except background)
  totalArea = 0;
  #pragma omp parallel for reduction(+:totalArea)
  for( int p = 0; p < newLabels->n; p++ ) {
    if(newLabels->val[p] > 0) totalArea++;
  }
    
  // For each label
  #pragma omp parallel for
  for( int i = 1; i <= numObj; i++ ) {
    // Variables
    int objArea;
    float objPerc;
    iftImage* objMask;

    // Init
    objArea = 0;

    // Assign
    objMask = iftThreshold(newLabels, i, i, 1);
    
    // Compute respective object area
    for( int p = 0; p < objMask->n; p++ ) {
      if(objMask->val[p] > 0) objArea++;
    }

    objPerc = objArea / (float)totalArea;

    // If the object is sufficiently large
    if( objPerc >= thresh ) {
      // Variables
      int index;
      iftVoxelArray *centers;
      
      // Assign
      centers = iftGeometricCentersFromLabelImage(objMask);
      index = iftGetVoxelIndex(objMask, centers->val[1]);

      // If it is permitted
      if( mask->val[index] != 0 ) {
        #pragma omp critical
        iftInsertSet(&seeds, index);
      } 

      // Free
      iftDestroyVoxelArray(&centers);
    }

    // Free
    iftDestroyImage(&objMask);
  }
  
  // Free
  iftDestroyImage(&newLabels);
  iftDestroyAdjRel(&A);
  
  return seeds;
}

iftImage *iftGrayObjMapGridSamplOnMaskByArea
(iftImage *objsm, iftImage *mask, int k, float thresh, float obj_perc )
{ 
  // Variables
  iftImage *seed_img;

  // Assign
  seed_img = iftCreateImage(objsm->xsize, objsm->ysize, objsm->zsize);

  // Since no pixel/voxel belongs to the same label (obj and bkg), it
  // is possible to parallelize it
  #pragma omp sections
  {
    #pragma omp section // For the background
    {
      // Variables
      int min, max, nseeds;
      iftImage *bin;
      iftSet* seeds, *S;

      // Init
      seeds = NULL; S = NULL;

      // Assign
      iftMinMaxValues(objsm, &min, &max);
      bin = iftThreshold(objsm, min, thresh*max, 1);
      
      // The number of seeds must be, always, > 0
      nseeds = iftMax( iftRound(k * (1 - obj_perc)) , 1 );
      seeds = iftMultiLabelGridSamplingOnMaskByArea(bin, mask, nseeds);

      // Write the sampled seeds
      S = seeds;  
      while (S != NULL) { seed_img->val[S->elem] = 1; S = S->next;}
      
      // Free
      iftDestroyImage(&bin);
      iftDestroySet(&seeds);
      iftDestroySet(&S);
    }    
    #pragma omp section // For the object
    {
      // Variables
      int max, nseeds;
      iftImage *bin;
      iftSet* seeds, *S;

      // Init
      seeds = NULL; S = NULL;

      // Assign
      max = iftMaximumValue(objsm);
      bin = iftThreshold(objsm, thresh*max, max, 1);
      
      // The number of seeds must be, always, > 0
      nseeds = iftMax( iftRound(k * obj_perc) , 1 );
      seeds = iftMultiLabelGridSamplingOnMaskByArea(bin, mask, nseeds);

      // Write the sampled seeds
      S = seeds;
      while (S != NULL) { seed_img->val[S->elem] = 1; S = S->next;}
      
      // Free
      iftDestroyImage(&bin);
      iftDestroySet(&seeds);
      iftDestroySet(&S);
    }
  }

  
  return seed_img;
}

iftImage *iftGrayObjMapCentrSamplOnMaskByArea
(iftImage *objsm, iftImage *mask, int k, float map_thr, float seed_thr )
{
  // Variable
  iftImage *seeds;

  // Assign
  seeds = iftCreateImage(objsm->xsize, objsm->ysize, objsm->zsize);

  // Since a pixel belonging to the background cannot belong to the object at the 
  // same time, it is possible to parallelize it
  #pragma omp parallel sections
  {
    #pragma omp section // For Background
    {
      // Variables
      int min, max;
      iftImage *bin;
      iftSet *sampled;

      // Assign
      iftMinMaxValues(objsm, &min, &max);

      // Find the components belonging to the background
      bin = iftThreshold(objsm, min, map_thr * max, 1);
      sampled = iftMultiLabelGridSamplingOnMaskByArea(bin, mask, k); 

      // Write seeds
      while (sampled != NULL) { 
        seeds->val[sampled->elem] = 1; 
        sampled = sampled->next;
      }
      
      // Free
      iftDestroyImage(&bin);
      iftDestroySet(&sampled); 
    }
    #pragma omp section // For Object
    {

      // Variables
      int max;
      iftImage *bin;
      iftSet *sampled;

      // Assign
      max = iftMaximumValue(objsm);

      // Find the components belonging to the object
      bin = iftThreshold(objsm, map_thr * max, max, 1);
      sampled = iftMultiLabelCentroidSamplingOnMaskByArea(bin, mask, seed_thr); 

      // Write seeds
      while (sampled != NULL) { 
        seeds->val[sampled->elem] = 1; 
        sampled = sampled->next;
      }
      
      // Free
      iftDestroyImage(&bin);
      iftDestroySet(&sampled); 
    }
  }
  
  return seeds;
}

iftIntArray *iftGridSamplingByMaxPathValues(const iftImage *img, const iftImage *bin_mask, long n_samples) {
    puts("iftGridSamplingByMaxPathValues");

    iftIntArray *grid = iftCreateIntArray(n_samples);

    iftAdjRel *A = (iftIs3DImage(img)) ? iftSpheric(1.0) : iftCircular(1.0);
    iftImage *pathvals = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftImage *label_img = iftCreateImageFromImage(img);

    int initial_seed;
    if (bin_mask == NULL) {
        iftVoxel cv = {(img->xsize + 1) / 2, (img->ysize + 1) / 2, (img->zsize + 1) / 2};
        initial_seed = iftGetVoxelIndex(img, cv);

        #pragma omp parallel for
        for (int p = 0; p < img->n; p++)
            pathvals->val[p] = IFT_INFINITY_INT;
    }
    else {
        initial_seed = iftGetVoxelIndex(img, iftGeometricCenterVoxel(bin_mask));

        #pragma omp parallel for
        for (int p = 0; p < bin_mask->n; p++)
            pathvals->val[p] = (bin_mask->val[p] == 0) ? IFT_INFINITY_INT_NEG : IFT_INFINITY_INT;
    }

    
    iftGQueue *Q = iftCreateGQueue(IFT_QSIZE, pathvals->n, pathvals->val);
    iftGQueue *F = iftCreateGQueue(IFT_QSIZE, pathvals->n, pathvals->val);
    iftSetRemovalPolicy(F, MAXVALUE);
    iftLabeledSet *seeds = NULL;

    pathvals->val[initial_seed] = 0;
    iftInsertGQueue(&F, initial_seed);


    for (int i = 0; i < n_samples; i++) {
        int r = iftRemoveGQueue(F);
        int label = i + 1;

        grid->val[i] = r;
        label_img->val[r] = label;
        iftInsertLabeledSet(&seeds, r, label);

        _iftPerformDIFTForGridSampling(img, seeds, A, Q, pathvals, label_img, F);

        int removed_label;
        iftRemoveLabeledSet(&seeds, &removed_label);
    }
    // iftWriteImageByExt(label_img, "label_img.nii.gz");

    iftDestroyAdjRel(&A);
    iftDestroyImage(&pathvals);
    iftDestroyImage(&label_img);
    iftDestroyGQueue(&Q);
    iftDestroyGQueue(&F);
    iftDestroyLabeledSet(&seeds);

    return grid;
}


iftIntArray *iftHybridGridSamplingOnMask(const iftImage *img, const iftImage *bin_mask, int n_samples) {
    puts("iftHybridGridSamplingOnMask");
    if (bin_mask == NULL)
        iftError("Binary mask is NULL", "iftHybridGridSamplingOnMask");

    int n_samples_per_grid = n_samples / 2;

    float radius = iftEstimateGridOnMaskSamplingRadius(bin_mask, -1, n_samples_per_grid);
    iftIntArray *first_grid = iftGridSamplingOnMask(bin_mask, radius, -1, n_samples_per_grid);
    n_samples_per_grid = first_grid->n;

    iftAdjRel *A = (iftIs3DImage(img)) ? iftSpheric(1.0) : iftCircular(1.0);
    iftImage *pathvals = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftImage *label_img = iftCreateImageFromImage(img);

    iftGQueue *Q = iftCreateGQueue(IFT_QSIZE, pathvals->n, pathvals->val);

    #pragma omp parallel for
    for (int p = 0; p < bin_mask->n; p++)
        pathvals->val[p] = (bin_mask->val[p] == 0) ? IFT_INFINITY_INT_NEG : IFT_INFINITY_INT;

    for (int i = 0; i < n_samples_per_grid; i++) {
        int r = first_grid->val[i];
        pathvals->val[r] = 0;
        iftInsertGQueue(&Q, r);
        label_img->val[r] = i + 1;
    }


    while (!iftEmptyGQueue(Q)) {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftGetVoxelCoord(img, p);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(img, v)) {
                int q = iftGetVoxelIndex(img, v);
                long tmp = pathvals->val[p] + abs(img->val[q] - img->val[p]);

                if (tmp < pathvals->val[q]) {
                    if (Q->L.elem[q].color == IFT_GRAY)
                        iftRemoveGQueueElem(Q, q);

                    pathvals->val[q] = tmp;
                    iftInsertGQueue(&Q, q);
                    label_img->val[q] = label_img->val[p];
                }
            }
        }
    }    


    // the number of object labels in label_img is n_samples_per_grid
    iftIntArray *second_grid = iftCreateIntArray(n_samples_per_grid);
    iftIntArray *min_costs = iftIntRepeat(IFT_INFINITY_INT_NEG, n_samples_per_grid);

    for (int p = 0; p < label_img->n; p++) {
        int label = label_img->val[p];

        if ((label > 0) && (pathvals->val[p] > min_costs->val[label - 1])) {
            min_costs->val[label - 1] = pathvals->val[p];
            second_grid->val[label - 1] = p;
        }
    }

    iftIntArray *grid = iftCreateIntArray(n_samples_per_grid * 2);

    for (int i = 0; i < n_samples_per_grid; i++) {
        grid->val[i] = first_grid->val[i];
        grid->val[i + n_samples_per_grid] = second_grid->val[i];
    }
    // iftWriteImageByExt(label_img, "label_img.nii.gz");

    iftDestroyIntArray(&first_grid);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&pathvals);
    iftDestroyImage(&label_img);
    iftDestroyGQueue(&Q);
    iftDestroyIntArray(&second_grid);
    iftDestroyIntArray(&min_costs);

    return grid;
}


iftIntArray *iftHybridGridSamplingOnMaskISFCostFunction(const iftImage *img, const iftImage *bin_mask, int n_samples,
                                                        float alpha, float beta) {
    puts("iftHybridGridSamplingOnMaskISFCostFunction");
    if (bin_mask == NULL)
        iftError("Binary mask is NULL", "iftHybridGridSamplingOnMaskISFCostFunction");

    int n_samples_per_grid = n_samples / 2;

    float radius = iftEstimateGridOnMaskSamplingRadius(bin_mask, -1, n_samples_per_grid);
    iftIntArray *first_grid = iftGridSamplingOnMask(bin_mask, radius, -1, n_samples_per_grid);
    n_samples_per_grid = first_grid->n;

    iftAdjRel *A = (iftIs3DImage(img)) ? iftSpheric(1.0) : iftCircular(1.0);
    iftFImage *pathvals = iftCreateFImage(img->xsize, img->ysize, img->zsize);
    iftImage *root = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftImage *label_img = iftCreateImageFromImage(img);

    iftFHeap *H = iftCreateFHeap(pathvals->n, pathvals->val);


    #pragma omp parallel for
    for (int p = 0; p < bin_mask->n; p++) {
        pathvals->val[p] = (bin_mask->val[p] == 0) ? IFT_INFINITY_INT_NEG : IFT_INFINITY_INT;
        root->val[p] = -1;
    }

    for (int i = 0; i < n_samples_per_grid; i++) {
        int r = first_grid->val[i];
        pathvals->val[r] = 0;
        root->val[r] = r;
        iftInsertFHeap(H, r);
        label_img->val[r] = i + 1;
    }


    while (!iftEmptyFHeap(H)) {
        int p = iftRemoveFHeap(H);
        iftVoxel u = iftGetVoxelCoord(img, p);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(img, v)) {
                int q = iftGetVoxelIndex(img, v);
                // long tmp = pathvals->val[p] + abs(img->val[q] - img->val[p]);
                float tmp = powf(0.08 * abs(img->val[q] - img->val[root->val[p]]), 3) + iftVoxelDistance(u, v);

                if (tmp < pathvals->val[q]) {
                    if (H->color[q] == IFT_GRAY)
                        iftRemoveFHeapElem(H, q);

                    pathvals->val[q] = tmp;
                    iftInsertFHeap(H, q);
                    root->val[q] = root->val[p];
                    label_img->val[q] = label_img->val[p];
                }
            }
        }
    }    


    // the number of object labels in label_img is n_samples_per_grid
    iftIntArray *second_grid = iftCreateIntArray(n_samples_per_grid);
    iftIntArray *min_costs = iftIntRepeat(IFT_INFINITY_INT_NEG, n_samples_per_grid);

    for (int p = 0; p < label_img->n; p++) {
        int label = label_img->val[p];

        if ((label > 0) && (pathvals->val[p] > min_costs->val[label - 1])) {
            min_costs->val[label - 1] = pathvals->val[p];
            second_grid->val[label - 1] = p;
        }
    }

    iftIntArray *grid = iftCreateIntArray(n_samples_per_grid * 2);

    for (int i = 0; i < n_samples_per_grid; i++) {
        grid->val[i] = first_grid->val[i];
        grid->val[i + n_samples_per_grid] = second_grid->val[i];
    }
    // iftWriteImageByExt(label_img, "label_img.nii.gz");

    iftDestroyIntArray(&first_grid);
    iftDestroyAdjRel(&A);
    iftDestroyFImage(&pathvals);
    iftDestroyImage(&root);
    iftDestroyImage(&label_img);
    iftDestroyFHeap(&H);
    iftDestroyIntArray(&second_grid);
    iftDestroyIntArray(&min_costs);

    return grid;
}


iftIntArray *iftGridSamplingForPatchExtraction(iftImageDomain img_dom, int patch_xsize, int patch_ysize,
                                               int patch_zsize, int stride_x, int stride_y, int stride_z) {
    // area/volume inside the image where the grid points can be selected so that the grids are in the image's domain 
    iftBoundingBox bb;
    bb.begin.x = patch_xsize / 2;
    bb.begin.y = patch_ysize / 2;
    bb.end.x = img_dom.xsize - (patch_xsize / 2);
    bb.end.y = img_dom.ysize - (patch_ysize / 2);

    int n_points_x = ceil((bb.end.x - bb.begin.x + 1) / ((float) stride_x));
    int n_points_y = ceil((bb.end.y - bb.begin.y + 1) / ((float) stride_y));
    int n_points_z;

    // if it is 3D
    if (img_dom.zsize > 1) {
        bb.begin.z = patch_zsize / 2;
        bb.end.z = img_dom.zsize - (patch_zsize / 2);
        n_points_z = ceil((bb.end.z - bb.begin.z + 1)/ ((float) stride_z));
    }
    else {
        stride_z = 1; // to guarantee that the main loop will work for 2D images
        bb.begin.z = 0;
        bb.end.z = 0;
        n_points_z = 1;
    }

    int n = n_points_x * n_points_y * n_points_z;
    iftIntArray *grid = iftCreateIntArray(n);

    long i = 0;

    iftVoxel v;
    for (v.z = bb.begin.z; v.z <= bb.end.z; v.z += stride_z)
        for (v.y = bb.begin.y; v.y <= bb.end.y; v.y += stride_y)
            for (v.x = bb.begin.x; v.x <= bb.end.x; v.x += stride_x)
                grid->val[i++] = v.x + (v.y *img_dom.xsize) + (v.z * img_dom.xsize * img_dom.ysize);

    return grid;
}


iftBoundingBoxArray *iftBoundingBoxesAroundVoxels(const iftImage *img, const iftIntArray *voxel_indices, int size) {
    iftBoundingBoxArray *patches = iftCreateBoundingBoxArray(voxel_indices->n);
    iftImageDomain domain = iftGetImageDomain(img);

    iftImageDomain bb_sizes = {size, size, size};
    if (!iftIs3DImage(img))
        bb_sizes.zsize = 1;

    #pragma omp parallel for
    for (int i = 0; i < voxel_indices->n; i++) {
        iftVoxel center = iftGetVoxelCoord(img, voxel_indices->val[i]);

        iftBoundingBox bb = {.begin = {0, 0, 0}, .end = {bb_sizes.xsize-1, bb_sizes.ysize-1, bb_sizes.zsize-1}};
        bb = iftCentralizeBoundingBox(bb, center);
        patches->val[i] = iftFitBoundingBoxOnImageDomain(bb, domain); // the bb might be out of the image domain
    }


    return patches;
}


iftImage *iftLabelComponentsBySeeds(iftImage *comp, iftLabeledSet *seeds, bool incr)
{
    iftImage *label = iftCreateImageFromImage(comp);

    iftFIFO *F = iftCreateFIFO(label->n);
    iftAdjRel *A = (iftIs3DImage(comp) ? iftSpheric(1.74f) : iftCircular(1.5f));

    for (iftLabeledSet *S = seeds; S != NULL; S = S->next)
    {
        int p = S->elem;
        int l = S->label + incr;

        if (F->color[p] == IFT_WHITE)
        {
            iftInsertFIFO(F, p);

            while (!iftEmptyFIFO(F))
            {
                int q = iftRemoveFIFO(F);
                label->val[q] = l;
                iftVoxel u = iftGetVoxelCoord(label, q);

                for (int i = 1; i < A->n; i++)
                {
                    iftVoxel v = iftGetAdjacentVoxel(A, u, i);
                    if (iftValidVoxel(label, v))
                    {
                        int t = iftGetVoxelIndex(label, v);

                        if (F->color[t] == IFT_WHITE && comp->val[t] == comp->val[p])
                        {
                            if (label->val[t] == 0) {
                                iftInsertFIFO(F, t);
                            } else { /* if a component have two different labels, all component is set to label zero */
                                for (int s = 0; s < label->n; s++) {
                                    if (label->val[s] == l) label->val[s] = 0;
                                }
                                l = 0;
                            }
                        }
                    }
                }
            }
        }
    }

    iftDestroyFIFO(&F);
    iftDestroyAdjRel(&A);

    return label;
}

iftLabeledSet *iftDataSetToLabeledSeeds(iftDataSet *Z, iftImage *comp)
{
  iftLabeledSet *S=NULL;

  if (Z->ntrainsamples != Z->nsamples)
      iftError("It requires that a full training set", "iftDataSetToLabeledSeeds");
  
  if (Z->ref_data_type != IFT_REF_DATA_IMAGE && Z->ref_data_type != IFT_REF_DATA_FIMAGE && Z->ref_data_type != IFT_REF_DATA_MIMAGE)
      iftError("Reference data must be an iftImage, iftMImage or iftFImage", "iftDataSetToLabeledSeeds");

  if (Z->nclasses == 0)
      iftError("There are no labeled seeds", "iftDataSetToLabeledSeeds");
    

  if (comp == NULL){ // voxel dataset
    for (int s=0; s < Z->nsamples; s++)
      if (Z->sample[s].truelabel > 0)
	iftInsertLabeledSet(&S,Z->sample[s].id,Z->sample[s].truelabel-1);
  } else { // region (component) dataset
    for (int p=0; p < comp->n; p++)
      if (Z->sample[comp->val[p]-1].truelabel > 0){
	iftInsertLabeledSet(&S,p,Z->sample[comp->val[p]-1].truelabel-1);
      }
  }

  return(S);
}

iftLabeledSet *iftConnectObjectSeeds(iftImage *img, iftLabeledSet *seeds)
{
    if (iftNumberOfLabels(seeds) > 2) {
        iftWarning("Function only works for binary image", "iftConnectObjectSeeds");
    }

    iftImage *pred    = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftImage *pathval = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftImage *dst_marker = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftAdjRel *A      = (iftIs3DImage(img) ? iftSpheric(1.74) : iftCircular(1.45));
    iftImage *basins  = iftImageBasins(img,A);

    
    /* It first runs one IFT with fmax from the background seeds on
       the image of basins to obtain a superior reconstruction image,
       whose complement is defined as the new image of basins (i.e., a
       background saliency map). Subsequently, another IFT is executed
       with flast on the complement of the new image of basins to
       connect the object seeds */

    iftSet *bkgSeed=NULL;
    for (iftLabeledSet *S = seeds; S != NULL; S = S->next){
      int p = S->elem;
      if (S->label == 0) // background seed
    	iftInsertSet(&bkgSeed,p);
    }
    iftImage   *suprec = iftCloseBasins(basins, bkgSeed, NULL);
    iftImage   *erode  = iftErode(suprec,A,NULL);    
    iftDestroySet(&bkgSeed);
    iftDestroyImage(&basins);
    iftDestroyImage(&suprec);    
    basins = iftComplement(erode);
    iftDestroyImage(&erode);    

    iftDestroyAdjRel(&A);
    A = (iftIs3DImage(img) ? iftSpheric(1.0f) : iftCircular(1.0f));
    
    int max_val  = iftMaximumValue(basins);
    iftGQueue *Q = iftCreateGQueue(max_val+1, pathval->n, pathval->val);

    for (int p = 0; p < img->n; p++)
    {
      pred->val[p]    = IFT_NIL;
      pathval->val[p] = IFT_INFINITY_INT;
    }

    iftLabeledSet *out = NULL;
    bool         first = false;
    int   first_marker = IFT_NIL;
    int   num_markers  = 0;
    
    for (iftLabeledSet *S = seeds; S != NULL; S = S->next)
    {
        int p = S->elem;

	if (S->label == 0)
	  pathval->val[p]=IFT_INFINITY_INT_NEG; // avoid to cross background seeds
	else {
	  
	  dst_marker->val[p] = S->marker; // create marker image from object seeds

	  if (S->marker > num_markers) // count the number of markers
	    num_markers = S->marker;
	  
	  if (first == false)
	    {   /* Insert only the first object seed in the queue */
	      pathval->val[p] = 0;
	      iftInsertGQueue(&Q, p);
	      first = true;
	      first_marker = S->marker;
	    }
	}
        /* insert all seeds in the output seed set */
        iftInsertLabeledSetMarkerAndHandicap(&out, p, S->label, S->marker, 0);
    }

    char   *visited_marker       = iftAllocCharArray(num_markers+1); 
    visited_marker[first_marker] = 1;
    visited_marker[0]            = 1; /* there is no marker with
				    number 0 and background markers
				    will be mapped to this bucket, so
				    it will avoid to add seeds from
				    paths that reach background
				    markers. */
    
    /* image foresting transform */
    while (!iftEmptyGQueue(Q))
    {
        int p      = iftRemoveGQueue(Q);
        iftVoxel u = iftGetVoxelCoord(img, p);

	if (visited_marker[dst_marker->val[p]]==0){ /* add optimum path to destination as seeds */
	  int q = p;
	  while (pred->val[q] != IFT_NIL) {
	    iftInsertLabeledSetMarkerAndHandicap(&out, q, 1, dst_marker->val[p], 0);
	    /* make it thicker */
	    iftVoxel v = iftGetVoxelCoord(img, q);
	    for (int i = 1; i < A->n; i++)
	      {
		iftVoxel w = iftGetAdjacentVoxel(A, v, i);		
		if (iftValidVoxel(img, w))
		  {
		    int r = iftGetVoxelIndex(img, v);
		    iftInsertLabeledSetMarkerAndHandicap(&out, r, 1, dst_marker->val[p], 0);		    
		  }
	      }
	    q = pred->val[q];
	  }
	  visited_marker[dst_marker->val[p]]=1; /* avoid to add seeds
						   from another path that
						   reaches the same
						   marker. */
	}
	  
        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(img, v))
            {
                int q = iftGetVoxelIndex(img, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
		  if (basins->val[q] < pathval->val[q])
                    {
		      if (Q->L.elem[q].color == IFT_GRAY)
			iftRemoveGQueueElem(Q, q);
		      
		      pred->val[q] = p;
		      pathval->val[q] = basins->val[q];
		      iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyImage(&basins);
    iftDestroyImage(&pathval);
    iftDestroyGQueue(&Q);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&pred);
    iftFree(visited_marker);
    iftDestroyImage(&dst_marker);

    return(out); 
}


iftLabeledSet *iftPropagateSeedsToCluster(const iftDataSet *Z, int truelabel, float purity)
{
    int *grp_size   = iftAllocIntArray(Z->ngroups + 1);
    float *grp_perc = iftAllocFloatArray(Z->ngroups + 1);

    for (int i = 0; i < Z->nsamples; i++)
    {
        int grp = Z->sample[i].group;
        grp_size[grp] += 1;
        if (Z->sample[i].truelabel == truelabel)
            grp_perc[grp] += 1;
    }

    for (int i = 0; i < Z->ngroups + 1; i++) {
        grp_perc[i] /= grp_size[i];
    }

    iftLabeledSet *out = NULL;

    for (int i = 0; i < Z->nsamples; i++)
    {
        int grp = Z->sample[i].group;
        if (grp_perc[grp] > purity) {
            iftInsertLabeledSet(&out, Z->sample[i].id, truelabel - 1); /* decrement because of diff
                                                                          between dataset and image labels */
        }
    }

    iftFree(grp_size);
    iftFree(grp_perc);

    return out;
}


iftLabeledSet *iftSelectSeedsForEnhancement(iftLabeledSet *seeds, iftImage *groups, int label, float threshold)
{
    int ngroups = iftMaximumValue(groups) + 1;
    int nclasses = iftNumberOfLabels(seeds);

    float **perc = iftAllocFloatMatrix(ngroups, nclasses);
    int *size = iftAllocIntArray(ngroups);

    for (iftLabeledSet *S = seeds; S != NULL; S = S->next)
    {
        int grp = groups->val[S->elem];
        perc[S->label][grp] += 1;
        size[grp] += 1;
    }

    for (int i = 1; i < ngroups; i++)
    {
        for (int lb = 0; lb < nclasses; lb++) {
            perc[lb][i] /= size[i];
        }
    }

    iftLabeledSet *out = NULL;
    for (iftLabeledSet *S = seeds; S != NULL; S = S->next)
    {
        int grp = groups->val[S->elem];
        int lb = S->label;

        if ( (lb == label && perc[lb][grp] >= (1 - threshold)) || perc[lb][grp] > threshold) {
            iftInsertLabeledSet(&out, S->elem, S->label);
        }
    }

    return out;
}

iftImage *iftSamplingByOSMOX(iftImage* objsm, iftImage *mask, int num_seeds, float obj_perc, float dist_penalty)
{
  if(objsm->n < num_seeds || num_seeds < 0) {
    iftError("Invalid number of seeds!", "iftSamplingByOSMOX");
  } 
  else if(obj_perc < 0.0 || obj_perc > 1.0) {
    iftError("Invalid object percentage!", "iftSamplingByOSMOX");
  }
  if(mask != NULL) iftVerifyImageDomains(objsm, mask, "iftSamplingByOSMOX");

  int obj_seeds, bkg_seeds, max_val, min_val;
  iftSet *obj_set, *bkg_set, *s;
  iftImage *seed_img, *mask_copy, *invsm;;

  obj_set = NULL;
  bkg_set = NULL;
  s = NULL;
  
  iftMinMaxValues(objsm, &min_val, &max_val);

  if( mask == NULL ) mask_copy = iftSelectImageDomain(objsm->xsize, objsm->ysize, objsm->zsize);
  else mask_copy = iftCopyImage(mask);

  seed_img = iftCreateImage(objsm->xsize, objsm->ysize, objsm->zsize);

  // Establish the number of seeds
  if(max_val == min_val)
  {
    obj_seeds = num_seeds;
    bkg_seeds = 0;
  }
  else
  {
    obj_seeds = iftRound(num_seeds * obj_perc);
    bkg_seeds = num_seeds - obj_seeds;
  }
  
  obj_set = iftObjSalMapSamplByHighestValue(objsm, mask_copy, obj_seeds, dist_penalty);

  invsm = iftComplement(objsm);

  s = obj_set;
  while( s != NULL ) {
    mask_copy->val[s->elem] = 0;
    seed_img->val[s->elem] = 1;
    s = s->next;
  }

  bkg_set = iftObjSalMapSamplByHighestValue(invsm, mask_copy, bkg_seeds, dist_penalty);

  s = bkg_set;
  while( s != NULL ) {
    seed_img->val[s->elem] = 1;
    s = s->next;
  }

  // Free
  iftDestroyImage(&invsm);
  iftDestroySet(&obj_set);
  iftDestroySet(&bkg_set);
  iftDestroySet(&s);
  iftDestroyImage(&mask_copy);

  return (seed_img);
}

iftSet *iftObjSalMapSamplByHighestValue(iftImage *objsm, iftImage *mask, int num_seeds, float dist_penalty)
{
  if(objsm->n < num_seeds || num_seeds < 0) {
    iftError("Invalid number of seeds!", "iftSamplingByOSMOX");
  }

  if(mask != NULL) iftVerifyImageDomains(objsm, mask, "iftObjSalMapSamplByHighestValue");

  int patch_width, seed_count, total_area;
  float stdev;
  double *pixel_val;
  iftImage *mask_copy;
  iftAdjRel *A, *B;
  iftKernel *gaussian;
  iftDHeap *heap;
  iftSet *seed;
  
  // Copy the mask values
  if( mask == NULL ) mask_copy = iftSelectImageDomain(objsm->xsize, objsm->ysize, objsm->zsize);
  else mask_copy = mask;
  
  total_area = 0;

  #pragma omp parallel for reduction(+:total_area)
  for(int p = 0; p < objsm->n; p++)
  {
    if(mask_copy->val[p] != 0)
    {
      total_area++;
    }
  }

  // Estimate the gaussian influence zone
  patch_width = iftRound(sqrtf(total_area/(float)(num_seeds)));
  stdev = patch_width/dist_penalty;

  if(iftIs3DImage(objsm)) { A = iftSpheric(patch_width); B = iftSpheric(sqrtf(patch_width));}
  else { A = iftCircular(patch_width); B = iftCircular(sqrtf(patch_width));}
  
  // Copy of the object saliency map values
  pixel_val = (double *)calloc(objsm->n, sizeof(double));
  heap = iftCreateDHeap(objsm->n, pixel_val);

  // Gaussian penalty
  gaussian = iftCreateKernel(A);

  #pragma omp parallel for
  for(int i = 0; i < A->n; i++) {
    float dist;

    dist = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];
    gaussian->weight[i] = exp(-dist/(2*stdev*stdev));
  }

  seed = NULL;

  // Removing the highest values
  iftSetRemovalPolicyDHeap(heap, MAXVALUE);
  #pragma omp parallel for
  for( int p = 0; p < objsm->n; p++ ) {
    if(mask_copy->val[p] != 0) {
      iftVoxel p_voxel;
      
      p_voxel = iftGetVoxelCoord(objsm, p);
      
      pixel_val[p] = objsm->val[p];
      for(int i = 1; i < B->n; i++)
      {
        iftVoxel q_voxel;
        
        q_voxel = iftGetAdjacentVoxel(B, p_voxel, i);
        
        if(iftValidVoxel(objsm, q_voxel))
        {
          int q;
          
          q = iftGetVoxelIndex(objsm, q_voxel);
          
          pixel_val[p] += objsm->val[q];  
        }
      }
    }
    else pixel_val[p] = IFT_NIL;
  }
  
  for( int p = 0; p < objsm->n; p++ ) {
    if(pixel_val[p] != IFT_NIL) iftInsertDHeap(heap, p);
  }

  seed_count = 0;
  while( seed_count < num_seeds && !iftEmptyDHeap(heap) ) {
    int p;
    iftVoxel voxel_p;

    p = iftRemoveDHeap(heap);
    voxel_p = iftGetVoxelCoord(objsm, p);

    iftInsertSet(&seed, p);

    // Mark as removed
    pixel_val[p] = IFT_NIL;

    // For every adjacent voxel in the influence zone
    for(int i = 1; i < gaussian->A->n; i++) {
      iftVoxel voxel_q;

      voxel_q = iftGetAdjacentVoxel(gaussian->A, voxel_p, i);

      if(iftValidVoxel(objsm, voxel_q)) {
        int q;

        q = iftGetVoxelIndex(objsm, voxel_q);
        
        // If it was not removed (yet)  
        if( pixel_val[q] != IFT_NIL ) {
          // Penalize
          pixel_val[q] = (1.0 - gaussian->weight[i]) * pixel_val[q];

          iftGoDownDHeap(heap, heap->pos[q]);
        }
      }
    }

    seed_count++;
  }

  // Free
  iftDestroyKernel(&gaussian);
  iftDestroyDHeap(&heap);
  iftDestroyAdjRel(&A);
  iftDestroyAdjRel(&B);
  if(mask == NULL) iftDestroyImage(&mask_copy);

  return (seed);
}


iftDoubleMatrix *iftSeedsFeatures(iftLabeledSet *seeds, const iftMImage *mimg, int label)
{
    iftLabeledSet *set = ((label < 0) ? iftCopyLabeledSet(seeds) : iftCopyLabels(seeds, label));

    int size = iftLabeledSetSize(set);

    iftDoubleMatrix *M = iftCreateDoubleMatrix(mimg->m, size);

    int i = 0;
    for (iftLabeledSet *s = set; s; s = s->next) {
        for (int j = 0; j < mimg->m; j++) {
            iftMatrixElem(M, j, i) = mimg->val[s->elem][j];
        }
        i++;
    }

    return M;
}


void iftSeedsSuperpixelBins(iftLabeledSet *seeds, const iftMImage *mimg, const iftImage *superpixel,
                            iftDoubleMatrix *data, iftIntArray *label)
{
    int n_spp = iftMaximumValue(superpixel) + 1;

    double *bins_data = iftAlloc(mimg->m * n_spp, sizeof *bins_data);
    int *bins_count = iftAlloc(n_spp, sizeof *bins_count);
    int *bins_label = iftAlloc(n_spp, sizeof *bins_label);

    for (iftLabeledSet *s = seeds; s; s = s->next) {
        int p = s->elem;
        int spp = superpixel->val[p];
        int row = spp * mimg->m;
        for (int b = 0; b < mimg->m; b++) {
            bins_data[row + b] += mimg->val[p][b];
            bins_count[spp] += 1;
            if (s->label > 0) bins_label[spp] += 1;
            else bins_label[spp] -= 1;
        }
    }

    int count = 0;
    for (int i = 0; i < n_spp; i++) {
        if (bins_count[i] > 0) count++;
    }

    label->val = iftRealloc(label->val, count * (sizeof *label->val));
    label->n = count;

    iftReallocDoubleMatrix(data, mimg->m, count);

    int ii = 0;
    for (int i = 0; i < n_spp; i++) {
        int row = i * mimg->m;
        if (bins_count[i] > 0) {
            for (int j = 0; mimg->m; j++) {
                iftMatrixElem(data, j, ii) = bins_data[row + j] / bins_count[i];
            }
            ii++;
            label->val[i] = ((bins_label[i] > 0) ? 1 : 0);
        }
    }

    iftFree(bins_count);
    iftFree(bins_data);
    iftFree(bins_label);
}

iftImage *iftDrawDilatedSeeds(const iftImage *image, const iftLabeledSet *seeds,
                              const iftAdjRel *A, const iftColorTable *ctb_rgb)
{
    if (iftNumberOfLabels(seeds) > ctb_rgb->ncolors)
        iftError("Number of labels is greater than color table size", "DrawLabeledSet");

    iftImage *out = iftCopyImage(image);

    for (const iftLabeledSet *s = seeds; s; s = s->next)
    {
        iftVoxel u = iftGetVoxelCoord(image, s->elem);
        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (!iftValidVoxel(image, v))
                continue;

            int q = iftGetVoxelIndex(image, v);
            iftSetYCbCr(out, q, ctb_rgb->color[s->label]);
        }
    }

    return out;
}

void iftRemoveFrameCompInPlace(iftImage *bin)
{
  iftGQueue *Q    = iftCreateGQueue(iftMaximumValue(bin)+1,bin->n,bin->val);
  iftAdjRel *A    = NULL;

  if (iftIs3DImage(bin)){
    A = iftSpheric(sqrt(3.0));
    for (int p = 0; p < bin->n; p++){
      if (bin->val[p] != 0) {
	iftVoxel u = iftGetVoxelCoord(bin,p);
	if ((u.x == 0)||(u.y == 0)||(u.z == 0))
	  if (Q->L.elem[p].color == IFT_WHITE)
	    iftInsertGQueue(&Q,p);
      }
    }
  } else {
    A = iftCircular(sqrt(2.0));
    for (int p = 0; p < bin->n; p++){
      if (bin->val[p] != 0) {
	iftVoxel u = iftGetVoxelCoord(bin,p);
	if ((u.x == 0)||(u.y == 0))
	  if (Q->L.elem[p].color == IFT_WHITE)
	    iftInsertGQueue(&Q,p);
      }
    }
  }

  while (!iftEmptyGQueue(Q)) {
    int p        = iftRemoveGQueue(Q);
    bin->val[p]  = 0;
    iftVoxel u   = iftGetVoxelCoord(bin,p);
    for (int i=1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(bin,v)){
	int q = iftGetVoxelIndex(bin,v);
	if ((bin->val[q] != 0)&&(Q->L.elem[p].color == IFT_WHITE)){
	  iftInsertGQueue(&Q,q);
	}
      }
    }
  }
  iftDestroyAdjRel(&A);
  iftDestroyGQueue(&Q);
  
}

iftImage *iftRemoveFrameComp(iftImage *bin)
{
  iftGQueue *Q    = iftCreateGQueue(iftMaximumValue(bin)+1,bin->n,bin->val);
  iftAdjRel *A    = NULL;
  iftImage *fbin  = iftCopyImage(bin);
  
  if (iftIs3DImage(fbin)){
    A = iftSpheric(sqrt(3.0));
    for (int p = 0; p < fbin->n; p++){
      if (fbin->val[p] != 0) {
	iftVoxel u = iftGetVoxelCoord(fbin,p);
	if ((u.x == 0)||(u.y == 0)||(u.z == 0))
	  if (Q->L.elem[p].color == IFT_WHITE)
	    iftInsertGQueue(&Q,p);
      }
    }
  } else {
    A = iftCircular(sqrt(2.0));
    for (int p = 0; p < fbin->n; p++){
      if (fbin->val[p] != 0) {
	iftVoxel u = iftGetVoxelCoord(fbin,p);
	if ((u.x == 0)||(u.y == 0))
	  if (Q->L.elem[p].color == IFT_WHITE)
	    iftInsertGQueue(&Q,p);
      }
    }
  }

  while (!iftEmptyGQueue(Q)) {
    int p        = iftRemoveGQueue(Q);
    fbin->val[p]  = 0;
    iftVoxel u   = iftGetVoxelCoord(fbin,p);
    for (int i=1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(fbin,v)){
	int q = iftGetVoxelIndex(fbin,v);
	if ((fbin->val[q] != 0)&&(Q->L.elem[p].color == IFT_WHITE)){
	  iftInsertGQueue(&Q,q);
	}
      }
    }
  }
  iftDestroyAdjRel(&A);
  iftDestroyGQueue(&Q);

  return(fbin);
}

