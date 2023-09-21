#include "iftSimilarity.h"

#include "ift/core/dtypes/GQueue.h"
#include "ift/core/dtypes/IntArray.h"
#include "ift/core/dtypes/IntMatrix.h"
#include "ift/core/io/Stream.h"


/********************** PRIVATE FUNCTIONS *************************/

/******************************************************************/


/********************** PUBLIC FUNCTIONS *************************/
double iftDiceSimilarity(const iftImage *bin_source, const iftImage *bin_target) {
    iftVerifyImageDomains(bin_source, bin_target, "iftDiceSimilarity");

    ulong vol_intersec = 0; // |intersection(volume(bin_source), volume(bin_target))|
    ulong vol_sum      = 0; // |volume(bin_source)| + |volume(bin_target)|

    iftVerifyImageDomains(bin_source, bin_target, "iftSimilarityByDice");

    for (size_t i = 0; i < bin_target->n; i++) {
      vol_intersec += ((bin_source->val[i] != 0) && (bin_target->val[i] != 0));
      vol_sum      += ((bin_source->val[i] != 0) + (bin_target->val[i] != 0));
    }

    double dice = (2.0 * vol_intersec) / (vol_sum);

    return dice;
}


iftDblArray *iftDiceSimilarityMultiLabel(const iftImage *label_source,
                                         const iftImage *label_target, int n_objects) {
    // CHECKERS
    iftVerifyImageDomains(label_source, label_target, "iftDiceSimilarityMultiLabel");

    iftIntArray *vol_intersec = iftCreateIntArray(n_objects + 1);
    iftIntArray *vol_sum      = iftCreateIntArray(n_objects + 1);
    iftDblArray *dices        = iftCreateDblArray(n_objects + 1);

    // binary images
    if (n_objects == 1) {
        dices->val[0] = dices->val[1] = iftDiceSimilarity(label_source, label_target);
    }
    // label images
    else {
        for (size_t p = 0; p < label_source->n; p++) {
            vol_intersec->val[label_source->val[p]] += (label_source->val[p] == label_target->val[p]);
            vol_sum->val[label_source->val[p]]++;
            vol_sum->val[label_target->val[p]]++;
        }

        for (int o = 1; o <= n_objects; o++) {
            dices->val[o] = (2.0 * vol_intersec->val[o]) / vol_sum->val[o];
            dices->val[0] += dices->val[o];
        }
        dices->val[0] /= n_objects;
    }

    iftDestroyIntArray(&vol_intersec);
    iftDestroyIntArray(&vol_sum);

    return dices;
}


/* For sake of efficiency, it assumes that the objects are centralized in the same image domain */
double iftASD(const iftImage *bin_source, const iftImage *bin_target) {
    iftVerifyImageDomains(bin_source, bin_target, "iftASD");

    // get the adjacency relations
    iftAdjRel *A = NULL, *B = NULL;
    if (iftIs3DImage(bin_source)) {
        A = iftSpheric(1.75); // ~= sqrtf(3.0) -- used in EDT
        B = iftSpheric(1.0); // used to find out the object borders
    }
    else {
        A = iftSpheric(1.42); // ~= sqrtf(2.0) -- used in EDT
        B = iftCircular(1.0); // used to find out the object borders
    }

    // get the border spels
    int n_source_border_spels;
    iftBMap *Sb = iftObjectBorderBitMap(bin_source, B, &n_source_border_spels); // source border spels
    iftSet *Tb  = iftObjectBorderSet(bin_target, B); // target border spels
    iftDestroyAdjRel(&B);

    // IFT initialization
    iftImage *dist = iftCreateImage(bin_source->xsize, bin_source->ysize, bin_source->zsize);
    iftImage *root = iftCreateImage(bin_source->xsize, bin_source->ysize, bin_source->zsize);
    iftGQueue *Q   = iftCreateGQueue(IFT_QSIZE, bin_source->n, dist->val);
    iftSetImage(dist, IFT_INFINITY_INT);

    // Initialization of Target Border spels (seeds) - all has distance 0
    while (Tb != NULL) {
        int p        = iftRemoveSet(&Tb);
        dist->val[p] = 0;
        root->val[p] = p;
        iftInsertGQueue(&Q, p);
    }
    iftDestroySet(&Tb);


    // Image Foresting Transform
    double mean_dist           = 0.0;
    int n_visited_border_spels = 0;
    // for each object border of the target binary image
    while (!iftEmptyGQueue(Q)) {
        int p = iftRemoveGQueue(Q); // the spel p of the Source Image already reached the minumum distance to the target border

        //////////////////////////////////////////////////////////////////////////////////////////////////////
        // If the voxel p is on the border of Source Binary Image, we increment the number of visited       //
        // source border spels.                                                                             //
        // The IFT-EDT stops when all border spels from the Source Image have figured out its minimum       //
        // distance to a border spel from the Target Image, which happens when the former leave the Queue.  //
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        if (iftBMapValue(Sb, p)) {
            mean_dist += sqrtf(dist->val[p]); // euclidean distance

            n_visited_border_spels++;
            // all border spels from the Source Image figured out the minimum distance to the border of the Target Image
            if (n_visited_border_spels == n_source_border_spels)
                break; // it terminates the IFT-EDT
        }

        // Gets the voxel and its root.
        iftVoxel u = iftGetVoxelCoord(bin_source, p);
        iftVoxel r = iftGetVoxelCoord(bin_source, root->val[p]);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftValidVoxel(bin_source, v)) {
                int q   = iftGetVoxelIndex(bin_source, v);
                int tmp = iftSquaredVoxelDistance(v, r);

                if (tmp < dist->val[q]) {
                    if (dist->val[q] != IFT_INFINITY_INT)
                        iftRemoveGQueueElem(Q, q);

                    dist->val[q] = tmp;
                    root->val[q] = root->val[p];
                    iftInsertGQueue(&Q, q);
                }
            }
        }
    }

    mean_dist /= (n_source_border_spels*1.0); // to ensure that will be double

    iftDestroyGQueue(&Q);
    iftDestroyImage(&root);
    iftDestroyImage(&dist);
    iftDestroyAdjRel(&A);
    iftDestroyBMap(&Sb);

    return mean_dist;
}


iftDblArray *iftASDMultiLabel(const iftImage *label_source, const iftImage *label_target, int n_objects) {
    iftVerifyImageDomains(label_source, label_target, "iftASDMultiLabel");

    iftDblArray *asd_array = iftCreateDblArray(n_objects + 1);

    // binary images
    if (n_objects == 1) {
        asd_array->val[0] = asd_array->val[1] = iftASD(label_source, label_target);
    }
    // label images
    else {
        #pragma omp parallel for
        for (int o = 1; o <= n_objects; o++) {
            // it creates a binary image for each object, where the object has value 1
            iftImage *tmp_bin_source = iftCreateImage(label_source->xsize, label_source->ysize, label_source->zsize);
            iftImage *tmp_bin_target = iftCreateImage(label_target->xsize, label_target->ysize, label_target->zsize);

            for (int p = 0; p < label_source->n; p++) {
                tmp_bin_source->val[p] = (label_source->val[p] == o);
                tmp_bin_target->val[p] = (label_target->val[p] == o);
            }

            asd_array->val[o] = iftASD(tmp_bin_source, tmp_bin_target);

            iftDestroyImage(&tmp_bin_source);
            iftDestroyImage(&tmp_bin_target);
        }

        for (int o = 1; o <= n_objects; o++)
            asd_array->val[0] += asd_array->val[o];
        asd_array->val[0] /= n_objects;
    }

    return asd_array;
}


double iftASSD(const iftImage *bin_source, const iftImage *bin_target) {
    // CHECKERS
    iftVerifyImageDomains(bin_source, bin_target, "iftASSD");

    double border_dist1, border_dist2;
    double result = 0.0;

    border_dist1 = iftASD(bin_source, bin_target);
    border_dist2 = iftASD(bin_target, bin_source);

    result = (border_dist1 + border_dist2) / 2.0;

    return (result * bin_source->dx);
}


iftDblArray *iftASSDMultiLabel(const iftImage *label_source, const iftImage *label_target, int n_objects) {
    iftVerifyImageDomains(label_source, label_target, "iftASSDMultiLabel");

    iftDblArray *assd_array = iftCreateDblArray(n_objects + 1);

    // binary images
    if (n_objects == 1) {
        assd_array->val[0] = assd_array->val[1] = iftASSD(label_source, label_target);
    }
    // label images
    else {
        #pragma omp parallel for
        for (int o = 1; o <= n_objects; o++) {
            // it creates a binary image for each object, where the object has value 1
            iftImage *tmp_bin_source = iftCreateImage(label_source->xsize, label_source->ysize, label_source->zsize);
            iftImage *tmp_bin_target = iftCreateImage(label_target->xsize, label_target->ysize, label_target->zsize);

            for (int p = 0; p < label_source->n; p++) {
                tmp_bin_source->val[p] = (label_source->val[p] == o);
                tmp_bin_target->val[p] = (label_target->val[p] == o);
            }

            assd_array->val[o] = iftASSD(tmp_bin_source, tmp_bin_target);

            iftDestroyImage(&tmp_bin_source);
            iftDestroyImage(&tmp_bin_target);
        }

        for (int o = 1; o <= n_objects; o++)
            assd_array->val[0] += assd_array->val[o];
        assd_array->val[0] /= n_objects;
    }

    return assd_array;
}


double iftBorderGradMatchingScore(iftSet *borders, iftImage *grad) {
    if (borders == NULL)
        iftError("Border Set is NULL", "iftBorderGradMatchingScore");
    if (grad == NULL)
        iftError("Gradient Image is NULL", "iftBorderGradMatchingScore");

    double score = 0.0;
    for(iftSet *S = borders; S != NULL; S = S->next) {
        int p = S->elem;
        score += grad->val[p];
    }

    return score;
}


double iftBorderRegionOfInfluenceGradMatchingScore(iftSet *borders, iftImage *grad, iftAdjRel *A, double max_dist) {
    iftImage *dist = NULL;
    iftImage *root = NULL;
    iftImage *region_of_inf_grad = NULL;
    iftSet *S = NULL;
    iftGQueue *Q = NULL;
    double score = 0.0;
    int p, q, tmp;
    iftVoxel u,v,r;

    if (borders == NULL)
        iftError("Border Set is NULL", "iftBorderRegionOfInfluenceGradMatchingScore");
    if (grad == NULL)
        iftError("Gradient Image is NULL", "iftBorderRegionOfInfluenceGradMatchingScore");
    if (grad == NULL)
        iftError("Gradient Image is NULL", "iftBorderRegionOfInfluenceGradMatchingScore");

    max_dist           *= max_dist; // Computing the maximum squared distance

    dist               = iftCreateImage(grad->xsize,grad->ysize,grad->zsize);
    root               = iftCreateImage(grad->xsize,grad->ysize,grad->zsize);
    region_of_inf_grad = iftCreateImage(grad->xsize,grad->ysize,grad->zsize);
    Q                  = iftCreateGQueue(IFT_QSIZE, grad->n, dist->val);

    // Setting the initial distance to all voxels to infinity, since we don't care about whether the maximum
    // gradient is inside or outside the object
    iftSetImage(dist, IFT_INFINITY_INT);
    // Setting the region of influence gradient values to negative infinity, since this map will only hold maximum
    // gradient values to voxels on the given object border
    iftSetImage(region_of_inf_grad, IFT_INFINITY_INT_NEG);

    // Initializing seeds
    for(S = borders; S != NULL; S = S->next) {
        p = S->elem;
        dist->val[p]=0;
        root->val[p]=p;
        // Initializing the region of influence gradient values using the original gradient for the seeds.
        // This is unecessary since max_dist will always be a non-negative value, but we leave it here anyway.
        // One may skip the next loop altogether if max_dist is 0.0, but we leave it for greater clarity
        region_of_inf_grad->val[p] = grad->val[p];
        iftInsertGQueue(&Q,p);
    }

    // Image Foresting Transform
    while(!iftEmptyGQueue(Q)) {
        p=iftRemoveGQueue(Q);

        // Stopping computation if the maximum distance to a voxel p is more than the accepted threshold
        if (dist->val[p] <= max_dist) {
            u = iftGetVoxelCoord(grad, p);
            r = iftGetVoxelCoord(grad, root->val[p]);

            // Computing the maximum gradient value for a border voxel (p's root) in its influence region.
            // We do this here since p's root is certainly defined at this point
            region_of_inf_grad->val[root->val[p]] = iftMax(grad->val[p], region_of_inf_grad->val[root->val[p]]);

            for (int i = 1; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(grad, v)) {
                    q = iftGetVoxelIndex(grad, v);
                    if (dist->val[q] > dist->val[p]) {
                        tmp = iftSquaredVoxelDistance(v, r);
                        if (tmp < dist->val[q]) {
                            if (dist->val[q] != IFT_INFINITY_INT)
                                iftRemoveGQueueElem(Q, q);
                            dist->val[q] = tmp;
                            root->val[q] = root->val[p];
                            iftInsertGQueue(&Q, q);
                        }
                    }
                }
            }
        }
    }


    // Computing the score for the border voxels using the maximum gradients of the region of influence
    score = iftBorderGradMatchingScore(borders, region_of_inf_grad);

    /* Cleaning up! */
    iftDestroyImage(&dist);
    iftDestroyImage(&root);
    iftDestroyImage(&region_of_inf_grad);
    iftDestroyGQueue(&Q);


    return score;
}



/* For sake of efficiency, it assumes that the objects are centralized
   in the same image domain */
float iftMeanDistBetweenBoundaries(iftImage *bin1, iftImage * bin2)
{
  iftImage  *dist=NULL,*root=NULL,*label=NULL;
  int        i,p,q,tmp;
  iftVoxel   u,v,r;
  iftGQueue *Q=NULL;
  iftSet    *S1=NULL, *S2=NULL;
  iftAdjRel *A=NULL,*B=NULL;
  int        meansize, *seed;
  float      distance=0.0, nvoxels;

  iftVerifyImageDomains(bin1,bin2,"iftMeanDistBetweenBoundaries");

  // Initialization

  if (iftIs3DImage(bin1)) {
    A = iftSpheric(sqrtf(3.0));
    B = iftSpheric(1.0);
  }
  else {
    A = iftSpheric(sqrtf(2.0));
    B = iftCircular(1.0);
  }

  S1          = iftObjectBorderSet(bin1, B);
  S2          = iftObjectBorderSet(bin2, B);
  iftDestroyAdjRel(&B);

  dist  = iftCreateImage(bin1->xsize,bin1->ysize,bin1->zsize);
  root  = iftCreateImage(bin1->xsize,bin1->ysize,bin1->zsize);
  label = iftCreateImage(bin1->xsize,bin1->ysize,bin1->zsize);
  Q  = iftCreateGQueue(IFT_QSIZE,bin1->n,dist->val);

  iftSetImage(dist, IFT_INFINITY_INT);

  /* Insert boundaries in queue in alternate order to minimize the
     differences between distance from A to B and distance from B to
     A */

  seed = iftAllocIntArray(bin1->n);

  meansize=0; i=0;
  while (S1 != NULL) {
    p = iftRemoveSet(&S1);
    dist->val[p]=0;
    root->val[p]=p;
    label->val[p]=1;
    meansize++;
    seed[i]=p; i+=2;
  }
  i=1;
  while (S2 != NULL) {
    p = iftRemoveSet(&S2);
    dist->val[p]=0;
    root->val[p]=p;
    label->val[p]=2;
    meansize++;
    seed[i]=p; i+=2;
  }

  for (i=0; i < meansize; i++)
    {
      p = seed[i];
      if (Q->L.elem[p].color == IFT_WHITE)
    iftInsertGQueue(&Q,p);
    }
  iftFree(seed);
  meansize /= 2;

  // Image Foresting Transform

  distance=0.0; nvoxels = meansize;

  while((!iftEmptyGQueue(Q))&&(nvoxels > 0)) {

    p=iftRemoveGQueue(Q);

    //Gets the voxel and its root.
    u = iftGetVoxelCoord(bin1,p);
    r = iftGetVoxelCoord(bin1,root->val[p]);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(bin1,v)){
    q = iftGetVoxelIndex(bin1,v);
    if (Q->L.elem[q].color != IFT_BLACK){
      tmp = iftSquaredVoxelDistance(v,r);
      if (tmp < dist->val[q]){
        if (dist->val[q] != IFT_INFINITY_INT)
          iftRemoveGQueueElem(Q, q);
        dist->val[q]  = tmp;
        root->val[q]  = root->val[p];
        label->val[q]  = label->val[p];
        iftInsertGQueue(&Q, q);
      }
    } else {
      if (label->val[q]!=label->val[p]){
        distance += sqrtf(dist->val[p])+sqrtf(dist->val[q])+(iftVoxelDistance(u,v));
        nvoxels--;
      }
    }
      }
    }
  }

  iftDestroyGQueue(&Q);
  iftDestroyImage(&root);
  iftDestroyImage(&dist);
  iftDestroyImage(&label);
  iftDestroyAdjRel(&A);

  return(distance/meansize);
}


double iftShannonEntropy(const iftImage *image)
{
   int min, max, i, size;
   double *frequencies, entropy, logn;

   entropy = 0.0;
   min = iftMinimumValue(image);
   max = iftMaximumValue(image);
   size = max - min + 1;
   frequencies = iftAllocDoubleArray(size);
   for (i = 0; i < image->n; i++)
      frequencies[image->val[i]-min]+= 1.0;
   logn = log(image->n);
   for (i = 0; i < size; i++)
      if (frequencies[i] > 0.0)
          entropy += frequencies[i] * (log(frequencies[i])-logn);
   entropy = (-1.0)*entropy / image->n;
   iftFree(frequencies);
   return entropy;
}

double iftJointEntropy(const iftImage *image1, const iftImage *image2)
{
   int min1, max1, min2, max2, i, size, size2;
   double *frequencies, entropy, logn;

   if ((image1->xsize != image2->xsize)||(image1->ysize != image2->ysize)||(image1->zsize != image2->zsize))
      iftWarning("Images do not have the same size. The behavior of this function is not determined.", "iftJointEntropy");

   entropy = 0.0;
   min1 = iftMinimumValue(image1);
   max1 = iftMaximumValue(image1);
   min2 = iftMinimumValue(image2);
   max2 = iftMaximumValue(image2);
   size2 = (max2-min2+1);
   size = (max1-min1+1) * (max2-min2+1);
   logn = log(image1->n);
   frequencies = iftAllocDoubleArray(size);
   for (i = 0; i < image1->n; i++)
      frequencies[(image1->val[i]-min1)*size2 + image2->val[i]-min2]+= 1.0;
   for (i = 0; i < size; i++)
      if (frequencies[i] > 0.0)
         entropy += frequencies[i] * (log(frequencies[i])-logn);

   //Obs: Divided by 2 in order to get a result between 0.0 and 1.0
   entropy = (-1.0) * entropy * 2 / image1->n;
   iftFree(frequencies);
   return entropy;
}

double iftNormalizedMutualInformation(const iftImage *image1, const iftImage *image2)
{
   double entropy1, entropy2, jointEntropy, nmi;
   entropy1 = iftShannonEntropy(image1);
   entropy2 = iftShannonEntropy(image2);
   jointEntropy = iftJointEntropy(image1, image2);
   nmi = (entropy1 + entropy2) / jointEntropy;
   return nmi;
}


iftMatrix *iftComputeScoreMatrixByDice(const iftFileSet *label_img_files, int n_objs) {
    iftMatrix *score_matrix = iftCreateMatrix(label_img_files->n, label_img_files->n);

    #pragma omp parallel for
    for (size_t i = 0; i < label_img_files->n; i++) {
        iftImage *source_img = iftReadImageByExt(label_img_files->files[i]->path);

        for (size_t j = i; j < label_img_files->n; j++) {
            iftImage *target_img = iftReadImageByExt(label_img_files->files[j]->path);
            iftDblArray *dices   = iftDiceSimilarityMultiLabel(source_img, target_img, n_objs);

            score_matrix->val[iftGetMatrixIndex(score_matrix, i, j)] = dices->val[0]; // average dice
            score_matrix->val[iftGetMatrixIndex(score_matrix, j, i)] = score_matrix->val[iftGetMatrixIndex(score_matrix, i, j)];

            iftDestroyDblArray(&dices);
            iftDestroyImage(&target_img);
        }
        iftDestroyImage(&source_img);
    }

    return score_matrix;
}



iftMatrix *iftComputeScoreMatrixByASSD(const iftFileSet *label_img_files, int n_objs) {
    iftMatrix *score_matrix = iftCreateMatrix(label_img_files->n, label_img_files->n);

    #pragma omp parallel for
    for (size_t i = 0; i < label_img_files->n; i++) {
        iftImage *source_img = iftReadImageByExt(label_img_files->files[i]->path);

        for (size_t j = i; j < label_img_files->n; j++) {
            iftImage *target_img = iftReadImageByExt(label_img_files->files[j]->path);
            iftDblArray *assds   = iftASSDMultiLabel(source_img, target_img, n_objs);

            score_matrix->val[iftGetMatrixIndex(score_matrix, i, j)] = assds->val[0];
            score_matrix->val[iftGetMatrixIndex(score_matrix, j, i)] = score_matrix->val[iftGetMatrixIndex(score_matrix, i, j)];

            iftDestroyDblArray(&assds);
            iftDestroyImage(&target_img);
        }
        iftDestroyImage(&source_img);
    }

    return score_matrix;
}


double iftImageEuclideanDistance(const iftImage *img1, const iftImage *img2) {
    iftVerifyImages(img1, img2, "iftImageEuclideanDistance");

    double imed = 0.0;

    #pragma omp parallel for reduction(+:imed)
    for (int p = 0; p < img1->n; p++)
        for (int q = 0; q < img1->n; q++) {
            iftVoxel u = iftGetVoxelCoord(img1, p);
            iftVoxel v = iftGetVoxelCoord(img1, q);

            double dist = iftSquaredVoxelDistance(u, v) / 2.0;
            imed += (exp(-dist) * (img1->val[p] - img2->val[p]) * (img1->val[q] - img2->val[q]));
        }

    imed = sqrtf(imed / (2 * IFT_PI));

    return imed;
}


double iftASA(const iftImage* super_map, const iftImage* gt) {
    iftVerifyImageDomains(super_map, gt, "iftASA");
    
    int min_superpixel;
    int n_superpixels;
    iftMinMaxValues(super_map, &min_superpixel, &n_superpixels);
    
    if (min_superpixel <= 0)
        iftError("Superpixel Map must only have values >= 1.\nFound: %d", "iftASA", min_superpixel);
    
    int min_obj;
    int n_objs;
    iftMinMaxValues(gt, &min_obj, &n_objs);
    
    if (min_obj < 0)
        iftError("Ground-Truth Label Image must only have values >= 0, where 0 is the background.\nFound: %d",
                 "iftASA", min_obj);
    
    
    // [i][j] = number of pixels with from the object j inside the superpixel i
    // 0 is the background
    iftIntMatrix* n_objs_per_superpixels_mat = iftCreateIntMatrix(n_objs + 1, n_superpixels + 1);
    
    for (int p = 0; p < super_map->n; p++) {
        int superpixel = super_map->val[p]; // row
        int label = gt->val[p]; // col
        
        iftMatrixElem(n_objs_per_superpixels_mat, label, superpixel)++;
    }
    
    double asa = 0.0;

    #pragma omp parallel for
    for (int superpixel = 1; superpixel < n_objs_per_superpixels_mat->nrows; superpixel++) {
        int max_intersec = -1;
        
        // finds the maximum intersection between the superpixel and the GT objects
        for (int label = 0; label < n_objs_per_superpixels_mat->ncols; label++) {
            if (max_intersec < iftMatrixElem(n_objs_per_superpixels_mat, label, superpixel)) {
                max_intersec = iftMatrixElem(n_objs_per_superpixels_mat, label, superpixel);
            }
        }

        #pragma omp atomic
        asa += max_intersec;
    }
    
    asa /= gt->n;
    
    iftDestroyIntMatrix(&n_objs_per_superpixels_mat);
    
    return asa;
}


double iftGeneralBalancedCoeff(const iftImage *bin_source, const iftImage *bin_target) {
    iftVerifyImageDomains(bin_source, bin_target, "iftGeneralBalancedCoeff(");

    long obj_size = 0, false_pos = 0, false_neg = 0;
    for (int i = 0; i < bin_target->n; i++) {
        obj_size += ( bin_target->val[i] != 0 );
        false_pos += ( bin_target->val[i] != 0 && ! bin_source->val[i] );
        false_neg += ( ! bin_target->val[i] &&  bin_source->val[i] != 0);
    }

    if (false_neg >= obj_size)
        return 0;

    double gbc = (obj_size - false_pos) * (obj_size - false_neg) / (double)( obj_size * obj_size);

    return gbc;
}


