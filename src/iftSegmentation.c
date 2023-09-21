#include "iftSegmentation.h"

#include "ift/core/dtypes/FIFO.h"
#include "ift/core/dtypes/IntMatrix.h"
#include "ift/core/dtypes/IntQueue.h"
#include "ift/core/dtypes/LIFO.h"
#include "ift/core/io/Stream.h"
#include "ift/imgproc/basic/Histogram.h"
#include "iftFunctions.h"


/*---------------------- Private functions ------------------------------- */
/**
 * @brief It defines the API of functions to compute the color distance between two voxels from an image.
 *
 * @author Samuka Martins
 * @date Dec 7, 2017
 *
 * @param img Input Image.
 * @param p   Voxel p.
 * @param q   Voxel q.
 * @return The computed distance.
 */
typedef float (*iftDistColorVoxels)(const iftImage *img, int p, int q);


float iftDistColorVoxelsColorImage(const iftImage *img, int p, int q) {
    // added a factor to color component influence, since the range of values from
    // the Y is much higher than the Cb and Cr
    return sqrtf(0.2 * iftPowerOfTwo(img->val[q] - img->val[p]) +
                 1.0 * iftPowerOfTwo(img->Cb[q] - img->Cb[p]) +
                 1.0 * iftPowerOfTwo(img->Cr[q] - img->Cr[p]));   
}


float iftDistColorVoxelsGreyImage(const iftImage *img, int p, int q) {
    return abs(img->val[q] - img->val[p]);   
}


int iftFindRelaxationRoot(iftImage *pred, iftImage *new_label, iftImage *old_label, int p, char *valid)
{
  int  q = pred->val[p], r = p;

  while ((q != IFT_NIL) && (new_label->val[q] != old_label->val[q])) {
    r     = q;
    q     = pred->val[q];
  }
  
  /* verify if the predecessor of the root candidate is connected to
     the root by a path that did not change the label */

  if (iftIsLabelRootConnected(pred,new_label,q)){
    *valid = 1;
  }else{
    *valid = 0;
  }

  return(r);
}

void iftFindRelaxationRoot1(iftImageForest *fst, iftImage *new_label, int p, int *root, char *color)
{
  /* Find the root of the path recursively, while color is IFT_WHITE, and
    identify the relaxation root when returning from recursion, as the
    first node whose predecessor has a distinct new label, or a root
    who changed label. */
  if (fst->pred->val[p] == IFT_NIL)
  {
    *color = IFT_GRAY;
    if (new_label->val[p] != fst->label->val[p])
    {
      *root  = p; *color = IFT_BLACK;
    }
  }
  else
  {
    iftFindRelaxationRoot1(fst, new_label, fst->pred->val[p], root, color);
    if (*color == IFT_GRAY)
    {
      if (new_label->val[p] != new_label->val[fst->pred->val[p]])
      {
        *root  = p; *color = IFT_BLACK;
      }
    }
  }
}

iftFImage *iftWeightNormFactor(const iftFImage *weight, iftAdjRel *A)
{
  iftFImage *norm_factor=iftCreateFImage(weight->xsize,weight->ysize,weight->zsize);

#pragma omp parallel for shared(weight,norm_factor)
  for (int p=0; p < weight->n; p++) {
    iftVoxel u = iftFGetVoxelCoord(weight,p);
    for (int i=1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftFValidVoxel(weight,v)){
	     int q = iftFGetVoxelIndex(weight,v);
	     norm_factor->val[p] += weight->val[q];
      }
    }
  }

  return(norm_factor);
}

/*
@brief This function is deprecated. See iftRelaxObjects.
*/
void iftSmoothFrontier(iftImageForest *fst, iftSet **Frontier, iftBMap *inFrontier,
                       iftFImage *border_weight, iftFImage *norm_factor,
                       int num_smooth_iterations) {
    iftImage  *prev_label, *next_label;
    iftFImage *prev_weight, *next_weight;
    float     *sum, max_membership;
    int       l, i, p, q, r, max_label, iter, tmp;
    iftSet    *prev_frontier = *Frontier, *next_frontier = NULL, *Subtree = NULL, *Processed = NULL;
    iftVoxel  u, v;
    iftAdjRel *A             = fst->A;

    /* Initialization */

    prev_label = iftCopyImage(fst->label);
    next_label = iftCopyImage(prev_label);
    int prev_label_max_val = iftMaximumValue(prev_label);
    sum         = iftAllocFloatArray(prev_label_max_val + 1);
    prev_weight = iftCreateFImage(prev_label->xsize, prev_label->ysize, prev_label->zsize);
    next_weight = iftCreateFImage(next_label->xsize, next_label->ysize, next_label->zsize);

    for (p = 0; p < prev_label->n; p++)
        prev_weight->val[p] = next_weight->val[p] = 1.0;

    /* Smooth frontier and reset its path values */

    for (iter = 0; iter < num_smooth_iterations; iter++) {
        while (prev_frontier != NULL) {
            p = iftRemoveSet(&prev_frontier);
            iftInsertSet(&next_frontier, p);
            u = iftGetVoxelCoord(prev_label, p);

            for (l = 0; l <= prev_label_max_val; l++) {
                sum[l] = 0.0;
            }

            for (i = 1; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(prev_label, v)) {
                    q = iftGetVoxelIndex(prev_label, v);
                    sum[prev_label->val[q]] += prev_weight->val[q] * border_weight->val[q];
                    if (iftBMapValue(inFrontier, q) == 0) /* expand frontier */
                    {
                        if (fst->pred->val[q] != IFT_NIL) {
                            iftInsertSet(&next_frontier, q);
                            iftBMapSet1(inFrontier, q);
                        }
                    }
                }
            }

            for (l = 0; l <= prev_label_max_val; l++)
                sum[l] = sum[l] / norm_factor->val[p];

            max_membership = IFT_INFINITY_FLT_NEG;
            max_label = IFT_NIL;
            for (l = 0; l <= prev_label_max_val; l++) {
                if (sum[l] > max_membership) {
                    max_membership = sum[l];
                    max_label      = l;
                }
            }
            next_label->val[p]  = max_label;
            next_weight->val[p] = sum[max_label];
        }

        prev_frontier = next_frontier;
        next_frontier = NULL;

        for (r = 0; r < prev_label->n; r++) {
            prev_weight->val[r] = next_weight->val[r];
            prev_label->val[r]  = next_label->val[r];
        }
    }


    iftFree(sum);
    iftDestroyFImage(&prev_weight);
    iftDestroyFImage(&next_weight);
    iftDestroyImage(&next_label);
    iftDestroySet(&prev_frontier);
    *Frontier = NULL;

    /* Fix the forest by first making available to be conquered all
       voxels whose labels have changed and their subtrees. */

    for (p = 0; p < fst->label->n; p++) {
        /* If the label has changed */

        if (fst->label->val[p] != prev_label->val[p]) {
            fst->pathval->val[p] = IFT_INFINITY_INT; /* Allow the voxel to be conquered */
            iftInsertSet(&Subtree, p);
            /* Make its subtree available to be conquered as well */
            while (Subtree != NULL) {
                r = iftRemoveSet(&Subtree);
                u = iftGetVoxelCoord(fst->pred, r);

                for (i = 1; i < A->n; i++) {
                    v = iftGetAdjacentVoxel(A, u, i);
                    if (iftValidVoxel(fst->pred, v)) {
                        q = iftGetVoxelIndex(fst->pred, v);
                        if (fst->pred->val[q] == r) /* q is in the subtree of r */
                        {
                            fst->pathval->val[q] = IFT_INFINITY_INT;
                            iftInsertSet(&Subtree, q);
                        }
                    }
                }
            }
        }
    }

    /* Insert in priority queue the seed voxels, which will be voxels in
       the frontier of the available region whose labels are the same of
       the neighbor in the region. The initial path value of such seeds
       must be the current path value. Also insert as seeds in the queue
       the voxels in that region, whose predecessor outside the region
       maintained the old label, as long as they are true seeds (i.e.,
       they also have a neighbor outside the region with the same new
       label of them). In this case, the initial cost is the one of the
       predecessor voxel to avoid backward path propagation. */

    for (p = 0; p < fst->label->n; p++) {
        if (fst->pathval->val[p] == IFT_INFINITY_INT) {
            /* p is in the available region */
            int pp = fst->pred->val[p]; /* predecessor of p */

            u = iftGetVoxelCoord(fst->pred, p);
            fst->pred->val[p] = IFT_NIL;

            for (i = 1; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(fst->pred, v)) {
                    q = iftGetVoxelIndex(fst->pred, v);

                    if (fst->pathval->val[q] !=
                        IFT_INFINITY_INT) /* q is outside the available region */
                    {
                        if (prev_label->val[q] ==
                            prev_label->val[p]) /* the label of q is the new label of p. Then p can be a true seed, if its predecessor is outside the available region. */
                        {
                            if (fst->Q->L.elem[q].color ==
                                IFT_WHITE)  /* q is not in the queue yet */
                            {
                                iftInsertGQueue(&fst->Q, q);
                            }
                            if (fst->label->val[pp] ==
                                prev_label->val[pp]) /* pp is outside the available region. => its label is the  old one, different from the new label of p. */
                            {
                                fst->pathval->val[p] = fst->pathval->val[pp];
                                if (fst->Q->L.elem[p].color == IFT_WHITE)
                                    iftInsertGQueue(&fst->Q, p);
                            }
                        }
                    }
                }
            }
        }
    }


    iftDestroyImage(&fst->label);
    fst->label = prev_label;


    /* execute the IFT to reconstruct the forest under the new labeling
       constraint. This forest is not optimum, since this is a relaxed
       IFT, but it maintains the connectivity between roots and voxels
       of the same label, respecting the filtering process. */

    while (!iftEmptyGQueue(fst->Q)) {
        p = iftRemoveGQueue(fst->Q);
        iftInsertSet(&Processed, p);

        u = iftGetVoxelCoord(fst->img, p);

        for (i = 1; i < A->n; i++) {
            v = iftGetAdjacentVoxel(A, u, i);
            if (iftValidVoxel(fst->img, v)) {
                q = iftGetVoxelIndex(fst->img, v);
                if (fst->Q->L.elem[q].color != IFT_BLACK) {
                    tmp = iftMax(fst->pathval->val[p], fst->img->val[q]);
                    if (tmp < fst->pathval->val[q] || ((fst->pred->val[q] == p))) {
                        if (fst->Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(fst->Q, q);

                        fst->root->val[q]    = fst->root->val[p];
                        fst->pred->val[q]    = p;
                        fst->label->val[q]   = fst->label->val[p];
                        fst->pathval->val[q] = tmp;
                        iftInsertGQueue(&fst->Q, q);
                    }
                }
            }
        }
    }

    while (Processed != NULL) {
        p = iftRemoveSet(&Processed);
        fst->Q->L.elem[p].color = IFT_WHITE;
    }

    /* Verify forest consistency: This can be removed when we decide for */
    /*    the final procedure. */

    for (p = 0; p < fst->label->n; p++) {
        r = p;
        while (fst->pred->val[r] != IFT_NIL) {
            if (fst->label->val[r] != fst->label->val[fst->pred->val[r]]) {
                iftWarning("Incorrect reconstruction of the label map", "iftSmoothFrontier");
                fprintf(stderr, "%d, %d:", iftGetXCoord(fst->label, r),
                        iftGetYCoord(fst->label, r));
                fprintf(stderr, "Label: %d, Pred Label: %d\n", fst->label->val[r],
                        fst->label->val[fst->pred->val[r]]);
                fprintf(stderr, "Cost: %d, Pred cost: %d\n", fst->pathval->val[r],
                        fst->pathval->val[fst->pred->val[r]]);
                fprintf(stderr, "Pred of predecessor: %d\n", fst->pred->val[fst->pred->val[r]]);
            }
            r = fst->pred->val[r];
        }
    }
}


void iftMeansAboveBelowT(const iftImage *img, int T, int *T1, int *T2) {
    long long mean1 = 0, mean2 = 0, nv1 = 0, nv2 = 0;

    for(int p = 0; p < img->n; p++) {
        if (img->val[p] < T){
            mean1 += img->val[p];
            nv1++;
        }
        if (img->val[p] > T){
            mean2 += img->val[p];
            nv2++;
        }
    }

    int delta = (int) ((mean2/nv2) - (mean1/nv1)) /2;
    *T1   = T+delta;
    *T2   = T-delta;

    //printf("T %d T1 %d T2 %d error %d\n",T,*T1,*T2,error);
}

/*---------------------- Public Functions ---------------------------------*/

iftMaskJoinOperation iftMaskJoinOperStrToMaskJoinOper(char *maskJoinOperStr) {
    iftMaskJoinOperation maskJoinOper = 0;

    if(!strcmp(maskJoinOperStr, "union"))
        maskJoinOper = IFT_MASK_JOIN_OPERATION_UNION;
    else if(!strcmp(maskJoinOperStr, "intersection"))
        maskJoinOper = IFT_MASK_JOIN_OPERATION_INTERSECTION;
    
    return maskJoinOper;
}

iftImage *iftBrainGrad(const iftImage *brain_img) {
    int a, b, c;

    b = iftOtsu(brain_img);
    iftMeansAboveBelowT(brain_img, b, &c, &a);
    iftImage *enha  = iftApplySShape(brain_img, a, b, c);
    iftImage *grad  = iftTextGradient(enha);
    iftDestroyImage(&enha);
    iftCopyVoxelSize(brain_img, grad);

    return grad;
}


iftImage *iftApplySShape(const iftImage *img, int a, int b, int c) {
    iftImage *enha = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftCopyVoxelSize(img, enha);

    double weight;
    for (int p = 0; p < img->n; p++) {
        if (img->val[p] <= a){
            weight = 0.0;
        } else {
            if ((img->val[p] > a) && (img->val[p] <= b)) {
                weight = (2.0*((float)img->val[p]-a)*((float)img->val[p]-a)/(((float)c-a)*((float)c-a)));
            } else {
                if ((img->val[p] > b) && (img->val[p] <= c)) {
                    weight = (1.0 - 2.0*((float)img->val[p]-c)*((float)img->val[p]-c)/(((float)c-a)*((float)c-a)));
                } else{
                    weight = 1.0;
                }
            }

        }
        enha->val[p]= (int) (img->val[p] * weight);
    }

    return enha;
}

iftFImage **iftGradientVectorField(iftImage *img, iftImage *mask)
{
  iftAdjRel *A;
  float   dist,gx,gy,gz, g, gmax;
  float   gxCb , gyCb , gzCb, gxCr , gyCr , gzCr;
  int     i,p,q;
  iftVoxel   u,v;
  float     _2sigma2;
  int        dx, dy, dz;
  iftFImage **grad = (iftFImage **)iftAlloc(3,sizeof(iftFImage *));
  char      no_mask=0;

  if (mask==NULL){
    mask = iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftSetImage(mask,1);
    no_mask=1;
  }

  if (iftIs3DImage(mask)){
    A = iftSpheric(sqrtf(3.0));
  }else{
    A = iftCircular(sqrtf(2.0));
  }

  float     *mag =iftAllocFloatArray(A->n), *weight = iftAllocFloatArray(A->n);

  for (int i=0; i < 3; i++) {
    grad[i] = iftCreateFImage(img->xsize,img->ysize,img->zsize);
    grad[i]->dx = img->dx;
    grad[i]->dy = img->dy;
    grad[i]->dz = img->dz;
  }

  iftMaxAdjShifts(A, &dx, &dy, &dz);
  _2sigma2 = 2.0*(dx*dx+dy*dy+dz*dz)/9.0;
  for (i=0; i < A->n; i++){
    mag[i]=sqrtf(A->dx[i]*A->dx[i]+A->dy[i]*A->dy[i]+A->dz[i]*A->dz[i]);
    weight[i]=exp(-mag[i]/_2sigma2)/mag[i];
  }

  if ( !iftIsColorImage(img) ) {
    for (u.z=0; u.z < img->zsize; u.z++)
      for (u.y=0; u.y < img->ysize; u.y++)
        for (u.x=0; u.x < img->xsize; u.x++) {
          p = iftGetVoxelIndex(img,u);
	  if (mask->val[p]!=0){
	    gx = gy = gz = 0.0;
	    for (i=1; i < A->n; i++) {
	      v.x = u.x + A->dx[i];
	      v.y = u.y + A->dy[i];
	      v.z = u.z + A->dz[i];
	      if (iftValidVoxel(img,v)){
		q = iftGetVoxelIndex(img,v);
		dist = img->val[q]-img->val[p];
		gx  += dist*A->dx[i]*weight[i];
		gy  += dist*A->dy[i]*weight[i];
		gz  += dist*A->dz[i]*weight[i];
	      }
	    }
	    grad[0]->val[p]=gx;
	    grad[1]->val[p]=gy;
	    grad[2]->val[p]=gz;
	  }
	}
  }else{ // colored image
    for (u.z=0; u.z < img->zsize; u.z++)
      for (u.y=0; u.y < img->ysize; u.y++)
        for (u.x=0; u.x < img->xsize; u.x++) {
          p = iftGetVoxelIndex(img,u);
	  if (mask->val[p]!=0){
	    gx = gy = gz = 0.0;
	    gxCb = gyCb = gzCb = 0.0;
	    gxCr = gyCr = gzCr = 0.0;
	    for (i=1; i < A->n; i++) {
	      v.x = u.x + A->dx[i];
	      v.y = u.y + A->dy[i];
	      v.z = u.z + A->dz[i];
	      if (iftValidVoxel(img,v)){
		q = iftGetVoxelIndex(img,v);
		dist = img->val[q]-img->val[p];
		gx  += dist*A->dx[i]*weight[i];
		gy  += dist*A->dy[i]*weight[i];
		gz  += dist*A->dz[i]*weight[i];
		dist = img->Cb[q]-img->Cb[p];
		gxCb  += dist*A->dx[i]*weight[i];
		gyCb  += dist*A->dy[i]*weight[i];
		gzCb  += dist*A->dz[i]*weight[i];
		dist = img->Cr[q]-img->Cr[p];
		gxCr  += dist*A->dx[i]*weight[i];
		gyCr  += dist*A->dy[i]*weight[i];
		gzCr  += dist*A->dz[i]*weight[i];
	      }
	    }
	    
	    gmax = gx;
	    g    = gxCb;
	    if (g > gmax)
	      gmax = g;
	    g    = gxCr;
	    if (g > gmax)
	      gmax = g;
	    grad[0]->val[p]=gmax;
	    
	    gmax = gy;
	    g    = gyCb;
	    if (g > gmax)
	      gmax = g;
	    g    = gyCr;
	    if (g > gmax)
	      gmax = g;
	    grad[1]->val[p]=gmax;

	    gmax = gz;
	    g    = gzCb;
	    if (g > gmax)
	      gmax = g;
	    g    = gzCr;
	    if (g > gmax)
	      gmax = g;
	    grad[2]->val[p]=gmax;
	  }
	}
  }
  iftFree(mag);
  iftFree(weight);
  if (no_mask)
    iftDestroyImage(&mask);

  return(grad);
}

iftImage *iftTextGradient(const iftImage *img) {
    float   dist,gx,gy,gz;
    int     i,p,q;
    iftVoxel   u,v;
    iftAdjRel *A=iftSpheric(1.0),*A6=iftSpheric(1.0);
    float   *mg=iftAllocFloatArray(A6->n);
    iftImage   *grad=iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftCopyVoxelSize(img, grad);

    typedef struct _features {
        float *f;
    } Features;

    Features *feat=(Features *)iftAlloc(img->n,sizeof(Features));
    for (p=0; p < img->n; p++)
        feat[p].f = iftAllocFloatArray(A->n);


    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++) {
                p = u.x + img->tby[u.y] + img->tbz[u.z];
                for (i=0; i < A->n; i++) {
                    v.x = u.x + A->dx[i];
                    v.y = u.y + A->dy[i];
                    v.z = u.z + A->dz[i];
                    if (iftValidVoxel(img, v)){
                        q = v.x + img->tby[v.y] + img->tbz[v.z];
                        feat[p].f[i]=(float)img->val[q];///(float)Imax;
                    }
                }
            }

    for (i=0; i < A6->n; i++)
        mg[i]=sqrt(A6->dx[i]*A6->dx[i]+A6->dy[i]*A6->dy[i]+A6->dz[i]*A6->dz[i]);

    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++) {
                p = u.x + img->tby[u.y] + img->tbz[u.z];
                gx = gy = gz = 0.0;
                for (i=1; i < A6->n; i++) {
                    v.x = u.x + A6->dx[i];
                    v.y = u.y + A6->dy[i];
                    v.z = u.z + A6->dz[i];
                    if (iftValidVoxel(img, v)){
                        q = v.x + img->tby[v.y] + img->tbz[v.z];
                        dist = iftSignedManhattanDistance(feat[p].f,feat[q].f,A->n);
                        gx  += dist*A6->dx[i]/mg[i];
                        gy  += dist*A6->dy[i]/mg[i];
                        gz  += dist*A6->dz[i]/mg[i];
                    }
                }
                grad->val[p]=(int)sqrt(gx*gx + gy*gy + gz*gz);//(100000.0*sqrt(gx*gx + gy*gy + gz*gz));
            }


    for (p=0; p < img->n; p++)
        iftFree(feat[p].f);
    iftFree(feat);

    iftFree(mg);
    iftDestroyAdjRel(&A);
    iftDestroyAdjRel(&A6);
    return(grad);
}


iftImage *iftComputeGradient(const iftImage *img, iftGradientAlg grad_alg) {
    iftImage *grad_img = NULL;
    if (grad_alg == IFT_NONE_GRAD) {
        grad_img = iftCopyImage(img);
    }
    else if (grad_alg == IFT_IMAGE_BASINS) {
        grad_img = iftImageBasins(img, NULL);
    }
    else if (grad_alg == IFT_BRAIN_GRAD) {
        grad_img = iftBrainGrad(img);
    }
    else if (grad_alg == IFT_IMAGE_GRAD_MAGNITUDE) {
        grad_img = iftImageGradientMagnitude(img, NULL);
    }

    return grad_img;
}


iftImage *iftCombineGradImages(iftImage *grad_img1, iftImage *grad_img2, float alpha, int max_img_range) {
    if (iftAlmostZero(alpha)) {
        return iftCopyImage(grad_img2);
    }
    else if (iftAlmostZero(1.0 - alpha)) {
        return iftCopyImage(grad_img1);
    }
    else {
        iftImage *norm_grad_img1 = iftNormalize(grad_img1, 0, max_img_range);
        iftImage *norm_grad_img2 = iftNormalize(grad_img2, 0, max_img_range);

        iftImage *combined_grad_img = iftImageLinearCombination(norm_grad_img1, norm_grad_img2, alpha);


        iftDestroyImage(&norm_grad_img1);
        iftDestroyImage(&norm_grad_img2);

        return combined_grad_img;
    }
}


iftImage *iftMEnhanceObject(iftMImage *img, iftLabeledSet *seed, int obj)
{
    iftImage      *objmap=NULL;
    iftDataSet    *Z=NULL;
    iftKnnGraph   *graph=NULL;
    iftLabeledSet *Sobj=NULL, *Sbkg=NULL;
    float          train_perc;

    /* Create dataset of markers */

    Z  = iftMImageSeedsToDataSet(img, seed);

    if (Z->nsamples > 500)
        train_perc   = (float) 500 /Z->nsamples;
    else
        train_perc   = 1.0f;

    iftSelectSupTrainSamples(Z,train_perc);
    graph = iftCreateKnnGraph(Z, (int) (0.3 * Z->ntrainsamples));
    iftUnsupTrain(graph, iftNormalizedCut);

    /* Find object and background roots */

    for(int u = 0; u < graph->nnodes; u++) {
        int s = graph->node[u].sample;
//        int r = graph->node[graph->node[u].root].sample;

        if(Z->sample[s].truelabel > 1) {
            /* needs to be sample.id because s sample index and not pixel index */
            iftInsertLabeledSet(&Sobj, Z->sample[s].id, 1);
        } else if(Z->sample[s].truelabel == 1) {
            iftInsertLabeledSet(&Sbkg, Z->sample[s].id, 0);
        }
    }

    if (Sobj == NULL){
        char msg[200];
        sprintf(msg,"Markers make no clusters for object %d",obj);
        iftWarning(msg,"iftEnhanceObject");
        iftDestroyKnnGraph(&graph);
        iftDestroyDataSet(&Z);
        iftDestroyLabeledSet(&Sbkg);
        return(objmap);
    }

    if (Sbkg == NULL){
        char msg[200];
        sprintf(msg,"Markers make no clusters for background");
        iftWarning(msg,"iftEnhanceObject");
        iftDestroyKnnGraph(&graph);
        iftDestroyDataSet(&Z);
        iftDestroyLabeledSet(&Sobj);
        return(objmap);
    }

    iftDestroyKnnGraph(&graph);

    iftLabeledSet *S = NULL;

    iftConcatLabeledSet(&S, &Sobj);
    iftConcatLabeledSet(&S, &Sbkg);

    iftDestroyLabeledSet(&Sobj);
    iftDestroyLabeledSet(&Sbkg);

    iftDestroyDataSet(&Z);
    Z = iftMImageSeedsToDataSet(img, S);

    iftSetStatus(Z, IFT_TRAIN);
    iftCplGraph *cgraph = iftCreateCplGraph(Z);
    iftSupTrain(cgraph);

    iftDataSet *Zimg = iftMImageToDataSet(img, NULL, 0);

    iftClassifyWithCertaintyValues(cgraph, Zimg);

    objmap = iftDataSetObjectMap(Zimg, NULL, iftNormalizationValue((int) iftMMaximumValue(img, -1)), obj + 1);

    iftDestroyDataSet(&Zimg);
    iftDestroyDataSet(&Z);
    iftDestroyCplGraph(&cgraph);
    iftDestroyLabeledSet(&S);

    return objmap;
}

iftImage *iftEnhanceObject(iftImage *img, iftLabeledSet *seed, int obj)
{
    iftMImage *aux = iftImageToMImage(img, GRAY_CSPACE);
    iftImage* objmap = iftMEnhanceObject(aux, seed, obj);
    iftDestroyMImage(&aux);
    return objmap;
}

iftImage *iftMEnhanceEdges(iftMImage *img, iftAdjRel *A, iftLabeledSet *seed, float alpha)
{
  iftImage      *objmap=NULL,*objbasins=NULL;
  iftImage      *basins=NULL;
  int            p;

  if ((alpha < 0.0)||(alpha > 1.0))
      iftError("Combination factor must be in [0,1]", "iftMEnhanceEdges");

  /* Compute edge enhancement based on image properties */

  basins = iftMImageBasins(img,A);

  if (seed == NULL) 
    return(basins);

  /* Compute edge enhancement based on object properties */

  iftLabeledSet *S=seed;
  int Lmax= IFT_INFINITY_INT_NEG;
  while (S != NULL) {
    if (S->label > Lmax) 
      Lmax = S->label;
    S = S->next;
  }
  Lmax += 1; /* To count with the fact that seed labels start at 0,
		but class labels start at 1 */

  objbasins = iftCreateImage(basins->xsize,basins->ysize,basins->zsize);
  char noobjmap = 1;
 
  for (int l=2; l <= Lmax; l++) {
    objmap = iftMEnhanceObject(img, seed, l);
    if (objmap != NULL) {
      iftImage *aux = iftImageBasins(objmap,A);
      for (p=0; p < aux->n; p++) 
	objbasins->val[p] = iftMax(objbasins->val[p], aux->val[p]);
      iftDestroyImage(&aux);
      iftDestroyImage(&objmap);
      noobjmap=0;
    }
  }
  
  if (noobjmap)
    return(basins);

  /* Combine edge enhancements */ 

  for (p=0; p < img->n; p++) {       
    basins->val[p] = (int)((1.0-alpha)*(float)basins->val[p] +  
			   alpha*(float)objbasins->val[p]);
  }   

  iftDestroyImage(&objbasins); 

  basins->dx = img->dx;
  basins->dy = img->dy;
  basins->dz = img->dz;

  return(basins);
}

iftImage *iftEnhanceEdges(iftImage *img, iftAdjRel *A, iftLabeledSet *seed, float alpha)
{
  iftImage      *objmap=NULL,*objbasins=NULL;
  iftImage      *basins=NULL;
  int            p;

  if ((alpha < 0.0)||(alpha > 1.0))
      iftError("Combination factor must be in [0,1]", "iftEnhanceEdges");

  /* Compute edge enhancement based on image properties */

  basins = iftImageBasins(img,A);

  if (seed == NULL) 
    return(basins);

  /* Compute edge enhancement based on object properties */

  iftLabeledSet *S=seed;
  int Lmax= IFT_INFINITY_INT_NEG;
  while (S != NULL) {
    if (S->label > Lmax) 
      Lmax = S->label;
    S = S->next;
  }
  Lmax += 1; /* To count with the fact that seed labels start at 0,
		but class labels start at 1 */

  objbasins = iftCreateImage(basins->xsize,basins->ysize,basins->zsize);
  char noobjmap = 1;
 
  for (int l=2; l <= Lmax; l++) {
    objmap = iftEnhanceObject(img, seed, l);
    if (objmap != NULL) {
      iftImage *aux = iftImageBasins(objmap,A);
      for (p=0; p < aux->n; p++) 
	objbasins->val[p] = iftMax(objbasins->val[p], aux->val[p]);
      iftDestroyImage(&aux);
      iftDestroyImage(&objmap);
      noobjmap=0;
    }
  }
  
  if (noobjmap)
    return(basins);

  /* Combine edge enhancements */ 

  for (p=0; p < img->n; p++) {       
    basins->val[p] = (int)((1.0-alpha)*(float)basins->val[p] +  
			   alpha*(float)objbasins->val[p]);
  }   

  iftDestroyImage(&objbasins); 

  iftCopyVoxelSize(img,basins);

  return(basins);
}

iftImage *iftWatershedWithPredMap(const iftImage *basins, iftAdjRel *Ain, iftLabeledSet *seeds, iftSet *forbidden, iftImage **pred)
{
    iftImage  *pathval = NULL, *label = NULL, *predecessor = NULL;
    iftGQueue  *Q = NULL;
    int      i, p, q, tmp;
    iftVoxel    u, v;
    iftLabeledSet *S = seeds;

    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(basins))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }
    else A = Ain;

    // Initialization
    pathval  = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
    label = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
    predecessor = iftCreateImage(basins->xsize,basins->ysize,basins->zsize);
    Q     = iftCreateGQueue(iftMaximumValue(basins) + 1, basins->n, pathval->val);
    for (p = 0; p < basins->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
    }

    // sets the forbidden region
    iftSet *F = forbidden;
    while (F != NULL) {
        p   =  F->elem;
        pathval->val[p] = IFT_INFINITY_INT_NEG;
        F = F->next;
    }

    while (S != NULL)
    {
        p = S->elem;
        label->val[p] = S->label;
        pathval->val[p] = 0;
        predecessor->val[p] = IFT_NIL;
        iftInsertGQueue(&Q, p);
        S = S->next;
    }


    // Image Foresting Transform

    while (!iftEmptyGQueue(Q))
    {
        p = iftRemoveGQueue(Q);
        u = iftGetVoxelCoord(basins, p);

        for (i = 1; i < A->n; i++)
        {
            v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(basins, v))
            {
                q = iftGetVoxelIndex(basins, v);
                if (pathval->val[q] > pathval->val[p])
                {
                    tmp = iftMax(pathval->val[p], basins->val[q]);
                    if (tmp < pathval->val[q])  // For this path-value function,
                    {
                        // this implies that q has never
                        // been inserted in Q.
                        label->val[q] = label->val[p];
                        predecessor->val[q] = p;
                        pathval->val[q]  = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }
    iftDestroyGQueue(&Q);
    iftDestroyImage(&pathval);

    iftCopyVoxelSize(basins, label);

    *pred = predecessor;

    return (label);
}
 
iftImage *iftWatershed(const iftImage *basins, iftAdjRel *Ain, iftLabeledSet *seeds, iftSet *forbidden) {
  iftImage  *pathval = NULL, *label = NULL;
  iftGQueue  *Q = NULL;
  int      i, p, q, tmp;
  iftVoxel    u, v;
  iftLabeledSet *S = seeds;

  iftAdjRel *A = NULL;
  if (Ain == NULL) {
    if (iftIs3DImage(basins))
      A = iftSpheric(1.0);
    else A = iftCircular(1.0);
  }
  else A = Ain;

  // Initialization
  pathval  = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  label = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
  Q     = iftCreateGQueue(iftMaximumValue(basins) + 1, basins->n, pathval->val);
  for (p = 0; p < basins->n; p++)
  {
    pathval->val[p] = IFT_INFINITY_INT;
  }

  // sets the forbidden region
  iftSet *F = forbidden;
  while (F != NULL) {
    p   =  F->elem;
    pathval->val[p] = IFT_INFINITY_INT_NEG;
    F = F->next;
  }

  while (S != NULL)
  {
    p = S->elem;
    label->val[p] = S->label;
    pathval->val[p] = 0;
    iftInsertGQueue(&Q, p);
    S = S->next;
  }


  // Image Foresting Transform

  while (!iftEmptyGQueue(Q))
  {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(basins, p);

    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftValidVoxel(basins, v))
      {
        q = iftGetVoxelIndex(basins, v);
        if (pathval->val[q] > pathval->val[p])
        {
          tmp = iftMax(pathval->val[p], basins->val[q]);
          if (tmp < pathval->val[q])  // For this path-value function,
          {
            // this implies that q has never
            // been inserted in Q.
            label->val[q] = label->val[p];
            pathval->val[q]  = tmp;
            iftInsertGQueue(&Q, q);
          }
        }
      }
    }
  }

  if (Ain == NULL) {
    iftDestroyAdjRel(&A);
  }
  iftDestroyGQueue(&Q);
  iftDestroyImage(&pathval);

  iftCopyVoxelSize(basins, label);

  return (label);
}


iftImage *iftWaterCut(iftMImage *mimg, iftAdjRel *Ain, iftLabeledSet *seeds, iftSet *forbidden)
{
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DMImage(mimg))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }
    else A = Ain;

    // Initialization
    float max = iftMMaximumValue(mimg, -1);
    iftImage *pathval  = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *label = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftGQueue *Q     = iftCreateGQueue(sqrt(max * max * mimg->m), mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
    }

    // sets the forbidden region
    iftSet *F = forbidden;
    while (F != NULL) {
        int p   =  F->elem;
        pathval->val[p] = IFT_INFINITY_INT_NEG;
        F = F->next;
    }

    iftLabeledSet *S = seeds;
    while (S != NULL)
    {
        int p = S->elem;
        label->val[p] = S->label;
        pathval->val[p] = 0;
        iftInsertGQueue(&Q, p);
        S = S->next;
    }

    // Image Foresting Transform
    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(mimg, v))
            {
                int q = iftGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    int tmp = iftMax(pathval->val[p], iftMImageDist(mimg, p, q));
//                    int tmp = iftMImageDist(mimg, p, q);

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q]  = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }
    iftDestroyGQueue(&Q);
    iftDestroyImage(&pathval);

    iftCopyVoxelSize(mimg, label);

    return (label);
}



iftImage *iftWatershedMeanCut(const iftImage *basins, iftAdjRel *Ain, const iftLabeledSet *seeds,
                              const iftSet *forbidden, float *result) {
    if (basins == NULL)
        iftError("Gradient Image is NULL", "iftWatershedMeanCut");

    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(basins))
            A = iftSpheric(1.0);
        else A = iftCircular(1.5);
    }
    else A = Ain;

    ////////////////// Initialization //////////////////
    iftImage *pathval = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
    iftImage *label   = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
    iftCopyVoxelSize(basins, label);
    iftGQueue *Q      = iftCreateGQueue(iftMaximumValue(basins)+1, basins->n, pathval->val);

    for (int p = 0; p < basins->n; p++) {
        pathval->val[p] = IFT_INFINITY_INT;
    }

    // sets the forbidden region
    const iftSet *F = forbidden;
    while (F != NULL) {
        int p           = F->elem;
        pathval->val[p] = IFT_INFINITY_INT_NEG;
        F               = F->next;
    }

    // sets the cost of the seeds
    const iftLabeledSet *S = seeds;
    while (S != NULL) {
        int p           = S->elem;
        label->val[p]   = S->label;
        pathval->val[p] = 0;
        iftInsertGQueue(&Q, p);
        S = S->next;
    }
    //////////////////////////////////////////////////////

    // checks if a given boundary voxel was already computed/visited in the Boundary's Mean Gradient computation
    iftBMap *visited = iftCreateBMap(label->n);
    int count        = 0; // counts the number of bourdary voxels

    // Image Foresting Transform
    while (!iftEmptyGQueue(Q)) {
        int p      = iftRemoveGQueue(Q);
        iftVoxel u = iftGetVoxelCoord(basins, p);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(basins, v)) {
                int q = iftGetVoxelIndex(basins, v);

                if (pathval->val[q] > pathval->val[p]) {
                    int tmp = iftMax(pathval->val[p], basins->val[q]);
                    // For this path-value function, this implies that q has never been inserted in Q.
                    if (tmp < pathval->val[q]) {
                        label->val[q] = label->val[p];
                        pathval->val[q]  = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
                // q can be in the Queue or out of it (q has already optimum cost) 
                else {
                    // if q is out of the Queue and does not have the same label from the adjacent p
                    if ((Q->L.elem[q].color == IFT_BLACK) && (label->val[p] != label->val[q])) {
                        // if q was not still visited in the Mean Gradient Computation
                        if (iftBMapValue(visited, q) == 0) {
                            (*result) += basins->val[q];
                            count++;
                            iftBMapSet1(visited, q);
                        }
                    }
                }
            }
        }
    }

    (*result) = (*result) / iftMax(count, 1);

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }
    iftDestroyGQueue(&Q);
    iftDestroyImage(&pathval);
    iftDestroyBMap(&visited);

    return label;
}


/**
@brief Computes a single shot watershed transform and applies relaxation.
 
@param  basins  iftImage    Gradient image created in the client.
@param  A iftAdjRel.     Adjacency Relation structure.
@param  seed iftLabeledSet. Set of labeled seeds from the client.
@param  num_smooth_iterations. Integer with the number of relaxation iterations.
@param  smooth_factor. Float number within [0, 1] with the relaxation factor.

@return void. All forest maps are updated by reference.
*/
iftImage *iftRelaxedWatershed(iftImage *basins, iftAdjRel *A, iftLabeledSet *seed, int num_smooth_iterations, float smooth_factor)
{
  iftImageForest  *fst = NULL;
  iftSmoothBorder *smooth = NULL;
  iftImage *label = NULL;

  fst = iftCreateImageForest(basins, A);
  smooth = iftCreateSmoothBorder(fst->img, fst->A, num_smooth_iterations, smooth_factor);

  iftDiffWatershed(fst, seed, NULL);
  iftRelaxObjects(fst, smooth);

  label = iftCopyImage(fst->label);

  iftDestroyImageForest(&fst);
  iftDestroySmoothBorder(&smooth);

  return label;
}

/**
@brief Computes the watershed transform in a forest structure.
 
Complexity: O(|A| * |V|).

@param  fst  iftImageForest.    Forest structure created with a gradient image.
@param  seed iftLabeledSet.     List of spels with image index and seed label.
@param  removal_markers iftSet. List of spels marked for removal. NULL if empty

@return void. All forest maps are updated by reference.
*/
void iftDiffWatershed(iftImageForest *fst, iftLabeledSet *seed, iftSet * removal_markers)
{
  iftAdjRel *A = fst->A;
  iftGQueue *Q = fst->Q;
  iftVoxel   u, v;
  int        i, p, q, tmp;
  iftSet    *Frontier = NULL;
  iftLabeledSet *S;
  iftBMap   *processed = fst->processed;
  iftImage  *pathval = fst->pathval, *pred = fst->pred, *label = fst->label;
  iftImage  *root = fst->root, *basins = fst->img, *marker = fst->marker;

  iftFillBMap(processed, 0);

  if (removal_markers != NULL)
  {
    Frontier = iftTreeRemoval(fst, removal_markers);
    while (Frontier != NULL)
    {
      p = iftRemoveSet(&Frontier);
      iftInsertGQueue(&Q, p);
    }
  }

  S = seed;
  while (S != NULL)
  {
    p = S->elem;

    if (Q->L.elem[p].color == IFT_GRAY)
    {
      /* p is also a frontier voxel, but the priority is it as a seed. */
      iftRemoveGQueueElem(Q, p);
    }

    label->val[p]   = S->label;
    pathval->val[p] = 0;
    root->val[p]    = p;
    pred->val[p]    = IFT_NIL;
    marker->val[p]  = S->marker;
    iftInsertGQueue(&Q, p);
    S = S->next;
  }

  /* Image Foresting Transform */
  while (!iftEmptyGQueue(Q))
  {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(basins, p);
    iftBMapSet1(processed, p);

    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);
      if (iftValidVoxel(basins, v))
      {
        q = iftGetVoxelIndex(basins, v);

	if (Q->L.elem[q].color != IFT_BLACK){
	  
	  tmp = iftMax(pathval->val[p], basins->val[q]);
	  
	  /* if pred[q]=p then p and q belong to a tie-zone */
	  if ((tmp < pathval->val[q]) || ((pred->val[q] == p)))
	    {
	      if (Q->L.elem[q].color == IFT_GRAY)
		{
		  iftRemoveGQueueElem(Q, q);
		}
	      pred->val[q]     = p;
	      root->val[q]     = root->val[p];
	      label->val[q]    = label->val[p];
	      marker->val[q]   = marker->val[p];
	      pathval->val[q]  = tmp;
	      iftInsertGQueue(&Q, q);
	    }
	}
      }
    }
  }
  
  iftResetGQueue(Q);

}

/**
@brief Applies a relaxation in the label map.
 
This function relies on the bitmap 'processed' inside the iftImageForest
structure. The bitmaps is automatically created with the forest but it
is only updated in iftDiffWatershed function. The object's border is computed
only in the region marked with 1 in the bitmap. To apply the relaxation
to the entire label map, fill the processed bitmap with 1 before invoking
this function.

Besides the object relaxation, this function also corrects any possible
inconsistencies in the forest that could be caused during the relaxation. More
details about forest inconsistency during post-processing filters 
are described in Nikolas Moya's masters dissertation.

@param  fst  iftImageForest updated using iftDiffWatershed.
@param  smooth iftSmoothBorder. Structure with the relaxation parameters.

@return void. Consistent forest with relaxation roots.
*/

void iftRelaxObjects(iftImageForest *fst, iftSmoothBorder *smooth) {
    iftImage  *prev_label  = smooth->prev_label, *next_label = smooth->next_label, *new_label = NULL;
    iftImage  *prev_marker = smooth->prev_marker, *next_marker = smooth->next_marker;
    iftImage  *new_marker  = NULL;
    iftFImage *prev_weight = smooth->prev_weight, *next_weight = smooth->next_weight;
    float     *sum, *mk, max_membership;
    int       l, i, p, q, r;
    int       borderSize   = 0, max_label, iter;
    iftBMap   *processed   = fst->processed, *inBorder = iftCreateBMap(fst->img->n);
    iftSet    *prev_border = NULL, *next_border = NULL, *dilated_border = NULL;
    iftSet    *tmp         = NULL;
    iftVoxel  u, v;
    iftAdjRel *A           = fst->A;
    char      color;

    /* Initialization */
    int fst_label_maxval = iftMaximumValue(fst->label);
    sum = iftAllocFloatArray(fst_label_maxval + 1);
    mk  = iftAllocFloatArray(fst_label_maxval + 1);

    /* Detect border of the processed regions */
    for (p = 0; p < processed->n; p++) {
        /* Reset data structures for object boundary relaxation */
        prev_weight->val[p] = next_weight->val[p] = 1.0;
        prev_label->val[p]  = next_label->val[p]  = fst->label->val[p];
        prev_marker->val[p] = next_marker->val[p] = fst->marker->val[p];

        u = iftGetVoxelCoord(fst->img, p);
        if (iftBMapValue(processed, p)) {
            for (i = 1; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(prev_label, v)) {
                    q = iftGetVoxelIndex(fst->img, v);
                    if (fst->label->val[p] != fst->label->val[q]) {
                        if (iftBMapValue(inBorder, p) == 0) {
                            iftBMapSet1(inBorder, p);
                            iftInsertSet(&prev_border, p);
                            borderSize++;
                        }
                        if (iftBMapValue(inBorder, q) == 0) {
                            iftBMapSet1(inBorder, q);
                            iftInsertSet(&prev_border, q);
                            borderSize++;
                        }
                    }
                }
            }
        }
    }
    int prev_label_max_val = fst_label_maxval;
    if (borderSize == 0) {
        iftFree(sum);
        iftDestroyBMap(&inBorder);
        iftWarning("There is no need for boundary smoothing", "iftRelaxObjects");
        return;
    }

    /* Relax object boundaries in border region while it is dilated */
    for (iter = 0; iter < smooth->smooth_iterations; iter++) {

        next_border = NULL;

        tmp = prev_border;
        while (prev_border != NULL) {
            //p = iftRemoveSet(&prev_border);
            p = prev_border->elem;
            iftInsertSet(&next_border, p);
            u = iftGetVoxelCoord(prev_label, p);

            for (l = 0; l <= prev_label_max_val; l++) {
                sum[l] = 0.0;
                mk[l] = IFT_NIL;
            }

            for (i = 1; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(prev_label, v)) {
                    q                      = iftGetVoxelIndex(prev_label, v);
                    sum[prev_label->val[q]] += prev_weight->val[q] * smooth->border_weight->val[q];
                    mk[prev_label->val[q]] = prev_marker->val[q];

                    if (iftBMapValue(inBorder, q) == 0) /* expand border */
                    {
                        iftInsertSet(&next_border, q);
                        iftBMapSet1(inBorder, q);
                    }
                }
            }

            for (l = 0; l <= prev_label_max_val; l++)
                sum[l] = sum[l] / smooth->norm_factor->val[p];

            max_membership = IFT_INFINITY_FLT_NEG;
            max_label = IFT_NIL;
            for (l = 0; l <= prev_label_max_val; l++) {
                if (sum[l] > max_membership) {
                    max_membership = sum[l];
                    max_label      = l;
                }
            }
            next_label->val[p]  = max_label;
            next_weight->val[p] = sum[max_label];
            next_marker->val[p] = mk[max_label];

            prev_border = prev_border->next;
        }

        prev_border = tmp;
        while (prev_border != NULL) {
            r = iftRemoveSet(&prev_border);
            prev_weight->val[r] = next_weight->val[r];
            prev_label->val[r]  = next_label->val[r];
            prev_marker->val[r] = next_marker->val[r];
        }
        prev_border = next_border;
    }

    iftFree(sum);
    iftFree(mk);
    iftDestroyBMap(&inBorder);
    new_label      = next_label;
    new_marker     = next_marker;
    dilated_border = next_border;

    /* Fix possible segmentation inconsistencies */
    while (dilated_border != NULL) {
        p = iftRemoveSet(&dilated_border);
        /* If this is the case, find the new root defined by relaxation */
        if (new_label->val[p] != fst->label->val[p]) {
            color = IFT_WHITE;
            iftFindRelaxationRoot1(fst, new_label, p, &r, &color);
            if (color == IFT_BLACK) {
                iftPropagateLabelMarkerAndRootToSubtree(fst, A, new_label, new_marker, r);
                fst->pred->val[r] = IFT_NIL;
            }
        }
    }

    /* Fix possible segmentation inconsistencies */
    // while (dilated_border != NULL)
    // {
    //   p = iftRemoveSet(&dilated_border);
    //   if (new_label->val[p] != old_label->val[p])
    //   {
    //     /* the label of p has changed */
    //     /* Find the new root defined by relaxation */
    //     r = iftFindRelaxationRoot(fst->pred, new_label, old_label, p, &valid);
    //     if (valid)
    //     {
    //       iftInsertSet(&new_roots, r);
    //       iftPropagateLabelMarkerAndRootToSubtree(fst, A, new_label->val[r], new_marker->val[r], r);
    //     }
    //   }
    // }

    // while (new_roots != NULL)
    // {
    //   p = iftRemoveSet(&new_roots);
    //   fst->pred->val[p] = IFT_NIL;
    // }

}


void iftWaterGrayForest(iftImageForest *fst, iftImage *marker)
{
    iftImage   *label=fst->label, *basins=fst->img, *pathval=fst->pathval;
    iftImage   *pred=fst->pred, *root=fst->root;
    iftGQueue  *Q=fst->Q;
    iftAdjRel  *A=fst->A;
    int         i,p,q,l=1,tmp;
    iftVoxel    u,v;

    // Initialization

    for (p=0; p < pathval->n; p++) {
        pathval->val[p] = marker->val[p] + 1;
        iftInsertGQueue(&Q,p);
    }

    // Image Foresting Transform

    while(!iftEmptyGQueue(Q)) {
        p=iftRemoveGQueue(Q);

        if (pred->val[p] == IFT_NIL) { // root voxel
            pathval->val[p]  -= 1;
            label->val[p]=l; l++;
        }

        u = iftGetVoxelCoord(basins,p);

        for (i=1; i < A->n; i++){
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(basins,v)){
                q = iftGetVoxelIndex(basins,v);
                if (pathval->val[q] > pathval->val[p]){
                    tmp = iftMax(pathval->val[p], basins->val[q]);
                    if (tmp < pathval->val[q]){
                        iftRemoveGQueueElem(Q,q);
                        label->val[q]      = label->val[p];
                        root->val[q]       = root->val[p];
                        pred->val[q]       = p;
                        pathval->val[q]    = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

}

iftImage *iftWaterGray(iftImage *basins, iftImage *pathval, iftAdjRel  *A)
{
  iftImage   *label=NULL;
  iftGQueue  *Q=NULL;
  int         i,p,q,l=1,tmp;
  iftVoxel    u,v;
 
  // Initialization 
  
  label   = iftCreateImage(basins->xsize,basins->ysize,basins->zsize);
  Q       = iftCreateGQueue(iftMaximumValue(pathval)+2,pathval->n,pathval->val);

  for (p=0; p < basins->n; p++) {
    pathval->val[p] += 1;
    label->val[p]= IFT_NIL;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    if (label->val[p] == IFT_NIL) { // root voxel
      pathval->val[p]  -= 1;
      label->val[p]=l; l++;
    }

    basins->val[p] = pathval->val[p]; // set the reconstruction value

    u = iftGetVoxelCoord(basins,p);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(basins,v)){	
	q = iftGetVoxelIndex(basins,v);
	if (Q->L.elem[q].color != IFT_BLACK){
	  tmp = iftMax(pathval->val[p], basins->val[q]);
	  if (tmp < pathval->val[q]){ 
	    iftRemoveGQueueElem(Q,q);
	    label->val[q]      = label->val[p];
	    pathval->val[q]    = tmp;
	    iftInsertGQueue(&Q, q);
	  }
	}
      }
    }
  }
  
  iftDestroyGQueue(&Q);
  iftCopyVoxelSize(basins,label);

  return(label);
}


iftImage *iftDualWaterGray(iftImage *domes, iftImage *pathval, iftAdjRel  *A)
{
  iftImage   *label=NULL;
  iftGQueue  *Q=NULL;
  int         i,p,q,l=1,tmp;
  iftVoxel    u,v;
 
  // Initialization 
  
  for (p=0; p < domes->n; p++){ /* make sure that the domes image and
				the path value map do not contain
				zeros, due to the reason below. */ 
    pathval->val[p] += 1;
    domes->val[p]   += 1;
  }

  label   = iftCreateImage(domes->xsize,domes->ysize,domes->zsize);
  Q       = iftCreateGQueue(iftMaximumValue(domes)+1,domes->n,pathval->val);
  iftSetRemovalPolicy(Q, MAXVALUE);

  for (p=0; p < domes->n; p++) {
    pathval->val[p] = pathval->val[p]-1; /* this will never be
					  negative due to the above
					  operation, so we can process
					  original regions with value zero */

    label->val[p]= IFT_NIL;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    if (label->val[p] == IFT_NIL) { // root voxel
      pathval->val[p]  += 1;
      label->val[p]=l; l++;
    }

    domes->val[p] = pathval->val[p]-1; // set the reconstruction value

    u  = iftGetVoxelCoord(domes,p);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(domes,v)){	
	q = iftGetVoxelIndex(domes,v);
	if (pathval->val[q] < pathval->val[p]){
	  tmp = iftMin(pathval->val[p], domes->val[q]);
	  if (tmp > pathval->val[q]){ 
	    iftRemoveGQueueElem(Q,q);
	    label->val[q]      = label->val[p];
	    pathval->val[q]    = tmp;
	    iftInsertGQueue(&Q, q);
	  }
	}
      }
    }
  }
  
  iftDestroyGQueue(&Q);
  iftCopyVoxelSize(domes,label);

  return(label);
}

iftImage *iftDualWaterGrayOnMask(iftImage *domes, iftImage *marker, iftImage *mask, iftAdjRel  *A, iftSet **roots)
{
    iftImage   *label=NULL;
    iftGQueue  *Q=NULL;
    int         i,p,q,l=1,tmp;
    iftVoxel    u,v;

    // Initialization
    for (p=0; p < domes->n; p++)
    {
        if(mask->val[p] != 0) {
            /* make sure that the domes image and
                the path value map do not contain
                zeros, due to the reason below. */
            marker->val[p] += 1;
            domes->val[p] += 1;
        }
    }

    label   = iftCreateImage(domes->xsize,domes->ysize,domes->zsize);
    Q       = iftCreateGQueue(iftMaximumValue(domes)+1, domes->n, marker->val);
    iftResetGQueue(Q);
    iftSetRemovalPolicy(Q, MAXVALUE);

    for (p=0; p < domes->n; p++)
    {
        if(mask->val[p] != 0)
        {
            marker->val[p] = marker->val[p] - 1; /* this will never be
			negative due to the above
			operation, so we can process
			original regions with value zero */

            label->val[p]=IFT_NIL;
            iftInsertGQueue(&Q,p);
        }
    }

    // Image Foresting Transform

    while(!iftEmptyGQueue(Q))
    {
        p=iftRemoveGQueue(Q);

        if (label->val[p]==IFT_NIL)   // root voxel
        {
            if(roots != NULL) iftInsertSet(roots, p);
            marker->val[p]  += 1;
            label->val[p]=l;
            l++;
        }

        domes->val[p] = marker->val[p]; // set the reconstruction value

        u  = iftGetVoxelCoord(domes,p);

        for (i=1; i < A->n; i++)
        {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(domes,v))
            {
                q = iftGetVoxelIndex(domes,v);
                if (mask->val[q] != 0 && marker->val[q] < marker->val[p])
                {
                    tmp = iftMin(marker->val[p], domes->val[q]);
                    if (tmp > marker->val[q])
                    {
                        iftRemoveGQueueElem(Q,q);
                        label->val[q]      = label->val[p];
                        marker->val[q]    = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftDestroyGQueue(&Q);
    iftCopyVoxelSize(domes,label);

    return(label);
}

iftImage *iftWaterDist(iftImage *dist, iftImage *label, int H, iftAdjRel  *A)
{
  iftImage   *nlabel=NULL, *pathval=NULL;
  iftGQueue  *Q=NULL;
  int         i,p,q,l=1,tmp;
  iftVoxel    u,v;
 
   
  // Initialization 
  
  nlabel    = iftCreateImage(dist->xsize,dist->ysize,dist->zsize);
  pathval   = iftCreateImage(dist->xsize,dist->ysize,dist->zsize);
  Q         = iftCreateGQueue(iftMaximumValue(dist)+H+1,dist->n,pathval->val);
  iftSetRemovalPolicy(Q, MAXVALUE);

  for (p=0; p < dist->n; p++){ 
    if (label->val[p]==0)
      pathval->val[p]= IFT_INFINITY_INT;
    else{      
      nlabel->val[p]  = IFT_NIL;
      dist->val[p]   += H;
      pathval->val[p] = dist->val[p]-H; 
      iftInsertGQueue(&Q,p);
    }
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    if (nlabel->val[p] == IFT_NIL) { // root voxel
      pathval->val[p]  += H;
      nlabel->val[p]=l; l++;
    }
    u = iftGetVoxelCoord(dist,p);
    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(dist,v)){	
	q = iftGetVoxelIndex(dist,v);
	if (label->val[q] == label->val[p]){
	  if (pathval->val[q] < pathval->val[p]){
	    tmp = iftMin(pathval->val[p], dist->val[q]);
	    if (tmp > pathval->val[q]){ 
	      iftRemoveGQueueElem(Q,q);
	      nlabel->val[q]     = nlabel->val[p];
	      pathval->val[q]    = tmp;
	      iftInsertGQueue(&Q, q);
	    }
	  }
	}
      }
    }
  }

  for (p=0; p < dist->n; p++){ // restore the original distance value
    if (label->val[p]>0)
      dist->val[p] -= H;
  }


  iftDestroyGQueue(&Q);
  iftDestroyImage(&pathval);

  iftCopyVoxelSize(dist,nlabel);

  return(nlabel);
}

int iftCumHistogramThres(iftImage *img, float perc) {
    int    img_max_val = iftMaximumValue(img);
    double hist[img_max_val + 1];
    int    p, i;

    if ((perc <= 0.0) && (perc >= 1.0))
        iftError("Invalid area percentage", "iftCumHistogramThres");


    for (i = 0; i <= img_max_val; i++)
        hist[i] = 0;

    for (p = 0; p < img->n; p++)
        hist[img->val[p]]++;

    for (i = 0; i <= img_max_val; i++)
        hist[i] /= img->n;

    for (i = 1; i <= img_max_val; i++) {
        hist[i] += hist[i - 1];
        if (hist[i] > perc) {
            return (i - 1);
        }
    }
    return (i - 1);
}


iftImage *iftThreshold(const iftImage *img, int lowest, int highest, int value) {
    iftImage *bin = iftCreateImageFromImage(img);

    for (int p = 0; p < img->n; p++) 
        if ((img->val[p] >= lowest) && (img->val[p] <= highest))
            bin->val[p] = value;
        else bin->val[p] = 0;

    return bin;
}


iftImage *iftFThreshold(const iftFImage *img, float lowest, float highest, int value) {
    iftImage *bin = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftCopyVoxelSize(img, bin);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++) {
        if ((img->val[p] >= lowest) && (img->val[p] <= highest))
            bin->val[p] = value;
        else
            bin->val[p] = 0;
  } 

  return bin;
}

// It is assuming only non-negative intensity values. 
int iftOtsu(const iftImage *img) {
    double *hist                           = NULL;
    double p1, p2, m1, m2, s1, s2, J, Jmax = -1.0;
    int    i, T, Topt                      = 0;

    int img_max_val = iftMaximumValue(img);
    hist = iftAllocDoubleArray(img_max_val + 1);

    for (i = 0; i < img->n; i++)
        hist[img->val[i]] += 1.0;

    for (i = 0; i <= img_max_val; i++) {
        hist[i] /= img->n;
    }

    for (T = 1; T < img_max_val; T++) {
        p1 = 0.0;
        for (i = 0; i <= T; i++)
            p1 += hist[i];
        p2     = 1.0 - p1;
        if ((p1 > 0.0) && (p2 > 0.0)) {
            m1 = 0.0;
            for (i = 0; i <= T; i++)
                m1 += hist[i] * i;
            m1 /= p1;
            m2     = 0.0;
            for (i = T + 1; i <= img_max_val; i++)
                m2 += hist[i] * i;
            m2 /= p2;
            s1     = 0.0;
            for (i = 0; i <= T; i++)
                s1 += hist[i] * (i - m1) * (i - m1);
            s1 /= p1;
            s2     = 0.0;
            for (i = T + 1; i <= img_max_val; i++)
                s2 += hist[i] * (i - m2) * (i - m2);
            s2 /= p2;
            J      = (p1 * p2 * (m1 - m2) * (m1 - m2)) / (p1 * s1 + p2 * s2);
        } else {
            J = 0.0;
        }
        if (J > Jmax) {
            Jmax = J;
            Topt = T;
        }
    }

    iftFree(hist);
    return (Topt);
}

int iftOtsuInRegion(iftImage *img, iftImage *mask)
{
  double *hist=NULL;
  double p1,p2,m1,m2,s1,s2,J,Jmax=-1.0;
  int p,i,T,Topt=0, Imax,nelems=0;
  
  Imax = IFT_INFINITY_INT_NEG;
  for (p=0; p < img->n; p++) 
    if (mask->val[p] != 0){
      nelems++;
      if (img->val[p] > Imax)
	Imax = img->val[p];
    }

  hist = iftAllocDoubleArray(Imax+1);

  for (p=0; p < img->n; p++) 
    if (mask->val[p] != 0)
      hist[img->val[p]]+=1.0;
  
  for (i=0; i <= Imax; i++) {
    hist[i] /= nelems;
  }

  for (T=1; T < Imax; T++){
    p1 = 0.0;
    for (i=0; i <= T; i++) 
      p1 += hist[i];
    p2 = 1.0 - p1;
    if ((p1 > 0.0)&&(p2 > 0.0)){
      m1 = 0.0;
      for (i=0; i <= T; i++) 
	m1 += hist[i]*i;
      m1 /= p1;
      m2 = 0.0;
      for (i=T+1; i <= Imax; i++) 
	m2 += hist[i]*i;
      m2 /= p2;
      s1 = 0.0;
      for (i=0; i <= T; i++) 
	s1 += hist[i]*(i-m1)*(i-m1);
      s1 /= p1;
      s2 = 0.0;
      for (i=T+1; i <= Imax; i++) 
	s2 += hist[i]*(i-m2)*(i-m2);
      s2 /= p2;
      J = (p1*p2*(m1-m2)*(m1-m2))/(p1*s1+p2*s2);
    }else{
      J = 0.0;      
    }
    if (J > Jmax){
      Jmax = J;
      Topt = T;
    }
  }

  iftFree(hist);
  return(Topt);
}


/* Returns a weight image for object smoothing (relaxed IFT). */

iftFImage *iftSmoothWeightImage(const iftImage *basins, float beta) {
    iftFImage *weight = iftCreateFImage(basins->xsize,basins->ysize,basins->zsize);

    #pragma omp parallel for
    for (int p = 0; p < basins->n; p++)
        weight->val[p] = 1.0 / (1.0 + (beta * basins->val[p])); // pseudo inverse
    

    return weight;
}


iftImage *iftFastSmoothObjects(const iftImage *labelIn, const iftFImage *weight, int n_iters) {
  iftAdjRel *A;
  iftImage  *label2;
  iftFImage *flabel1, *flabel2;
  float     *sum, max_membership;
  int       l, i, p, q, max_label, iter;
  iftVoxel  u, v;
  iftSet    *Fprev = NULL, *Fnext = NULL;
  iftBMap   *inFrontier;
  iftFImage *norm_factor;
  
  if (n_iters < 1)
    iftError("Invalid number of iterations", "iftFastSmoothObjects");
  
  if (labelIn->zsize == 1) // 2D
    A = iftCircular(sqrtf(2.0));
  else // 3D
    A = iftSpheric(sqrtf(3.0));
  
  norm_factor = iftWeightNormFactor(weight, A);
  int label1_max_val = iftMaximumValue(labelIn);
  sum        = iftAllocFloatArray(label1_max_val + 1);
  label2     = iftCopyImage(labelIn);
  flabel1    = iftCreateFImage(labelIn->xsize, labelIn->ysize, labelIn->zsize);
  flabel2    = iftCreateFImage(label2->xsize, label2->ysize, label2->zsize);
  inFrontier = iftCreateBMap(labelIn->n);
  
  for (p = 0; p < labelIn->n; p++)
    flabel1->val[p] = flabel2->val[p] = 1.0;
  
  /* Find initial frontier set: this implementation is more
     efficient when the background is large. */
  
  for (int z = 0; z < labelIn->zsize; z++)
    for (int y = 0; y < labelIn->ysize; y++)
      //#pragma omp parallel for shared(labelIn,A) private(p,i,v,u)
      for (int x = 0; x < labelIn->xsize; x++) {
	u.x=x; u.y=y; u.z=z;
	p = iftGetVoxelIndex(labelIn, u);
	if (labelIn->val[p] > 0) {
	  for (i = 1; i < A->n; i++) {
	    v.x = u.x + A->dx[i];
	    v.y = u.y + A->dy[i];
	    v.z = u.z + A->dz[i];
	    if (iftValidVoxel(labelIn, v)) {
	      q = iftGetVoxelIndex(labelIn, v);
	      if (labelIn->val[q] != labelIn->val[p]) {
		
		//                              #pragma omp critical
		//                              {
		if (iftBMapValue(inFrontier, p) == 0) {
		  iftInsertSet(&Fprev, p);
		  iftBMapSet1(inFrontier, p);
		}
		if (iftBMapValue(inFrontier, q) == 0) {
		  iftInsertSet(&Fprev, q);
		  iftBMapSet1(inFrontier, q);
		}
		//                              }
	      }
	    }
	  }
	}
      }
  
  /* Smooth objects */
  
  for (iter = 0; iter < n_iters; iter++) {
    //      printf("Processing iteration %d\n",iter+1);
    
    while (Fprev != NULL) {
      
      p = iftRemoveSet(&Fprev);
      
      iftInsertSet(&Fnext, p);
      
      u.x = iftGetXCoord(labelIn, p);
      u.y = iftGetYCoord(labelIn, p);
      u.z = iftGetZCoord(labelIn, p);
      
      for (l = 0; l <= label1_max_val; l++)
	sum[l] = 0.0;
      
      for (i = 1; i < A->n; i++) {
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	v.z = u.z + A->dz[i];
	if (iftValidVoxel(labelIn, v)) {
	  q = iftGetVoxelIndex(labelIn, v);
	  sum[labelIn->val[q]] += flabel1->val[q] * weight->val[q];
	  if (iftBMapValue(inFrontier, q) == 0) { // expand frontier set
	    iftInsertSet(&Fnext, q);
	    iftBMapSet1(inFrontier, q);
	  }
	}
      }
      
      for (l = 0; l <= label1_max_val; l++) {
	sum[l] = sum[l] / norm_factor->val[p];
      }
      
      max_membership = IFT_INFINITY_FLT_NEG;
      max_label = IFT_NIL;
      for (l = 0; l <= label1_max_val; l++)
	if (sum[l] > max_membership) {
	  max_membership = sum[l];
	  max_label      = l;
	}
      
      label2->val[p]  = max_label;
      flabel2->val[p] = sum[max_label];
      
    }
    
    Fprev = Fnext;
    Fnext = NULL;
    
    for (p = 0; p < flabel1->n; p++) {
      flabel1->val[p] = flabel2->val[p];
      labelIn->val[p]  = label2->val[p];
    }
  }

  
  iftFree(sum);
  iftDestroyFImage(&flabel1);
  iftDestroyFImage(&flabel2);
  iftDestroyFImage(&norm_factor);
  iftDestroyAdjRel(&A);
  iftDestroySet(&Fprev);
  iftDestroyBMap(&inFrontier);
  
  iftCopyVoxelSize(labelIn, label2);
  
  return (label2);
}



iftImage *iftBinarizeByOPF(iftImage *orig, iftImage *enha, int init_thres, float train_perc)
{
  iftCplGraph *graph;
  iftDataSet  *Z, *Z1;
  int i, j, o, *value, *index, N=(int)(train_perc*enha->n);
  iftImage *bin; 

  /* Sort enhanced image values from the highest to the lowest
     values */
 
  value = iftAllocIntArray(enha->n);
  index = iftAllocIntArray(enha->n);

  for (i=0; i < enha->n; i++) {
    value[i] = enha->val[i];
    index[i] = i;
  }

  iftBucketSort(value, index, enha->n, IFT_INCREASING);

  /* Create image dataset for classification and select the training
     samples. The initial classes are determined by the initial
     threshold. Then select N training samples around it in order to
     reclassify samples by OPF. */

  Z = iftImageToDataSet(orig);
  if (Z->nfeats==3){ // for color images only
    Z->alpha[0] = 0.20;
    Z->alpha[1] = 1.00;
    Z->alpha[2] = 1.00;
  }


  Z->nclasses = 2;
  for (i=0; i < Z->nsamples; i++) {
    if (enha->val[Z->sample[i].id] <= init_thres)
      Z->sample[i].truelabel = 1;
    else
      Z->sample[i].truelabel = 2;
  }

  o = 0;
  for (i=0; i < enha->n-1; i++) 
    if ((value[i] <= init_thres)&&(value[i+1]>init_thres)){
      o = i; 
      break;
    }

  i = o; j = 0;
  while ((i >= 0)&&(j < N/2)){
    iftHasSampleStatus(Z->sample[index[i]], IFT_TRAIN);
    Z->ntrainsamples++;
    j++; i--;
  }
  i = o+1; 
  while ((i < enha->n)&&(j < N)){
    iftSetSampleStatus(&Z->sample[index[i]], IFT_TRAIN);
    Z->ntrainsamples++;
    j++; i++;
  }
  iftFree(value);
  iftFree(index);
  Z1  = iftExtractSamples(Z,IFT_TRAIN);  

  /* Train and execute the OPF classifier */

  iftSelectSupTrainSamples(Z1,1.0);
  graph = iftCreateCplGraph(Z1);    
  iftSupTrain(graph);               
  iftClassify(graph,Z);             
  iftDestroyCplGraph(&graph);
  iftDestroyDataSet(&Z1);
  
  /* Copy object voxels */
    
  bin = iftCreateImage(enha->xsize,enha->ysize,enha->zsize);
  for (i=0; i < enha->n; i++) {
    if (Z->sample[i].label == 1)
      bin->val[i] = 255;
  }

  iftDestroyDataSet(&Z);
  iftCopyVoxelSize(orig,bin);

  return(bin);
    
}

iftImage *iftEnhanceWhenCrIsGreaterThanCb(iftImage *img)
{
  iftImage *enha=iftCreateImage(img->xsize, img->ysize, img->zsize);
  int p;

  if (img->Cb == NULL)
      iftError("Image must be colored", "iftEnhanceWhenCrIsGreaterThanCb");

  for(p=0; p < img->n; p++) 
    if (img->Cr[p] > img->Cb[p]) 
      enha->val[p] = (int)(255.0*((float)img->Cr[p]-(float)img->Cb[p])/
			   ((float)img->Cr[p]+(float)img->Cb[p]));

  if (iftMaximumValue(enha)==0){
    iftDestroyImage(&enha);
    iftWarning("No differences between color channels","iftEnhanceWhenCrIsGreaterThanCb");
    enha = iftCopyImage(img);
  }
  
  return(enha);
}

iftImage  *iftWatershedOnVoxelDist(iftDataSet *dataset, iftAdjRel *A, iftLabeledSet *seed){
  iftImage *img = (iftImage *)dataset->ref_data;

	float *pathval = iftAllocFloatArray(img->n);
	iftFHeap *H = iftCreateFHeap(img->n, pathval);
	iftSetRemovalPolicyFHeap(H,MINVALUE);

	iftImage *label = iftCreateImage(img->xsize,img->ysize,img->zsize);

	int p;
	for(p = 0; p < img->n; p++){
		pathval[p] = IFT_INFINITY_FLT;
	}

	iftLabeledSet *S = seed;
	while(S != NULL){
		p = S->elem;
		label->val[p] = S->label;
		pathval[p] = 0;

		iftInsertFHeap(H,p);

		S = S->next;
	}

	while(!iftEmptyFHeap(H)){
		p = iftRemoveFHeap(H);

		iftVoxel u = iftGetVoxelCoord(img,p);

		int i;
		for(i = 1; i < A->n; i++){
			iftVoxel v = iftGetAdjacentVoxel(A,u,i);

			if(iftValidVoxel(img,v)){
				int q = iftGetVoxelIndex(img,v);

				if(H->color[q] != IFT_BLACK){
					//Rounded squared distance
					float weight = dataset->iftArcWeight(dataset->sample[p].feat,dataset->sample[q].feat, dataset->alpha,dataset->nfeats);

					float tmp = iftMax(pathval[p], weight);
					if (tmp < pathval[q]){
						label->val[q] = label->val[p];
						pathval[q] = tmp;

						if(H->color[q] == IFT_WHITE)
							iftInsertFHeap(H,q);
						else
							iftGoUpFHeap(H,H->pos[q]);
					}
				}
			}
		}
	}

	iftFree(pathval);
	iftDestroyFHeap(&H);

	return label;
}

/**
@brief It transforms an image of borders into an image of regions. The
regions are the interior of the closed contours and the lines
disappear.
@author Alexandre Falcao. 
@date   4/12/2015
@param  border: Binary image with value 0 for background and value >= 1 for border voxels.
@return An image of regions labeled from 1 to N. 
*/

iftImage* iftBorderImageToLabelImage(iftImage* border)
{
  iftAdjRel *A; 
  iftImage  *marker = iftAddValue(border,1);
  iftImage  *label;
  
  if (iftIs3DImage(border))
    A = iftSpheric(sqrtf(3.0));
  else
    A = iftCircular(1);

  label = iftWaterGray(border, marker, A);

  iftDestroyAdjRel(&A);
  iftDestroyImage(&marker);

  return(label);
}

iftImage *iftAboveAdaptiveThreshold(iftImage *img, iftImage *mask, iftAdjRel *A, float perc, int niters, int value)
{
  iftVoxel u,v;
  int i,it,p,q,n;
  float mean;
  iftImage  *prev  = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage  *next  = iftCreateImage(img->xsize,img->ysize,img->zsize);


  /* Repeat the following process during a certain number of
     iterations (niters), excluding previously selected voxels: For
     each voxel p, compute the mean value among its adjacent
     voxels. If the intensity of p is above a percentile of the mean
     value, then p is selected (i.e., assign value to the output
     intensity of p). */
  
  if (mask != NULL){
    for (it=1; it <= niters; it++) {
      for (p=0; p < img->n; p++) {
	if (mask->val[p]!=0){
	  u = iftGetVoxelCoord(img,p);
	  if (prev->val[p] == 0) {
	    mean = 0.0; n=0;
	    for (i=1; i < A->n; i++) {
	      v = iftGetAdjacentVoxel(A,u,i);
	      if (iftValidVoxel(img,v)){
		q     = iftGetVoxelIndex(img,v);
		if ((prev->val[q]==0)&&(mask->val[q]!=0)){
		  mean += img->val[q];
		  n++;
		}
	      }
	    }
	    if (n > 0) mean /= n;	
	    if (img->val[p] > perc*mean){
	      next->val[p] = value;
	    }
	  }
	}
      }

      for (p=0; p < img->n; p++) {
	prev->val[p]  = iftMax(prev->val[p], next->val[p]);
	next->val[p]  = 0;
      }
    }
  } else {
    for (it=1; it <= niters; it++) {
      for (p=0; p < img->n; p++) {
	u = iftGetVoxelCoord(img,p);
	if (prev->val[p] == 0) {
	  mean = 0.0; n=0;
	  for (i=1; i < A->n; i++) {
	    v = iftGetAdjacentVoxel(A,u,i);
	    if (iftValidVoxel(img,v)){
	      q     = iftGetVoxelIndex(img,v);
	      if (prev->val[q]==0){
		mean += img->val[q];
		n++;
		}
	      }
	    }
	  if (n > 0) mean /= n;	
	  if (img->val[p] > perc*mean){
	    next->val[p] = value;
	  }
	}
      }

      for (p=0; p < img->n; p++) {
	prev->val[p]  = iftMax(prev->val[p], next->val[p]);
	next->val[p]  = 0;
      }
    }
  }
 
  iftDestroyImage(&next);

  return(prev);
}

iftImage *iftBelowAdaptiveThreshold(iftImage *img, iftImage *mask, iftAdjRel *A, float perc, int niters, int value)
{
  iftVoxel u,v;
  int i,it,p,q,n;
  float mean;
  iftImage  *prev  = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage  *next  = iftCreateImage(img->xsize,img->ysize,img->zsize);


  /* Repeat the following process during a certain number of
     iterations (niters), excluding previously selected voxels: For
     each voxel p, compute the mean value among its adjacent
     voxels. If the intensity of p is below a percentile of the mean
     value, then p is selected (i.e., assign value to the output
     intensity of p). */

  if (mask != NULL){
    for (it=1; it <= niters; it++) {
      for (p=0; p < img->n; p++) {
	if (mask->val[p]!=0){
	  u = iftGetVoxelCoord(img,p);
	  if (prev->val[p] == 0) {
	    mean = 0.0; n=0;
	    for (i=1; i < A->n; i++) {
	      v = iftGetAdjacentVoxel(A,u,i);
	      if (iftValidVoxel(img,v)){
		q     = iftGetVoxelIndex(img,v);
		if ((prev->val[q]==0)&&(mask->val[q]!=0)){
		  mean += img->val[q];
		  n++;
		}
	      }
	    }
	    if (n > 0) mean /= n;	
	    if (img->val[p] < perc*mean){
	      next->val[p] = value;
	    }
	  }
	}
      }

      for (p=0; p < img->n; p++) {
	prev->val[p]  = iftMax(prev->val[p], next->val[p]);
	next->val[p]  = 0;
      }
    }
  } else {
    for (it=1; it <= niters; it++) {
      for (p=0; p < img->n; p++) {
	u = iftGetVoxelCoord(img,p);
	if (prev->val[p] == 0) {
	  mean = 0.0; n=0;
	  for (i=1; i < A->n; i++) {
	    v = iftGetAdjacentVoxel(A,u,i);
	    if (iftValidVoxel(img,v)){
	      q     = iftGetVoxelIndex(img,v);
	      if (prev->val[q]==0){
		mean += img->val[q];
		n++;
	      }
	    }
	  }
	  if (n > 0) mean /= n;	
	  if (img->val[p] < perc*mean){
	    next->val[p] = value;
	  }
	}
      }

      for (p=0; p < img->n; p++) {
	prev->val[p]  = iftMax(prev->val[p], next->val[p]);
	next->val[p]  = 0;
      }
    }
  }

  iftDestroyImage(&next);

  return(prev);
}

/* For binary segmentation only */
float iftFBetaScore(iftImage *bin, iftImage *gt, float beta)
{
  float TP=0.0,FP=0.0,FN=0.0,FBeta=0.0,beta2=beta*beta;
  int p; 
  
  iftVerifyImageDomains(bin,gt,"iftFBetaScore");

  for (p=0; p < bin->n; p++) {
    if ((bin->val[p]!=0)&&(gt->val[p]!=0))
      TP++;
    if ((bin->val[p]!=0)&&(gt->val[p]==0))
      FP++;
    if ((bin->val[p]==0)&&(gt->val[p]!=0))
      FN++;
  }

  FBeta = ((1.0+beta2)*TP)/((1.0+beta2)*TP+beta2*FN+FP);

  return(FBeta);
}


/**
@brief F1-Score for binary segmentation.\n

Given a groundtruth and binary label image, compute the F1-Score.\n
 
@param  bin                 iftImage.    Binary label image. 0 for background and 1 for object.
@param  gt                  iftImage.    Correct segmentation. 0 for background and 1 for object.

@return A float with the FScore [0, 1].
*/
float iftFScoreError(iftImage *bin, iftImage *gt)
{
  float fscore = 0.0;
  iftErrorClassification errors;

  errors = iftSegmentationErrors(gt, bin);
  fscore = iftFScoreGivenErrors(&errors);

  return fscore;
}



/**
@brief F1-Score extension for multilabels.\n

Given a multilabel image and a multilabel ground truth, it computes the Fscore measure indivudally. 
The average is stored at index 0. The other i indexes store the fscore error for the ith object.
The perfect segmentation produces error 1.\n
 
@param  mlabel                iftImage.    Multi label image resultant from a segmentation.
@param  mlgt                  iftImage.    Multi label groundtruth read from disk.
@param  number_of_objects     Integer.     Number of objets. For binary segmentation, this value must be 1.

@return A float array of size (number_of_objects + 1) with the fscore of each object. Index 0 is the average between all objects.
*/
float *iftFScoreMultiLabel (iftImage *mlabel, iftImage *mlgt, int number_of_objects)
{
    int p = 0, k = 0;
    iftImage *tmp_label, *tmp_gt;

    float *output = iftAllocFloatArray(number_of_objects + 1);
    float accu_sum = 0.0;
    iftVerifyImageDomains(mlabel, mlgt, "iftFScoreMultiLabel");

    for (k = 1; k <= number_of_objects; k++)
    {
        tmp_label = iftCreateImage(mlabel->xsize, mlabel->ysize, mlabel->zsize);
        tmp_gt    = iftCreateImage(mlgt->xsize, mlgt->ysize, mlgt->zsize);
        for (p = 0; p < mlabel->n; p++)
        {
            if (mlabel->val[p] == k)
                tmp_label->val[p] = 1;
            if (mlgt->val[p] == k)
                tmp_gt->val[p] = 1;
        }
        output[k] = iftFScoreError(tmp_gt, tmp_label);
        accu_sum += output[k];
        iftDestroyImage(&tmp_label);
        iftDestroyImage(&tmp_gt);
    }
    output[0] = (float) accu_sum / (number_of_objects * 1.0);
    return output;
}


/* IMPORTANT: This method only works for binary gtruth and segmentation image */
iftErrorClassification iftSegmentationErrors(iftImage* gt_image, iftImage* cl_image)
{
  iftErrorClassification errors;
  errors.tp = 0; errors.fp = 0; errors.tn = 0; errors.fn = 0;

  int i;
  for(i = 0; i < cl_image->n; i++){
    //Object Voxel
    if(gt_image->val[i] != 0){
      if(cl_image->val[i] == gt_image->val[i])
        errors.tp++;
      else
        errors.fn++;
    }
    //Background Voxel
    else{
      if(cl_image->val[i] == 0)
        errors.tn++;
      else
        errors.fp++;
    }
  }
  return errors;
}

/**

@brief Preserves the input label of border voxels inside and outside
the object and sets the label of the remaining voxels to zero.

@author Alexandre Falcao

@param label: it can be an image with labeled objects or labeled
supervoxels. If the label image contains a single object, the function
will preserve only the internal border voxels.
@param get_margins: set if the margins of the original image are needed. 1 for yes and 0 for not. Parameter added by Adan Echemendia on May, 2017.

@date November, 30th 2015

@return border: A labeled image of border voxels. Each border voxel
preserves the orginal label of its supervoxel/object in the input
image.

*/

iftImage *iftBorderImage(const iftImage *label, bool get_margins)
{
 iftAdjRel *A;
 iftImage  *border = iftCreateImage(label->xsize,label->ysize,label->zsize);
 int        p,q,i; 
 iftVoxel   u, v;
    
  if (iftIs3DImage(label))
    A = iftSpheric(1.0);
  else
    A = iftCircular(1.0);

  if (get_margins){
    for(p=0; p < label->n; p++){
      u = iftGetVoxelCoord(label, p);
      for(i=1; i < A->n; i++){
        v = iftGetAdjacentVoxel(A,u,i);
        if (iftValidVoxel(label, v)){
          q = iftGetVoxelIndex(label, v);
          if (label->val[p] != label->val[q]){
            border->val[p] = label->val[p];
            break;
          }
        } else {
          border->val[p] = label->val[p];
        }
      }
    }
  }
  else{
    for(p=0; p < label->n; p++) {
      u = iftGetVoxelCoord(label, p);
      for (i = 1; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A, u, i);
        if (iftValidVoxel(label, v)) {
          q = iftGetVoxelIndex(label, v);
          if (label->val[p] != label->val[q]) {
            border->val[p] = label->val[p];
            break;
          }
        }
      }
    }
  }

    iftDestroyAdjRel(&A);
    return(border);
}


/**

@brief Computes the recall between two border images, the ground truth and the resulting image from segmentation 
@author Alexandre Falcao
@param  gt: the ground truth image of borders 
@param  border: the border image resulting from segmentation
@param  tolerance_dist: a distance of tolerance from the ground truth, since
the borders of ground-truth masks created by manual tracing do not
follow the border of the object in the image.
@date   November, 30th 2015
@return a number within [0,1] that represents the boundary recall.

*/

float iftBoundaryRecall(iftImage *gt, iftImage *border, float tolerance_dist)
{
  iftVerifyImageDomains(gt, border,"iftBoundaryRecall");
  float number_of_gt_voxels=0,number_of_matchings=0; 
  int   p, q, i;
  iftVoxel u, v;
  iftAdjRel *A;

  if (iftIs3DImage(gt))
    A = iftSpheric(tolerance_dist);
  else
    A = iftCircular(tolerance_dist);

  for(p = 0; p < gt->n; p++)
    if (gt->val[p]!=0){      
      number_of_gt_voxels++;
      u = iftGetVoxelCoord(gt,p);
        for (i=0; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A,u,i);
        if (iftValidVoxel(gt,v)){
          q = iftGetVoxelIndex(gt,v);
          if (border->val[q]!=0){
            number_of_matchings++;
            break;
          }
        }
      }
    }

  if (number_of_gt_voxels == 0.0)
      iftError("Empty ground-truth image", "iftBoundaryRecall");

  iftDestroyAdjRel(&A);

  return (number_of_matchings/number_of_gt_voxels);
}


/**

@brief Computes the under segmentation error --- total of leaking voxels across the boundaries of a ground truth image.
@author Alexandre Falcao and John Vargas
@param  gt: the ground truth image of regions 
@param  label: the label image resulting from segmentation
@param  perc_of_intersection: Percentage of intersection between region and ground truth segment to consider the region as part of the object (true intersection).
@date   December, 7th 2015
@return A number that represents the under segmentation error.

*/

float iftUnderSegmentationSLIC(iftImage *gt, iftImage *label, float perc_of_intersection)
{
  iftVerifyImageDomains(gt, label, "iftUnderSegmentation");
  int i,j, p, number_of_gt_segments, number_of_regions;

  number_of_gt_segments = iftMaximumValue(gt)    + 1;
  number_of_regions     = iftMaximumValue(label); 

  int *number_of_voxels_per_region = iftAllocIntArray(number_of_regions);
  for(p = 0; p < label->n; p++)
    number_of_voxels_per_region[label->val[p] - 1]++;

  int *number_of_voxels_from_intersecting_regions = iftAllocIntArray(number_of_regions);
  int  total_number_of_voxels_from_intersecting_regions = 0;

  for(j = 0; j < number_of_gt_segments; j++){ /* for each segment of the gt image */

    for(p = 0; p < label->n; p++){ /* compute the number of voxels from each intersecting region */
      if(gt->val[p] == j){
	number_of_voxels_from_intersecting_regions[label->val[p] - 1]++;
      }
    }

    for(i = 0; i < number_of_regions; i++){
      /* percentage of intersection between region and segment */
      float perc = (float)number_of_voxels_from_intersecting_regions[i]/(float)number_of_voxels_per_region[i];

      if(perc >= perc_of_intersection){ /* the region is considered a true intersection */
	total_number_of_voxels_from_intersecting_regions += number_of_voxels_per_region[i];
      }

      number_of_voxels_from_intersecting_regions[i] = 0;

    }
  }

  float ue = (1.0/(float)label->n)*((float)(total_number_of_voxels_from_intersecting_regions - label->n));
  iftFree(number_of_voxels_per_region);
  iftFree(number_of_voxels_from_intersecting_regions);

  return(ue);
}

float iftUnderSegmentation(iftImage *gt_image, iftImage *label_image)
{
  iftVerifyImageDomains(gt_image, label_image, "iftUnderSegmentationMin");
  int i, j, p, num_obj, num_regions;
  float area, total_err;
  num_obj = iftMaximumValue(gt_image) + 1;
  num_regions = iftMaximumValue(label_image);
  
  int *num_pix_total = iftAllocIntArray(num_regions);
  for(p = 0; p < label_image->n; p++)
    num_pix_total[label_image->val[p] - 1]++;
  
  int *num_pix_obj = iftAllocIntArray(num_regions);
  
  int sum_err = 0;
  for(j = 0; j < num_obj; j++){
    for(p = 0; p < label_image->n; p++){
      if(gt_image->val[p] == j){
        num_pix_obj[label_image->val[p] - 1]++;
      }
    }
    for(i = 0; i < num_regions; i++){
      area =  num_pix_obj[i];
      if((num_pix_total[i] - num_pix_obj[i]) < area){
        area = (num_pix_total[i] - num_pix_obj[i]);
      }
      sum_err = sum_err + area;
      num_pix_obj[i] = 0;
    }
  }
  
  total_err = sum_err/(float)label_image->n;
  
  free(num_pix_total);
  free(num_pix_obj);
  
  return total_err;
}


float iftCompactness2D(iftImage *label) {
  iftImage  *border = iftBorderImage(label,1);
  iftAdjRel *A=NULL;
  int      i, p , q, *area, number_of_regions;
  iftVoxel u, v; 
  float    *perim, co;
  
  number_of_regions = iftMaximumValue(label);
  perim = iftAllocFloatArray(number_of_regions);
  area  = iftAllocIntArray(number_of_regions);

  if (iftIs3DImage(label))
      iftError("Label image must be 2D", "iftCompactness2D");
  else
    A    = iftCircular(sqrtf(2.0));

  for (p=0; p < label->n; p++) {
    area[label->val[p] - 1]++; 
    if (border->val[p]) {
      u = iftGetVoxelCoord(label,p);
      for (i=1; i < A->n; i++) {
	v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(label,v)){
	  q = iftGetVoxelIndex(label,v);
	  if (border->val[p]==border->val[q]){
	    perim[border->val[p] - 1] += iftSmoothEuclideanDistance((float) iftSquaredVoxelDistance(u, v));
	  }
	}
      }
    }
  }
  for (i=0; i < number_of_regions; i++) 
    perim[i] /= 2.0;

  iftDestroyImage(&border);

  co = 0; 
  for (i=0; i < number_of_regions; i++) 
    if (perim[i] > 0.0) 
      co += (4.0 * IFT_PI * area[i] / (perim[i] * perim[i])) * ((float)area[i] / label->n);

  iftFree(area);
  iftFree(perim);
  iftDestroyAdjRel(&A);

  return(co);
}


float iftTopologyMeasure(iftImage *label) {
  iftMatrix *adjMatrix;
  int        number_of_regions, number_of_valid_regions;
  float     *average_size, topology, average_number_of_adjacent_segments;
  char      *valid_region;
  int        p, q, i, index, r, c, nadjs;
  iftAdjRel *A;
  iftVoxel   u, v;

  if (iftIs3DImage(label))
    A    = iftSpheric(1.0);
  else
    A    = iftCircular(1.0);

  number_of_regions = iftMaximumValue(label);
  adjMatrix         = iftCreateMatrix(number_of_regions, number_of_regions);
  valid_region      = iftAllocCharArray(number_of_regions);
  average_size      = iftAllocFloatArray(number_of_regions);

  /* Find the valid regions as those with no adjacent voxels outside
     the image domain */

  for (r=0; r < number_of_regions; r++)
    valid_region[r]=1;

  number_of_valid_regions = number_of_regions;

  for (p=0; p < label->n; p++) {
    if (valid_region[label->val[p] - 1]){
      u = iftGetVoxelCoord(label,p);
      for (i=1; i < A->n; i++) {
	v = iftGetAdjacentVoxel(A,u,i);
	if (!iftValidVoxel(label,v)){
	  valid_region[label->val[p] - 1]=0;
	  number_of_valid_regions--;
	  break;
	}
      }
    }
  }

  if (number_of_valid_regions==0){
    iftWarning("Infinity topology","iftTopologyMeasure");
    return(IFT_INFINITY_FLT);
  }

  /* For each valid region, compute in the adjacency matrix the size
     of its border segments with each adjacent region. */

  for (p=0; p < label->n; p++) {
    if (valid_region[label->val[p] - 1]){
      u = iftGetVoxelCoord(label,p);
      r = label->val[p] - 1;
      for (i=1; i < A->n; i++) {
	v = iftGetAdjacentVoxel(A,u,i);
	q = iftGetVoxelIndex(label,v);
	if (label->val[p] != label->val[q]) {
	  c     = label->val[q] - 1;
	  index = iftGetMatrixIndex(adjMatrix, r, c);
	  adjMatrix->val[index]++;
	}
      }
    }
  }

  /* For each valid region, compute the average size of its border
     segments between regions. Compute also the average number of adjacent segments. */

  average_number_of_adjacent_segments=0.0;
 
  for (r=0; r < adjMatrix->nrows; r++) {
    if (valid_region[r]){
      nadjs = 0;
      for (c=0; c < adjMatrix->ncols; c++) {
	index = iftGetMatrixIndex(adjMatrix, r, c);
	if (adjMatrix->val[index]>0.0){
	  average_size[r] += adjMatrix->val[index];
	  nadjs++;
	}
      }
      average_number_of_adjacent_segments += nadjs;
      average_size[r] /= nadjs;
    }
  }

  average_number_of_adjacent_segments /= number_of_valid_regions;

  /* compute the topology measure */

  topology = 0;
  for (r=0; r < adjMatrix->nrows; r++) {
    if (valid_region[r]){
      for (c=0; c < adjMatrix->ncols; c++) {
	index = iftGetMatrixIndex(adjMatrix, r, c);
	if (adjMatrix->val[index]>0.0){
	  topology += fabs(adjMatrix->val[index]-average_size[r]);
	}
      }
    }
  }

  if (average_number_of_adjacent_segments > IFT_EPSILON)
    topology /= (number_of_valid_regions * average_number_of_adjacent_segments);
  else{
    iftWarning("Infinity topology","iftTopologyMeasure");
    topology = IFT_INFINITY_FLT;
  }

  iftFree(average_size);
  iftDestroyMatrix(&adjMatrix);
  iftFree(valid_region);

  return(topology);
}


iftImage *iftOrientedWatershed(const iftImage *img, const iftImage *grad, iftAdjRel *Ain,
                               iftFloatArray *alpha, iftLabeledSet *seeds, iftSet *forbidden)
{
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(img))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }
    else A = Ain;

    if (alpha->n != iftNumberOfLabels(seeds)) {
        iftWarning("Number of labels and alphas are different", "iftOrientedByIntensityWatershed");
    }

    // Initialization
    iftImage *pathval  = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftImage *label = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftGQueue *Q     = iftCreateGQueue(4 * iftMaximumValue(grad) + 1, img->n, pathval->val);

    for (int p = 0; p < img->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
    }

    // sets the forbidden region
    iftSet *F = forbidden;
    while (F != NULL) {
        int p   =  F->elem;
        pathval->val[p] = IFT_INFINITY_INT_NEG;
        F = F->next;
    }

    iftLabeledSet *S = seeds;
    while (S != NULL)
    {
        int p = S->elem;
        label->val[p] = S->label;
        pathval->val[p] = 0;
        iftInsertGQueue(&Q, p);
        S = S->next;
    }

    // Image Foresting Transform
    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftGetVoxelCoord(img, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(img, v))
            {
                int q = iftGetVoxelIndex(img, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    int w_p_q;
                    /* Orientation */
                    if (img->val[p] > img->val[q]) {
                        /* Bright to dark */
                        w_p_q = grad->val[q] * (1 + alpha->val[label->val[p]]);
                    } else {
                        /* Dark to bright */
                        w_p_q = grad->val[q] * (1 - alpha->val[label->val[p]]);
                    }

                    int tmp;
                    if (label->val[p] != 0) {
                        tmp = iftMax(pathval->val[p], 2 * w_p_q + 1);
                    } else {
                        tmp = iftMax(pathval->val[p], 2 * w_p_q);
                    }

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q]  = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }
    iftDestroyGQueue(&Q);
    iftDestroyImage(&pathval);

    iftCopyVoxelSize(img, label);

    return (label);
}


iftImage *iftSelectAndPropagateRegionsAboveArea(iftImage *label, int area)
{
  iftHist *hist;
  iftAdjRel *A;
  iftImage  *rlabel,*nlabel;
  iftFImage *pathval=NULL;
  iftFHeap  *Q=NULL;
  int        i,p,q;
  iftVoxel   u,v;
  
  if (iftIs3DImage(label))
    A = iftSpheric(sqrtf(3.0));
  else
    A = iftCircular(sqrtf(2.0));

  if (iftMinimumValue(label)==0)
    for (p=0; p < label->n; p++) 
      label->val[p]++;

  rlabel     = iftRelabelRegions(label,A);

  int max_val = iftMaximumValue(rlabel);
  int nbins = max_val + 1;
  hist       = iftCalcGrayImageHist(rlabel, NULL, nbins, max_val, 0);

  int max_area = IFT_NIL;
  for (i=0; i < hist->nbins; i++) {
    if (hist->val[i]>max_area)
      max_area = hist->val[i];
  }

  for (p=0; p < rlabel->n; p++){
    if (hist->val[rlabel->val[p]] < area)
      rlabel->val[p]= IFT_NIL; 
  }

  // Initialization 

  pathval = iftCreateFImage(rlabel->xsize,rlabel->ysize,rlabel->zsize);
  nlabel  = iftCopyImage(rlabel);
  Q       = iftCreateFHeap(rlabel->n,pathval->val);
  for (p=0; p < rlabel->n; p++) {
    pathval->val[p]= IFT_INFINITY_INT;
    if (rlabel->val[p] != IFT_NIL){
      pathval->val[p] = max_area - hist->val[rlabel->val[p]];
      iftInsertFHeap(Q,p);
    }
  }

  // Image Foresting Transform

  while(!iftEmptyFHeap(Q)) {
    p=iftRemoveFHeap(Q);
    u = iftGetVoxelCoord(rlabel,p);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);

      if (iftValidVoxel(rlabel,v)){	
        q = iftGetVoxelIndex(rlabel,v);
        if (rlabel->val[q] == IFT_NIL){
          if (pathval->val[q] > pathval->val[p]){
            nlabel->val[q]   = nlabel->val[p];
            pathval->val[q]  = pathval->val[p];
            if(Q->color[q] == IFT_WHITE)
              iftInsertFHeap(Q, q);
            else
              iftGoUpFHeap(Q, Q->pos[q]);
          }
        }
      }
    }
  }


  iftDestroyFHeap(&Q);
  iftDestroyFImage(&pathval);
    iftDestroyHist(&hist);
  iftDestroyImage(&rlabel); 
  rlabel = iftRelabelRegions(nlabel,A);
  iftDestroyImage(&nlabel); 
  iftDestroyAdjRel(&A); 

  return(rlabel);
}

iftImage *iftSelectAndPropagateRegionsAboveAreaByColor(iftImage *img, iftImage *label, int area)
{
  iftHist *hist;
  iftAdjRel *A;
  iftImage  *rlabel,*nlabel;
  iftFImage *pathval=NULL;
  iftFHeap  *Q=NULL;
  int        i,p,q;
  iftVoxel   u,v;

  if (iftIs3DImage(label))
    A = iftSpheric(sqrtf(3.0));
  else
    A = iftCircular(sqrtf(2.0));

  if (iftMinimumValue(label)==0)
    for (p=0; p < label->n; p++)
      label->val[p]++;

  rlabel     = iftRelabelRegions(label,A);

  int max_val = iftMaximumValue(rlabel);
  int nbins = max_val + 1;
  hist       = iftCalcGrayImageHist(rlabel, NULL, nbins, max_val, 0);
  nlabel  = iftCopyImage(rlabel);

  for (p=0; p < rlabel->n; p++){
    if (hist->val[rlabel->val[p]] < area)
      rlabel->val[p]= IFT_NIL;
  }

  // Initialization
  pathval = iftCreateFImage(rlabel->xsize,rlabel->ysize,rlabel->zsize);

  Q = iftCreateFHeap(rlabel->n,pathval->val);
  for (p=0; p < rlabel->n; p++) {
    pathval->val[p]= IFT_INFINITY_INT;
    if (rlabel->val[p] != IFT_NIL){
      pathval->val[p] =0.0;
      iftInsertFHeap(Q,p);
    }
  }

  // Image Foresting Transform
  if (iftIsColorImage(img)){
    float tmp;
    float alpha[3]={1.0,1.0,1.0};
    float f1[3],f2[3];

    while(!iftEmptyFHeap(Q)) {
      p=iftRemoveFHeap(Q);
      u = iftGetVoxelCoord(rlabel,p);
      f1[0]=img->val[p];f1[1]=img->Cb[p];f1[2]=img->Cr[p];
      for (i=1; i < A->n; i++){
        v = iftGetAdjacentVoxel(A,u,i);

        if (iftValidVoxel(rlabel,v)){
          q = iftGetVoxelIndex(rlabel,v);
          if (rlabel->val[q] == IFT_NIL){
            f2[0]=img->val[q];f2[1]=img->Cb[q];f2[2]=img->Cr[q];
            tmp=iftMax(pathval->val[p],iftDistance1(f1,f2,alpha,3));
            if (pathval->val[q] > tmp){
              nlabel->val[q]   = nlabel->val[p];
              pathval->val[q]  = tmp;
              if(Q->color[q] == IFT_WHITE)
                iftInsertFHeap(Q, q);
              else
                iftGoUpFHeap(Q, Q->pos[q]);
            }
          }
        }
      }
    }
  }
  else{
    float tmp;
    float alpha[1]={1.0};
    float f1[1],f2[1];

    while(!iftEmptyFHeap(Q)) {
      p=iftRemoveFHeap(Q);
      u = iftGetVoxelCoord(rlabel,p);
      f1[0]=img->val[p];
      for (i=1; i < A->n; i++){
        v = iftGetAdjacentVoxel(A,u,i);

        if (iftValidVoxel(rlabel,v)){
          q = iftGetVoxelIndex(rlabel,v);
          if (rlabel->val[q] == IFT_NIL){
            f2[0]=img->val[q];
            tmp=iftMax(pathval->val[p],iftDistance1(f1,f2,alpha,1));
            if (pathval->val[q] > tmp){
              nlabel->val[q]   = nlabel->val[p];
              pathval->val[q]  = tmp;
              if(Q->color[q] == IFT_WHITE)
                iftInsertFHeap(Q, q);
              else
                iftGoUpFHeap(Q, Q->pos[q]);
            }
          }
        }
      }
    }
  }

  iftDestroyFHeap(&Q);
  iftDestroyFImage(&pathval);
    iftDestroyHist(&hist);
  iftDestroyImage(&rlabel);
  rlabel = iftRelabelRegions(nlabel,A);
  iftDestroyImage(&nlabel);
  iftDestroyAdjRel(&A);

  return(rlabel);
}

iftImage *iftSmoothRegionsByDiffusion(const iftImage *label_img, const iftImage *orig_img,
                                      float smooth_factor, int n_iters) {
    iftAdjRel *A = NULL;

    if (iftIs3DImage(label_img))
        A = iftSpheric(sqrtf(3.0));
    else
        A = iftCircular(sqrtf(2.0));

    iftImage  *grad   = iftImageGradientMagnitude(orig_img, A);
    iftFImage *weight = iftSmoothWeightImage(grad,smooth_factor);
    iftImage  *nlabel = iftFastSmoothObjects(label_img,weight,n_iters);

    iftDestroyImage(&grad);
    iftDestroyFImage(&weight);
    iftDestroyAdjRel(&A);

return(nlabel);
}

char iftIsSegmentationConsistent(iftImageForest *fst)
{
  int p, q;

  for (p = 0; p < fst->label->n; p++)
  {
    q = p;
    while (fst->pred->val[q] != IFT_NIL)
    {
      if ((fst->label->val[q] != fst->label->val[fst->pred->val[q]]) ||
      (fst->root->val[q]  != fst->root->val[fst->pred->val[q]])    )
      {
	printf("Segmentation inconsistency: A voxel has label %d and root %d,  \n and its predecessor has label %d and root %d\n", fst->label->val[q], fst->root->val[q], fst->label->val[fst->pred->val[q]], fst->root->val[fst->pred->val[q]]);
        return (0);
      }

      q = fst->pred->val[q];
    }

  }

  return (1);
}

iftSmoothBorder * iftCreateSmoothBorder(iftImage *basins, iftAdjRel *A, int smooth_iterations, float smooth_factor)
{
  iftSmoothBorder *smooth = (iftSmoothBorder*)iftAlloc(1, sizeof(iftSmoothBorder));

  smooth->border_weight       = iftSmoothWeightImage(basins, smooth_factor);
  smooth->norm_factor         = iftWeightNormFactor(smooth->border_weight, A);
  smooth->smooth_factor       = smooth_factor;
  smooth->smooth_iterations   = smooth_iterations;
  smooth->prev_label          = iftCreateImage(basins->xsize,basins->ysize,basins->zsize);
  smooth->next_label          = iftCreateImage(basins->xsize,basins->ysize,basins->zsize);
  smooth->prev_weight         = iftCreateFImage(basins->xsize, basins->ysize, basins->zsize);
  smooth->next_weight         = iftCreateFImage(basins->xsize, basins->ysize, basins->zsize);
  smooth->prev_marker         = iftCreateImage(basins->xsize,basins->ysize,basins->zsize);
  smooth->next_marker         = iftCreateImage(basins->xsize,basins->ysize,basins->zsize);
  return smooth;
}


void iftDestroySmoothBorder (iftSmoothBorder **smooth)
{
  iftSmoothBorder *aux;

  aux = *smooth;
  if(aux != NULL){
    if (aux->border_weight != NULL)  iftDestroyFImage(&(aux->border_weight));
    if (aux->norm_factor   != NULL)  iftDestroyFImage(&(aux->norm_factor));
    if (aux->prev_label    != NULL)  iftDestroyImage(&(aux->prev_label));
    if (aux->next_label    != NULL)  iftDestroyImage(&(aux->next_label));
    if (aux->prev_marker   != NULL)  iftDestroyImage(&(aux->prev_marker));
    if (aux->next_marker   != NULL)  iftDestroyImage(&(aux->next_marker));
    if (aux->prev_weight   !=NULL) iftDestroyFImage(&(aux->prev_weight)); 
    if (aux->next_weight   !=NULL) iftDestroyFImage(&(aux->next_weight)); 
    iftFree(aux);
    *smooth = NULL;
  }
}

int iftIsLabelRootConnected (iftImage *pred, iftImage *label, int p)
{
  while(pred->val[p] != IFT_NIL)
  {
    if (label->val[p] != label->val[pred->val[p]])
      return 0;
    p = pred->val[p];
  }
  return 1;
}

iftLabeledSet *iftLabelToForestGeodesicRobot(iftImage *gradient,
    iftImage *label,
    int seeds_per_iteration,
    int number_of_objects,
    int min_distance_border,
    int max_marker_size,
    int min_marker_size)
{
  int p, j, seeds_added, number_of_seeds = 0;
  iftLabeledSet  *available_seeds = NULL, *current_seeds = NULL;
  iftBMap *seeds_bmap            = NULL;
  iftImage  *seeds_image         = NULL;
  iftImage *seeds_image_copy     = NULL;
  iftImage *current_segmentation = NULL;
  iftImageForest *forest         = NULL;
  iftImage *gt_image             = NULL;
  iftLabeledSet *segmentation_seeds = NULL;
  iftAdjRel *A = NULL;
  
  // Creating adjacency relation
  if (iftIs3DImage(gradient))
    A = iftSpheric(sqrt(3));
  else
    A = iftCircular(sqrt(2));

  gt_image = label;

  forest = iftCreateImageForest(gradient, A);
  seeds_bmap  = iftCreateBMap(gradient->n);
  seeds_image = iftCreateImage(gradient->xsize, gradient->ysize, gradient->zsize);
  iftSetImage(seeds_image, -1);

  seeds_added = -1;
  for (j = 0; j < 50 && seeds_added != 0; j++)
  { 
    printf("Rebuilding %.2fp\n", (float) j*2);
    available_seeds  = iftGeodesicMarkersForSegmentation(gt_image, current_segmentation);
    seeds_image_copy = iftCopyImage(seeds_image);
    seeds_added      = iftMarkersFromMisclassifiedSeeds(seeds_image, available_seeds, seeds_bmap, seeds_per_iteration,
                                                        number_of_objects + 1, gt_image, current_segmentation,
                                                        min_distance_border, max_marker_size, min_marker_size);
    
    //This produces only the new seeds added this iteration
    for (p = 0; p < seeds_image_copy->n; p++)
    {
      if (seeds_image_copy->val[p] == seeds_image->val[p])
        seeds_image_copy->val[p] = IFT_NIL;
      else
        seeds_image_copy->val[p] = seeds_image->val[p];
    }

    current_seeds = iftLabeledSetFromSeedImage(seeds_image_copy, false);
    iftConcatLabeledSet(&segmentation_seeds, &current_seeds);

    iftDiffWatershed(forest, current_seeds, NULL);
    current_segmentation = forest->label;

    number_of_seeds += seeds_added;

    iftDestroyImage(&seeds_image_copy);
    iftDestroyLabeledSet(&available_seeds);
    iftDestroyLabeledSet(&current_seeds);
  }
  iftDestroyBMap(&seeds_bmap);
  iftDestroyImage(&seeds_image);
  printf("Rebuilding 100p\n");

  iftDestroyImageForest(&forest);
  iftDestroyAdjRel(&A);

  return segmentation_seeds;
}


iftLabeledSet *iftShiftCoordsFromLabeledSet(const iftLabeledSet *S, const iftImage *img, int dx, int dy, int dz) {
    iftLabeledSet *shift     = NULL;

    if (S != NULL) {
        if (img == NULL)
            iftError("Image is NULL", "iftShiftCoordsFromLabeledSet");

        const iftLabeledSet *aux = S;
        while (aux != NULL) {
            int label  = aux->label;
            int p      = aux->elem;
            iftVoxel v = iftGetVoxelCoord(img, p);

            // shifting
            v.x += dx;
            v.y += dy;
            v.z += dz;

            if (iftValidVoxel(img, v)) {
                int q = iftGetVoxelIndex(img, v);
                iftInsertLabeledSet(&shift, q, label);
            }

            aux = aux->next;
        }
    }

    return shift;
}


iftSet *iftShiftCoordsFromSet(const iftSet *S, const iftImage *img, int dx, int dy, int dz) {
    iftSet *shift = NULL;

    if (S != NULL) {
        if (img == NULL)
            iftError("Image is NULL", "iftShiftCoordsFromSet");

        const iftSet *aux = S;
        while (aux != NULL) {
            int p      = aux->elem;
            iftVoxel v = iftGetVoxelCoord(img, p);

            // shifting
            v.x += dx;
            v.y += dy;
            v.z += dz;

            if (iftValidVoxel(img, v)) {
                int q = iftGetVoxelIndex(img, v);
                iftInsertSet(&shift, q);
            }

            aux = aux->next;
        }
    }

    return shift;
}

iftImage *iftConnectedThreshold(const iftImage *img, const iftAdjRel *A, iftImage *marker, const int minthres, const int maxthres)
{
  iftImage *extended_marker = iftCopyImage(marker);
  iftFIFO  *Q = iftCreateFIFO(img->n);
  int       i, p, q;
  iftVoxel  u, v;

  iftVerifyImageDomains(img,marker,"iftConnectedThreshold");
  
  for (p=0; p < img->n; p++) {
    if (extended_marker->val[p]!=0)
      iftInsertFIFO(Q,p);
  }
  
  while (!iftEmptyFIFO(Q)){
    p = iftRemoveFIFO(Q);
    u = iftGetVoxelCoord(img,p);

    for (i=1; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){
	q = iftGetVoxelIndex(img,v);
	if (Q->color[q]==IFT_WHITE){
	  if ((img->val[q]>=minthres)&&
	      (img->val[q]<=maxthres)){
	    extended_marker->val[q]=extended_marker->val[p];
	    iftInsertFIFO(Q,q);
	  }	    
	}
      }
    }
  }

  iftDestroyFIFO(&Q);

  return(extended_marker);
}


iftImage *iftMajorityVoting(const iftImage *img, const iftIntMatrix *label_occurs) {
    if (img->n != label_occurs->nrows)
        iftError("Number of Voxels is different in image and in matrix (%d, %d)", "iftMajorityVoting",
                 img->n, label_occurs->nrows);

    iftImage *seg_img = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftCopyVoxelSize(img, seg_img);

    #pragma omp parallel for
    for (int p = 0; p < label_occurs->nrows; p++) {
        int maj_label = 0;
        int max_occ   = -1;

        for (int label = 0; label < label_occurs->ncols; label++) {
            if (iftMatrixElem(label_occurs, label, p) > max_occ) {
                maj_label = label;
                max_occ   = iftMatrixElem(label_occurs, label, p);
            }
        }
        seg_img->val[p] = maj_label;
    }

    return seg_img;
}


iftImage* iftDescendantMap(iftImage *pred, iftAdjRel *A, iftBMap *set_of_interest) {
    iftImage *descendants = NULL;
    iftFIFO *Q = NULL;
    iftLIFO *S = NULL;
    bool count_on_set = (set_of_interest != NULL);

    descendants = iftCreateImage(pred->xsize, pred->ysize, pred->zsize);
    Q = iftCreateFIFO(pred->n);
    S = iftCreateLIFO(pred->n);

    for(int p = 0; p < pred->n; p++) {
        if(pred->val[p] == IFT_NIL) {
            iftInsertFIFO(Q, p);
        }
    }

    while(!iftEmptyFIFO(Q)) {
        int p = iftRemoveFIFO(Q);
        iftVoxel u = iftGetVoxelCoord(pred, p);
        iftInsertLIFO(S, p);

        for(int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if(iftValidVoxel(pred, v) && iftImgVoxelVal(pred, v) == p) {
                iftInsertFIFO(Q, iftGetVoxelIndex(pred, v));
            }
        }
    }

    while(!iftEmptyLIFO(S)) {
        int p = iftRemoveLIFO(S);

        if(pred->val[p] != IFT_NIL) {
            // Counting all leaves/tree nodes as descendants
            if(!count_on_set) {
                descendants->val[pred->val[p]] += descendants->val[p] + 1;
            } else {
                // Counting only descendants in the set of interest voxels
                descendants->val[pred->val[p]] += descendants->val[p];

                if(iftBMapValue(set_of_interest, p)) {
                    descendants->val[pred->val[p]]++;
                }
            }
        }
    }

    return descendants;
}


iftImage *iftConnCompLabeling(const iftImage *label_img, int *n_objs) {
    iftImage *out_label_img = iftCreateImageFromImage(label_img);

    iftIntQueue *Q = iftCreateIntQueue(label_img->n);
    iftAdjRel *A = iftSpheric(sqrt(3.0));

    int label = 1;

    for (int t = 0; t < label_img->n; t++) {
        if ((label_img->val[t] != 0) && (out_label_img->val[t] == 0)) {
            out_label_img->val[t] = label;
            iftInsertIntQueue(Q, t);
            
            while (!iftIsIntQueueEmpty(Q)) {
                int p;
                iftRemoveIntQueue(Q, &p);
                iftVoxel u = iftGetVoxelCoord(label_img, p);

                for (int i = 1; i < A->n; i++) {
                    iftVoxel v = {.x = u.x + A->dx[i], .y = u.y + A->dy[i], .z = u.z + A->dz[i]};

                    if (iftValidVoxel(label_img, v)) {
                        int q = iftGetVoxelIndex(label_img, v);

                        if ((out_label_img->val[q] == 0) && (label_img->val[p] == label_img->val[q])) {
                            out_label_img->val[q] = label;
                            iftInsertIntQueue(Q, q);
                        }
                    }
                }
            }
            label++;
        }
    }

    if (n_objs != NULL)
        *n_objs = label - 1;

    iftDestroyIntQueue(&Q);

    return out_label_img;
}


iftImage *iftOrientedWaterCut(const iftImage *img, iftAdjRel *Ain, iftFloatArray *alpha,
                                          iftLabeledSet *seeds, iftSet *forbidden)
{
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(img))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }
    else A = Ain;

    if (alpha->n != iftNumberOfLabels(seeds)) {
        iftWarning("Number of labels and alphas are different", "iftOrientedByIntensityWatershed");
    }

    // Initialization
    iftImage *pathval  = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftImage *label = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftGQueue *Q     = iftCreateGQueue(4 * iftMaximumValue(img) + 1, img->n, pathval->val);

    for (int p = 0; p < img->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
    }

    // sets the forbidden region
    iftSet *F = forbidden;
    while (F != NULL) {
        int p   =  F->elem;
        pathval->val[p] = IFT_INFINITY_INT_NEG;
        F = F->next;
    }

    iftLabeledSet *S = seeds;
    while (S != NULL)
    {
        int p = S->elem;
        label->val[p] = S->label;
        pathval->val[p] = 0;
        iftInsertGQueue(&Q, p);
        S = S->next;
    }

    // Image Foresting Transform
    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftGetVoxelCoord(img, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(img, v))
            {
                int q = iftGetVoxelIndex(img, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    int w_p_q;
                    /* Orientation */
                    if (img->val[p] > img->val[q]) {
                        /* Bright to dark */
                        w_p_q = abs(img->val[p] - img->val[q]) * (1 + alpha->val[label->val[p]]);
                    } else {
                        /* Dark to bright */
                        w_p_q = abs(img->val[p] - img->val[q]) * (1 - alpha->val[label->val[p]]);
                    }

                    int tmp;
                    if (label->val[p] != 0) {
                        tmp = iftMax(pathval->val[p], 2 * w_p_q + 1);
                    } else {
                        tmp = iftMax(pathval->val[p], 2 * w_p_q);
                    }

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q]  = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }
    iftDestroyGQueue(&Q);
    iftDestroyImage(&pathval);

    iftCopyVoxelSize(img, label);

    return (label);
}


iftImage *iftOrientedColorWaterCut(iftMImage *mimg, iftImage *orient, iftAdjRel *Ain, float beta, iftLabeledSet *seeds)
{
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(orient))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }
    else A = Ain;

	if (beta > 1 || beta < -1) {
        iftError("Beta must be between (-1, 1)", "iftOrientedColorWaterCut");
	}

    // Initialization
    iftImage *pathval  = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *label = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftGQueue *Q     = iftCreateGQueue(iftMMaximumValue(mimg, -1) * 3, mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
    }

    iftLabeledSet *S = seeds;
    while (S != NULL)
    {
        int p = S->elem;
        label->val[p] = S->label;
        pathval->val[p] = 0;
        iftInsertGQueue(&Q, p);
        S = S->next;
    }

    // Image Foresting Transform
    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftMValidVoxel(mimg, v))
            {
                int q = iftMGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    int w_p_q;
                    /* Orientation */
                    if (orient->val[p] > orient->val[q]) {
                        /* Bright to dark */
                        w_p_q = iftMImageDist(mimg, p, q) * (1 + beta);
                    } else if (orient->val[p] < orient->val[q]) {
                        /* Dark to bright */
                        w_p_q = iftMImageDist(mimg, p, q) * (1 - beta); 
                    } else {
						w_p_q = iftMImageDist(mimg, p, q);
					}

                    int tmp;
                    if (label->val[p] != 0) {
                        tmp = iftMax(pathval->val[p], 2 * w_p_q + 1);
                    } else {
                        tmp = iftMax(pathval->val[p], 2 * w_p_q);
                    }

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q]  = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }
    iftDestroyGQueue(&Q);
    iftDestroyImage(&pathval);

    iftCopyVoxelSize(mimg, label);

    return (label);
}


iftImage *iftEnhancedWaterCut(iftMImage *mimg, iftImage *objmap, iftAdjRel *Ain, iftLabeledSet *seeds, float alpha)
{
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DMImage(mimg))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }
    else A = Ain;

    // Initialization
    float max = iftMMaximumValue(mimg, -1);
    iftImage *pathval  = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    iftImage *label = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    iftGQueue *Q     = iftCreateGQueue(sqrt(max * max * mimg->m), mimg->n, pathval->val);

    for (int p = 0; p < mimg->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
    }

    iftLabeledSet *S = seeds;
    while (S != NULL)
    {
        int p = S->elem;
        label->val[p] = S->label;
        pathval->val[p] = 0;
        iftInsertGQueue(&Q, p);
        S = S->next;
    }

    // Image Foresting Transform
    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(mimg, v))
            {
                int q = iftGetVoxelIndex(mimg, v);

                if (Q->L.elem[q].color != IFT_BLACK)
                {
					int arc_weight = iftMImageDist(mimg, p, q) * (1-alpha) + abs(objmap->val[p] - objmap->val[q]) * alpha;
					int tmp = iftMax(pathval->val[p], arc_weight);
//                    int tmp = arc_weight;

                    if (tmp < pathval->val[q])
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);

                        label->val[q] = label->val[p];
                        pathval->val[q]  = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }
    iftDestroyGQueue(&Q);
    iftDestroyImage(&pathval);

    iftCopyVoxelSize(mimg, label);

    return (label);
}


iftImage *iftEnhancedWatershed(iftImage *basins, iftImage *objmap, iftAdjRel *Ain, iftLabeledSet *seeds, float alpha)
{
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(basins))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }
    else A = Ain;

    // Initialization
    iftImage *pathval = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
    iftImage *label = iftCreateImage(basins->xsize, basins->ysize, basins->zsize);
    iftImage *objbasins = iftImageBasins(objmap, A);
    iftGQueue *Q = iftCreateGQueue(iftMaximumValue(basins) + 1, basins->n, pathval->val);

    for (int p = 0; p < basins->n; p++)
    {
        pathval->val[p] = IFT_INFINITY_INT;
    }


    for (iftLabeledSet *S = seeds; S != NULL; S = S->next){
        int p = S->elem;
        label->val[p] = S->label;
        pathval->val[p] = 0;
        iftInsertGQueue(&Q, p);
    }


    // Image Foresting Transform

    while (!iftEmptyGQueue(Q))
    {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftGetVoxelCoord(basins, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(basins, v))
            {
                int q = iftGetVoxelIndex(basins, v);
                if (Q->L.elem[q].color != IFT_BLACK)
                {
                    int arcw = basins->val[q] * alpha + (1- alpha) * objbasins->val[q];

                    int tmp = iftMax(pathval->val[p], arcw);

                    if (tmp < pathval->val[q])  // For this path-value function,
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);
                        // this implies that q has never
                        // been inserted in Q.
                        label->val[q] = label->val[p];
                        pathval->val[q]  = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }
    iftDestroyGQueue(&Q);
    iftDestroyImage(&pathval);
    iftDestroyImage(&objbasins);

    iftCopyVoxelSize(basins, label);

    return (label);
}


iftImage *iftSuperPixelMajorityVote(iftImage *comp, iftImage* objmap, float threshold)
{
    iftImage *label = iftCreateImage(comp->xsize, comp->ysize, comp->zsize);

    int n_comp = iftMaximumValue(comp) + 1;
    float max = iftMaximumValue(objmap);

    /* positive are foreground, otherwise background */
    int *votes = iftAllocIntArray(n_comp);

    for (int p = 0; p < comp->n; p++)
    {
        float prob = objmap->val[p] / max;
        int suppxl = comp->val[p];
        if (prob > threshold) votes[suppxl] += 1;
        else votes[suppxl] -= 1;
    }

    for (int p = 0; p < comp->n; p++)
    {
        int suppxl = comp->val[p];
        if (votes[suppxl] > 0) label->val[p] = 1;
        else label->val[p] = 0;
    }

    iftFree(votes);

    return label;
}

iftImage *iftCreateJointMaskFromFileDir(char *maskDir, iftMaskJoinOperation operation, int maskVal)
{
    iftFileSet *fileset = iftLoadFileSetFromDirOrCSV(maskDir, 1, true);
    iftImage *mask = iftCreateJointMaskFromFileset(fileset, operation, maskVal);
    iftDestroyFileSet(&fileset);

    return mask;
}

iftImage *iftCreateJointMaskFromFileset(iftFileSet *maskFileset, iftMaskJoinOperation operation, int maskVal)
{
    iftImage *img0 = iftReadImageByExt(maskFileset->files[0]->path);
    int maxVal = iftMaximumValue(img0);
    iftImage *img0Bin = iftThreshold(img0, floor((float)maxVal/2.0), maxVal, maskVal);
    iftImage *mask = iftCopyImage(img0Bin);

    for(int i = 0; i < maskFileset->n; i++) {
        iftImage *img = iftReadImageByExt(maskFileset->files[i]->path);
        iftImage *imgBin = iftThreshold(img, floor((float)maxVal/2.0), maxVal, maskVal);

        /* apply the join operation */
        if(operation == IFT_MASK_JOIN_OPERATION_UNION)
        {
            #pragma omp parallel for
            for(int p = 0; p < img->n; p++) {
                mask->val[p] = iftMax(mask->val[p], imgBin->val[p]);
            }
        }
        else if (operation == IFT_MASK_JOIN_OPERATION_INTERSECTION)
        {
            #pragma omp parallel for
            for(int p = 0; p < img->n; p++) {
                mask->val[p] = iftMin(mask->val[p], imgBin->val[p]);
            }
        }

        iftDestroyImage(&img);
        iftDestroyImage(&imgBin);
    }

    iftDestroyImage(&img0);
    iftDestroyImage(&img0Bin);

    return mask;
}

iftImage *iftCreateSuperpixelBoundariesMask(iftImage *spixLabels, int maskVal)
{
    iftAdjRel *A = iftCircular(1.0);
    iftImage *mask = iftCreateImage(spixLabels->xsize, spixLabels->ysize, spixLabels->zsize);

    for (int p = 0; p < spixLabels->n; p++) {
        iftVoxel u = iftGetVoxelCoord(spixLabels,p);
        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftValidVoxel(spixLabels, v)) {
                int q = iftGetVoxelIndex(spixLabels, v);
                if (spixLabels->val[p] != spixLabels->val[q] && mask->val[q] != 255) {
                    mask->val[p] = maskVal;
                    break;
                }
            }
        }
    }

    return mask;
}


iftImage *iftBoundingBoxArrayToLabel(const iftImage *img, const iftBoundingBoxArray *bb_ary)
{
    iftImage *label = iftCreateImageFromImage(img);

#pragma omp parallel for
    for (int i = 0; i < bb_ary->n; i++) {
        iftBoundingBox *bb = &bb_ary->val[i];
        iftVoxel v;
        for (v.z = bb->begin.z; v.z < bb->end.z; v.z++) {
            for (v.y = bb->begin.y; v.y < bb->end.y; v.y++) {
                for (v.x = bb->begin.x; v.x < bb->end.x; v.x++) {
                    int p = iftGetVoxelIndex(img, v);
                    label->val[p] = i;
                }
            }
        }
    }

    return label;
}


iftSet *iftFindOptPath(const iftMImage *mimg, int src, int dst, float sigma, iftFImage *pathval, iftImage *pred)
{
    iftAdjRel *A = iftCircular(sqrtf(2.0f));
    iftAdjRel *L = iftLeftSide(A, 1.0f), *R = iftRightSide(A, 1.0f);

    pathval->val[src] = 0;
    pathval->val[dst] = IFT_INFINITY_FLT;

    iftFHeap *Q = iftCreateFHeap(pathval->n, pathval->val);
    iftInsertFHeap(Q, src);
    iftSet *colored = NULL;
    iftInsertSet(&colored, src);

    while (!iftEmptyFHeap(Q))
    {
        int p = iftRemoveFHeap(Q);

        if (p == dst)
            break;

        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (!iftValidVoxel(pathval, v))
                continue;
            int q = iftGetVoxelIndex(pathval, v);
            if (Q->color[q] == IFT_BLACK)
                continue;

            v = iftGetAdjacentVoxel(L, u, i);
            int left = (iftValidVoxel(pathval, v) ? iftGetVoxelIndex(pathval, v) : p);

            v = iftGetAdjacentVoxel(R, u, i);
            int right = (iftValidVoxel(pathval, v) ? iftGetVoxelIndex(pathval, v) : p);

            double fdist = iftMImageDist(mimg, left, right);

            double arc_w = pathval->val[p] + iftApproxExp( - fdist / sigma);

            if (arc_w < pathval->val[q])
            {
                pathval->val[q] = (float) arc_w;
                pred->val[q] = p;
                if (Q->color[q] == IFT_GRAY) {
                    iftGoUpFHeap(Q, Q->pos[q]);
                } else {
                    iftInsertFHeap(Q, q);
                    iftInsertSet(&colored, q);
                }
            }
        }
    }

    iftSet *path = NULL;
    for (int p = dst; p != IFT_NIL; p = pred->val[p])
        iftInsertSet(&path, p);

    // reseting pathval and pred
    for (iftSet *S = colored; S; S = S->next) {
        pathval->val[S->elem] = IFT_INFINITY_FLT;
        pred->val[S->elem] = IFT_NIL;
    }

    // marking optimum path (blocking for no crossing)
    for (iftSet *S = path; S; S = S->next) {
        pathval->val[S->elem] = IFT_INFINITY_FLT_NEG;
    }

    iftDestroySet(&colored);
    iftDestroyFHeap(&Q);
    iftDestroyAdjRel(&A);
    iftDestroyAdjRel(&L);
    iftDestroyAdjRel(&R);

    return path;
}


iftSet *iftFindOptPathWithProb(const iftMImage *mimg, const iftFImage *obj, const iftFImage *bkg, int src, int dst,
                               float sigma, float gamma, iftFImage *pathval, iftImage *pred)
{
    iftAdjRel *A = iftCircular(sqrtf(2.0f));
    iftAdjRel *L = iftLeftSide(A, 1.0f), *R = iftRightSide(A, 1.0f);

    pathval->val[src] = 0;
    pathval->val[dst] = IFT_INFINITY_FLT;

    iftFHeap *Q = iftCreateFHeap(pathval->n, pathval->val);
    iftInsertFHeap(Q, src);
    iftSet *colored = NULL;
    iftInsertSet(&colored, src);

    while (!iftEmptyFHeap(Q))
    {
        int p = iftRemoveFHeap(Q);

        if (p == dst)
            break;

        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (!iftValidVoxel(pathval, v))
                continue;
            int q = iftGetVoxelIndex(pathval, v);
            if (Q->color[q] == IFT_BLACK)
                continue;

            v = iftGetAdjacentVoxel(L, u, i);
            int left = (iftValidVoxel(pathval, v) ? iftGetVoxelIndex(pathval, v) : p);

            v = iftGetAdjacentVoxel(R, u, i);
            int right = (iftValidVoxel(pathval, v) ? iftGetVoxelIndex(pathval, v) : p);

            double fdist = iftMImageDist(mimg, left, right);
            float  odist = fabsf(obj->val[right] - obj->val[left]);
            float  bdist = fabsf(bkg->val[right] - bkg->val[left]);

            double arc_w = pathval->val[p] + iftApproxExp( - fdist / sigma - odist / gamma - bdist / gamma);

            if (arc_w < pathval->val[q])
            {
                pathval->val[q] = (float) arc_w;
                pred->val[q] = p;
                if (Q->color[q] == IFT_GRAY) {
                    iftGoUpFHeap(Q, Q->pos[q]);
                } else {
                    iftInsertFHeap(Q, q);
                    iftInsertSet(&colored, q);
                }
            }
        }
    }

    iftSet *path = NULL;
    for (int p = dst; p != IFT_NIL; p = pred->val[p])
        iftInsertSet(&path, p);

    // reseting pathval and pred
    for (iftSet *S = colored; S; S = S->next) {
        pathval->val[S->elem] = IFT_INFINITY_FLT;
        pred->val[S->elem] = IFT_NIL;
    }

    // marking optimum path (blocking for no crossing)
    for (iftSet *S = path; S; S = S->next) {
        pathval->val[S->elem] = IFT_INFINITY_FLT_NEG;
    }

    iftDestroySet(&colored);
    iftDestroyFHeap(&Q);
    iftDestroyAdjRel(&A);
    iftDestroyAdjRel(&L);
    iftDestroyAdjRel(&R);

    return path;
}


// private
static int _get_direction(const iftImage *obj, const iftAdjRel *A, int src)
{
    iftVoxel u = iftGetVoxelCoord(obj, src);

    // find object pixel
    int j = 1;
    for (; j < A->n; j++) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, j);
        if (iftValidVoxel(obj, v)) {
            int q = iftGetVoxelIndex(obj, v);
            if (obj->val[q])
                break;
        }
    }

    for (int i = 0; i < A->n; i++, j = (j + 1) % A->n) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, j);
        if (iftValidVoxel(obj, v)) {
            int q = iftGetVoxelIndex(obj, v);
            if (!obj->val[q])
                break;
        } else break;
    }

    return j;
}


iftImage *iftSinglePathContour(const iftImage *label)
{
    // pre processing
    iftAdjRel *disk = iftCircular(2.0);
    iftAdjRel *A = iftClockCircular(sqrtf(2.0f));
    iftImage *opened = iftOpen(label, disk, NULL);
    iftImage *single_comp = iftSelectLargestComp(opened, A);
    iftImage *contour = iftBorderImage(single_comp, true);

    iftDestroyImage(&single_comp);
    iftDestroyImage(&opened);
    iftDestroyAdjRel(&disk);

    iftImage *pred = iftCreateImage(contour->xsize, contour->ysize, contour->zsize);
    iftImage *dist = iftCreateImage(contour->xsize, contour->ysize, contour->zsize);
    iftLIFO *Q = iftCreateLIFO(contour->n);

    for (int i = 0; i < contour->n; i++)
        pred->val[i] = IFT_NIL;

    int src = IFT_NIL;
    for (int i = 0; i < contour->n; i++) {
        if (contour->val[i]) {
            src = i;
            break;
        }
    }

    iftImage *out = iftCreateImage(contour->xsize, contour->ysize, contour->zsize);
    if (src == IFT_NIL) {
        iftWarning("No contour found.", "iftReduceContour");
        goto free_memory;
    }

    int ref_orientation = _get_direction(label, A, src);
    iftInsertLIFO(Q, src);
    while (!iftEmptyLIFO(Q))
    {
        int p = iftRemoveLIFO(Q);
        iftVoxel u = iftGetVoxelCoord(contour, p);

        int j = ref_orientation; // because of the selection of the starting point (no predecessor)
        if (pred->val[p] != IFT_NIL) {
            iftVoxel t = iftGetVoxelCoord(pred, pred->val[p]);
            j = ift6NeighborsClockFromLeft(&t, &u);
        }
        iftSet *adj = NULL;
        for (int i = 0; i < A->n; i++, j = (j + 1) % A->n)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, j);
            if (!iftValidVoxel(contour, v))
                continue;
            int q = iftGetVoxelIndex(contour, v);
            if (contour->val[q] != 0 && Q->color[q] != IFT_BLACK &&
                dist->val[q] <= dist->val[p])
            {
                dist->val[q] = dist->val[p] + 1;
                pred->val[q] = p;
                if (Q->color[q] == IFT_WHITE)
                    iftInsertSet(&adj, q);
            }
        }

        for (iftSet *s = adj; s; )
        {
            iftInsertLIFO(Q, s->elem);
            iftSet *aux = s;
            s = s->next;
            iftFree(aux);
        }
    }

    int p = src;
    iftVoxel u = iftGetVoxelCoord(contour, src);
    for (int i = 1; i < A->n; i++) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, i);
        if (iftValidVoxel(contour, v)) {
            int q = iftGetVoxelIndex(contour, v);
            if (dist->val[p] < dist->val[q])
                p = q;
        }
    }

    for ( ; p != IFT_NIL; p = pred->val[p]) {
        out->val[p] = contour->val[src];
    }

    free_memory:
    iftDestroyImage(&contour);
    iftDestroyLIFO(&Q);
    iftDestroyImage(&dist);
    iftDestroyImage(&pred);
    iftDestroyAdjRel(&A);

    return out;
}


iftImage *iftFillContour(const iftImage *contour)
{
    iftImage *label = iftCopyImage(contour);

    iftSetImage(label, 1);

    iftFIFO *Q = iftCreateFIFO(label->n);

    iftAdjRel *A = NULL;
    if (iftIs3DImage(contour)) {
        A = iftSpheric(1.0);
        for (int z = 0, p = 0; z < contour->zsize; z++) {
            for (int y = 0; y < contour->ysize; y++) {
                for (int x = 0; x < contour->xsize; x++, p++) {
                    if (!contour->val[p] && (y == 0 || x == 0 || z == 0 ||
                                             x == (contour->xsize - 1) ||
                                             y == (contour->ysize - 1) ||
                                             z == (contour->zsize - 1))) {
                        iftInsertFIFO(Q, p);
                    }
                }
            }
        }
    } else {
        A = iftCircular(1.0);
        for (int y = 0, p = 0; y < contour->ysize; y++) {
            for (int x = 0; x < contour->xsize; x++, p++) {
                if (!contour->val[p] && (y == 0 || x == 0 ||
                                         x == (contour->xsize - 1) ||
                                         y == (contour->ysize - 1))) {
                    iftInsertFIFO(Q, p);
                }
            }
        }
    }

    while (!iftEmptyFIFO(Q))
    {
        int p = iftRemoveFIFO(Q);
        label->val[p] = 0;
        iftVoxel u = iftGetVoxelCoord(label, p);

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (!iftValidVoxel(label, v))
                continue;
            int q = iftGetVoxelIndex(label, v);
            if (Q->color[q] == IFT_WHITE && contour->val[p] == contour->val[q])
                iftInsertFIFO(Q, q);
        }
    }

    iftDestroyFIFO(&Q);
    iftDestroyAdjRel(&A);

    return label;
}

iftImage *iftDISF(iftMImage *mimg, iftAdjRel *A, int num_init_seeds, int num_superpixels, iftImage *mask)
{    
  //    int num_sampl_seeds;
    int *root_map;
    double *cost_map;
    iftImage *label_img, *mask_img;
    iftDHeap *heap;
    iftSet *seeds;
    float decrement_rate = (float) num_superpixels / (float) num_init_seeds;
    

    // INIT =======================================================================
    label_img = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    root_map = (int*)calloc(mimg->n, sizeof(int));
    cost_map = (double*)calloc(mimg->n, sizeof(double));
    heap = iftCreateDHeap(mimg->n, cost_map);

    if(mask == NULL) mask_img = iftSelectImageDomain(mimg->xsize, mimg->ysize, mimg->zsize);
    else mask_img = mask;

    // GRID SAMPLING ==============================================================
    iftImage *seed_img;

    seeds = NULL;
    seed_img = iftGridSampling(mimg, mask_img, num_init_seeds);

    for(int p_index = 0; p_index < mimg->n; ++p_index)
    {
        if(seed_img->val[p_index] > 0) iftInsertSet(&seeds, p_index); 
        if(mask_img->val[p_index] == 0)
          cost_map[p_index] = IFT_INFINITY_DBL_NEG;
    }

    //    num_sampl_seeds = iftSetSize(seeds);
    
    iftDestroyImage(&seed_img);
    if(mask == NULL) iftDestroyImage(&mask_img);

    // DISF ========================================================================
    int it, num_seeds;

    it = 0;
    do
    {
        int num_maintain;
        int *superpixel_area, *seed_index_lut;
        bool **is_tree_adj;
        float **superpixel_feats;
        iftSet **superpixel_adj;

        num_seeds = iftSetSize(seeds);

        // Init --------------------------------------------------------------------
        seed_index_lut = (int*)calloc(num_seeds, sizeof(int));
        superpixel_area = (int *)calloc(num_seeds, sizeof(int));
        is_tree_adj = (bool**)calloc(num_seeds, sizeof(bool*));
        superpixel_feats = (float **)calloc(num_seeds, sizeof(float *));
        superpixel_adj = (iftSet **)calloc(num_seeds, sizeof(iftSet *));
        
        #pragma omp parallel for
        for(int i = 0; i < num_seeds; i++)
        {   
            superpixel_feats[i] = (float *)calloc(mimg->m, sizeof(float));
            is_tree_adj[i] = (bool *)calloc(num_seeds, sizeof(bool));

            for(int j = 0; j < mimg->m; ++j)
                superpixel_feats[i][j] = 0;

            for(int j = 0; j < num_seeds; ++j)
                is_tree_adj[i][j] = false;

            superpixel_adj[i] = NULL;
        }

        #pragma omp parallel for
        for(int p_index = 0; p_index < mimg->n; p_index++)
        {
            root_map[p_index] = IFT_NIL;
            label_img->val[p_index] = 0;
            
            if(cost_map[p_index] != IFT_INFINITY_DBL_NEG)
                cost_map[p_index] = IFT_INFINITY_DBL;
        }

        int seed_count;
        seed_count = 0;
        for(iftSet *s = seeds; s != NULL; s = s->next)
        {
            root_map[s->elem] = s->elem;
            cost_map[s->elem] = 0.0;
            seed_index_lut[seed_count] = s->elem;
            label_img->val[s->elem] = seed_count + 1;
            seed_count++;

            iftInsertDHeap(heap, s->elem);
        }

        // IFT ---------------------------------------------------------------------
        while(!iftEmptyDHeap(heap))
        {
            int p_index, p_label;
            iftVoxel p_voxel;

            p_index = iftRemoveDHeap(heap);
            p_voxel = iftMGetVoxelCoord(mimg, p_index);
            p_label = label_img->val[p_index];

            // Add spel p to its tree
            superpixel_area[p_label - 1]++;
            for(int j = 0; j < mimg->m; j++) superpixel_feats[p_label - 1][j] += mimg->val[p_index][j];

            // For each adjacent
            for(int i = 1; i < A->n; i++)
            {
                iftVoxel q_voxel;

                q_voxel = iftGetAdjacentVoxel(A, p_voxel, i);

                if(iftMValidVoxel(mimg, q_voxel))
                {
                    int q_index, q_label;

                    q_index = iftMGetVoxelIndex(mimg, q_voxel);
                    q_label = label_img->val[q_index];

                    if(heap->color[q_index] != IFT_BLACK)
                    {
                        double cost, color_dist;

                        color_dist = 0.0;
                        for(int j = 0; j < mimg->m; j++)
                        {
                            double tmp;

                            tmp = mimg->val[q_index][j] - (superpixel_feats[p_label - 1][j]/(float)superpixel_area[p_label-1]);

                            color_dist += tmp*tmp;
                        }
                        color_dist = sqrtf(color_dist);

                        cost = iftMax(cost_map[p_index], color_dist);

                        if(cost < cost_map[q_index])
                        {
                            cost_map[q_index] = cost;
                            root_map[q_index] = root_map[p_index];
                            label_img->val[q_index] = p_label;

                            if(heap->color[q_index] == IFT_GRAY) iftGoUpDHeap(heap, heap->pos[q_index]);
                            else iftInsertDHeap(heap, q_index);
                        }
                    }
                    else if(q_label != p_label)
                    {
                        if(!is_tree_adj[p_label-1][q_label - 1])
                        {
                            iftInsertSet(&superpixel_adj[p_label-1], q_label-1);
                            iftInsertSet(&superpixel_adj[q_label-1], p_label-1);
                            is_tree_adj[p_label-1][q_label-1] = true;
                            is_tree_adj[q_label-1][p_label-1] = true;
                        }
                    }
                }
            }
        }

        // Seed selection ----------------------------------------------------------
        iftDHeap *prio_heap;
        double *superpixel_prior;

        superpixel_prior = (double *)calloc(num_seeds, sizeof(double));
        prio_heap = iftCreateDHeap(num_seeds, superpixel_prior);

        iftSetRemovalPolicyDHeap(prio_heap, MAXVALUE);
        
        // Estimate the number of seeds to be removed 
        //num_maintain = iftMax(num_sampl_seeds * exp(-(it + 1)), num_superpixels);
	num_maintain  = iftMax(num_seeds*(1.0-decrement_rate),num_superpixels);
	  
        // Compute the medoids
        for(iftSet *s = seeds; s != NULL; s = s->next)
        {
            int s_index, s_label;

            s_index = s->elem;
            s_label = label_img->val[s_index];

            for(int j = 0; j < mimg->m; j++) 
                superpixel_feats[s_label - 1][j] /= superpixel_area[s_label - 1];
        }
        
        // Compute the relevance of each seed
        for(iftSet *s = seeds; s != NULL; s = s->next)
        {
            int s_index, s_label;
            double area_prio, grad_prio;

            s_index = s->elem;
            s_label = label_img->val[s_index];

            area_prio = superpixel_area[s_label - 1]/(double)mimg->n;

            grad_prio = IFT_INFINITY_DBL;
            for(iftSet *t = superpixel_adj[s_label - 1]; t != NULL; t = t->next)
            {
                int t_label;
                double grad;

                t_label = t->elem;

                grad = iftFeatDistance(superpixel_feats[s_label - 1], superpixel_feats[t_label], mimg->m);

                if(grad_prio > grad) grad_prio = grad;
            }

            superpixel_prior[s_label - 1] = area_prio * grad_prio;

            iftInsertDHeap(prio_heap, s_label - 1);
        }

        // Select the most relevant seeds
        iftDestroySet(&seeds);

        for(int i = 0; i < num_maintain; i++) 
        {
            int s_index, s_label;

            s_label = iftRemoveDHeap(prio_heap);
            s_index = seed_index_lut[s_label];

            iftInsertSet(&seeds, s_index);
        }
        
        ++it;

        // Clear -------------------------------------------------------------------
        for(int i = 0; i < num_seeds; i++) 
        {
            free(superpixel_feats[i]);
            iftDestroySet(&superpixel_adj[i]);
            free(is_tree_adj[i]);
        }
        free(is_tree_adj);
        free(superpixel_feats);
        free(superpixel_adj);
        free(seed_index_lut);
        free(superpixel_area);
        free(superpixel_prior);
        iftDestroyDHeap(&prio_heap);
        iftResetDHeap(heap);
    }while(num_seeds > num_superpixels);

    // CLEAR =======================================================================
    free(root_map);
    free(cost_map);
    iftDestroyDHeap(&heap);
    iftDestroySet(&seeds);

    return label_img;
}
