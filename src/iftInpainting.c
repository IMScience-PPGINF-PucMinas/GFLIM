#include "iftInpainting.h"

#include "ift/core/dtypes/GQueue.h"
#include "ift/core/io/Stream.h"
#include <xmmintrin.h>
#include <emmintrin.h>


/*----------------- Private structures and methods -------------------- */

typedef struct ift_texture_patch {
    float *feat;
    float *weight;
} iftTexturePatch;

typedef struct ift_fast_texture_comparison_data {
    iftMatrix *valid_feat;
    float *features;
} iftFastTextureComparisonData;

/* Texture patch storage functions */

iftTexturePatch *iftCreateTexturePatch(int n);
void iftDestroyTexturePatch(iftTexturePatch **patch);

/* Texture sampling functions */
iftImage          *iftCreateExternalShell(iftImage *mask, float adj_radius1, float adj_radius2);
iftImage 			*iftCreateInternalShell(iftImage *mask, float adj_radius1, float adj_radius2);
iftTextureSamples *iftTextureSamplesAroundMaskBoundaryOutside(iftImage *img, iftImage *mask, float adj_radius);
iftTextureSamples *iftTextureSamplesAroundMaskBoundaryInside(iftImage *img, iftImage *mask, float adj_radius);

/* Inpainting functions */
iftFImage         *iftConfidenceMap(iftImage *mask, iftBMap *filled);
iftSet            *iftBoundaryOfFillingRegion(iftFImage *confidence, iftBMap *filled);
void             iftVoxelTexture(int p, iftImage *img, iftAdjRel *A, iftTexturePatch *patch, iftBMap *filled);
void 				   iftBoundaryVoxelTexture(int p, iftImage *img, iftImage *mask, iftAdjRel *A, iftTexturePatch *patch);
void 					iftFindBestTextureMatch(iftTexturePatch *patch, iftTexturePatch *best_patch, iftFastTextureComparisonData *texture);
void 					iftMapTexture(iftTexturePatch *patch, int p, iftAdjRel *A, iftImage *img, iftBMap *filled);
float            iftVoxelConfidence(int p, iftFImage *confidence, iftAdjRel *A, iftBMap *filled);
float            iftImageTerm(int p, iftImage *img, iftAdjRel *A, iftAdjRel *Abdr, iftAdjRel *Aimg, iftBMap *filled);
void 					iftSetPriorityMapAndHeap(iftFHeap *H, iftFImage *priority, iftSet *B, iftImage *img, iftFImage *confidence, iftAdjRel *A, iftBMap *filled);
void 					iftUpdatePatchBoundaryPriority(int p, iftFHeap *H, iftFImage *priority, iftImage *img, iftFImage *confidence, iftAdjRel *A, iftBMap *filled);
void             iftUpdateConfidenceAndFilledMaps(int p, iftFImage *confidence, iftAdjRel *A, iftBMap *filled);

/* Optimization functions */
iftFastTextureComparisonData *iftCreateFastTextureComparisonData(iftTextureSamples *samples);
void 									iftDestroyFastTextureComparisonData(iftFastTextureComparisonData **samples);
iftMatrix 			*iftValidTextureSamplesToMatrix(iftTextureSamples *samples);
float* 				 iftAlignPatchFeatures(iftMatrix *texture);
int 					  iftFindBestTextureMatchIndex(iftTexturePatch * __restrict__ patch, float * __restrict__ _features, int nsamples, int nfeats);


iftTexturePatch *iftCreateTexturePatch(int n)
{
    iftTexturePatch *patch = (iftTexturePatch*)iftAlloc(1,sizeof(iftTexturePatch));

    patch->feat = iftAllocFloatArray(n);
    patch->weight = iftAllocFloatArray(n);

    return patch;
}

void iftDestroyTexturePatch(iftTexturePatch **patch)
{
    if(patch != NULL && *patch != NULL)
    {
        iftFree((*patch)->feat);
        iftFree((*patch)->weight);

        iftFree(*patch);

        *patch = NULL;
    }
}


/* Create list for texture samples */

iftTextureSamples *iftCreateTextureSamples(int nsamples, int nfeats)
{
    iftTextureSamples *texture = (iftTextureSamples *) iftAlloc(1,sizeof(iftTextureSamples));
    int s;

    texture->sample   = (iftTexture *) iftAlloc(nsamples,sizeof(iftTexture));
    texture->nsamples = nsamples;

    for (s=0; s < nsamples; s++) {
        texture->sample[s].feat  = iftAllocIntArray(nfeats);
        texture->sample[s].valid = 1;
    }
    texture->nfeats   = nfeats;

    return(texture);
}

/* Destroy list of texture samples */

void iftDestroyTextureSamples(iftTextureSamples **texture)
{
    if (texture != NULL && *texture != NULL) {
        for (int s=0; s < (*texture)->nsamples; s++)
            iftFree((*texture)->sample[s].feat);
        iftFree((*texture)->sample);
        iftFree((*texture));
        *texture = NULL;
    }
}

/* Define an external shell around a given mask */

iftImage *iftCreateExternalShell(iftImage *mask, float adj_radius1, float adj_radius2)
{
    iftSet   *S=NULL;
    iftImage *dil_mask1, *dil_mask2, *shell;

    if ((iftMaximumValue(mask)!=1)||(iftMinimumValue(mask)!=0))
        iftError("It requires a binary mask", "iftCreateExternalShell");

    dil_mask1 = iftDilateBin(mask,&S,adj_radius1);
    iftDestroySet(&S);
    dil_mask2 = iftDilateBin(mask,&S,adj_radius2);
    iftDestroySet(&S);
    shell     = iftSub(dil_mask2,dil_mask1);
    iftDestroyImage(&dil_mask1);
    iftDestroyImage(&dil_mask2);

    return(shell);
}

iftImage *iftCreateInternalShell(iftImage *mask, float adj_radius1, float adj_radius2)
{
    iftSet   *S=NULL;
    iftImage *ero_mask1, *ero_mask2, *shell;

    if ((iftMaximumValue(mask)!=1)||(iftMinimumValue(mask)!=0))
        iftError("It requires a binary mask", "iftCreateInternalShell");

    ero_mask1 = iftErodeBin(mask,&S,adj_radius1);
    iftDestroySet(&S);
    ero_mask2 = iftErodeBin(mask,&S,adj_radius2);
    iftDestroySet(&S);
    shell     = iftSub(ero_mask1,ero_mask2);
    iftDestroyImage(&ero_mask1);
    iftDestroyImage(&ero_mask2);

    return(shell);
}

/* Extract texture samples outside the binary mask from a shell around it ranging
	from adj_radius to 5*adj_radius pixels. We start from adj_radius pixels
	to ensure that textures touching the filling area's border be added for
	inpainting.
 */
iftTextureSamples *iftTextureSamplesAroundMaskBoundaryOutside(iftImage *img, iftImage *mask, float adj_radius)
{
    iftAdjRel         *A=NULL;
    iftImage          *shell;
    iftTextureSamples *texture=NULL;

    /* Create adjacency relation for the texture template */

    if (iftIs3DImage(img))
        A = iftSpheric(adj_radius);
    else
        A = iftCircular(adj_radius);

    /* Create a shell around mask to extract external texture samples */

    /* We start from adj_radius because iftTextureSamplesOnMask does not
        select as samples patches that overlap the area not selected
        for texture sampling (i.e., patches centered at pixels with
        0-valued neighbors).
     */
    shell = iftCreateExternalShell(mask, adj_radius, 5*adj_radius);

    texture = iftTextureSamplesOnMask(img, shell, A);

    iftDestroyAdjRel(&A);
    iftDestroyImage(&shell);

    return texture;
}

/* This function samples texture inside the mask's boundary */

iftTextureSamples *iftTextureSamplesAroundMaskBoundaryInside(iftImage *img, iftImage *mask, float adj_radius)
{
    iftTextureSamples *texture=NULL;
    iftImage          *shell;
    iftAdjRel         *A=NULL;

    /* Create adjacency relation for the texture template */

    if (iftIs3DImage(img))
        A = iftSpheric(adj_radius);
    else
        A = iftCircular(adj_radius);

    /* Create a shell inside the mask to extract internal texture samples */

    shell = iftCreateInternalShell(mask, adj_radius, 5*adj_radius);
    texture = iftTextureSamplesOnMask(img, shell, A);

    iftDestroyAdjRel(&A);
    iftDestroyImage(&shell);

    return(texture);

}

/* This function samples texture on the foreground portion of the mask. A sample is
	only valid if it is *entirely* contained within the mask, hence, those near
	the mask's boundary are disconsidered because they may overlap with the filling
	area. This function can be used to determine that the entire background should
	be used for sampling, for example.
 */

iftTextureSamples *iftTextureSamplesOnMask(iftImage *img, iftImage *mask, iftAdjRel *A)
{
    int                s,i,j,p,q,nsamples,nfeats;
    iftVoxel           u,v;
    iftTextureSamples *texture=NULL;

    /* Compute the size of the mask */

    nsamples = 0;
    for (p=0; p < mask->n; p++)
        if (mask->val[p]!=0)
            nsamples++;


    /* Extract texture samples for filled/grayscale images */

    if (iftIsColorImage(img)) { /* filled Image */

        nfeats   = 3*A->n;
        texture  = iftCreateTextureSamples(nsamples,nfeats);

        s = 0;
        for (p=0; p < mask->n; p++) {
            if (mask->val[p]!=0){
                u = iftGetVoxelCoord(img,p);
                for (i=0,j=0; i < A->n; i++,j+=3) {
                    v = iftGetAdjacentVoxel(A,u,i);
                    if (iftValidVoxel(img,v)){
                        q = iftGetVoxelIndex(img,v);
                        // Only patches entirely contained in the marked region
                        // are used for sampling
                        if(mask->val[q] != 0)
                        {
                            texture->sample[s].feat[j]   = img->val[q];
                            texture->sample[s].feat[j+1] = img->Cb[q];
                            texture->sample[s].feat[j+2] = img->Cr[q];
                        }
                        else
                        {
                            texture->sample[s].valid = 0;
                            break;
                        }
                    }else{
                        texture->sample[s].valid = 0;
                        break;
                    }
                }
                s++;
            }
        }

    } else { /* GrayScale Image */

        nfeats      = A->n;
        texture     = iftCreateTextureSamples(nsamples,nfeats);

        s  = 0;
        for (p=0; p < mask->n; p++) {
            if (mask->val[p]!=0){
                u = iftGetVoxelCoord(img,p);
                for (i=0; i < A->n; i++) {
                    v = iftGetAdjacentVoxel(A,u,i);
                    if (iftValidVoxel(img,v)){
                        q = iftGetVoxelIndex(img,v);
                        // Only patches entirely contained in the marked region
                        // are used for sampling
                        if(mask->val[q] != 0)
                        {
                            texture->sample[s].feat[i]   = img->val[q];
                        }
                        else
                        {
                            texture->sample[s].valid = 0;
                            break;
                        }
                    }else{
                        texture->sample[s].valid = 0;
                        break;
                    }
                }
                s++;
            }
        }
    }


    return(texture);

}


/* Create confidence map */

iftFImage *iftConfidenceMap(iftImage *mask, iftBMap *filled)
{
    iftFImage *confidence = iftCreateFImage(mask->xsize, mask->ysize, mask->zsize);

    if ((iftMaximumValue(mask)!=1)||(iftMinimumValue(mask)!=0))
        iftError("It requires a binary mask", "iftConfidenceMap");

    for (int p=0; p < mask->n; p++)
    {
        confidence->val[p] = 1.0 - mask->val[p];

        if(confidence->val[p] > IFT_EPSILON)
            iftBMapSet1(filled, p);
    }


    return(confidence);
}

/* Compute boundary outside region under filling */

iftSet *iftBoundaryOfFillingRegion(iftFImage *confidence, iftBMap *filled)
{
    int        i,p,q;
    iftVoxel   u,v;
    iftAdjRel *A;
    iftSet    *B=NULL;

    if (iftIs3DFImage(confidence))
        A = iftSpheric(1.0);
    else
        A = iftCircular(1.0);

    for (p=0; p < confidence->n; p++) {
        if (iftBMapValue(filled, p)){
            u = iftFGetVoxelCoord(confidence,p);
            for (i=1; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A,u,i);
                if (iftFValidVoxel(confidence,v)){
                    q = iftFGetVoxelIndex(confidence,v);
                    if (!iftBMapValue(filled, q)){
                        iftInsertSet(&B,p);
                        break;
                    }
                }
            }
        }
    }

    iftDestroyAdjRel(&A);
    return(B);
}

/* Copies the filled portions of a patch into a feature vector */

void iftVoxelTexture(int p, iftImage *img, iftAdjRel *A, iftTexturePatch *patch, iftBMap *filled)
{
    int i, j, q;
    iftVoxel u, v;

    u = iftGetVoxelCoord(img,p);
    if (iftIsColorImage(img)){ /* filled Image */
        for (i=0,j=0; i < A->n; i++,j+=3) {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(img,v)){
                q = iftGetVoxelIndex(img,v);
                if (iftBMapValue(filled,q)){
                    patch->feat[j]   = img->val[q];
                    patch->feat[j+1] = img->Cb[q];
                    patch->feat[j+2] = img->Cr[q];

                    patch->weight[j] = patch->weight[j+1] = patch->weight[j+2] = 1.0;
                }else{
                    patch->feat[j]   = patch->feat[j+1] = patch->feat[j+2] = 0.0;
                    patch->weight[j] = patch->weight[j+1] = patch->weight[j+2] = 0.0;
                }
            }
            else  // Ensuring that pixels outside the image do not influence
                // the euclidean distance computation
            {
                patch->feat[j]   = patch->feat[j+1] = patch->feat[j+2] = 0.0;
                patch->weight[j] = patch->weight[j+1] = patch->weight[j+2] = 0.0;
            }
        }
    } else { /* GrayScale Image */
        for (i=0; i < A->n; i++) {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(img,v)){
                q = iftGetVoxelIndex(img,v);
                if (iftBMapValue(filled,q)){
                    patch->feat[i] = img->val[q];
                    patch->weight[i] = 1.0;
                }else{
                    patch->feat[i] = 0.0;
                    patch->weight[i] = 0.0;
                }
            }
            else  // Ensuring that pixels outside the image do not influence
                // the euclidean distance computation
            {
                patch->feat[i] = 0.0;
                patch->weight[i] = 0.0;
            }
        }
    }
}

/* Copies the filled portions of a boundary patch into a feature vector */

void iftBoundaryVoxelTexture(int p, iftImage *img, iftImage *mask, iftAdjRel *A, iftTexturePatch *patch)
{
    int i, j, q;
    iftVoxel u, v;

    u = iftGetVoxelCoord(img,p);
    if (iftIsColorImage(img)){ /* filled Image */
        for (i=0,j=0; i < A->n; i++,j+=3) {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(img,v)){
                q = iftGetVoxelIndex(img,v);
                if (mask->val[q]>0){
                    patch->feat[j]   = img->val[q];
                    patch->feat[j+1] = img->Cb[q];
                    patch->feat[j+2] = img->Cr[q];

                    patch->weight[j] = patch->weight[j+1] = patch->weight[j+2] = 1.0;
                }else{
                    patch->feat[j]   = patch->feat[j+1] = patch->feat[j+2] = 0.0;
                    patch->weight[j] = patch->weight[j+1] = patch->weight[j+2] = 0.0;
                }
            }
            else  // Ensuring that pixels outside the image do not influence
                // the euclidean distance computation
            {
                patch->feat[j]   = patch->feat[j+1] = patch->feat[j+2] = 0.0;
                patch->weight[j] = patch->weight[j+1] = patch->weight[j+2] = 0.0;
            }
        }
    } else { /* GrayScale Image */
        for (i=0; i < A->n; i++) {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(img,v)){
                q = iftGetVoxelIndex(img,v);
                if (mask->val[q]>0){
                    patch->feat[i] = img->val[q];
                    patch->weight[i] = 1.0;
                }else{
                    patch->feat[i] = 0.0;
                    patch->weight[i] = 0.0;
                }
            }
            else  // Ensuring that pixels outside the image do not influence
                // the euclidean distance computation
            {
                patch->feat[i] = 0.0;
                patch->weight[i] = 0.0;
            }
        }
    }
}

void iftFindBestTextureMatch(iftTexturePatch *patch, iftTexturePatch *best_patch, iftFastTextureComparisonData *texture)
{
    int i, best_s = IFT_NIL, nfeats, nsamples;
    iftMatrix *valid_feat = texture->valid_feat;
    float * __restrict__ features = texture->features;

    nfeats 	= valid_feat->ncols;
    nsamples = valid_feat->nrows;

    best_s = iftFindBestTextureMatchIndex(patch, features, nsamples, nfeats);

    for (i=0; i < nfeats; i++) {
        best_patch->feat[i]=valid_feat->val[iftGetMatrixIndex(valid_feat, i, best_s)];
    }
}

/* Maps the texture from the source template to the current unfilled regions
	of the patch
*/
void iftMapTexture(iftTexturePatch *patch, int p, iftAdjRel *A, iftImage *img, iftBMap *filled)
{
    int i, j, q;
    iftVoxel u, v;

    u = iftGetVoxelCoord(img,p);
    if (iftIsColorImage(img)){ /* filled Image */
        for (i=0,j=0; i < A->n; i++,j+=3) {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(img,v)){
                q = iftGetVoxelIndex(img,v);
                if(!iftBMapValue(filled, q))
                {
                    img->val[q]	= patch->feat[j];
                    img->Cb[q] 	= patch->feat[j+1];
                    img->Cr[q] 	= patch->feat[j+2];
                }
            }
        }
    } else { /* GrayScale Image */
        for (i=0; i < A->n; i++) {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(img,v)){
                q = iftGetVoxelIndex(img,v);

                if(!iftBMapValue(filled, q))
                {
                    img->val[q] = patch->feat[i];
                }
            }
        }
    }
}

iftVector iftComputeImageGradientOnFilledPixel(int q, iftImage *img, iftBMap *filled,
                                               iftAdjRel *B, float *mag)
{
    int j;
    iftVector i;
    iftVoxel v;

    v = iftGetVoxelCoord(img, q);

    i.x = i.y = i.z = 0.0;

    if(iftBMapValue(filled, q))
    {
        for(j = 1; j < B->n; j++)
        {
            iftVoxel w;

            w = iftGetAdjacentVoxel(B, v, j);

            if(iftValidVoxel(img, w))
            {
                int r = iftMGetVoxelIndex(img, w);

                if(iftBMapValue(filled, r))
                {
                    float diff = (img->val[r]-img->val[q]);
                    i.x += (diff*B->dx[j]/mag[j]);
                    i.y += (diff*B->dy[j]/mag[j]);
                    i.z += (diff*B->dz[j]/mag[j]);
                }
            }
        }
    }

    return i;
}

/* Estimate the confidence value for a voxel */

float iftVoxelConfidence(int p, iftFImage *confidence, iftAdjRel *A, iftBMap *filled)
{
    int        i,q;
    iftVoxel   u,v;
    float      val=0.0;

    u = iftFGetVoxelCoord(confidence,p);
    for (i=0; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A,u,i);
        if (iftFValidVoxel(confidence,v)){
            q = iftFGetVoxelIndex(confidence,v);

            if(iftBMapValue(filled, q))
                val += confidence->val[q];
        }
    }
    val /= A->n;

    return(val);
}

/* Estimate the image term for a voxel */

float iftImageTerm(int p, iftImage *img, iftAdjRel *A, iftAdjRel *Abdr, iftAdjRel *Aimg, iftBMap *filled)
{
    int        i,q;
    iftVoxel   u,v;
    iftVector  n,cur,Imax;
    float      mag_bdr[Abdr->n], mag_img[Aimg->n], diff, m, max_grad_mag = IFT_INFINITY_FLT_NEG;

    for (i=0; i < Abdr->n; i++)
        mag_bdr[i]=sqrtf(Abdr->dx[i]*Abdr->dx[i]+Abdr->dy[i]*Abdr->dy[i]+Abdr->dz[i]*Abdr->dz[i]);

    for (i=0; i < Aimg->n; i++)
        mag_img[i]=sqrtf(Aimg->dx[i]*Aimg->dx[i]+Aimg->dy[i]*Aimg->dy[i]+Aimg->dz[i]*Aimg->dz[i]);


    n.x = n.y = n.z = 0.0;
    cur.x = cur.y = cur.z = 0.0;
    Imax.x = Imax.y = Imax.z = 0.0;

    u = iftGetVoxelCoord(img,p);

    // Computing the direction of the local vector on the filling region's boundary
    for (i=1; i < Abdr->n; i++) {
        v = iftGetAdjacentVoxel(Abdr,u,i);
        if (iftValidVoxel(img,v)){ /* estimate normal on the boundary of the filling region */
            q    = iftGetVoxelIndex(img,v);
            diff = (iftBMapValue(filled, q)-iftBMapValue(filled, p));
            n.x += (diff*Abdr->dx[i]/mag_bdr[i]);
            n.y += (diff*Abdr->dy[i]/mag_bdr[i]);
            n.z += (diff*Abdr->dz[i]/mag_bdr[i]);
        }
    }

    // Computing the maximum gradient inside the patch among all valid pixels
    for (i=0; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A,u,i);
        if (iftValidVoxel(img,v)){ /* estimate normal on the boundary of the filling region */
            q    = iftGetVoxelIndex(img,v);
            if (iftBMapValue(filled, q)){ /* estimate image gradient outside filling region */
                m = 0.0;

                cur = iftComputeImageGradientOnFilledPixel(q, img, filled, Aimg, mag_img);
                m = iftVectorMagnitude(cur);

                if(m > max_grad_mag)
                {
                    max_grad_mag = m;
                    Imax = cur;
                }
            }
        }
    }

    m = iftVectorMagnitude(n);
    if (m > IFT_EPSILON) {
        n.x /= m;    n.y /= m;    n.z /= m;
    }

    /* the image term must be higher as more orthogonal is the image gradient to the normal vector */

    m = iftVectorMagnitude(Imax);
    if (m > IFT_EPSILON) {
        Imax.x /= m;    Imax.y /= m;    Imax.z /= m;
        return(1.0 - fabs(n.x*Imax.x + n.y*Imax.y + n.z*Imax.z));
    } else {
        return(IFT_EPSILON);
    }

}


/* Sets the priority and confidence maps for all the pixels in set B. The heap
	is updated accordingly.
 */

void iftSetPriorityMapAndHeap(iftFHeap *H, iftFImage *priority, iftSet *B, iftImage *img,
                              iftFImage *confidence, iftAdjRel *A, iftBMap *filled)
{
    iftSet  *Baux=B;
    int      p;
    float    C = 0.0, D = 0.0;
    iftAdjRel *Abdr = NULL, *Aimg = NULL;

    if(iftIs3DImage(img))
        Abdr = iftSpheric(3.0);
    else
        Abdr = iftCircular(3.0);

    if(iftIs3DImage(img))
        Aimg = iftSpheric(1.5);
    else
        Aimg = iftCircular(1.5);

    while (Baux != NULL) {
        p = Baux->elem;

        C = iftVoxelConfidence(p,confidence,A,filled);
        D = iftImageTerm(p,img,A,Abdr,Aimg,filled);

        priority->val[p] = C*D;
        confidence->val[p] = C;

        if(H->color[p] == IFT_GRAY)
            iftGoUpFHeap(H, p);
        else
            iftInsertFHeap(H, p);

        Baux = Baux->next;
    }

    iftDestroyAdjRel(&Abdr);
    iftDestroyAdjRel(&Aimg);
}

/* Updates the priority and confidence for patches locally, around the currently
	filled template
*/

void iftUpdatePatchBoundaryPriority(int p,iftFHeap *H, iftFImage *priority,
                                    iftImage *img, iftFImage *confidence, iftAdjRel *A,
                                    iftBMap *filled)
{
    int        i,j,q,r;
    iftVoxel   u,v,w;
    iftAdjRel *C = NULL;
    iftSet    *B=NULL;
    char 		is_boundary;

    if (iftIs3DFImage(confidence))
        C = iftSpheric(1.0);
    else
        C = iftCircular(1.0);

    u = iftFGetVoxelCoord(confidence,p);
    for (i=0; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A, u, i);
        if(iftFValidVoxel(confidence, v)) {
            is_boundary = 0;
            q = iftFGetVoxelIndex(confidence, v);
            if (iftBMapValue(filled,q)){
                for (j=1; j < C->n && !is_boundary; j++) {
                    w = iftGetAdjacentVoxel(C,v,j);
                    if (iftFValidVoxel(confidence,w)){
                        r = iftFGetVoxelIndex(confidence,w);
                        if (!iftBMapValue(filled, r)){
                            is_boundary = 1;
                        }
                    }
                }
            }

            if(is_boundary)
            {
                iftInsertSet(&B,q);
            }
            else
            {
                if(H->color[q] == IFT_GRAY)
                    iftRemoveFHeapElem(H, q);
            }
        }
    }

    iftSetPriorityMapAndHeap(H, priority, B, img, confidence, A, filled);
    iftDestroyAdjRel(&C);
    iftDestroySet(&B);
}


/* Update confidence values inside template  */

void iftUpdateConfidenceAndFilledMaps(int p, iftFImage *confidence, iftAdjRel *A, iftBMap *filled)
{
    int        i,q;
    iftVoxel   u,v;

    u = iftFGetVoxelCoord(confidence,p);
    for (i=0; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A,u,i);
        if (iftFValidVoxel(confidence,v)){
            q = iftFGetVoxelIndex(confidence,v);

            if(!iftBMapValue(filled, q))
            {
                confidence->val[q] = confidence->val[p];
                iftBMapSet1(filled, q); // Marking pixel as filled
            }
        }
    }
}

/** SSE optimization **/

/* Allocation/Deallocation of auxiliary data for fast texture sample comparison */

iftFastTextureComparisonData *iftCreateFastTextureComparisonData(iftTextureSamples *samples)
{
    iftFastTextureComparisonData *comparison_data = NULL;
    iftMatrix *valid_feat = NULL;

    valid_feat = iftValidTextureSamplesToMatrix(samples);

    comparison_data = (iftFastTextureComparisonData*)iftAlloc(1, sizeof(iftFastTextureComparisonData));
    comparison_data->valid_feat = valid_feat;
    comparison_data->features = iftAlignPatchFeatures(valid_feat);

    return comparison_data;
}

void iftDestroyFastTextureComparisonData(iftFastTextureComparisonData **samples)
{
    if(samples != NULL && *samples != NULL)
    {
        _mm_free((*samples)->features);

        iftDestroyMatrix(&(*samples)->valid_feat);

        iftFree(*samples);

        *samples = NULL;
    }
}

/* This function filters out the sampled texture selecting only the valid ones
	to compose a Samples x Features matrix (i.e., the features of a sample are
	set as columns).
 */
iftMatrix *iftValidTextureSamplesToMatrix(iftTextureSamples *samples)
{
    int i, i0, j, nfeats, nvalid = 0;
    iftMatrix *feat = NULL;

    nfeats = samples-> nfeats;

    for(i = 0; i < samples->nsamples; i++)
    {
        if(samples->sample[i].valid)
            nvalid++;
    }

    feat = iftCreateMatrix(nfeats, nvalid);

    i0 = 0;
    for(i = 0; i < samples->nsamples; i++)
    {
        if(samples->sample[i].valid)
        {
            for(j = 0; j < nfeats; j++)
            {
                feat->val[iftGetMatrixIndex(feat, j, i0)] = samples->sample[i].feat[j];
            }

            i0++;
        }
    }

    return feat;
}

/* Aligning feature patches in the following way:
	P0_0, P1_0, P2_0, P3_0, P0_1, P1_1, P2_1, P3_1,...,Pn-3_m, Pn-2_m, Pn-1_m, Pn_m,

	where P1 represents Patch 1 and _0 the 0th feature of the corresponding path.

	This alignment allows the computation of the euclidean distance to 4 patches
	simultaneously using SSE, while avoiding cache misses because if we were to
	store the features in a Patch by Feature matrix (using features as rows) the
	SSE computation would often request extra rows from memory.

	We pad the array with zeros to ensure its boundary is aligned to IFT_MEMORY_ALIGNMENT
*/
float* iftAlignPatchFeatures(iftMatrix *texture)
{
    int i, j, k, i0;
    float *features = NULL;
    const int nsamples = texture->nrows;
    const int nfeats = texture->ncols;
    const int nsamples_padded = nsamples + 4 - nsamples % 4;
    const int n = nsamples_padded*nfeats;

    features = iftAllocAlignedFloatArray(n, IFT_MEMORY_ALIGNMENT);

    i0 = 0;
    for(k = 0; k < nsamples_padded/4; k++)
    {
        for(i = 0; i < nfeats; i++)
        {
            for(j = 0; j < 4; j++)
            {
                int sample = k*4+j;

                // We only copy valid samples, i.e., the extra padding
                // is set to zero
                if(sample < nsamples)
                {
                    int matrix_id = iftGetMatrixIndex(texture, i, sample);

                    features[i0] = texture->val[matrix_id];
                }

                i0++;
            }
        }
    }

    return features;
}

int iftFindBestTextureMatchIndex(iftTexturePatch * __restrict__ patch, float * __restrict__ _features,
                                 int nsamples, int nfeats)
{
    int i, j, best_i = 0;
    // total number of nodes + extra padding
    const int n = nsamples + 4 - nsamples % 4;
    const int nloop = n / 4;

    float * __restrict__ features = (float * __restrict__)_features;
    float * __restrict__ feat = patch->feat;
    float * __restrict__ weight = patch->weight;

    // Variables for storing the feature vector of a given pixel p
    // using format F0,F0,F0,F0,F1,F1,F1,F1,...,Fn,Fn,Fn,Fn
    // for both the features and weights. Note that __m128 stores
    // 4 floats at once
    __m128 * __restrict__ f = (__m128*)iftAlloc(nfeats, sizeof(__m128));
    __m128 * __restrict__ w = (__m128*)iftAlloc(nfeats, sizeof(__m128));

    // Variables for extracting the patch's index of minimum distance
    int min_index_arr[4] __attribute__((aligned(IFT_MEMORY_ALIGNMENT)));
    float min_dist_arr[4] __attribute__((aligned(IFT_MEMORY_ALIGNMENT)));

    for(i = 0; i < nfeats; i++)
    {
        // Setting features with format F0,F0,F0,F0,F1,F1,F1,F1,...,Fn,Fn,Fn,Fn
        f[i] = _mm_set_ps(feat[i], feat[i], feat[i], feat[i]);
        // Setting weight with the same format as above W0,W0,W0,W0,W1,W1,W1,W1,...,Wn,Wn,Wn,Wn
        w[i] = _mm_set_ps(weight[i], weight[i], weight[i], weight[i]);
    }


    __m128 min_dist = _mm_set_ps(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX);
    __m128i min_dist_i = _mm_set_epi32(0,0,0,0);

    __m128 *p2 = (__m128*)features;

    // Computing the euclidean distance to all patches, computing for 4 patches
    // simultaneously using SSE
    for(j = 0; j < nloop; j++)
    {
        const int index = j * 4;
        __m128 dist = _mm_setzero_ps();
        // Numbers are set in reverse indexing order using _mm_set_epi32
        __m128i dist_i = _mm_set_epi32(index+3, index+2, index+1, index);

        // Euclidean distance
        for(i = 0; i < nfeats; i++)
        {
            __m128 sub = _mm_sub_ps(*p2, f[i]);

            dist = _mm_add_ps(dist, _mm_mul_ps(_mm_mul_ps(sub, sub), w[i]));
            p2++;
        }

        // Computing the minimum distance and the index of the corresponding path
        __m128 mask = _mm_cmplt_ps(dist, min_dist);
        __m128i mask_i = _mm_castps_si128(mask);

        min_dist = _mm_or_ps(_mm_and_ps(mask, dist), _mm_andnot_ps(mask, min_dist));

        min_dist_i = _mm_or_si128(_mm_and_si128(mask_i, dist_i), _mm_andnot_si128(mask_i, min_dist_i));
    }

    // Obtaining the index of the patch with minimum distance
    // among the four minimum distances encountered.
    // IMPORTANT: the distance to the padding patches (i.e., the last ones
    // with feature values set to 0.0) will be "maximum" and is only going
    // to be "considered" if the minimum distance to all patches is
    // maximum. In other words, we compute the following distance:
    // d = (x-0.0)*(x-0.0) = x^2 to the padding patches, which is
    // certainly not minimum unless all patches have feature values 0.0
    // as well.
    _mm_store_si128((__m128i *)min_index_arr, min_dist_i);
    _mm_store_ps(min_dist_arr, min_dist);

    best_i = min_index_arr[iftFArgmin(min_dist_arr, 4)];

    iftFree(f);
    iftFree(w);

    return best_i;
}



/* ------------------- Public method --------------------------------------*/

/* This function samples texture from a band around the mask's boundary, ranging
	from 2*adj_radius to 4*adj_radius outside it, and uses iftInpainting for
	foreground completion.
 */

iftImage *iftInpaintingFromSamplesAroundMaskBoundaryOutside(iftImage *img, iftImage *mask, float adj_radius)
{
    iftAdjRel  *A      = NULL;
    iftImage *inpainted = NULL;
    iftTextureSamples *texture = NULL;

    iftVerifyImageDomains(img,mask,"iftInpainting");

    if (adj_radius <= 1.0)
        iftError("Adjacency radius must be greater than 1", "iftInpainting");

    if (iftIs3DImage(img))
        A = iftSpheric(adj_radius);
    else
        A = iftCircular(adj_radius);

    texture    = iftTextureSamplesAroundMaskBoundaryOutside(img, mask, adj_radius);

    inpainted = iftInpainting(img, mask, texture, A);

    iftDestroyTextureSamples(&texture);
    iftDestroyAdjRel(&A);

    return inpainted;
}

iftImage *iftInpainting(iftImage *img, iftImage *mask, iftTextureSamples *texture,
                        iftAdjRel *A)
{
    int         p;
    iftImage   *nimg   = iftCopyImage(img);
    iftFImage  *confidence, *priority;

    iftSet     *B      = NULL;
    iftFastTextureComparisonData *fast_texture_comp = NULL;
    iftFHeap 	*H = NULL;
    iftBMap 		*filled = NULL;
    iftTexturePatch       *patch, *best_patch;

    iftVerifyImageDomains(img,mask,"iftInpainting");

    if (iftIsColorImage(img)){
        patch = iftCreateTexturePatch(3*A->n);
        best_patch = iftCreateTexturePatch(3*A->n);
    } else {
        patch = iftCreateTexturePatch(A->n);
        best_patch = iftCreateTexturePatch(A->n);
    }

    filled	  = iftCreateBMap(nimg->n);
    confidence = iftConfidenceMap(mask, filled);
    priority   = iftCreateFImage(mask->xsize,mask->ysize,mask->zsize);
    B          = iftBoundaryOfFillingRegion(confidence, filled);
    fast_texture_comp = iftCreateFastTextureComparisonData(texture);

    H 			  = iftCreateFHeap(priority->n, priority->val);
    H->removal_policy = MAXVALUE;

    iftSetPriorityMapAndHeap(H, priority, B, nimg, confidence, A, filled);

    while (!iftEmptyFHeap(H)) {
        p = iftRemoveFHeap(H);
        iftVoxelTexture(p,nimg,A,patch, filled);
        iftFindBestTextureMatch(patch, best_patch, fast_texture_comp);
        iftMapTexture(best_patch,p,A,nimg, filled);
        iftUpdateConfidenceAndFilledMaps(p, confidence, A, filled);
        iftUpdatePatchBoundaryPriority(p, H, priority, nimg, confidence, A, filled);
    }

    iftDestroyTexturePatch(&patch);
    iftDestroyTexturePatch(&best_patch);
    iftDestroyFImage(&priority);
    iftDestroyFImage(&confidence);
    iftDestroyFHeap(&H);
    iftDestroyBMap(&filled);
    iftDestroySet(&B);
    iftDestroyFastTextureComparisonData(&fast_texture_comp);

    return(nimg);
}

iftImage *iftInpaintBoundary(iftImage *img, iftImage *mask)
{
    iftImage   *nimg   = NULL;
    iftAdjRel  *A      = NULL;
    iftSet     *B      = NULL;
    int         p;
    iftTextureSamples *texture;
    iftFastTextureComparisonData *fast_texture_comp = NULL;
    iftBMap 	  *filled = NULL;
    float      adj_radius;
    iftTexturePatch       *patch, *best_patch;

    iftVerifyImageDomains(img,mask,"iftInpaintBoundary");

    adj_radius = 1.0;

    nimg = iftMask(img, mask);

    if (iftIs3DImage(img)){
        A = iftSpheric(adj_radius);
    }else{
        A = iftCircular(adj_radius);
    }

    B = iftObjectBorderSet(mask,A);

    if (iftIsColorImage(img)){
        patch = iftCreateTexturePatch(3*A->n);
        best_patch = iftCreateTexturePatch(3*A->n);
    } else {
        patch = iftCreateTexturePatch(A->n);
        best_patch = iftCreateTexturePatch(A->n);
    }

    texture	= iftTextureSamplesAroundMaskBoundaryInside(nimg, mask, adj_radius);
    fast_texture_comp = iftCreateFastTextureComparisonData(texture);
    filled 	= iftCreateBMap(nimg->n);

    while (B != NULL) {
        p = iftRemoveSet(&B);
        iftBoundaryVoxelTexture(p,nimg,mask,A,patch);
        iftFindBestTextureMatch(patch,best_patch,fast_texture_comp);
        iftMapTexture(best_patch,p,A,nimg,filled);
    }

    iftDestroyTexturePatch(&patch);
    iftDestroyTexturePatch(&best_patch);
    iftDestroyAdjRel(&A);
    iftDestroyBMap(&filled);
    iftDestroyTextureSamples(&texture);
    iftDestroyFastTextureComparisonData(&fast_texture_comp);

    return(nimg);
}































