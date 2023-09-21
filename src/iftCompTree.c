#include "iftCompTree.h"

#include "ift/core/dtypes/GQueue.h"
#include "ift/core/io/Stream.h"



void iftCumSize(iftCompTree *ctree, int i)
{
    int s,j;

    for (j=0; j < ctree->node[i].numsons; j++)
    {
        s = ctree->node[i].son[j];
        iftCumSize(ctree,s);
        ctree->node[i].size = ctree->node[i].size + ctree->node[s].size;
    }
    return;
}

int iftAreaLevel(iftCompTree *ctree, int *level, int i, int thres)
{
    if ((ctree->node[i].size > thres)||(i==ctree->root))
        return(ctree->node[i].level);
    else
        return(level[i]=iftAreaLevel(ctree,level,ctree->node[i].dad,thres));
}

int iftVolumeLevel(iftCompTree *ctree, int *level, int i, int thres, int cumvol)
{
    int dad,vol=cumvol;

    dad = ctree->node[i].dad;
    if (dad != IFT_NIL)
        vol = cumvol+
              abs(ctree->node[i].level-ctree->node[dad].level)*ctree->node[i].size;

    if ((vol > thres)||(i==ctree->root))
        return(ctree->node[i].level);
    else
        return(level[i]=iftVolumeLevel(ctree,level,dad,thres,vol));
}


int iftRepresentative(iftImage *cmap, int p)
{
    if (cmap->val[p]==p)
        return(p);
    else
        return(cmap->val[p]=iftRepresentative(cmap,cmap->val[p]));
}

int iftAncestor(iftImage *dad, iftImage *cmap, int rq)
{
    int r,ro;

    ro = r  = dad->val[rq];
    while (r != IFT_NIL)
    {
        ro = r = iftRepresentative(cmap,r);
        r  = dad->val[r];
    }
    return(ro);
}

iftCompTree *iftCreateMaxTree(const iftImage *img)
{
    iftCompTree *ctree=(iftCompTree *)iftAlloc(1,sizeof(iftCompTree));
    iftImage    *dad,*cmap,*tmp;
    int i,r,p,q,rp,rq,Imax=iftMaximumValue(img);
    iftGQueue *Q;
    iftVoxel u,v;
    int *nsons=NULL;
    int *size=NULL;
    iftAdjRel *A;

    if (img->zsize==1)
        A = iftCircular(1.0);
    else
        A = iftSpheric(1.0);


    ctree->cmap = iftCreateImage(img->xsize,img->ysize,img->zsize);
    if (iftIsColorImage(img))
      iftCopyCbCr(img,ctree->cmap);
    iftCopyVoxelSize(img,ctree->cmap);
    cmap        = ctree->cmap;
    ctree->root = IFT_NIL; /* Tree is empty */
    dad         = iftCreateImage(img->xsize,img->ysize,img->zsize);
    size        = iftAllocIntArray(img->n);
    Q           = iftCreateGQueue(Imax+1,img->n,img->val);
    iftSetRemovalPolicy(Q,MAXVALUE);
    iftSetTieBreak(Q,LIFOBREAK);

    for (p=0; p < cmap->n; p++)
    {
        dad->val[p]  = IFT_NIL;
        cmap->val[p] =p;
        size[p]=1;
        iftInsertGQueue(&Q,p);
    }

    while(!iftEmptyGQueue(Q))
    {
        p  = iftRemoveGQueue(Q);
        rp = iftRepresentative(cmap,p);
        u  = iftGetVoxelCoord(cmap,p);
        for (i=1; i < A->n; i++)
        {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(cmap,v))
            {
                q = iftGetVoxelIndex(cmap,v);
                if (img->val[q]==img->val[p])   /* propagate component */
                {
                    if (Q->L.elem[q].color == IFT_GRAY)
                    {
                        iftRemoveGQueueElem(Q,q);
                        cmap->val[q]=rp;
                        size[rp]=size[rp]+1;
                        iftInsertGQueue(&Q,q);
                    }
                }
                else
                {
                    if (img->val[p] < img->val[q]) /* find current dad of rq */
                    {
                        rq = iftRepresentative(cmap,q);
                        r  = iftAncestor(dad,cmap,rq);
                        if (r == IFT_NIL)   /* rp is dad of the rq */
                        {
                            dad->val[rq]=rp;
                        }
                        else
                        {
                            if (img->val[r]==img->val[rp])  /* merge components */
                            {
                                if (r != rp)
                                {
                                    if (size[rp] <= size[r])
                                    {
                                        cmap->val[rp] = r;
                                        size[r] = size[r] + size[rp];
                                        rp = r;
                                    }
                                    else
                                    {
                                        cmap->val[r] = rp;
                                        size[rp] = size[rp] + size[r];
                                    }
                                }
                            }
                            else     /* img->val[r] > img->val[rp] */
                            {
                                dad->val[r] = rp; /* rp is dad of r */
                            }
                        }
                    }
                }
            }
        }
    }
    iftFree(size);
    iftDestroyGQueue(&Q);
    iftDestroyAdjRel(&A);

    /* Compress cmap map and count number of nodes */

    ctree->numnodes = 0;
    for (p=0; p < cmap->n; p++)
    {
        if (dad->val[cmap->val[p]] != IFT_NIL)
            r = cmap->val[p];
        cmap->val[p] = iftRepresentative(cmap,p);

        if (cmap->val[p]==p)
            ctree->numnodes++;
    }

    /* Create and initialize nodes of the MaxTree. */

    ctree->node = (iftCompTreeNode *)iftAlloc(ctree->numnodes,sizeof(iftCompTreeNode));
    tmp         = iftCreateImage(img->xsize,img->ysize,img->zsize);
    for (p=0; p < cmap->n; p++)
    {
        tmp->val[p]= IFT_NIL;
    }

    i = 0;
    for (p=0; p < cmap->n; p++)
    {
        if (cmap->val[p]==p)
        {
            ctree->node[i].level = img->val[p];
            ctree->node[i].comp  = p;
            tmp->val[p]          = i;
            ctree->node[i].dad   = IFT_NIL;
            ctree->node[i].son   = NULL;
            ctree->node[i].numsons = 0;
            ctree->node[i].size  = 0;
            i++;
        }
    }

    /* Make the component map to point back to the maxtree. */

    for (p=0; p < tmp->n; p++)
    {
        if (tmp->val[p] == IFT_NIL)
            tmp->val[p] = tmp->val[cmap->val[p]];
    }

    for (p=0; p < cmap->n; p++)
    {
        cmap->val[p] = tmp->val[p];
    }
    iftDestroyImage(&tmp);

    /* Copy dad information to the maxtree and find its root */

    for (i=0; i < ctree->numnodes; i++)
    {
        if (dad->val[ctree->node[i].comp] != IFT_NIL)
            ctree->node[i].dad = cmap->val[dad->val[ctree->node[i].comp]];
        else
        {
            ctree->node[i].dad = IFT_NIL;
            ctree->root = i;
        }
    }
    iftDestroyImage(&dad);

    /* Copy son information to the maxtree */

    nsons = iftAllocIntArray(ctree->numnodes);
    for (i=0; i < ctree->numnodes; i++)
    {
        p = ctree->node[i].dad;
        if (p != IFT_NIL)
        {
            nsons[p]++;
        }
    }
    for (i=0; i < ctree->numnodes; i++)
    {
        if (nsons[i] != 0)
        {
            ctree->node[i].son = iftAllocIntArray(nsons[i]);
        }
    }
    iftFree(nsons);

    for (i=0; i < ctree->numnodes; i++)
    {
        p = ctree->node[i].dad;
        if (p != IFT_NIL)
        {
            ctree->node[p].son[ctree->node[p].numsons]=i;
            ctree->node[p].numsons++;
        }
    }

    /* Compute size of each node */

    for (p=0; p < cmap->n; p++)
        ctree->node[cmap->val[p]].size++;


    return(ctree);
}


iftCompTree *iftCreateMinTree(const iftImage *img)
{
    iftCompTree *ctree=(iftCompTree *)iftAlloc(1,sizeof(iftCompTree));
    iftImage    *dad,*cmap,*tmp;
    int i,r,p,q,rp,rq,Imax=iftMaximumValue(img);
    iftGQueue *Q;
    iftVoxel u,v;
    int *nsons=NULL;
    int *size=NULL;
    iftAdjRel *A;

    if (img->zsize==1)
        A = iftCircular(1.0);
    else
        A = iftSpheric(1.0);

    ctree->cmap = iftCreateImage(img->xsize,img->ysize,img->zsize);
    if (iftIsColorImage(img))
      iftCopyCbCr(img,ctree->cmap);
    iftCopyVoxelSize(img,ctree->cmap);
    cmap        = ctree->cmap;
    ctree->root = IFT_NIL; /* Tree is empty */
    dad         = iftCreateImage(img->xsize,img->ysize,img->zsize);
    size        = iftAllocIntArray(img->n);
    Q           = iftCreateGQueue(Imax+1,img->n,img->val);
    iftSetTieBreak(Q,LIFOBREAK);

    for (p=0; p < cmap->n; p++)
    {
        dad->val[p]  = IFT_NIL;
        cmap->val[p] =p;
        size[p]=1;
        iftInsertGQueue(&Q,p);
    }

    while(!iftEmptyGQueue(Q))
    {
        p  = iftRemoveGQueue(Q);
        rp = iftRepresentative(cmap,p);
        u  = iftGetVoxelCoord(cmap,p);
        for (i=1; i < A->n; i++)
        {
            v = iftGetAdjacentVoxel(A,u,i);
            if (iftValidVoxel(cmap,v))
            {
                q = iftGetVoxelIndex(cmap,v);
                if (img->val[q]==img->val[p])   /* propagate component */
                {
                    if (Q->L.elem[q].color == IFT_GRAY)
                    {
                        iftRemoveGQueueElem(Q,q);
                        cmap->val[q]=rp;
                        size[rp]=size[rp]+1;
                        iftInsertGQueue(&Q,q);
                    }
                }
                else
                {

                    if (img->val[p] > img->val[q]) /* find current dad of rq */
                    {
                        rq = iftRepresentative(cmap,q);
                        r  = iftAncestor(dad,cmap,rq);
                        if (r == IFT_NIL)   /* rp is dad of the rq */
                        {
                            dad->val[rq]=rp;
                        }
                        else
                        {
                            if (img->val[r]==img->val[rp])  /* merge components */
                            {
                                if (r != rp)
                                {
                                    if (size[rp] <= size[r])
                                    {
                                        cmap->val[rp] = r;
                                        size[r] = size[r] + size[rp];
                                        rp = r;
                                    }
                                    else
                                    {
                                        cmap->val[r] = rp;
                                        size[rp] = size[rp] + size[r];
                                    }
                                }
                            }
                            else     /* img->val[r] < img->val[rp] */
                            {
                                dad->val[r] = rp; /* rp is dad of r */
                            }
                        }
                    }
                }
            }
        }
    }
    iftFree(size);
    iftDestroyGQueue(&Q);
    iftDestroyAdjRel(&A);

    /* Compress cmap map and count number of nodes */

    ctree->numnodes = 0;
    for (p=0; p < cmap->n; p++)
    {
        if (dad->val[cmap->val[p]] != IFT_NIL)
            r = cmap->val[p];
        cmap->val[p] = iftRepresentative(cmap,p);

        if (cmap->val[p]==p)
            ctree->numnodes++;
    }

    /* Create and initialize nodes of the MinTree. */

    ctree->node = (iftCompTreeNode *)iftAlloc(ctree->numnodes,sizeof(iftCompTreeNode));
    tmp         = iftCreateImage(img->xsize,img->ysize,img->zsize);
    for (p=0; p < cmap->n; p++)
    {
        tmp->val[p]= IFT_NIL;
    }

    i = 0;
    for (p=0; p < cmap->n; p++)
    {
        if (cmap->val[p]==p)
        {
            ctree->node[i].level = img->val[p];
            ctree->node[i].comp  = p;
            tmp->val[p]          = i;
            ctree->node[i].dad   = IFT_NIL;
            ctree->node[i].son  = NULL;
            ctree->node[i].numsons = 0;
            ctree->node[i].size  = 0;
            i++;
        }
    }


    /* Make the component map to point back to the mintree. */

    for (p=0; p < tmp->n; p++)
    {
        if (tmp->val[p] == IFT_NIL)
            tmp->val[p] = tmp->val[cmap->val[p]];
    }

    for (p=0; p < cmap->n; p++)
    {
        cmap->val[p] = tmp->val[p];
    }
    iftDestroyImage(&tmp);

    /* Copy dad information to the mintree and find its root */

    for (i=0; i < ctree->numnodes; i++)
    {
        if (dad->val[ctree->node[i].comp] != IFT_NIL)
            ctree->node[i].dad = cmap->val[dad->val[ctree->node[i].comp]];
        else
        {
            ctree->node[i].dad = IFT_NIL;
            ctree->root = i;
        }
    }
    iftDestroyImage(&dad);

    /* Copy son information to the mintree */

    nsons = iftAllocIntArray(ctree->numnodes);
    for (i=0; i < ctree->numnodes; i++)
    {
        p = ctree->node[i].dad;
        if (p != IFT_NIL)
        {
            nsons[p]++;
        }
    }
    for (i=0; i < ctree->numnodes; i++)
    {
        if (nsons[i] != 0)
        {
            ctree->node[i].son = iftAllocIntArray(nsons[i]);
        }
    }
    iftFree(nsons);

    for (i=0; i < ctree->numnodes; i++)
    {
        p = ctree->node[i].dad;
        if (p != IFT_NIL)
        {
            ctree->node[p].son[ctree->node[p].numsons]=i;
            ctree->node[p].numsons++;
        }
    }

    /* Compute size of each node */

    for (p=0; p < cmap->n; p++)
        ctree->node[cmap->val[p]].size++;

    return(ctree);
}


void iftDestroyCompTree(iftCompTree **ctree)
{
    int i;
    iftCompTree *aux=*ctree;

    if (aux != NULL)
    {
      if (aux->cmap != NULL) iftDestroyImage(&aux->cmap);
      for (i=0; i < aux->numnodes; i++)
        {
	  if (aux->node[i].son!=NULL)	    
	    iftFree(aux->node[i].son);
        }
      if (aux->node != NULL) iftFree(aux->node);
      iftFree(aux);
      *ctree = NULL;
    }
}
