#include "iftFImageForest.h"

#include "ift/core/io/Stream.h"


iftFImageForest  *iftCreateFImageForest(iftImage *img)
{
  iftFImageForest *fst=(iftFImageForest *)iftAlloc(1,sizeof(iftFImageForest));

  if (fst != NULL) {
    fst->pathval = iftCreateFImage(img->xsize,img->ysize,img->zsize);
    fst->root    = iftCreateImage(img->xsize,img->ysize,img->zsize);
    fst->label   = iftCreateImage(img->xsize,img->ysize,img->zsize);
    fst->pred    = iftCreateImage(img->xsize,img->ysize,img->zsize);
  }else{
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateFImageForest");
  }

  return(fst);
}

void     iftDestroyFImageForest(iftFImageForest **fst)
{
  iftFImageForest *aux;

  aux = *fst;
  if(aux != NULL){
    if (aux->pathval  != NULL)  iftDestroyFImage(&(aux->pathval));
    if (aux->root     != NULL)  iftDestroyImage(&(aux->root));
    if (aux->label    != NULL)  iftDestroyImage(&(aux->label));
    if (aux->pred     != NULL)  iftDestroyImage(&(aux->pred));
    iftFree(aux);
    *fst = NULL;
  }
}
