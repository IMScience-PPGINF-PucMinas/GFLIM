#ifndef IFT_FIMAGEFOREST_H_
#define IFT_FIMAGEFOREST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftFImage.h"

typedef struct ift_fimageforest {
  iftFImage  *pathval;
  iftImage  *root;
  iftImage  *label;
  iftImage  *pred;
} iftFImageForest;

iftFImageForest  *iftCreateFImageForest(iftImage *img);
void            iftDestroyFImageForest(iftFImageForest **img);

#ifdef __cplusplus
}
#endif

#endif


