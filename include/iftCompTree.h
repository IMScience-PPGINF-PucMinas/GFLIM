#ifndef IFT_COMPTREE_H_
#define IFT_COMPTREE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftImage.h"
#include "iftAdjacency.h"

typedef struct ift_comptreenode {
  int  level;   /* gray level */
  int  comp;    /* representative pixel of this node */
  int  dad;     /* dad node in the maxtree */
  int *son;     /* son nodes in the maxtree */
  int  numsons; /* number of sons in the maxtree */
  int  size;    /* number of pixels of the node */
} iftCompTreeNode;

typedef struct ift_comptree {
  iftCompTreeNode *node;     /* nodes of the mtree */
  iftImage        *cmap;     /* component map */
  int              root;     /* root node of the mtree */
  int              numnodes; /* number of nodes of the maxtree */
} iftCompTree;

int iftRepresentative(iftImage *cmap, int p);
int iftAncestor(iftImage *dad, iftImage *cmap, int rq);

//! swig(newobject)
iftCompTree *iftCreateMaxTree(const iftImage *img);

//! swig(newobject)
iftCompTree *iftCreateMinTree(const iftImage *img);
void         iftDestroyCompTree(iftCompTree **ctree);
void         iftCumSize(iftCompTree *ctree, int i);
int          iftAreaLevel(iftCompTree *ctree, int *level, int i, int thres);
int          iftVolumeLevel(iftCompTree *ctree, int *level, int i, int thres, int cumvol);

#ifdef __cplusplus
}
#endif

#endif
