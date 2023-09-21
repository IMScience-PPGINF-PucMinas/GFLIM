//
// Created by ilan on 14/02/2022.
//

#ifndef IFT_IMAGESEQUENCE_H
#define IFT_IMAGESEQUENCE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftMImage.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "ift/core/dtypes/GQueue.h"
#include "ift/core/tools/String.h"

typedef struct ift_image_sequence {
    float          *****val;
    int            xsize, ysize, zsize, tsize;
    float          dx, dy, dz, dt;
    int            m;
    long           n;
    iftColorSpace  cspace;

    bool           *allocated;
    iftFileSet     *files;
} iftImageSequence;

#define IMGSEQ_FORLOOP(img)       \
    for(int t=0;t<img->tsize;t++) \
    for(int z=0;z<img->zsize;z++) \
    for(int y=0;y<img->ysize;y++) \
    for(int x=0;x<img->xsize;x++)

//#define IMGSEQ_ITE_VAL(img) (img->val[t][z][y][x])
#define IMGSEQ_ITE_VAL(img, b) (img->val[t][z][y][x][b])
#define IMGSEQ_VOXEL_VAL(img, v, b) (img->val[v.t][v.z][v.y][v.x][b])

#define IMGSEQ_FORLOOP_PARALLEL(img)      \
    _Pragma("omp parallel for collapse(4)") \
    for(int t=0;t<img->tsize;t++)         \
    for(int z=0;z<img->zsize;z++)         \
    for(int y=0;y<img->ysize;y++)         \
    for(int x=0;x<img->xsize;x++)

iftImageSequence  *iftCreateImageSequence(int xsize,int ysize,int zsize, int tsize, int nbands);
void  iftImageSequenceAllocVolume(iftImageSequence *img, int t);
iftImageSequence *iftLoadImageSequence(iftFileSet *files);
void  iftDestroyImageSequence(iftImageSequence **img);
iftImageSequence *iftCopyImageSequence(const iftImageSequence *img);
iftImageSequence *iftCreateImageSequenceFromImageSequence(const iftImageSequence *img);

void iftSetImageSequence(iftImageSequence *img, float value);

iftImageSequence *ift4DFindAndLabelObjectBorders(const iftImageSequence *label_img,
                                                 const iftAdjRel *Ain);



double ift4DDiceSimilarity(const iftImageSequence *label_source,
                                           const iftImageSequence *label_target);

iftDblArray *ift4DDiceSimilarityMultiLabel(const iftImageSequence *label_source,
                                           const iftImageSequence *label_target);

iftImageSequence *ift4DEuclDistTrans(const iftImageSequence *label_img,
                                     const iftAdjRel *Ain, iftSide side,
                                     iftImageSequence **root_out,
                                     iftImageSequence **edt_label_out,
                                     iftImageSequence **pred_out);

char ift4DValidVoxel(const iftImageSequence *img, iftVoxel v);
long ift4DGetVoxelIndex(const iftImageSequence *img, iftVoxel u);
iftVoxel ift4DGetVoxelCoord(const iftImageSequence *img, long p);
float * ift4DGetVal(const iftImageSequence *img, long p);

iftGQueue *iftCreateGQueueWithImageSequence(int nbuckets, int nelems, iftImageSequence *img);

iftMImage *iftImageSequenceToMImage(const iftImageSequence *img, int time);
iftImageSequence *iftMImageToImageSequence(const iftMImage *mimg);

void iftWriteImageSequence(const iftImageSequence *img, const char *filename);
iftImageSequence * iftReadImageSequence(const char *filename);
//TODO rename
iftImageSequence * iftReadImageAsImageSequence(const char *filename, const iftImageSequence *imgseq_template, int current_time);

iftImage *iftExtractVolumeFromImageSequence(const iftImageSequence *img, int time);
iftImageSequence *iftVolumeToImageSequence(const iftImage *img);
void iftAppendVolumeToImageSequence(iftImageSequence *img, iftImage *vol);
void iftSetVolumeInImageSequence(iftImageSequence *img, iftImage *vol, int time);
void iftSetMImageInImageSequence(iftImageSequence *img, iftMImage *mimg, int time);
//void iftAppendMImageToImageSequence(iftImageSequence *img, iftMImage *mimg);
bool iftIsImageCompatibleToImageSequence(const iftImageSequence *img, iftImage *vol);
bool iftIs3DImageSequence(const iftImageSequence *img);
bool iftIs4DImageSequence(const iftImageSequence *img);
bool ift4DIsColorImage(const iftImageSequence *img);

iftMImage *iftGetXYSliceFromTime(const iftImageSequence *img, int zcoord, int time);
iftMImage *iftGetZXSliceFromTime(const iftImageSequence *img, int ycoord, int time);
iftMImage *iftGetYZSliceFromTime(const iftImageSequence *img, int xcoord, int time);

iftImageSequence *ift4DNormalize(iftImageSequence *img, float min, float max);
void ift4DMinMaxValues(const iftImageSequence *img, int band, float *min, float *max);
float ift4DMaximumValue(const iftImageSequence *img, int band);
float ift4DMinimumValue(const iftImageSequence *img, int band);
uchar ift4DImageDepth(const iftImageSequence *img);

iftImageSequence *ift4DGradient(const iftImageSequence *img, iftAdjRel *A, int Imax);
iftImageSequence *ift4DBorderImage(const iftImageSequence *label, bool get_margins);
long iftCountObjectSpelsFromTime(const iftImageSequence *img, int obj_label, int time);
long iftCountObjectSpelsImageSequence(const iftImageSequence *img, int obj_label);

iftLabeledSet * ift4DReadSeeds(const iftImageSequence *img, const char *_filename, ...);
void ift4DWriteSeeds(iftLabeledSet* seed, const iftImageSequence* img, const char *_filename, ...);


iftImageSequence *ift4DThreshold(const iftImageSequence *img, float lowest, float highest, float value);
iftImageSequence *ift4DThresholdInTime(const iftImageSequence *img, float lowest, float highest, float value, int time);

iftImageSequence *ift4DBinarize(const iftImageSequence *img);

typedef struct ift4d_imageforest {
    iftImageSequence  *pathval;
    iftImageSequence  *label;
    iftImageSequence  *root;
    iftImageSequence  *pred;
    iftGQueue *Q;
    iftImageSequence  *img;
    iftAdjRel *A;
} ift4DImageForest;

ift4DImageForest  *iftCreate4DImageForest(iftImageSequence *img, iftAdjRel *A);
void             iftReset4DImageForest(ift4DImageForest *fst);
void             iftDestroy4DImageForest(ift4DImageForest **fst);
void             ift4DDiffWatershed(ift4DImageForest *fst,
                                    iftLabeledSet *seed,
                                    iftSet * removal_markers);
iftSet *ift4DTreeRemoval(ift4DImageForest *fst, iftSet *trees_for_removal);

#define CMAX     65535

typedef struct ift4d_dyntrees {
    iftImageSequence *img;          /* original image in some color space */
    iftAdjRel        *A;            /* adjacency relation */
    iftImageSequence *label;        /* label map */
    iftImageSequence *root;         /* root map */
    iftImageSequence *cumfeat;      /* cumulative feature vector map */
    int              *treesize;     /* tree size map */
    iftImageSequence *cost;         /* path cost map */
    iftImageSequence *pred;         /* predecessor map */
    iftGQueue        *Q;            /* priority queue with integer costs */
    float             maxfeatdist;  /* maximum feature distance */
} ift4DDynTrees;

ift4DDynTrees *iftCreate4DDynTrees(iftImageSequence *img, iftAdjRel *A);
void                      iftDestroy4DDynTrees(ift4DDynTrees **dt);
void                      iftExecDiff4DDynTrees(ift4DDynTrees *dt, iftLabeledSet **S, iftSet **M);
iftSet                   *iftDiff4DDynTreeRemoval(ift4DDynTrees *dt, iftSet **M);
void                      iftDiff4DDynSubtreeUpdate(ift4DDynTrees *dt, long q);
#ifdef __cplusplus
}
#endif

#endif //IFT_IMAGESEQUENCE_H
