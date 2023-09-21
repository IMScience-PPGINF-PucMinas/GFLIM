//
// Created by ilan on 14/02/2022.
//

#include "iftImageSequence.h"
#include "ift/core/dtypes/IntArray.h"
#include "iftSimilarity.h"
#include "iftRepresentation.h"
#include "../include/iftImageSequence.h"

//TODO remove this workaround
iftVoxel U_NIL = {IFT_NIL,IFT_NIL,IFT_NIL,IFT_NIL};
long getPixelCoordinate(iftImageSequence *img, int p) {
    float *val = ift4DGetVal(img, p);
    int x = val[0];
    int y = val[1];
    int z = val[2];
    int t = val[3];
    iftVoxel u = {x,y,z,t};
    long q = ift4DGetVoxelIndex(img, u);
    return q;
}

void setPixelCoordinate(iftImageSequence *img, int p, iftVoxel u) {
    float *val = ift4DGetVal(img, p);
    val[0] = u.x;
    val[1] = u.y;
    val[2] = u.z;
    val[3] = u.t;
}

void setPixelCoordinateFromArray(iftImageSequence *img, int p, float *val) {
    iftVoxel u = {val[0], val[1], val[2], val[3]};
    setPixelCoordinate(img, p, u);
}

bool compPixelCoordinates(iftImageSequence *img, int p, int q) {
    float *valp = ift4DGetVal(img, p);
    float *valq = ift4DGetVal(img, q);
    return valp[0] == valq[0] && valp[1] == valq[1] && valp[2] == valq[2] && valp[3] == valq[3];
}
/////////////////////////////

iftImageSequence *
iftCreateImageSequence(int xsize, int ysize, int zsize, int tsize, int nbands) {
    iftImageSequence *img = (iftImageSequence *) iftAlloc(1, sizeof(iftImageSequence));
    if (img == NULL){
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateImageSequence");
    }

    img->m       = nbands;
    img->xsize   = xsize;
    img->ysize   = ysize;
    img->zsize   = zsize;
    img->tsize   = tsize;
    img->dx      = 1.0;
    img->dy      = 1.0;
    img->dz      = 1.0;
    img->dt      = 1.0;
    img->n       = xsize * ysize * zsize * tsize;
    img->cspace  = UNDEFINED_CSPACE;
    img->files   = NULL;

    img->allocated = iftAllocBoolArray(tsize);
    img->val = (float*****) iftAlloc(tsize, sizeof(float****));
    if (img->val == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateImageSequence");
    for (int t = 0; t < tsize; t++) {
        //TODO dont
        iftImageSequenceAllocVolume(img, t);
    }

    return img;
}

void  iftImageSequenceAllocVolume(iftImageSequence *img, int t) {
    if (!img || t < 0 || t >= img->tsize) {
        iftError("Invalid parameters.", "iftImageSequenceAllocVolume");
    }

    img->allocated[t] = true;
    img->val[t] = (float****) iftAlloc(img->zsize, sizeof(float***));
    if (img->val[t] == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftImageSequenceAllocVolume");
    for (int z = 0; z < img->zsize; z++) {
        img->val[t][z] = (float***) iftAlloc(img->ysize, sizeof(float**));
        if (img->val[t][z] == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftImageSequenceAllocVolume");
        for (int y = 0; y < img->ysize; y++) {
            img->val[t][z][y] = (float**) iftAlloc(img->xsize, sizeof(float*));
            if (img->val[t][z][y] == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftImageSequenceAllocVolume");
            for (int x = 0; x < img->xsize; x++) {
                img->val[t][z][y][x] = (float*) iftAlloc(img->m, sizeof(float));
                if (img->val[t][z][y][x] == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftImageSequenceAllocVolume");
            }
        }
    }
}

iftImageSequence *iftLoadImageSequence(iftFileSet *fileset) {
    iftImageSequence *img = NULL;
    int bands;
    for (int i = 0; i < fileset->n; i++) {
        iftImage *vol = iftReadImageByExt(fileset->files[i]->path);
        if (img == NULL) {
            bands = iftIsColorImage(vol)? 3 : 1;
            img = iftCreateImageSequence(vol->xsize, vol->ysize, vol->zsize, fileset->n, bands);
        } else {
            iftSetVolumeInImageSequence(img, vol, i);
        }
    }
    img->files = iftCopyFileSet(fileset);
    return img;
}

void iftDestroyImageSequence(iftImageSequence **img) {
    iftImageSequence *aux;

    aux = *img;
    if (aux != NULL) {
        for (int t = 0; t < aux->tsize; t++) {
            if (aux->allocated[t] || aux->val[t]) {
                for (int z = 0; z < aux->zsize; z++) {
                    for (int y = 0; y < aux->ysize; y++) {
                        for (int x = 0; x < aux->xsize; x++) {
                            iftFree(aux->val[t][z][y][x]);
                        }
                        iftFree(aux->val[t][z][y]);
                    }
                    iftFree(aux->val[t][z]);
                }
                iftFree(aux->val[t]);
            }
        }
        iftFree(aux->allocated);
        iftFree(aux->val);
        iftDestroyFileSet(&aux->files);

        iftFree(aux);
        *img = NULL;
    }
}

inline char ift4DValidVoxel(const iftImageSequence *img, iftVoxel v) {
    if ((v.x >= 0)&&(v.x < img->xsize)&&
        (v.y >= 0)&&(v.y < img->ysize)&&
        (v.z >= 0)&&(v.z < img->zsize)&&
        (v.t >= 0)&&(v.t < img->tsize)&&
        img->allocated[v.t])
        return(1);
    else
        return(0);
}

long ift4DGetVoxelIndex(const iftImageSequence *img, iftVoxel u) {
    long x = u.x;
    long y = u.y;
    long z = u.z;
    long t = u.t;

    long p = x + y * img->xsize + z * img->xsize * img->ysize +
             t * img->xsize * img->zsize * img->ysize;

    return p;
}

iftVoxel ift4DGetVoxelCoord(const iftImageSequence *img, long p) {
    iftVoxel u;
    long xsize = img->xsize;
    long ysize = img->ysize;
    long zsize = img->zsize;
    ldiv_t res1 = ldiv(p, xsize * ysize * zsize);
    ldiv_t res2 = ldiv(res1.rem, xsize * ysize);
    ldiv_t res3 = ldiv(res2.rem, xsize);

    u.x = res3.rem;
    u.y = res3.quot;
    u.z = res2.quot;
    u.t = res1.quot;

    return(u);
}

float *ift4DGetVal(const iftImageSequence *img, long p) {
    iftVoxel u = ift4DGetVoxelCoord(img, p);
    //if (ift4DValidVoxel(img, u)) {
        return img->val[u.t][u.z][u.y][u.x];
    //} else
    //    return NULL;
}

int image_sequence_cost_function(iftGQueue *Q, int p) {
    iftImageSequence *img = (iftImageSequence *) Q->L.value_data;
    return ift4DGetVal(img, p)[0];
}

//TODO nelems must be long
iftGQueue *iftCreateGQueueWithImageSequence(int nbuckets, int nelems,
                                            iftImageSequence *img) {
    return iftCreateGQueueWithCostFunction(nbuckets, nelems, img,
                                           image_sequence_cost_function);
}

iftMImage *iftImageSequenceToMImage(const iftImageSequence *img, int time) {
    iftMImage *mimg = iftCreateMImage(img->xsize, img->ysize,
                                      img->zsize, img->m);

    #pragma omp parallel for
    for (int z = 0; z < img->zsize; z++) {
        for (int y = 0; y < img->ysize; y++) {
            for (int x = 0; x < img->xsize; x++) {
                iftVoxel u = {x, y, z};
                long p = iftMGetVoxelIndex(mimg, u);
                for (int b = 0; b < img->m; b++)
                    mimg->val[p][b] = img->val[time][z][y][x][b];
            }
        }
    }

    return mimg;
}

iftImageSequence *iftMImageToImageSequence(const iftMImage *mimg) {
    iftImageSequence *img = iftCreateImageSequence(mimg->xsize, mimg->ysize,
                                                   mimg->zsize, 1, mimg->m);

#pragma omp parallel for
    for(long p = 0; p < img->n; p++) {
        for (long b = 0; b < img->m; b++) {
            ift4DGetVal(img, p)[b] = mimg->val[p][b];
        }
    }

    return img;
}


void iftWriteImageSequence(const iftImageSequence *img, const char *filename) {
    FILE *fp = fopen(filename,"wb");

    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageSequence", filename);

    fprintf(fp,"IMGSEQ\n");
    // Write X Y Z T nbands
    fprintf(fp,"%d %d %d %d %d\n",img->xsize,img->ysize,img->zsize,img->tsize,img->m);
    // Write dx dy dz dt
    fprintf(fp,"%f %f %f %f\n",img->dx,img->dy,img->dz, img->dt);
    // Write color space
    fprintf(fp, "%d\n", img->cspace);
    //TODO fwrite(&img->cspace, sizeof (iftColorSpace), 1, fp);
    // Write pixels
    for (int t = 0; t < img->tsize; t++) {
        for (int z = 0; z < img->zsize; z++) {
            for (int y = 0; y < img->ysize; y++) {
                for (int x = 0; x < img->xsize; x++) {
                        fwrite(img->val[t][z][y][x], sizeof(float), img->m, fp);
                }
            }
        }
    }
    fclose(fp);
}

//TODO fileset
iftImageSequence *iftReadImageSequence(const char *filename) {

    iftImageSequence *img=NULL;
    int xsize,ysize,zsize,tsize;
    int m;
    char msg[200], type[10];
    FILE *fp = fopen(filename,"rb");

    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageSequence", filename);
    else {
        if (fscanf(fp,"%s\n",type) != 1)
            iftError("Reading error", "iftReadImageSequence");

        if((strcmp(type,"IMGSEQ")==0)){
            // Read X Y Z T nbands
            if(fscanf(fp,"%d %d %d %d %d\n",&xsize,&ysize,&zsize,&tsize, &m)!=5)
                iftError("Reading error", "iftReadImageSequence");
            img = iftCreateImageSequence(xsize,ysize,zsize,tsize,m);
            // Read dx dy dz dt
            if(fscanf(fp,"%f %f %f %f\n",&img->dx,&img->dy,&img->dz,&img->dt)!=4)
                iftError("Reading error", "iftReadImageSequence");

            // Read color space
            if(fscanf(fp,"%d\n",(int*)&img->cspace)!=1)
                iftError("Reading error", "iftReadImageSequence");
            // todo if (fread(&img->cspace, sizeof (iftColorSpace), 1, fp) != 1)
                //iftError("Reading error", "iftReadImageSequence");

            for (int t = 0; t < img->tsize; t++) {
                img->allocated[t] = true;
                for (int z = 0; z < img->zsize; z++) {
                    for (int y = 0; y < img->ysize; y++) {
                        for (int x = 0; x < img->xsize; x++) {
                            if(fread(img->val[t][z][y][x], sizeof(float), img->m, fp) != img->m) {
                                sprintf(msg, "Could not read pixels of the ImageSequence %d %d %d %d", x, y, z, t);
                                iftWarning(msg, "iftReadImageSequence");
                            }
                        }
                    }
                }
            }
        } else {
            iftError("Input image must be a ImageSequence image", "iftReadImageSequence");
        }
        fclose(fp);
    }

    return(img);
}

iftImageSequence *
iftReadImageAsImageSequence(const char *filename, const iftImageSequence *imgseq_template,
                            int current_time) {
    char *ext     = iftLowerString(iftFileExt(filename));
    iftImageSequence *out = NULL;

    if (iftCompareStrings(ext, ".mimg")) {
        iftMImage *mimg = iftReadMImage(filename);
        out = iftCreateImageSequenceFromImageSequence(imgseq_template);
        iftSetMImageInImageSequence(out, mimg, current_time);
        iftDestroyMImage(&mimg);
    } else if (iftCompareStrings(ext, ".imsq")) {
        out = iftReadImageSequence(filename);
    } else {
        iftImage *vol = iftReadImageByExt(filename);
        out = iftCreateImageSequenceFromImageSequence(imgseq_template);
        iftSetVolumeInImageSequence(out, vol, current_time);
        iftDestroyImage(&vol);
    }

    iftFree(ext);
    return out;
}


iftImage *iftExtractVolumeFromImageSequence(const iftImageSequence *img, int time)  {
    iftImage *out = iftCreateImage(img->xsize, img->ysize, img->zsize);

    if (img->cspace == YCbCr_CSPACE) {
        iftSetCbCr(out, 0);
    }

    #pragma omp parallel for collapse(3)
    for (int z = 0; z < img->zsize; z++) {
        for (int y = 0; y < img->ysize; y++) {
            for (int x = 0; x < img->xsize; x++) {
                iftVoxel u = {x, y, z};
                long p = iftGetVoxelIndex(out, u);
                out->val[p] = img->val[time][z][y][x][0];
                if(ift4DIsColorImage(img)) {
                    out->val[p] = img->val[time][z][y][x][0];
                    out->Cb[p]  = img->val[time][z][y][x][1];
                    out->Cr[p]  = img->val[time][z][y][x][2];
                } else //if (img->m == 1) {
                    out->val[p] = img->val[time][z][y][x][0];
//                }else {
//                    iftError("Color space not recognized.", "iftExtractVolumeFromImageSequence");
//                }

            }
        }
    }

    return out;
}

iftImageSequence *iftVolumeToImageSequence(const iftImage *img) {
    int bands = iftIsColorImage(img)? 3 : 1;
    iftImageSequence *out = iftCreateImageSequence(img->xsize, img->ysize,
                                                   img->zsize, 1, bands);
    if (iftIsColorImage(img)) {
        out->cspace = YCbCr_CSPACE;
    } else {
        out->cspace = GRAY_CSPACE;
    }

    #pragma omp parallel for
    for(long p = 0; p < img->n; p++) {
        ift4DGetVal(out, p)[0] = img->val[p];
        if (bands == 3) {
            ift4DGetVal(out, p)[1] = img->Cb[p];
            ift4DGetVal(out, p)[2] = img->Cr[p];
        }
    }
    return out;
}

void iftAppendVolumeToImageSequence(iftImageSequence *img, iftImage *vol) {
    img->val = (float*****) iftRealloc(img->val, sizeof *img->val * (img->tsize + 1));
    img->allocated = (bool*) iftRealloc(img->allocated, sizeof *img->allocated * (img->tsize + 1));

    int t = img->tsize;
    img->tsize++;
    int xsize = img->xsize;
    int ysize = img->ysize;
    int zsize = img->zsize;
    img->n += xsize*zsize*ysize;
    int nbands = img->m;

    if ((iftIsColorImage(vol) && img->cspace == YCbCr_CSPACE) ||
        !(iftIsColorImage(vol) || ift4DIsColorImage(img))) {

        img->val[t] = (float****) iftAlloc(zsize, sizeof(float***));
        img->allocated[t] = true;
        if (img->val[t] == NULL) iftError(MSG_MEMORY_ALLOC_ERROR,
                                          "iftAppendVolumetoImageSequence");
        for (int z = 0; z < zsize; z++) {
            img->val[t][z] = (float ***) iftAlloc(ysize, sizeof(float **));
            if (img->val[t][z] == NULL)
                iftError(MSG_MEMORY_ALLOC_ERROR, "iftAppendVolumetoImageSequence");
            for (int y = 0; y < ysize; y++) {
                img->val[t][z][y] = (float **) iftAlloc(xsize, sizeof(float *));
                if (img->val[t][z][y] == NULL)
                    iftError(MSG_MEMORY_ALLOC_ERROR,
                             "iftAppendVolumetoImageSequence");
                for (int x = 0; x < xsize; x++) {
                    img->val[t][z][y][x] = (float *) iftAlloc(nbands,
                                                              sizeof(float));
                    if (img->val[t][z][y][x] == NULL)
                        iftError(MSG_MEMORY_ALLOC_ERROR,
                                 "iftAppendVolumetoImageSequence");

                    iftVoxel u = {x, y, z};
                    long p = iftGetVoxelIndex(vol, u);

                    img->val[t][z][y][x][0] = vol->val[p];
                    if (iftIsColorImage(vol)) {
                        img->val[t][z][y][x][1] = vol->Cb[p];
                        img->val[t][z][y][x][2] = vol->Cr[p];
                    }
                }
            }
        }
    } else {
        iftError("Color space of iftImage does not match with the color "
                 "space of iftImageSequence", "iftAppendVolumeToImageSequence");
    }



}

bool iftIs3DImageSequence(const iftImageSequence *img) {
    return img->zsize > 1;
}

bool iftIs4DImageSequence(const iftImageSequence *img) {
    return img->tsize > 1;
}

bool ift4DIsColorImage(const iftImageSequence *img) {
    return img->cspace != GRAY_CSPACE && img->cspace != GRAYNorm_CSPACE && img->cspace != UNDEFINED_CSPACE;
}

void
iftSetVolumeInImageSequence(iftImageSequence *img, iftImage *vol, int t) {
    int xsize = img->xsize;
    int ysize = img->ysize;
    int zsize = img->zsize;

    if (!img->allocated[t]) {
        iftImageSequenceAllocVolume(img, t);
    }

    if (img->xsize != vol->xsize || img->ysize != vol->ysize || img->zsize != vol->zsize || t < 0 || t >= img->tsize) {
        iftError("Dimensions of volumetric image does not match with the temporal image.", "iftSetVolumeInImageSequence");
    }

    #pragma omp parallel for collapse(3)
    for (int z = 0; z < zsize; z++) {
        for (int y = 0; y < ysize; y++) {
            for (int x = 0; x < xsize; x++) {
                iftVoxel u = {x, y, z};
                long p = iftGetVoxelIndex(vol, u);
                img->val[t][z][y][x][0] = vol->val[p];
                if (iftIsColorImage(vol) && ift4DIsColorImage(img)) {
                    img->val[t][z][y][x][1] = vol->Cb[p];
                    img->val[t][z][y][x][2] = vol->Cr[p];
                }
            }
        }
    }
}

void
iftSetMImageInImageSequence(iftImageSequence *img, iftMImage *mimg, int time) {
    int xsize = img->xsize;
    int ysize = img->ysize;
    int zsize = img->zsize;

    if (mimg->m != img->m || mimg->xsize != xsize ||
        mimg->ysize != ysize || mimg->zsize != zsize) {
        iftError("MImage does not match Image Sequence dimensions.",
                 "iftSetMImageInImageSequence");
    }

    #pragma omp parallel for collapse(3)
    for (int z = 0; z < zsize; z++) {
        for (int y = 0; y < ysize; y++) {
            for (int x = 0; x < xsize; x++) {
                iftVoxel u = {x, y, z};
                long p = iftMGetVoxelIndex(mimg, u);
                for (int b = 0; b < mimg->m; b++)
                    img->val[time][z][y][x][b] = mimg->val[p][b];
            }
        }
    }
}

inline bool iftIsImageCompatibleToImageSequence(const iftImageSequence *img,
                                         iftImage *vol) {
    //TODO check color space
    if (img->zsize == vol->zsize && img->ysize == vol->ysize &&
        img->xsize == vol->xsize /*&&  img->dx == vol->dx &&
        img->dy == vol->dy && img->dz == vol->dz*/)
        return true;
    else
        return false;
}

iftMImage *iftGetXYSliceFromTime(const iftImageSequence *img, int zcoord, int time) {
    iftMImage *slice;

    if ( (zcoord < 0) || (zcoord >= img->zsize))
        iftError("Invalid z coordinate", "iftGetXYSliceFromTime");

    slice = iftCreateMImage(img->xsize,img->ysize, 1, img->m);

    int z = zcoord;
    int t = time;
    int q = 0;
    for (int y = 0; y < img->ysize; y++)
        for (int x = 0; x < img->xsize; x++)
        {
            for (int b = 0; b < img->m; b++) {
                slice->val[q][b] = img->val[t][z][y][x][b];
            }
            q++;
        }
    //TODO iftMCopyVoxelSize(mimg,slice);

    return(slice);
}

iftMImage *iftGetZXSliceFromTime(const iftImageSequence *img, int ycoord, int time) {
    iftMImage *slice;

    if ( (ycoord < 0) || (ycoord >= img->ysize))
        iftError("Invalid y coordinate", "iftGetZXSliceFromTime");

    slice = iftCreateMImage(img->zsize,img->xsize, 1, img->m);

    int y = ycoord;
    int t = time;
    int q = 0;
    for (int x = 0; x < img->xsize; x++)
        for (int z = 0; z < img->zsize; z++)
        {
            for (int b = 0; b < img->m; b++) {
                slice->val[q][b] = img->val[t][z][y][x][b];
            }
            q++;
        }
    //TODO iftMCopyVoxelSize(mimg,slice);

    return(slice);
}

iftMImage *iftGetYZSliceFromTime(const iftImageSequence *img, int xcoord, int time) {
    iftMImage *slice;

    if ( (xcoord < 0) || (xcoord >= img->xsize))
        iftError("Invalid x coordinate", "iftGetYZSliceFromTime");

    slice = iftCreateMImage(img->ysize,img->zsize, 1, img->m);

    int x = xcoord;
    int t = time;
    int q = 0;
    for (int z = 0; z < img->zsize; z++)
        for (int y = 0; y < img->ysize; y++)
        {
            for (int b = 0; b < img->m; b++) {
                slice->val[q][b] = img->val[t][z][y][x][b];
            }
            q++;
        }
    //TODO iftMCopyVoxelSize(mimg,slice);

    return(slice);
}

iftImageSequence *ift4DNormalize(iftImageSequence *img, float min, float max) {
    iftImageSequence *nimg = iftCreateImageSequenceFromImageSequence(img);

    float img_min_val;
    float img_max_val;
    ift4DMinMaxValues(img, 0, &img_min_val, &img_max_val);


    if (img_min_val < img_max_val) {
        IMGSEQ_FORLOOP_PARALLEL(img) {
            for (int b = 0; b < img->m; b++) {
                IMGSEQ_ITE_VAL(nimg, b) = (int) ((max - min) * ((double) IMGSEQ_ITE_VAL(img, b) - (double) img_min_val) /
                                                 ((double) img_max_val - (double) img_min_val) + min);
            }
        }
    } else iftError("Cannot normalize empty image", "ift4DNormalize");

    return nimg;
}

float ift4DMaximumValue(const iftImageSequence *img, int band) {
    float max = IFT_INFINITY_FLT_NEG;

    IMGSEQ_FORLOOP(img) {
        if (IMGSEQ_ITE_VAL(img, band) > max)
            max = IMGSEQ_ITE_VAL(img, band);
    }

    return max;
}

float ift4DMinimumValue(const iftImageSequence *img, int band) {
    float min = IFT_INFINITY_FLT;

    IMGSEQ_FORLOOP(img) {
        if (IMGSEQ_ITE_VAL(img, band) < min)
            min = IMGSEQ_ITE_VAL(img, band);
    }

    return min;
}

iftImageSequence *iftCopyImageSequence(const iftImageSequence *img) {
    iftImageSequence *out = iftCreateImageSequenceFromImageSequence(img);

    IMGSEQ_FORLOOP_PARALLEL(img) {
        for (int b = 0; b < img->m; b++) {
            IMGSEQ_ITE_VAL(out, b) = IMGSEQ_ITE_VAL(img, b);
        }
    }

    for (int t = 0; t < img->tsize; t++) {
        out->allocated[t] = img->allocated[t];
    }

    return out;
}

iftImageSequence *
iftCreateImageSequenceFromImageSequence(const iftImageSequence *img) {
    iftImageSequence *out = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, img->m);
    int depth = ift4DImageDepth(img);
    int value = (iftMaxImageRange(depth)+1)/2;
    if (img->cspace == YCbCr_CSPACE) {
        IMGSEQ_FORLOOP_PARALLEL(out) {
            IMGSEQ_ITE_VAL(out, 1) = value;
            IMGSEQ_ITE_VAL(out, 2) = value;
        }
    }
    return out;
}

void iftSetImageSequence(iftImageSequence *img, float value) {
    IMGSEQ_FORLOOP_PARALLEL(img) {
        for (int b = 0; b < img->m; b++)
            IMGSEQ_ITE_VAL(img, b) = value;
    }
}

double ift4DDiceSimilarity(const iftImageSequence *label_source,
                           const iftImageSequence *label_target) {
    //TODO iftVerifyImageDomains(bin_source, bin_target, "iftDiceSimilarity");

    ulong vol_intersec = 0; // |intersection(volume(bin_source), volume(bin_target))|
    ulong vol_sum      = 0; // |volume(bin_source)| + |volume(bin_target)|

    IMGSEQ_FORLOOP(label_source) {
        iftVoxel u = {x,y,z,t};
        if (ift4DValidVoxel(label_source, u)) {
            vol_intersec += ((IMGSEQ_ITE_VAL(label_source, 0) != 0) &&
                             (IMGSEQ_ITE_VAL(label_target, 0) != 0));
            vol_sum += ((IMGSEQ_ITE_VAL(label_source, 0) != 0) +
                        (IMGSEQ_ITE_VAL(label_target, 0) != 0));
        }
    }

    double dice = (2.0 * vol_intersec) / (vol_sum);

    return dice;
}

iftDblArray *ift4DDiceSimilarityMultiLabel(const iftImageSequence *label_source,
                                           const iftImageSequence *label_target) {
    int n_objects = ift4DMaximumValue(label_target, 0);

    iftIntArray *vol_intersec = iftCreateIntArray(n_objects + 1);
    iftIntArray *vol_sum      = iftCreateIntArray(n_objects + 1);
    iftDblArray *dices        = iftCreateDblArray(n_objects + 1);

    IMGSEQ_FORLOOP(label_source) {
        iftVoxel u = {x,y,z,t};
        if (ift4DValidVoxel(label_source, u)) {
            int label_source_val = IMGSEQ_ITE_VAL(label_source, 0);
            int label_target_val = IMGSEQ_ITE_VAL(label_target, 0);
            vol_intersec->val[label_source_val] += (label_source_val == label_target_val);
            vol_sum->val[label_source_val]++;
            vol_sum->val[label_target_val]++;
        }
    }

    for (int o = 1; o <= n_objects; o++) {
        dices->val[o] = (2.0 * vol_intersec->val[o]) / vol_sum->val[o];
        dices->val[0] += dices->val[o];
    }
    dices->val[0] /= n_objects;

    iftDestroyIntArray(&vol_intersec);
    iftDestroyIntArray(&vol_sum);

    return dices;
}

//TODO continue refactoring from here
iftImageSequence *ift4DObjectBorders(const iftImageSequence *label_img, const iftAdjRel *Ain, bool keep_border_labels,
                                     bool include_image_frame) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = iftIs4DImageSequence(label_img) ? iftHyperSpheric(1.0) :
                (iftIs3DImageSequence(label_img)) ? iftSpheric(1.0) : iftCircular(1.0);
    else A = iftCopyAdjacency(Ain);

    iftImageSequence *border_img = iftCreateImageSequenceFromImageSequence(label_img);

    for (int p = 0; p < label_img->n; p++) {
        if (ift4DGetVal(label_img, p)[0] != 0) {
            iftVoxel u = ift4DGetVoxelCoord(label_img, p);
            if (ift4DValidVoxel(label_img, u)) {

                for (int i = 1; i < A->n; i++) {
                    iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                    if (ift4DValidVoxel(label_img, v)) {
                        if (label_img->val[u.t][u.z][u.y][u.x][0] != label_img->val[v.t][v.z][v.y][v.x][0]) {
                            border_img->val[u.t][u.z][u.y][u.x][0] = label_img->val[u.t][u.z][u.y][u.x][0];
                            break;
                        }
                    }
                        // voxel u belongs to the image's frame
                    else if (include_image_frame && v.t < label_img->tsize && v.t >= 0 && v.z < label_img->zsize &&
                             v.z >= 0 && label_img->allocated[v.t]) { // TODO
                        border_img->val[u.t][u.z][u.y][u.x][0] = ift4DValidVoxel(label_img, u);
                        break;
                    }
                }
            }
        }
    }

    if (!keep_border_labels) {
        //binarize
        IMGSEQ_FORLOOP_PARALLEL(border_img) {
            IMGSEQ_ITE_VAL(border_img, 0) = IMGSEQ_ITE_VAL(border_img, 0) != 0;
        }
    }

    if (Ain == NULL)
        iftDestroyAdjRel(&A);

    return border_img;
}

iftImageSequence *ift4DFastLabelComp(const iftImageSequence *bin, const iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        A = iftHyperCuboid(3,3,3,3);
    }
    else A = iftCopyAdjacency(Ain);

    iftImageSequence *label=NULL;
    int i,p,q,l=1, *cost;
    iftVoxel u,v;
    iftGQueue *Q;

    label  = iftCreateImageSequenceFromImageSequence(bin);
    cost   = iftAllocIntArray(bin->n);
    Q      = iftCreateGQueue(2,bin->n,cost);

    for (p=0; p < bin->n; p++){
        if ((ift4DGetVal(bin, p)[0] != 0) && (ift4DGetVal(label, p)[0] == 0)){
            cost[p] = 1;
            int nbuckets = Q->C.nbuckets;
            iftInsertGQueue(&Q,p);
            if (nbuckets != Q->C.nbuckets) {
                printf("HUMM1\n");
            }
        }
    }

    while(!iftEmptyGQueue(Q)){
        p = iftRemoveGQueue(Q);
        if (cost[p]==1){ /* p is the first in a component */
            cost[p]=0;
            ift4DGetVal(label, p)[0] = l; l++;
        }
        u = ift4DGetVoxelCoord(bin, p);
        for (i=1; i < A->n; i++){
            v = iftGetAdjacentVoxel(A,u,i);
            if (ift4DValidVoxel(bin, v)){
                q = ift4DGetVoxelIndex(bin, v);
                if ((ift4DGetVal(bin, p)[0] == ift4DGetVal(bin, q)[0]) && (ift4DGetVal(label, q)[0] == 0)){
                    ift4DGetVal(label, q)[0] = ift4DGetVal(label, p)[0];
                    iftRemoveGQueueElem(Q,q);
                    cost[q]=0;
                    int nbuckets = Q->C.nbuckets;
                    iftInsertGQueue(&Q,q);
                    if (nbuckets != Q->C.nbuckets) {
                        printf("HUMM2\n");
                    }
                }
            }
        }
    }

    iftDestroyAdjRel(&A);

    iftDestroyGQueue(&Q);
    iftFree(cost);
    return(label);
}

iftImageSequence *ift4DFindAndLabelObjectBorders(const iftImageSequence *label_img,
                                                 const iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = iftIs4DImageSequence(label_img) ? iftHyperCuboid(3, 3, 3, 3) :
                (iftIs3DImageSequence(label_img)) ? iftSpheric(sqrtf(3)) : iftCircular(sqrtf(2));
    else A = iftCopyAdjacency(Ain);

    iftImageSequence *border_img = ift4DObjectBorders(label_img, NULL, true, true);
    iftImageSequence *border_label_img = ift4DFastLabelComp(border_img, A);

    iftDestroyImageSequence(&border_img);

    iftDestroyAdjRel(&A);

    return border_label_img;
}

iftImageSequence *ift4DEuclDistTrans(const iftImageSequence *label_img,
                                     const iftAdjRel *Ain, iftSide side,
                                     iftImageSequence **root_out,
                                     iftImageSequence **edt_label_out,
                                     iftImageSequence **pred_out) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = iftHyperSpheric(1.0);
    else A = iftCopyAdjacency(Ain);

    // Initialization
    iftImageSequence *dist      = iftCreateImageSequenceFromImageSequence(label_img);
    iftImageSequence *root      = iftCreateImageSequence(label_img->xsize, label_img->ysize, label_img->zsize, label_img->tsize, 4);
    iftImageSequence *pred      = iftCreateImageSequence(label_img->xsize, label_img->ysize, label_img->zsize, label_img->tsize, 4);
    iftImageSequence *edt_label = ift4DFindAndLabelObjectBorders(label_img, iftSpheric(1.73));

    iftGQueue *Q = iftCreateGQueueWithImageSequence(IFT_QSIZE, label_img->n, dist);
    iftIntializeDistTransCost initializationCostFunction = iftGetIntializeDistTransCost(side);

    IMGSEQ_FORLOOP(label_img) {
        iftVoxel v = {x,y,z,t};
        if (ift4DValidVoxel(label_img, v)) {
            IMGSEQ_ITE_VAL(dist, 0) = initializationCostFunction((int) IMGSEQ_VOXEL_VAL(label_img, v, 0));
            if (IMGSEQ_VOXEL_VAL(edt_label, v, 0) != 0) {
                IMGSEQ_VOXEL_VAL(dist, v, 0) = 0;
                long p = ift4DGetVoxelIndex(label_img, v);
                setPixelCoordinate(root, p, v);
                setPixelCoordinate(pred, p, U_NIL);
                iftInsertGQueue(&Q, p);
            }
        }
    }
//    for (int p = 0; p < label_img->n; p++) {
//        ift4DGetVal(dist, p)[0] = initializationCostFunction((int) ift4DGetVal(label_img, p)[0]);
//
//        if (ift4DGetVal(edt_label, p)[0] != 0) {
//            iftVoxel u = ift4DGetVoxelCoord(root, p);
//            ift4DGetVal(dist, p)[0] = 0;
//            setPixelCoordinate(root, p, u);
//            // -> ift4DGetVal(root, p)[0] = p;
//            setPixelCoordinate(pred, p, U_NIL);
//            // -> ift4DGetVal(pred, p)[0] = IFT_NIL;
//            iftInsertGQueue(&Q,p);
//        }
//    }


    // Image Foresting Transform
    while(!iftEmptyGQueue(Q)) {
        long p = iftRemoveGQueue(Q);

        iftVoxel u = ift4DGetVoxelCoord(label_img, p);
        int _r = getPixelCoordinate(root, p);
        iftVoxel r = ift4DGetVoxelCoord(label_img,_r);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (ift4DValidVoxel(label_img, v)) {
                int q = ift4DGetVoxelIndex(label_img, v);

                if (ift4DGetVal(dist, q)[0] > ift4DGetVal(dist, p)[0]) {
                    float tmp = (v.x-r.x)*(v.x-r.x) + (v.y-r.y)*(v.y-r.y) + (v.z-r.z)*(v.z-r.z) + (v.t-r.t)*(v.t-r.t);

                    if (tmp < ift4DGetVal(dist, q)[0]) {
                        if (ift4DGetVal(dist, q)[0] != IFT_INFINITY_INT)
                            iftRemoveGQueueElem(Q, q);
                        ift4DGetVal(dist, q)[0]      = (int) tmp;
                        setPixelCoordinateFromArray(root, q, ift4DGetVal(root, p));
                        // -> ift4DGetVal(root, q)[0]      = ift4DGetVal(root, p)[0];
                        ift4DGetVal(edt_label, q)[0] = ift4DGetVal(edt_label, p)[0];
                        setPixelCoordinate(pred, q, u);
                        // -> ift4DGetVal(pred, q)[0]      = p;
                        iftInsertGQueue(&Q,q);
                    }
                }
            }
        }
    }

    iftDestroyGQueue(&Q);


    if (root_out == NULL)
        iftDestroyImageSequence(&root);
    else *root_out = root;

    if (edt_label_out == NULL)
        iftDestroyImageSequence(&edt_label);
    else *edt_label_out = edt_label;

    if (pred_out == NULL)
        iftDestroyImageSequence(&pred);
    else *pred_out = pred;

    iftDestroyAdjRel(&A);

    return dist;
}

void ift4DMinMaxValues(const iftImageSequence *img, int band, float *min,
                       float *max) {
    *min = IFT_INFINITY_FLT;
    *max = IFT_INFINITY_FLT_NEG;

    for (int t = 0; t < img->tsize; t++)
        for (int z = 0; z < img->zsize; z++)
            for (int y = 0; y < img->ysize; y++)
                for (int x = 0; x < img->xsize; x++)
                    if (img->val[t][z][y][x][band] < *min)
                        *min = img->val[t][z][y][x][band];
                    else if (img->val[t][z][y][x][band] > *max)
                        *max = img->val[t][z][y][x][band];
}

iftImageSequence *ift4DGradient(const iftImageSequence *img, iftAdjRel *A, int Imax) {
    iftImageSequence *grad = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, 1);

    float Gmin = IFT_INFINITY_FLT;
    float Gmax = IFT_INFINITY_FLT_NEG;

    IMGSEQ_FORLOOP(img) {
        iftVoxel u = {x, y, z, t};
        float value = 0;
        for (int i=1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);
            if (ift4DValidVoxel(img, v)){
                for (int b=0; b < img->m; b++)
                    value += pow(img->val[v.t][v.z][v.y][v.x][b]-img->val[t][z][y][x][b], 2);
            }
        }
        value = sqrtf(value)/A->n;
        if (value > Gmax) Gmax = value;
        if (value < Gmin) Gmin = value;
        grad->val[t][z][y][x][0] = value;
    }

    IMGSEQ_FORLOOP(img) {
        grad->val[t][z][y][x][0] = (Imax*(grad->val[t][z][y][x][0]-Gmin)/(Gmax-Gmin));
    }

    return grad;
}



iftImageSequence *
ift4DBorderImage(const iftImageSequence *label, bool get_margins) {
    iftAdjRel *A;
    iftImageSequence  *border = iftCreateImageSequenceFromImageSequence(label);
    iftVoxel   v;

    if (iftIs4DImageSequence(label))
        A = iftHyperSpheric(1.0);
    if (iftIs3DImageSequence(label))
        A = iftSpheric(1.0);
    else
        A = iftCircular(1.0);

    IMGSEQ_FORLOOP(label) {
        iftVoxel u = {x, y, z, t};
        for (long i = 1; i < A->n; i++) {
            v = iftGetAdjacentVoxel(A, u, i);
            if (ift4DValidVoxel(label, v)) {
                if (label->val[t][z][y][x][0] !=
                    label->val[v.t][v.z][v.y][v.x][0]) {
                    border->val[t][z][y][x][0] = label->val[t][z][y][x][0];
                    break;
                }
            } else if (get_margins) {
                border->val[t][z][y][x][0] = label->val[t][z][y][x][0];
            }
        }
    }

    iftDestroyAdjRel(&A);
    return(border);
}

long iftCountObjectSpelsFromTime(const iftImageSequence *img, int obj_label, int time) {
    int count = 0;
    for (int z = 0; z < img->zsize; z++)
        for (int y = 0; y < img->ysize; y++)
            for (int x = 0; x < img->xsize; x++)
                if (img->val[time][z][y][x][0] == obj_label)
                    count++;
    return count;
}

iftLabeledSet * ift4DReadSeeds(const iftImageSequence *img, const char *_filename, ...)
{
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];
    va_start(args, _filename);
    vsprintf(filename, _filename, args);
    va_end(args);


    FILE *fp=fopen(filename,"r");
    int i,l,mk,nseeds,handicap;
    int w,h,d,t;
    iftVoxel v;
    iftLabeledSet *seed=NULL;

    if (fscanf(fp,"%d %d %d %d %d",&nseeds, &w, &h, &d, &t) != 5) {
        iftError("Reading error", "iftReadSeeds");
    }

    for (i=0; i < nseeds; i++){
        if (fscanf(fp,"%d %d %d %d %d %d %d",&v.x,&v.y,&v.z,&v.t,&mk,&l, &handicap)!=7) {
            iftError("Reading error", "iftReadSeeds");
        }
        if (ift4DValidVoxel(img, v)){
            iftInsertLabeledSetMarkerAndHandicap(&seed, ift4DGetVoxelIndex(img, v), l, mk, handicap);
        }
    }

    fclose(fp);
    return(seed);
}


void ift4DWriteSeeds(iftLabeledSet *seed, const iftImageSequence *img, const char *_filename, ...)
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

    fprintf(file,"%d %d %d %d %d\n", nseeds, img->xsize, img->ysize, img->zsize, img->tsize);
    s = seed;
    while(s != NULL){
        iftVoxel voxel = ift4DGetVoxelCoord(img, s->elem);
        fprintf(file, "%d %d %d %d %d %d %d\n", voxel.x, voxel.y, voxel.z, voxel.t, s->marker, s->label, s->handicap);
        s = s->next;
    }

    fclose(file);
}

iftImageSequence *
ift4DThreshold(const iftImageSequence *img, float lowest,
               float highest, float value) {
    iftImageSequence *bin = iftCreateImageSequenceFromImageSequence(img);

    IMGSEQ_FORLOOP(img) {
        if ((img->val[t][z][y][x][0] >= lowest) && (img->val[t][z][y][x][0] <= highest))
            bin->val[t][z][y][x][0] = value;
        else bin->val[t][z][y][x][0] = 0;
    }
    return bin;
}

iftImageSequence *
ift4DThresholdInTime(const iftImageSequence *img, float lowest,
                     float highest, float value, int time) {
    iftImageSequence *bin = iftCreateImageSequenceFromImageSequence(img);


    for (int z = 0; z < bin->zsize; z++)
        for (int y = 0; y < bin->ysize; y++)
            for (int x = 0; x < bin->xsize; x++)
                if ((img->val[time][z][y][x][0] >= lowest) && (img->val[time][z][y][x][0] <= highest))
                    bin->val[time][z][y][x][0] = value;
                else bin->val[time][z][y][x][0] = 0;
    return bin;
}

uchar ift4DImageDepth(const iftImageSequence *img) {
    float img_min, img_max;
    ift4DMinMaxValues(img, 0, &img_min, &img_max);

    long max_range;

    if (img_min >= 0)
        max_range = iftNormalizationValue(img_max) + 1;
    else
        max_range = iftNormalizationValue(img_max - img_min) + 1;

    return (uchar) iftLog(max_range, 2);
}

ift4DImageForest *
iftCreate4DImageForest(iftImageSequence *img, iftAdjRel *A) {
    ift4DImageForest *fst=(ift4DImageForest *)iftAlloc(1, sizeof(ift4DImageForest));

    if (fst != NULL) {
        fst->pathval        = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, 1);
        fst->label          = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, 1);
        fst->root           = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, 1);
        fst->pred           = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, 1);
        fst->Q              = iftCreateGQueueWithImageSequence(iftMax(IFT_QSIZE, ift4DMaximumValue(img, 0) + 1),
                                                               img->n,
                                                               fst->pathval);
        fst->img            = img;
        fst->A              = A;
    } else {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreate4DImageForest");
    }

    // Default
    iftReset4DImageForest(fst);

//    TODO iftCopyVoxelSize(img, fst->pathval);
//    iftCopyVoxelSize(img, fst->label);
//    iftCopyVoxelSize(img, fst->marker);
//    iftCopyVoxelSize(img, fst->root);
//    iftCopyVoxelSize(img, fst->pred);

    return(fst);
}

void iftReset4DImageForest(ift4DImageForest *fst) {
    iftResetGQueue(fst->Q);

    if (fst->Q->C.removal_policy == MINVALUE){
        IMGSEQ_FORLOOP(fst->img) {
            iftVoxel u = {x,y,z,t};
            fst->pathval->val[t][z][y][x][0] = IFT_INFINITY_INT;
            fst->label->val[t][z][y][x][0]  = 0;
            fst->pred->val[t][z][y][x][0]  = IFT_NIL;
            fst->root->val[t][z][y][x][0] = ift4DGetVoxelIndex(fst->img, u);
        }
    }
}

void iftDestroy4DImageForest(ift4DImageForest **fst) {
    ift4DImageForest *aux = *fst;

    if (aux != NULL) {
//        if (aux->processed != NULL)
//            iftDestroyBMap(&aux->processed);
        if (aux->pathval != NULL)
            iftDestroyImageSequence(&aux->pathval);
        if (aux->label != NULL)
            iftDestroyImageSequence(&aux->label);
        if (aux->pred != NULL)
            iftDestroyImageSequence(&aux->pred);
        if (aux->root != NULL)
            iftDestroyImageSequence(&aux->root);
        if (aux->Q != NULL)
            iftDestroyGQueue(&aux->Q);

        iftFree(aux);
        *fst = NULL;
    }
}

void
ift4DDiffWatershed(ift4DImageForest *fst, iftLabeledSet *seed,
                   iftSet *removal_markers) {
    iftAdjRel *A = fst->A;
    iftGQueue *Q = fst->Q;
    iftVoxel   u, v;
    long       i, p, q;
    float      tmp;
    iftSet    *Frontier = NULL;
    iftLabeledSet *S;
    iftImageSequence  *pathval = fst->pathval, *pred = fst->pred, *label = fst->label;
    iftImageSequence  *root = fst->root, *basins = fst->img;

    if (removal_markers != NULL)
    {
        Frontier = ift4DTreeRemoval(fst, removal_markers);
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

        ift4DGetVal(label, p)[0]   = S->label;
        ift4DGetVal(pathval, p)[0] = 0;
        ift4DGetVal(root, p)[0]    = p;
        ift4DGetVal(pred, p)[0]    = IFT_NIL;
        iftInsertGQueue(&Q, p);
        S = S->next;
    }

    /* Image Foresting Transform */
    while (!iftEmptyGQueue(Q))
    {
        p = iftRemoveGQueue(Q);
        u = ift4DGetVoxelCoord(basins, p);

        for (i = 1; i < A->n; i++)
        {
            v = iftGetAdjacentVoxel(A, u, i);
            if (ift4DValidVoxel(basins, v))
            {
                q = ift4DGetVoxelIndex(basins, v);

                if (Q->L.elem[q].color != IFT_BLACK){

                    tmp = iftMax(ift4DGetVal(pathval, p)[0],
                                 ift4DGetVal(basins, q)[0]);

                    /* if pred[q]=p then p and q belong to a tie-zone */
                    if ((tmp < ift4DGetVal(pathval, q)[0]) || ((ift4DGetVal(pred, q)[0] == p)))
                    {
                        if (Q->L.elem[q].color == IFT_GRAY)
                        {
                            iftRemoveGQueueElem(Q, q);
                        }
                        ift4DGetVal(pred, q)[0]    = p;
                        ift4DGetVal(root, q)[0]    = ift4DGetVal(root, p)[0];
                        ift4DGetVal(label, q)[0]   = ift4DGetVal(label, p)[0];
                        ift4DGetVal(pathval, q)[0] = tmp;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
    }

    iftResetGQueue(Q);

}


iftSet *ift4DTreeRemoval(ift4DImageForest *fst,
                         iftSet *trees_for_removal) {
    long        i, p, q, r, n = fst->img->n;
    iftVoxel   u, v;
    iftAdjRel *A = iftHyperSpheric(1.0);//fst->A;
    iftSet   *Frontier = NULL;
    iftBMap  *inFrontier = iftCreateBMap(n);
    iftImageSequence *pathval = fst->pathval, *pred = fst->pred, *root = fst->root;
    iftImageSequence *img = fst->img;
    iftSet   *T1 = NULL, *T2 = NULL;

    /* Remove all marked trees and find the frontier voxels
       afterwards. */

    while (trees_for_removal != NULL)
    {
        p = trees_for_removal->elem;
        r = ift4DGetVal(root, p)[0];
        if (ift4DGetVal(pathval, r)[0] != IFT_INFINITY_INT) //tree not marked yet
        {
            ift4DGetVal(pathval, r)[0] = IFT_INFINITY_INT; // mark removed node
            ift4DGetVal(pred, r)[0] = IFT_NIL;
            iftInsertSet(&T1, r);
            while (T1 != NULL)
            {
                p = iftRemoveSet(&T1);
                iftInsertSet(&T2, p);

                u = ift4DGetVoxelCoord(img, p);
                for (i = 1; i < A->n; i++)
                {
                    v = iftGetAdjacentVoxel(A, u, i);
                    if (ift4DValidVoxel(img, v))
                    {
                        q   = ift4DGetVoxelIndex(img, v);
                        if (ift4DGetVal(pathval, q)[0] != IFT_INFINITY_INT)
                        {
                            if (ift4DGetVal(pred, q)[0] == p)
                            {
                                iftInsertSet(&T1, q);

                                ift4DGetVal(pathval, q)[0] = IFT_INFINITY_INT; // mark removed node
                                ift4DGetVal(pred, q)[0]    = IFT_NIL;
                            }
                        }
                    }
                }
            }
        }
        trees_for_removal = trees_for_removal->next;
    }

    /* Find the frontier voxels of non-removed trees */

    while (T2 != NULL)
    {
        p = iftRemoveSet(&T2);
        u = ift4DGetVoxelCoord(img, p);
        for (i = 1; i < A->n; i++)
        {
            v = iftGetAdjacentVoxel(A, u, i);
            if (ift4DValidVoxel(img, v))
            {
                q   = ift4DGetVoxelIndex(img, v);
                if (ift4DGetVal(pathval, q)[0] != IFT_INFINITY_INT)
                {
                    if (iftBMapValue(inFrontier, q) == 0)
                    {
                        iftInsertSet(&Frontier, q);
                        iftBMapSet1(inFrontier, q);
                    }
                }
            }
        }
    }
    iftDestroyBMap(&inFrontier);

    return (Frontier);
}

ift4DDynTrees *
iftCreate4DDynTrees(iftImageSequence *img, iftAdjRel *A) {
    ift4DDynTrees *dt = (ift4DDynTrees *)calloc(1, sizeof(ift4DDynTrees));

    dt->img      = img; dt->A = A;
    dt->label    = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, 1);
    dt->root     = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, 4);
    dt->cost     = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, 1);
    dt->pred     = iftCreateImageSequence(img->xsize, img->ysize, img->zsize, img->tsize, 4);
    dt->cumfeat  = iftCreateImageSequenceFromImageSequence(img);
    dt->treesize = iftAllocIntArray(img->n);
    dt->Q        = iftCreateGQueueWithImageSequence(CMAX+1, img->n, dt->cost);

    /* compute the maximum feature distance */

    dt->maxfeatdist = 0.0;
    for (long p=0; p < img->n; p++){
        iftVoxel u = ift4DGetVoxelCoord(img, p);
        for (int i=1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);
            if (ift4DValidVoxel(img, v)){
                float dist = 0.0;
                for (int b=0; b < img->m; b++){
                    dist += (img->val[v.t][v.z][v.y][v.x][b]-img->val[u.t][u.z][u.y][u.x][b])*
                            (img->val[v.t][v.z][v.y][v.x][b]-img->val[u.t][u.z][u.y][u.x][b]);
                }
                if (dist > dt->maxfeatdist)
                    dt->maxfeatdist = dist;
            }
        }
    }

    /* initialize maps */

    for (long p=0; p < img->n; p++){
        ift4DGetVal(dt->cost, p)[0]  = IFT_INFINITY_INT;
        setPixelCoordinate(dt->pred, p, U_NIL);
        //ift4DGetVal(dt->pred, p)[0]  = IFT_NIL;
    }

    return(dt);
}

void iftDestroy4DDynTrees(ift4DDynTrees **dt) {
    ift4DDynTrees *aux = *dt;

    if (aux != NULL) {
        iftDestroyGQueue(&aux->Q);
        iftDestroyImageSequence(&aux->cost);
        iftDestroyImageSequence(&aux->pred);
        iftDestroyImageSequence(&aux->label);
        iftDestroyImageSequence(&aux->root);
        iftDestroyImageSequence(&aux->cumfeat);
        iftFree(aux->treesize);
        iftFree(aux);
        *dt = NULL;
    }
}

void iftExecDiff4DDynTrees(ift4DDynTrees *dt, iftLabeledSet **S,
                           iftSet **M) {
    iftSet *F = iftDiff4DDynTreeRemoval(dt, M); /* remove trees and return their
                       frontier set, which is empty
                       when the marking set M is
                       empty */

    iftImageSequence *img      = dt->img;
    iftAdjRel *A               = dt->A;
    iftImageSequence  *cost    = dt->cost;
    iftImageSequence  *label   = dt->label;
    iftImageSequence  *pred    = dt->pred;
    iftImageSequence  *root    = dt->root;
    iftGQueue *Q               = dt->Q;
    iftImageSequence *cumfeat  = dt->cumfeat;
    int              *treesize = dt->treesize;

    /* Remove seeds from F, if it is the case since seeds have higher
       priority */

    /* if (F != NULL) { /\* non-empty set *\/ */
    /*   iftLabeledSet *seed = *S; */
    /*   while (seed != NULL){ */
    /*     int p = seed->elem; */
    /*     iftRemoveSetElem(&F,p); */
    /*     seed  = seed->next; */
    /*   } */
    /* } */

    /* Reinitialize maps for seed voxels and insert seeds and frontier
       in the queue */

    while (*S != NULL) {
        int lambda;
        long p         = iftRemoveLabeledSet(S,&lambda);
        iftVoxel u     = ift4DGetVoxelCoord(dt->img, p);
        ift4DGetVal(cost, p)[0]  = 0;
        ift4DGetVal(label, p)[0] = lambda;
        setPixelCoordinate(pred, p, U_NIL);
        setPixelCoordinate(root, p, u);
        for (int b=0; b < img->m; b++)
            ift4DGetVal(cumfeat, p)[b] = 0;
        treesize[p] = 0;
        iftUnionSetElem(&F,p);
    }

    while (F != NULL){
        long p = iftRemoveSet(&F);
        int nbuckets = Q->C.nbuckets;
        iftInsertGQueue(&Q,p);
        if (nbuckets != Q->C.nbuckets) {
            printf("HUMM5\n");
        }
    }

    /* Execute the Image Foresting Transform of DDT */

    int k = 0;
    while (!iftEmptyGQueue(Q)){
        long p      = iftRemoveGQueue(Q);
        iftVoxel u = ift4DGetVoxelCoord(img, p);

        if (k++ % 100000 == 0) {
            printf("%d\n", k);
        }
        /* set / update dynamic tree of p */

        long r = getPixelCoordinate(root, p);
        /* if (pred->val[p] == IFT_NIL) { /\* p is a root voxel *\/ */
        /*   for (int b=0; b < img->m; b++) */
        /* 	cumfeat->val[r][b] = img->val[r][b]; */
        /*   treesize[r] = 1; */
        /* } else { */
        for (int b=0; b < img->m; b++)
            ift4DGetVal(cumfeat, r)[b] += ift4DGetVal(img, p)[b]; // ift4DGetVal(img, p)[b];
        treesize[r] += 1;
        /* } */

        /* visit the adjacent voxels for possible path extension */

        for (int i=1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);

            if (ift4DValidVoxel(img, v)){
                long q   = ift4DGetVoxelIndex(img, v);
                if (Q->L.elem[q].color != IFT_BLACK){
                    float arcw = 0;
                    for (int b=0; b < img->m; b++){
                        arcw += iftPowerOfTwo(ift4DGetVal(img, q)[b] -
                                                      ift4DGetVal(cumfeat, r)[b] / treesize[r]);
                    }
                    arcw = CMAX * (arcw / dt->maxfeatdist);
                    int tmp = iftMax(ift4DGetVal(cost, p)[0], iftMin(arcw, CMAX));
                    if (tmp < ift4DGetVal(cost, q)[0])  {
                        if (Q->L.elem[q].color == IFT_GRAY)
                            iftRemoveGQueueElem(Q, q);
                        ift4DGetVal(cost, q)[0]  = tmp;
                        ift4DGetVal(label, q)[0] = ift4DGetVal(label, p)[0];
                        setPixelCoordinateFromArray(root, q, ift4DGetVal(root, p));
                        setPixelCoordinate(pred, q, u);
                        //ift4DGetVal(pred, q)[0]  = p;
                        iftInsertGQueue(&Q,q);
                    } else {
                        if ((ift4DGetVal(label, q)[0] != ift4DGetVal(label, p)[0]) &&
                            (getPixelCoordinate(pred, q) == p) ) {
                            iftDiff4DDynSubtreeUpdate(dt, q);
                        }
                    }
                }
            }
        }
    }

    iftResetGQueue(Q);
}

iftSet *iftDiff4DDynTreeRemoval(ift4DDynTrees *dt, iftSet **M) {
    iftSet /**F = NULL, *R = NULL,*/ *F = NULL, *T = NULL, *Q = NULL;
    iftAdjRel *B = NULL;
    iftBMap *_F = iftCreateBMap(dt->img->n);
    iftBMap *R = iftCreateBMap(dt->img->n);

    iftImageSequence *cost = dt->cost;
    iftImageSequence *root = dt->root;
    iftImageSequence *pred = dt->pred;
    iftImageSequence *img  = dt->img;

    if (iftIs4DImageSequence(dt->img))
        B = iftHyperCuboid(3,3,3,3);
    else if (iftIs3DImageSequence(dt->img))
        B = iftSpheric(1.74);
    else
        B = iftCircular(1.42);

    /* Find the roots of the trees that contain elements in the marking
       set. If the root has not been inserted yet in a root set R, do it
       for its entire marker and reset their cost and predecessor
       information. Set R must contain the roots of all trees marked for
       removal. Set T starts being equal to R and next it stores
       elements from the trees rooted in R. */

    while(*M != NULL) {
        long p = iftRemoveSet(M);
        long r = getPixelCoordinate(root, p);
        /* insert in R and T all roots in the marker of r for tree
           removal */
        iftInsertSet(&Q,r);
        while (Q != NULL) {
            long r = iftRemoveSet(&Q);
            if (ift4DGetVal(cost, r)[0] != IFT_INFINITY_INT){ /* r is not in R and T yet */
                ift4DGetVal(cost, r)[0] = IFT_INFINITY_INT;
                //iftInsertSet(&R,r); iftInsertSet(&T,r);
                //TODO bitset change from int to long
                iftBMapSet1(R, r); iftInsertSet(&T,r);
            }
            iftVoxel u = ift4DGetVoxelCoord(img, r);
            for (int i=1; i < B->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(B,u,i);
                if (ift4DValidVoxel(img, v)){
                    long s = ift4DGetVoxelIndex(img, v);
                    /* s belongs to the same marker of r and it has not been
                       inserted in R and T yet. */
                    if ((getPixelCoordinate(root, s) == s) &&
                        (ift4DGetVal(cost, s)[0] != IFT_INFINITY_INT))
                        iftInsertSet(&Q,s);
                }
            }
        }
    }

    /* Visit all nodes in each tree with root in R, while removing the
       tree. It also identifies the frontier voxels of trees that have
       not been marked for removal. */

    while (T != NULL) {
        long     p = iftRemoveSet(&T);
        iftVoxel u = ift4DGetVoxelCoord(img, p);
        for (int i=1; i < dt->A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(dt->A,u,i);
            if (ift4DValidVoxel(img, v)){
                long q = ift4DGetVoxelIndex(img, v);
                /* If q belongs to a tree for removal, then reset its cost and
                   predecessor information, inserting it in T to pursue tree
                   removal. */
                if (getPixelCoordinate(pred, q) == p){
                    ift4DGetVal(cost, q)[0]=IFT_INFINITY_INT;
                    setPixelCoordinate(pred, q, U_NIL);
                    // -> ift4DGetVal(pred, q)[0]=IFT_NIL;
                    iftInsertSet(&T,q);
                } else {
                    /* If q belongs to a tree that was not marked for removal,
                       then q is a frontier voxel. */
                    //if (!iftSetHasElement(R, dt->root->val[q])){
                    if (!iftBMapValue(R, getPixelCoordinate(root, q))){
                        if (!iftBMapValue(_F, q)) {
                            iftBMapSet1(_F,q);
                            iftInsertSet(&F, q);
                        }
                    }
                }
            }
        }
    }

    iftDestroyBMap(&_F);
    iftDestroyBMap(&R);
    iftDestroyAdjRel(&B);

    return(F);
}

void iftDiff4DDynSubtreeUpdate(ift4DDynTrees *dt, long q) {
    iftSet *T = NULL; /* Set to visit the nodes of the subtree rooted in
               q */
    /* If the subtree contains a single node, its root q, then update
       its cost, label and root information. */
    /* if (dt->Q->L.elem[q].color == IFT_GRAY){  */
    /*   iftRemoveGQueueElem(dt->Q, q); */
    /*   int p      = dt->pred->val[q]; */
    /*   int r      = dt->root->val[p]; */
    /*   float arcw = 0; */
    /*   for (int b=0; b < dt->img->m; b++){ */
    /*     arcw += (dt->img->val[q][b] - dt->cumfeat->val[r][b]/dt->treesize[r])* */
    /* 	      (dt->img->val[q][b] - dt->cumfeat->val[r][b]/dt->treesize[r]); */
    /*   } */
    /*   arcw              = CMAX * (arcw / dt->maxfeatdist);  */
    /*   dt->cost->val[q]  = iftMax(dt->cost->val[p], iftMin(arcw,CMAX)); */
    /*   dt->label->val[q] = dt->label->val[p]; */
    /*   dt->root->val[q]  = dt->root->val[p]; */
    /* } else { */
    /* the subtree contains one or more elements from a previous execution */

    iftInsertSet(&T,q);
    while (T != NULL) {
        long q = iftRemoveSet(&T);

        /* update properties of the previous tree of q */
        /* int r = dt->root->val[q]; */
        /* for (int b=0; b < dt->img->m; b++) */
        /* 	dt->cumfeat->val[r][b] -= dt->img->val[q][b]; */
        /* dt->treesize[r] -= 1; */

        /* update properties of the new tree of q */
        int p = getPixelCoordinate(dt->pred, q);
        int r = getPixelCoordinate(dt->root, p);
        for (int b=0; b < dt->img->m; b++)
            ift4DGetVal(dt->cumfeat, r)[b] += ift4DGetVal(dt->img, q)[b];
        dt->treesize[r] += 1;
        float arcw = 0;
        for (int b=0; b < dt->img->m; b++){
            arcw += iftPowerOfTwo(ift4DGetVal(dt->img, q)[b] -
                                          ift4DGetVal(dt->cumfeat, r)[b] / dt->treesize[r]);
        }
        arcw              = CMAX * (arcw / dt->maxfeatdist);
        ift4DGetVal(dt->cost, q)[0]  = iftMax(ift4DGetVal(dt->cost, p)[0], iftMin(arcw, CMAX));
        ift4DGetVal(dt->label, q)[0] = ift4DGetVal(dt->label, p)[0];
        setPixelCoordinateFromArray(dt->root, q, ift4DGetVal(dt->root, p));
        // -> ift4DGetVal(dt->root, q)[0]  = ift4DGetVal(dt->root, p)[0];

        /* Insert the childree of q in T to pursue the tree update
           process */
        p          = q;
        iftVoxel u = ift4DGetVoxelCoord(dt->img, p);
        for (int i=1; i < dt->A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(dt->A,u,i);
            if (ift4DValidVoxel(dt->img, v)){
                long q = ift4DGetVoxelIndex(dt->img, v);
                if (getPixelCoordinate(dt->pred, q) == p){
                    // -> ift4DGetVal(dt->pred, q)[0] == p){
                    iftInsertSet(&T,q);
                }
            }
        }
    }
/* } */
}

long
iftCountObjectSpelsImageSequence(const iftImageSequence *img, int obj_label) {
    long sum = 0;
    for (int t = 0; t < img->tsize; t++) {
        sum += iftCountObjectSpelsFromTime(img, obj_label, t);
    }
    return sum;
}

iftImageSequence *ift4DBinarize(const iftImageSequence *img) {
    iftImageSequence *bin = iftCreateImageSequenceFromImageSequence(img);
    IMGSEQ_FORLOOP_PARALLEL(img) {
        IMGSEQ_ITE_VAL(bin,0) = IMGSEQ_ITE_VAL(img,0) != 0;
    }
    return bin;
}
