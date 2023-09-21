#include "ift/imgproc/dtypes/BoundingBox.h"

#include "ift/core/io/CSV.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftImage.h"
#include "iftMatrix.h"


iftBoundingBoxArray *iftCreateBoundingBoxArray(long n) {
    iftBoundingBoxArray *bbs = iftAlloc(1, sizeof(iftBoundingBoxArray));
    bbs->val = iftAlloc(n, sizeof(iftBoundingBox));
    bbs->n = n;
    
    return bbs;
}

void iftDestroyBoundingBoxArray(iftBoundingBoxArray **bbs) {
    iftBoundingBoxArray *aux = *bbs;
    
    if (aux != NULL) {
        iftFree(aux->val);
        iftFree(aux);
    }
    
    *bbs = NULL;
}


iftBoundingBoxArray *iftReadBoundingBoxArray(const char *format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);
    
    if (!iftEndsWith(filename, ".csv"))
        iftError("Filename %s is not a CSV", "iftReadBoundingBoxArray", filename);
    
    iftCSV *csv = iftReadCSV(filename, ',');
    
    if (csv->ncols != 6)
        iftError("Number of columns %d in CSV is != 6", "iftReadBoundingBoxArray", csv->ncols);
    
    iftBoundingBoxArray *bbs = iftCreateBoundingBoxArray(csv->nrows);
#pragma omp parallel for
    for (int i = 0; i < csv->nrows; i++) {
        bbs->val[i].begin.x = atoi(csv->data[i][0]);
        bbs->val[i].begin.y = atoi(csv->data[i][1]);
        bbs->val[i].begin.z = atoi(csv->data[i][2]);
        bbs->val[i].end.x = atoi(csv->data[i][3]);
        bbs->val[i].end.y = atoi(csv->data[i][4]);
        bbs->val[i].end.z = atoi(csv->data[i][5]);
    }
    iftDestroyCSV(&csv);
    
    return bbs;
}


void iftWriteBoundingBoxArray(const iftBoundingBoxArray *bbs, const char *format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);
    
    if (!iftEndsWith(filename, ".csv"))
        iftError("Filename %s is not a CSV", "iftWriteBoundingBoxArray", filename);
    
    iftCSV *csv = iftCreateCSV(bbs->n, 6);
#pragma omp parallel for
    for (int i = 0; i < bbs->n; i++) {
        sprintf(csv->data[i][0], "%d", bbs->val[i].begin.x);
        sprintf(csv->data[i][1], "%d", bbs->val[i].begin.y);
        sprintf(csv->data[i][2], "%d", bbs->val[i].begin.z);
        sprintf(csv->data[i][3], "%d", bbs->val[i].end.x);
        sprintf(csv->data[i][4], "%d", bbs->val[i].end.y);
        sprintf(csv->data[i][5], "%d", bbs->val[i].end.z);
    }
    
    iftWriteCSV(csv, filename, ',');
    iftDestroyCSV(&csv);
}


iftVoxel iftBoundingBoxCenterVoxel(iftBoundingBox bb) {
    iftVoxel cv; // central voxel of the min bounding box
    cv.x = (bb.begin.x + bb.end.x + 1) / 2;
    cv.y = (bb.begin.y + bb.end.y + 1) / 2;
    cv.z = (bb.begin.z + bb.end.z + 1) / 2;
    
    return cv;
}


iftBoundingBox iftScaleBoundingBox(iftBoundingBox bb, float alpha) {
    iftMatrix *S = iftScaleMatrix(alpha, alpha, alpha);
    
    
    // centralize the bounding box on the origin (0,0,0)
    iftVoxel disp = iftBoundingBoxCenterVoxel(bb);
    bb.begin.x -= disp.x;
    bb.end.x   -= disp.x;
    bb.begin.y -= disp.y;
    bb.end.y   -= disp.y;
    bb.begin.z -= disp.z;
    bb.end.z   -= disp.z;
    
    iftBoundingBox sbb;
    iftPoint begin = {.x = bb.begin.x, .y = bb.begin.y, .z = bb.begin.z};
    begin          = iftTransformPoint(S, begin);
    
    iftPoint end = {.x = bb.end.x,   .y = bb.end.y,   .z = bb.end.z};
    end          = iftTransformPoint(S, end);
    
    // centralize the scaled bounding box to its original center
    sbb.begin.x = begin.x + disp.x;
    sbb.begin.y = begin.y + disp.y;
    sbb.begin.z = begin.z + disp.z;
    sbb.end.x   = end.x + disp.x;
    sbb.end.y   = end.y + disp.y;
    sbb.end.z   = end.z + disp.z;
    
    iftDestroyMatrix(&S);
    
    if (sbb.begin.x < 0) {
        printf("* Scaled Bounding Box: sbb.begin.x: %d < 0... it is gonna be 0\n", sbb.begin.x);
        sbb.begin.x = 0;
    }
    if (sbb.begin.y < 0) {
        printf("* Scaled Bounding Box: sbb.begin.y: %d < 0... it is gonna be 0\n", sbb.begin.y);
        sbb.begin.y = 0;
    }
    if (sbb.begin.z < 0) {
        printf("* Scaled Bounding Box: sbb.begin.z: %d < 0... it is gonna be 0\n", sbb.begin.z);
        sbb.begin.z = 0;
    }
    
    return sbb;
}


long iftBoundingBoxBoxVolume(iftBoundingBox bb) {
    int bb_xsize = bb.end.x - bb.begin.x + 1;
    int bb_ysize = bb.end.y - bb.begin.y + 1;
    int bb_zsize = bb.end.z - bb.begin.z + 1;
    
    return (bb_xsize * bb_ysize * bb_zsize);
}


void iftPrintBoundingBox(iftBoundingBox bb) {
    printf("begin: (%d, %d, %d)\n", bb.begin.x, bb.begin.y, bb.begin.z);
    printf("end: (%d, %d, %d)\n", bb.end.x, bb.end.y, bb.end.z);
}


iftBoundingBox iftCentralizeBoundingBox(iftBoundingBox bb, iftVoxel new_center) {
    iftVoxel bb_center = iftBoundingBoxCenterVoxel(bb);
    
    iftVoxel disp = iftVectorSub(new_center, bb_center);
    
    iftBoundingBox centralized_bb;
    centralized_bb.begin.x = bb.begin.x + disp.x;
    centralized_bb.begin.y = bb.begin.y + disp.y;
    centralized_bb.begin.z = bb.begin.z + disp.z;
    centralized_bb.end.x = bb.end.x + disp.x;
    centralized_bb.end.y = bb.end.y + disp.y;
    centralized_bb.end.z = bb.end.z + disp.z;
    
    return centralized_bb;
}


iftImageDomain iftBoundingBoxSize(iftBoundingBox bb) {
    iftImageDomain bb_sizes;
    if ((bb.begin.x > bb.end.x) || (bb.begin.y > bb.end.y) || (bb.begin.z > bb.end.z))
        iftError("Invalid Bounding Box: Some Begin coordinate > Ending coordinate", "iftBoundingBoxSize");
    
    bb_sizes.xsize = bb.end.x - bb.begin.x + 1;
    bb_sizes.ysize = bb.end.y - bb.begin.y + 1;
    bb_sizes.zsize = bb.end.z - bb.begin.z + 1;
    
    return bb_sizes;
}



iftBoundingBox iftMinBoundingBoxOfBoundingBoxes(iftBoundingBox bb1, iftBoundingBox bb2) {
    iftBoundingBox mbb;
    mbb.begin.x = iftMin(bb1.begin.x, bb2.begin.x);
    mbb.begin.y = iftMin(bb1.begin.y, bb2.begin.y);
    mbb.begin.z = iftMin(bb1.begin.z, bb2.begin.z);
    mbb.end.x   = iftMax(bb1.end.x, bb2.end.x);
    mbb.end.y   = iftMax(bb1.end.y, bb2.end.y);
    mbb.end.z   = iftMax(bb1.end.z, bb2.end.z);
    
    return mbb;
}


iftBoundingBox iftMinBoundingBoxOfVoxels(const iftVoxel *voxels, size_t n) {
    iftBoundingBox mbb = {.begin = voxels[0], .end = voxels[0]};
    
    for (size_t i = 1; i < n; i++) {
        mbb.begin.x = iftMin(mbb.begin.x, voxels[i].x);
        mbb.begin.y = iftMin(mbb.begin.y, voxels[i].y);
        mbb.begin.z = iftMin(mbb.begin.z, voxels[i].z);
        mbb.end.x   = iftMax(mbb.end.x, voxels[i].x);
        mbb.end.y   = iftMax(mbb.end.y, voxels[i].y);
        mbb.end.z   = iftMax(mbb.end.z, voxels[i].z);
    }
    
    return mbb;
}


iftBoundingBox iftReadBoundingBox(const char *format, ...) {
    va_list args;
    char path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(path, format, args);
    va_end(args);
    
    if (!iftEndsWith(path, ".csv"))
        iftError("Invalid Extension for Bounding Box: %s\nTry *.csv", "iftReadBoundingBox", path);
    
    iftBoundingBox bb;
    
    iftCSV *csv = iftReadCSV(path, ',');
    bb.begin.x = atoi(csv->data[0][0]);
    bb.begin.y = atoi(csv->data[0][1]);
    bb.begin.z = atoi(csv->data[0][2]);
    bb.end.x = atoi(csv->data[0][3]);
    bb.end.y = atoi(csv->data[0][4]);
    bb.end.z = atoi(csv->data[0][5]);
    
    iftDestroyCSV(&csv);
    
    return bb;
}

void iftWriteBoundingBox(iftBoundingBox bb, const char *format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);
    
    if (!iftEndsWith(filename, ".csv"))
        iftError("Invalid Extension for Bounding Box: %s\nTry *.csv", "iftWriteBoundingBox", filename);
    
    iftCSV *csv = iftCreateCSV(1, 6);
    sprintf(csv->data[0][0], "%d", bb.begin.x);
    sprintf(csv->data[0][1], "%d", bb.begin.y);
    sprintf(csv->data[0][2], "%d", bb.begin.z);
    sprintf(csv->data[0][3], "%d", bb.end.x);
    sprintf(csv->data[0][4], "%d", bb.end.y);
    sprintf(csv->data[0][5], "%d", bb.end.z);
    
    iftWriteCSV(csv, filename, ',');
    
    iftDestroyCSV(&csv);
}


iftBoundingBox iftFitBoundingBoxOnImageDomain(iftBoundingBox bb, iftImageDomain domain) {
    iftBoundingBox bb_fit = bb;
    
    iftVoxel upper_left_voxel = bb.begin;
    iftVoxel bottom_right_voxel = bb.end;
    
    if (!iftValidVoxel(&domain, upper_left_voxel)) {
        //iftWarning("Fitting Initial Point (%d, %d, %d) of the Bound. Box onto the image domain (%d, %d, %d)",
        //        "iftFitBoundingBoxOnImageDomain", upper_left_voxel.x, upper_left_voxel.y, upper_left_voxel.z, domain.xsize,
        //        domain.ysize, domain.zsize);

        bb_fit.begin.x = (upper_left_voxel.x >= 0) ? upper_left_voxel.x : 0;
        bb_fit.begin.y = (upper_left_voxel.y >= 0) ? upper_left_voxel.y : 0;
        bb_fit.begin.z = (upper_left_voxel.z >= 0) ? upper_left_voxel.z : 0;
    }
    
    if (!iftValidVoxel(&domain, bottom_right_voxel)) {
        //iftWarning("Fitting Ending Point (%d, %d, %d) of the Bound. Box onto the image domain (%d, %d, %d)",
        //        "iftFitBoundingBoxOnImageDomain", bottom_right_voxel.x, bottom_right_voxel.y, bottom_right_voxel.z, domain.xsize,
        //        domain.ysize, domain.zsize);
        
        bb_fit.end.x = (bottom_right_voxel.x >= domain.xsize) ? domain.xsize - 1 : bottom_right_voxel.x;
        bb_fit.end.y = (bottom_right_voxel.y >= domain.ysize) ? domain.ysize - 1 : bottom_right_voxel.y;
        bb_fit.end.z = (bottom_right_voxel.z >= domain.zsize) ? domain.zsize - 1 : bottom_right_voxel.z;
    }
    
    return bb_fit;
}







