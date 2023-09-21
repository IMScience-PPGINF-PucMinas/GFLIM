#include "ift/imgproc/dtypes/Roi.h"

iftRoi *iftCreateRoi(long n)
{
    if(n <= 0)
        iftError("The number of voxels in the ROI must be >= 0 (n=%lu)", "iftCreateRoi", n);

    iftRoi *roi = iftAlloc(1, sizeof(iftRoi));

    roi->val = iftAlloc(n, sizeof(iftVoxel));
    roi->n = n;
    roi->label = IFT_NIL;

    return roi;
}


void iftDestroyRoi(iftRoi **roi)
{
    iftRoi *aux = *roi;

    if (aux != NULL) {
        iftFree(aux->val);
        iftFree(aux);
        *roi = NULL;
    }
}

iftRoi *iftCopyRoi(iftRoi *roi)
{
    iftRoi *copy = iftCreateRoi(roi->n);

    for(int i = 0; i < roi->n; i++)
        copy->val[i] = roi->val[i];

    return copy;
}

iftRoiArray *iftCreateRoiArray(long n)
{
    if(n <= 0)
        iftError("The number of ROIs in the array must be >= 0 (n=%lu)", "iftCreateRoiArray", n);

    iftRoiArray *arr = iftAlloc(1, sizeof(iftRoiArray));

    arr->val = iftAlloc(n, sizeof(iftRoi*));
    arr->n = n;

    return arr;
}

void iftDestroyRoiArray(iftRoiArray **arr)
{
    iftRoiArray *aux = *arr;

    if (aux != NULL) {
        for(int i = 0; i < aux->n; i++) {
            if(aux->val[i] != NULL)
                iftDestroyRoi(&aux->val[i]);
        }
        iftFree(aux->val);
        iftFree(aux);
        *arr = NULL;
    }
}

iftRoiArray *iftCopyRoiArray(iftRoiArray *arr)
{
    iftRoiArray *copy = iftCreateRoiArray(arr->n);

    for(int i = 0; i < arr->n; i++)
        copy->val[i] = iftCopyRoi(arr->val[i]);

    return copy;
}

void iftWriteRoiArray(iftRoiArray *arr, char *filename)
{
    FILE *fp = fopen(filename,"wb");

    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteRoiArray", filename);

    fprintf(fp, "ROI_ARRAY\n");
    fprintf(fp, "%lu\n", arr->n);
    for(int i = 0; i < arr->n; i++) {
        iftRoi *roi = arr->val[i];
        fprintf(fp, "%lu\n", roi->n);
        for(int v = 0; v < roi->n; v++)
            fprintf(fp, "%d %d %d %d ", roi->val[v].x, roi->val[v].y, roi->val[v].z, roi->label);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

iftRoiArray *iftReadRoiArray(char *filename)
{
    FILE *fp = fopen(filename, "r");

    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftReadRoiArray", filename);

    char type[20];
    if (fscanf(fp,"%s\n",type) != 1)
        iftError("Reading error", "iftReadRoiArray");

    iftRoiArray *arr = NULL;
    if((strcmp(type,"ROI_ARRAY") == 0)){
        int n, n1, x, y, z, label;
        if(fscanf(fp,"%d\n", &n) != 1)
            iftError("Reading error", "iftReadRoiArray");
        arr = iftCreateRoiArray(n);
        for(int i = 0; i < arr->n; i++) {
            if(fscanf(fp,"%d\n", &n1) != 1)
                iftError("Reading error", "iftReadRoiArray");
            arr->val[i] = iftCreateRoi(n1);
            iftRoi *roi = arr->val[i];
            for(int v = 0; v < roi->n; v++) {
                if(fscanf(fp,"%d %d %d %d ", &x, &y, &z, &label) != 4)
                    iftError("Reading error", "iftReadRoiArray");
                roi->val[v].x = x;
                roi->val[v].y = y;
                roi->val[v].z = z;
            }
        }
    } else {
        iftError("Input file must be a ROI array", "iftReadRoiArray");
    }
    fclose(fp);

    return arr;
}

int iftBiggestRoiInArray(iftRoiArray *arr)
{
    int maxVal = IFT_INFINITY_INT_NEG, maxIdx = IFT_NIL;
    for(int i = 0; i < arr->n; i++) {
        if(arr->val[i]->n > maxVal) {
            maxIdx = i;
            maxVal = arr->val[i]->n;
        }
    }

    return maxIdx;
}

iftRoi *iftRoiFromBoundingBox(iftBoundingBox bb)
{
    int n = (bb.end.x - bb.begin.x + 1)*(bb.end.y - bb.begin.y + 1)*(bb.end.z - bb.begin.z + 1);
    iftRoi *roi = iftCreateRoi(n);

    iftVoxel v;
    int s = 0;
    for (v.z = bb.begin.z; v.z <= bb.end.z; v.z++) {
        for (v.y = bb.begin.y; v.y <= bb.end.y; v.y++) {
            for (v.x = bb.begin.x; v.x <= bb.end.x; v.x++) {
                roi->val[s] = v;
                s++;
            }
        }
    }

    return roi;
}

iftImage* iftExtractRoiNoBkgd(iftImage *img, iftRoi *roi)
{
    iftSize maskSize = {.x = img->xsize, .y = img->ysize, .z = img->zsize};
    iftImage *mask = iftMaskFromRoi(roi, maskSize);
    iftImage *obj = iftMask(img, mask);
    iftBoundingBox bb = iftMinBoundingBox(obj, NULL);
    iftImage* imgRoi = iftExtractROI(obj, bb);

    iftDestroyImage(&mask);
    iftDestroyImage(&obj);

    return imgRoi;
}

iftMImage* iftMExtractRoiNoBkgd(iftMImage *mimg, iftRoi *roi)
{
    iftSize maskSize = {.x = mimg->xsize, .y = mimg->ysize, .z = mimg->zsize};
    iftImage *mask = iftMaskFromRoi(roi, maskSize);
    iftBoundingBox bb = iftMinBoundingBox(mask, NULL);
    iftMImage *mimgRoi = iftCreateMImage(bb.end.x-bb.begin.x+1, bb.end.y-bb.begin.y+1, bb.end.z-bb.begin.z+1, mimg->m);

    #pragma omp parallel for
    for(int i = 0; i < roi->n; i++) {
        for(int b = 0; b < mimg->m; b++) {
            int p = iftMGetVoxelIndex(mimg, roi->val[i]);
            mimgRoi->val[i][b] = mimg->val[p][b];
        }
    }

    iftDestroyImage(&mask);

    return mimgRoi;
}

iftImage *iftMaskFromRoi(iftRoi *roi, iftSize maskSize)
{
    iftImage *mask = iftCreateImage(maskSize.x, maskSize.y, maskSize.z);

    for(int i = 0; i < roi->n; i++)
        iftImgVoxelVal(mask, roi->val[i]) = 1;

    return mask;
}

iftVoxel iftRoiCenterVoxel(iftImage *img, iftRoi *roi)
{
    iftSize maskSize = {.x = img->xsize, .y = img->ysize, .z = img->zsize};
    iftImage *mask = iftMaskFromRoi(roi, maskSize);
    iftVoxel gc = iftGeometricCenterVoxel(mask);
    iftDestroyImage(&mask);
    
    return gc;
}

iftRoi *iftRoiFromSuperpixel(iftImage *spixLabels, int label)
{
    int spixSize = iftGetSuperpixelSize(spixLabels, label);
    iftRoi *roi = iftCreateRoi(spixSize);

    int s = 0;
    for(int p = 0; p < spixLabels->n; p++) {
        if(spixLabels->val[p] == label) {
            roi->val[s] = iftGetVoxelCoord(spixLabels, p);
            s++;
        }
    }

    return roi;
}
