#include "ift/imgproc/dtypes/VoxelArray.h"
#include "iftCommon.h"

iftVoxelArray *iftCreateVoxelArray(long n)
{
    iftVoxelArray *varr = iftAlloc(1, sizeof(iftVoxelArray));

    varr->val = iftAlloc(n, sizeof(iftVoxel));
    varr->n = n;

    return varr;
}


void iftDestroyVoxelArray(iftVoxelArray **varr)
{
    iftVoxelArray *aux = *varr;

    if (aux != NULL) {
        iftFree(aux->val);
        iftFree(aux);
        *varr = NULL;
    }
}


int iftVoxelArrayFurthestPair(const iftVoxelArray *a, const iftVoxelArray *b)
{
    if (a->n != b->n)
        iftError("arrays must have the same length", "iftVoxelArrayFurthestPair");

    int dist = 0, index = -1;
    for (int i = 0; i < a->n; i++) {
        int current = iftSquaredVoxelDistance(a->val[i], b->val[i]);
        if (current > dist) {
            index = i;
            dist = current;
        }
    }

    return index;
}


void iftInsertVoxel(iftVoxelArray *arr, int index, iftVoxel *u)
{
    arr->n += 1;
    arr->val = iftRealloc(arr->val, arr->n * (sizeof *arr->val));
    for (int i = arr->n - 1; i > index; i--) {
        iftCopyVoxel(&arr->val[i - 1], &arr->val[i]);
    }
    iftCopyVoxel(u, &arr->val[index]);
}
