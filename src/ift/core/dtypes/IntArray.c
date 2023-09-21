#include "ift/core/dtypes/IntArray.h"

#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftCommon.h"


iftIntArray *iftCreateIntArray(long n) {
    iftIntArray *iarr = (iftIntArray*) iftAlloc(1, sizeof(iftIntArray));
    
    iarr->n = n;
    iarr->val = iftAllocIntArray(n);
    
    return iarr;
}


void iftDestroyIntArray(iftIntArray **iarr) {
    if (iarr != NULL && *iarr != NULL) {
        iftIntArray *iarr_aux = *iarr;
        
        if (iarr_aux->val != NULL)
            iftFree(iarr_aux->val);
        iftFree(iarr_aux);
        *iarr = NULL;
    }
}


iftIntArray *iftReadIntArray(const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);

    iftNumPyHeader *header = NULL;
    void *data = iftReadNumPy(npy_path, &header);
    iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);

    if (header->n_dims != 1)
        iftError("Number of dimensions %ld is > 1", "iftReadIntArray", header->n_dims);

    iftIntArray *arr = NULL;
    if (cdtype == IFT_INT_TYPE) {
        arr = iftCreateIntArray(header->size);
        iftFree(arr->val);
        arr->val = data;
    }
    else {
        arr = iftConvertToIntArray(data, header->size, cdtype);
        iftFree(data);
    }
    
    iftDestroyNumPyHeader(&header);

    return arr;
}


void iftWriteIntArray(const iftIntArray *arr, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    long shape[1] = {arr->n};
    iftNumPyHeader *header = iftCreateNumPyHeader(IFT_INT_TYPE, shape, 1);
    iftWriteNumPy(header, arr->val, npy_path);
    iftDestroyNumPyHeader(&header);
}


void iftResizeIntArray(iftIntArray **iarr, long n) {
    if (n == 0)
        iftError("The new size must be greater than 0", "iftResizeIntArray");
    
    if ((*iarr) == NULL) {
        (*iarr) = iftCreateIntArray(n);
    } else {
        if(n != (*iarr)->n) {
            iftIntArray *iarr1 = iftCreateIntArray((*iarr)->n);
            iftCopyData(iarr1->val, (*iarr)->val, (*iarr)->n, sizeof(int));
            iftDestroyIntArray(iarr);
            (*iarr) = iftCreateIntArray(n);
            iftCopyData((*iarr)->val, iarr1->val, iftMin(iarr1->n, (*iarr)->n), sizeof(int));
            iftDestroyIntArray(&iarr1);
        }
    }
}


void iftShuffleIntArray(int* array, int n) {
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    int j;
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        j = rand() % (i+1);
        // Swap arr[i] with the element at random index
        iftSwap(array[i],array[j]);
    }
}


iftIntArray *iftVoxelToIntArray(iftVoxel v) {
    iftIntArray *arr = iftCreateIntArray(3);
    arr->val[0] = v.x;
    arr->val[1] = v.y;
    arr->val[2] = v.z;
    
    return arr;
}


iftVoxel iftIntArrayToVoxel(const iftIntArray *arr) {
    if (arr->n != 3)
        iftError("Array should have size 3, instead of %d", "iftIntArrayToVoxel", arr->n);
    
    iftVoxel v = {.x = arr->val[0], .y = arr->val[1], .z = arr->val[2]};
    return v;
}


char *iftIntArrayAsString(const iftIntArray *iarr) {
    // 2 chars for square brackets + 10 chars per number (at most) + n-1 commas
    char *str = iftAlloc(2 + (10 * iarr->n) + (iarr->n-1), sizeof(char));
    
    sprintf(str, "[");
    for (long i = 0; i <= (long) (iarr->n-2); i++) {
        sprintf(str, "%s%d, ", str, iarr->val[i]);
    }
    sprintf(str, "%s%d]", str, iarr->val[iarr->n-1]);
    
    return str;
}


iftIntArray *iftConvertToIntArray(void *data, size_t length, iftCDataType cdtype) {
    if (data == NULL)
        iftError("Input array is NULL", "iftConvertToIntArray");
    if (cdtype != IFT_INT_TYPE)
        iftWarning("Array datatype is not INTEGER. A conversion will be done and can lost the origin array values.", "iftConvertToIntArray");
    
    iftIntArray *arr = iftCreateIntArray(length);
    
    if (cdtype == IFT_CHAR_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (int) ((char *) data)[i];
    else if (cdtype == IFT_UCHAR_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (int) ((uchar *) data)[i];
    else if (cdtype == IFT_SHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (int) ((short *) data)[i];
    else if (cdtype == IFT_USHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (int) ((ushort *) data)[i];
    else if (cdtype == IFT_INT_TYPE)
        memcpy(arr->val, data, length);
    else if (cdtype == IFT_UINT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (int) ((uint *) data)[i];
    else if (cdtype == IFT_LONG_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (int) ((long *) data)[i];
    else if (cdtype == IFT_ULONG_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (int) ((ulong *) data)[i];
    else if (cdtype == IFT_FLT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (int) ((float *) data)[i];
    else if (cdtype == IFT_DBL_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (int) ((double *) data)[i];
    else iftError("Invalid datatype during conversion: %s", "iftConvertToIntArray", iftCDataTypeToString(cdtype));
    
    return arr;
}


