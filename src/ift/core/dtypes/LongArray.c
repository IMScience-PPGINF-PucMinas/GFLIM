#include "ift/core/dtypes/LongArray.h"

#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "iftCommon.h"
#include "iftMemory.h"



iftLongArray *iftCreateLongArray(long n) {
    iftLongArray *iarr = (iftLongArray*) iftAlloc(1, sizeof(iftLongArray));
    
    iarr->n   = n;
    iarr->val = iftAllocLongIntArray(n);
    
    return iarr;
}


void iftDestroyLongArray(iftLongArray **iarr) {
    if (iarr != NULL && *iarr != NULL) {
        iftLongArray *iarr_aux = *iarr;
        
        if (iarr_aux->val != NULL)
            iftFree(iarr_aux->val);
        iftFree(iarr_aux);
        *iarr = NULL;
    }
}


iftLongArray *iftReadLongArray(const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    iftNumPyHeader *header = NULL;
    void *data = iftReadNumPy(npy_path, &header);
    iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);
    
    if (header->n_dims != 1)
        iftError("Number of dimensions %ld is > 1", "iftReadLongArray", header->n_dims);
    
    iftLongArray *arr = NULL;
    if (cdtype == IFT_LONG_TYPE) {
        arr = iftCreateLongArray(header->size);
        iftFree(arr->val);
        arr->val = data;
    }
    else {
        arr = iftConvertToLongArray(data, header->size, cdtype);
        iftFree(data);
    }
    
    iftDestroyNumPyHeader(&header);
    
    return arr;
}


void iftWriteLongArray(const iftLongArray *arr, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    long shape[1] = {arr->n};
    iftNumPyHeader *header = iftCreateNumPyHeader(IFT_LONG_TYPE, shape, 1);
    iftWriteNumPy(header, arr->val, npy_path);
    iftDestroyNumPyHeader(&header);
}


void iftResizeLongArray(iftLongArray **iarr, long n) {
    if (n == 0)
        iftError("The new size must be greater than 0", "iftResizeLongArray");
    
    if ((*iarr) == NULL) {
        (*iarr) = iftCreateLongArray(n);
    } else {
        if(n != (*iarr)->n) {
            iftLongArray *iarr1 = iftCreateLongArray((*iarr)->n);
            iftCopyData(iarr1->val, (*iarr)->val, (*iarr)->n, sizeof(long));
            iftDestroyLongArray(iarr);
            (*iarr) = iftCreateLongArray(n);
            iftCopyData((*iarr)->val, iarr1->val, iftMin(iarr1->n, (*iarr)->n), sizeof(long));
            iftDestroyLongArray(&iarr1);
        }
    }
}


iftLongArray *iftConvertToLongArray(void *data, size_t length, iftCDataType cdtype) {
    if (data == NULL)
        iftError("Input array is NULL", "iftConvertToLongArray");
    if (cdtype != IFT_LONG_TYPE)
        iftWarning("Array datatype is not LONG INTEGER. A conversion will be done and can lost the origin array values.", "iftConvertToLongArray");
    
    iftLongArray *arr = iftCreateLongArray(length);
    
    if (cdtype == IFT_CHAR_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (long) ((char *) data)[i];
    else if (cdtype == IFT_UCHAR_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (long) ((uchar *) data)[i];
    else if (cdtype == IFT_SHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (long) ((short *) data)[i];
    else if (cdtype == IFT_USHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (long) ((ushort *) data)[i];
    else if (cdtype == IFT_INT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (long) ((int *) data)[i];
    else if (cdtype == IFT_UINT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (long) ((uint *) data)[i];
    else if (cdtype == IFT_LONG_TYPE)
        memcpy(arr->val, data, length);
    else if (cdtype == IFT_ULONG_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (long) ((ulong *) data)[i];
    else if (cdtype == IFT_FLT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (long) ((float *) data)[i];
    else if (cdtype == IFT_DBL_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (long) ((double *) data)[i];
    else iftError("Invalid datatype during conversion: %s", "iftConvertToLongArray", iftCDataTypeToString(cdtype));
    
    return arr;
}





