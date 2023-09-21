#include "ift/core/dtypes/UCharArray.h"

#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "iftCommon.h"
#include "iftMemory.h"


iftUCharArray *iftCreateUCharArray(long n) {
    iftUCharArray *iarr = (iftUCharArray*) iftAlloc(1, sizeof(iftUCharArray));
    
    iarr->n = n;
    iarr->val = iftAllocUCharArray(n);
    
    return iarr;
}


void iftDestroyUCharArray(iftUCharArray **iarr) {
    if (iarr != NULL && *iarr != NULL) {
        iftUCharArray *iarr_aux = *iarr;
        
        if (iarr_aux->val != NULL)
            iftFree(iarr_aux->val);
        iftFree(iarr_aux);
        *iarr = NULL;
    }
}


iftUCharArray *iftReadUCharArray(const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);

    iftNumPyHeader *header = NULL;
    void *data = iftReadNumPy(npy_path, &header);
    iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);
    
    if (header->n_dims != 1)
        iftError("Number of dimensions %ld is > 1", "iftReadUCharArray", header->n_dims);

    iftUCharArray *arr = NULL;
    if (cdtype == IFT_UCHAR_TYPE) {
        arr = iftCreateUCharArray(header->size);
        iftFree(arr->val);
        arr->val = data;
    }
    else {
        arr = iftConvertToUCharArray(data, header->size, cdtype);
        iftFree(data);
    }
    
    iftDestroyNumPyHeader(&header);

    return arr;
}


void iftWriteUCharArray(const iftUCharArray *arr, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    long shape[1] = {arr->n};
    iftNumPyHeader *header = iftCreateNumPyHeader(IFT_UCHAR_TYPE, shape, 1);
    iftWriteNumPy(header, arr->val, npy_path);
    iftDestroyNumPyHeader(&header);
}


void iftResizeUCharArray(iftUCharArray **iarr, long n) {
    if (n == 0)
        iftError("The new size must be greater than 0", "iftResizeUCharArray");
    
    if ((*iarr) == NULL) {
        (*iarr) = iftCreateUCharArray(n);
    } else {
        if(n != (*iarr)->n) {
            iftUCharArray *iarr1 = iftCreateUCharArray((*iarr)->n);
            iftCopyData(iarr1->val, (*iarr)->val, (*iarr)->n, sizeof(uchar));
            iftDestroyUCharArray(iarr);
            (*iarr) = iftCreateUCharArray(n);
            iftCopyData((*iarr)->val, iarr1->val, iftMin(iarr1->n, (*iarr)->n), sizeof(uchar));
            iftDestroyUCharArray(&iarr1);
        }
    }
}


iftUCharArray *iftConvertToUCharArray(void *data, size_t length, iftCDataType cdtype) {
    if (data == NULL)
        iftError("Input array is NULL", "iftConvertToUCharArray");
    if (cdtype != IFT_UCHAR_TYPE)
        iftWarning("Array datatype is not UNSIGNED CHAR. A conversion will be done and can lost the origin array values.", "iftConvertToUCharArray");
    
    iftUCharArray *arr = iftCreateUCharArray(length);
    
    if (cdtype == IFT_CHAR_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (uchar) ((char *) data)[i];
    else if (cdtype == IFT_UCHAR_TYPE)
        memcpy(arr->val, data, length);
    else if (cdtype == IFT_SHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (uchar) ((short *) data)[i];
    else if (cdtype == IFT_USHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (uchar) ((ushort *) data)[i];
    else if (cdtype == IFT_INT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (uchar) ((int *) data)[i];
    else if (cdtype == IFT_UINT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (uchar) ((uint *) data)[i];
    else if (cdtype == IFT_LONG_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (uchar) ((long *) data)[i];
    else if (cdtype == IFT_ULONG_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (uchar) ((ulong *) data)[i];
    else if (cdtype == IFT_FLT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (uchar) ((float *) data)[i];
    else if (cdtype == IFT_DBL_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (uchar) ((double *) data)[i];
    else iftError("Invalid datatype during conversion: %s", "iftConvertToUCharArray", iftCDataTypeToString(cdtype));
    
    return arr;
}


