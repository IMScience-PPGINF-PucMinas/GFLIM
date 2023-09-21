#include "ift/core/dtypes/FloatArray.h"

#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "iftCommon.h"


iftFloatArray *iftCreateFloatArray(long n) {
    iftFloatArray *darr = (iftFloatArray *) iftAlloc(1, sizeof(iftFloatArray));
    
    darr->n   = n;
    darr->val = iftAllocFloatArray(n);
    
    return darr;
}


void iftDestroyFloatArray(iftFloatArray **darr) {
    if (darr != NULL && *darr != NULL) {
        iftFloatArray *darr_aux = *darr;
        
        if (darr_aux->val != NULL)
            iftFree(darr_aux->val);
        iftFree(darr_aux);
        *darr = NULL;
    }
}


iftFloatArray *iftReadFloatArray(const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    iftNumPyHeader *header = NULL;
    void *data = iftReadNumPy(npy_path, &header);
    iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);
    
    if (header->n_dims != 1)
        iftError("Number of dimensions %ld is > 1", "iftReadFloatArray", header->n_dims);
    
    iftFloatArray *arr = NULL;
    if (cdtype == IFT_FLT_TYPE) {
        arr = iftCreateFloatArray(header->size);
        iftFree(arr->val);
        arr->val = data;
    }
    else {
        arr = iftConvertToFloatArray(data, header->size, cdtype);
        iftFree(data);
    }
    
    iftDestroyNumPyHeader(&header);
    
    return arr;
}


void iftWriteFloatArray(const iftFloatArray *arr, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    long shape[1] = {arr->n};
    iftNumPyHeader *header = iftCreateNumPyHeader(IFT_FLT_TYPE, shape, 1);
    iftWriteNumPy(header, arr->val, npy_path);
    iftDestroyNumPyHeader(&header);
}


void iftResizeFloatArray(iftFloatArray **iarr, long n) {
    if (n == 0)
        iftError("The new size must be greater than 0", "iftResizeFloatArray");
    
    if ((*iarr) == NULL) {
        (*iarr) = iftCreateFloatArray(n);
    } else {
        if(n != (*iarr)->n) {
            iftFloatArray *iarr1 = iftCreateFloatArray((*iarr)->n);
            iftCopyData(iarr1->val, (*iarr)->val, (*iarr)->n, sizeof(float));
            iftDestroyFloatArray(iarr);
            (*iarr) = iftCreateFloatArray(n);
            iftCopyData((*iarr)->val, iarr1->val, iftMin(iarr1->n, (*iarr)->n), sizeof(float));
            iftDestroyFloatArray(&iarr1);
        }
    }
}


char *iftFloatArrayAsString(const iftFloatArray *farr) {
    // 2 chars for square brackets + 40 chars per number (at most) + n-1 commas
    char *str = iftAlloc(2 + (40 * farr->n) + (farr->n-1), sizeof(char));
    
    sprintf(str, "[");
    for (long i = 0; i <= (farr->n-2); i++)
        sprintf(str, "%s%f, ", str, farr->val[i]);
    sprintf(str, "%s%f]", str, farr->val[farr->n-1]);
    
    return str;
}


iftFloatArray *iftConvertToFloatArray(void *data, size_t length, iftCDataType cdtype) {
    if (data == NULL)
        iftError("Input array is NULL", "iftConvertToFloatArray");
    if (cdtype != IFT_FLT_TYPE)
        iftWarning("Array datatype is not FLOAT. A conversion will be done and can lost the origin array values.", "iftConvertToFloatArray");
    
    iftFloatArray *arr = iftCreateFloatArray(length);
    
    if (cdtype == IFT_CHAR_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (float) ((char *) data)[i];
    else if (cdtype == IFT_UCHAR_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (float) ((uchar *) data)[i];
    else if (cdtype == IFT_SHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (float) ((short *) data)[i];
    else if (cdtype == IFT_USHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (float) ((ushort *) data)[i];
    else if (cdtype == IFT_INT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (float) ((int *) data)[i];
    else if (cdtype == IFT_UINT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (float) ((uint *) data)[i];
    else if (cdtype == IFT_LONG_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (float) ((long *) data)[i];
    else if (cdtype == IFT_ULONG_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (float) ((ulong *) data)[i];
    else if (cdtype == IFT_FLT_TYPE)
        memcpy(arr->val, data, length);
    else if (cdtype == IFT_DBL_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (float) ((double *) data)[i];
    else iftError("Invalid datatype during conversion: %s", "iftConvertToFloatArray", iftCDataTypeToString(cdtype));
    
    return arr;
}

iftDblArray *iftFloatArrayToDblArray(iftFloatArray *array_src) {
    iftDblArray *array_dst = iftCreateDblArray(array_src->n);
    iftCopyFloatArrayToDblArray(array_dst->val, array_src->val, array_src->n);

    return array_dst;
}

iftFloatArray *iftDblArrayToFloatArray(iftDblArray *array_src) {
    iftFloatArray *array_dst = iftCreateFloatArray(array_src->n);
    iftCopyDblArrayToFloatArray(array_dst->val, array_src->val, array_src->n);

    return array_dst;
}
