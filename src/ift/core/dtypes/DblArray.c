#include "ift/core/dtypes/DblArray.h"

#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "iftCommon.h"
#include "iftMemory.h"


iftDblArray *iftCreateDblArray(long n) {
    iftDblArray *darr = (iftDblArray *) iftAlloc(1, sizeof(iftDblArray));
    
    darr->n   = n;
    darr->val = iftAllocDoubleArray(n);
    
    return darr;
}

iftDblArray *iftCopyDblArray(const double* array, long n) {
    iftDblArray* out = iftCreateDblArray(n);
    
    memcpy(out->val, array, n*sizeof(double));
    
    return out;
}


void iftDestroyDblArray(iftDblArray **darr) {
    
    if(darr == NULL){
        return;
    }
    if(*darr == NULL){
        return;
    }
    
    
    if (darr != NULL && *darr != NULL) {
        iftDblArray *darr_aux = *darr;
        
        if (darr_aux->val != NULL)
            iftFree(darr_aux->val);
        iftFree(darr_aux);
        *darr = NULL;
        
    }
}


iftDblArray *iftReadDblArray(const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    iftNumPyHeader *header = NULL;
    void *data = iftReadNumPy(npy_path, &header);
    iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);
    
    if (header->n_dims != 1)
        iftError("Number of dimensions %ld is > 1", "iftReadDblArray", header->n_dims);
    
    iftDblArray *arr = NULL;
    if (cdtype == IFT_DBL_TYPE) {
        arr = iftCreateDblArray(header->size);
        iftFree(arr->val);
        arr->val = data;
    }
    else {
        arr = iftConvertToDblArray(data, header->size, cdtype);
        iftFree(data);
    }
    
    iftDestroyNumPyHeader(&header);
    
    return arr;
}


void iftWriteDblArray(const iftDblArray *arr, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    long shape[1] = {arr->n};
    iftNumPyHeader *header = iftCreateNumPyHeader(IFT_DBL_TYPE, shape, 1);
    iftWriteNumPy(header, arr->val, npy_path);
    iftDestroyNumPyHeader(&header);
}


void iftResizeDblArray(iftDblArray **iarr, long n) {
    if (n == 0)
        iftError("The new size must be greater than 0", "iftResizeDblArray");
    
    if ((*iarr) == NULL) {
        (*iarr) = iftCreateDblArray(n);
    } else {
        if(n != (*iarr)->n) {
            iftDblArray *iarr1 = iftCreateDblArray((*iarr)->n);
            iftCopyData(iarr1->val, (*iarr)->val, (*iarr)->n, sizeof(double));
            iftDestroyDblArray(iarr);
            (*iarr) = iftCreateDblArray(n);
            iftCopyData((*iarr)->val, iarr1->val, iftMin(iarr1->n, (*iarr)->n), sizeof(double));
            iftDestroyDblArray(&iarr1);
        }
    }
}

char *iftDblArrayAsString(const iftDblArray *darr) {
    // 2 chars for square brackets + 40 chars per number (at most) + n-1 commas
    char *str = iftAlloc(2 + (40 * darr->n) + (darr->n-1), sizeof(char));
    
    sprintf(str, "[");
    for (long i = 0; i <= (darr->n-2); i++)
        sprintf(str, "%s%.10lf, ", str, darr->val[i]);
    sprintf(str, "%s%.10lf]", str, darr->val[darr->n-1]);
    
    return str;
}


iftDblArray *iftConvertToDblArray(void *data, size_t length, iftCDataType cdtype) {
    if (data == NULL)
        iftError("Input array is NULL", "iftConvertToDblArray");
    if (cdtype != IFT_DBL_TYPE)
        iftWarning("Array datatype is not DOUBLE. A conversion will be done and can lost the origin array values.", "iftConvertToDblArray");
    
    iftDblArray *arr = iftCreateDblArray(length);
    
    if (cdtype == IFT_CHAR_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (double) ((char *) data)[i];
    else if (cdtype == IFT_UCHAR_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (double) ((uchar *) data)[i];
    else if (cdtype == IFT_SHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (double) ((short *) data)[i];
    else if (cdtype == IFT_USHORT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (double) ((ushort *) data)[i];
    else if (cdtype == IFT_INT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (double) ((int *) data)[i];
    else if (cdtype == IFT_UINT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (double) ((uint *) data)[i];
    else if (cdtype == IFT_LONG_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (double) ((long *) data)[i];
    else if (cdtype == IFT_ULONG_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (double) ((ulong *) data)[i];
    else if (cdtype == IFT_FLT_TYPE)
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++)
            arr->val[i] = (double) ((float *) data)[i];
    else if (cdtype == IFT_DBL_TYPE)
        memcpy(arr->val, data, length);
    else iftError("Invalid datatype during conversion: %s", "iftConvertToDblArray", iftCDataTypeToString(cdtype));
    
    return arr;
}

