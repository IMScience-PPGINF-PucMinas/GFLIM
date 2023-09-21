#include "ift/core/dtypes/CharArray.h"

#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftMemory.h"


iftCharArray *iftCreateCharArray(long n) {
    iftCharArray *carr = (iftCharArray *) iftAlloc(1, sizeof(iftCharArray));
    carr->n = n;
    carr->val = iftAlloc(n, sizeof(char));
    
    return carr;
}


void iftDestroyCharArray(iftCharArray **carr) {
    if (carr) {
        iftCharArray *aux = *carr;
        iftFree(aux->val);
        iftFree(aux);
        *carr = NULL;
    }
}


iftCharArray *iftReadCharArray(const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    iftNumPyHeader *header = NULL;
    void *data = iftReadNumPy(npy_path, &header);
    iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);
    
    if (header->n_dims != 1)
        iftError("Number of dimensions %ld is > 1", "iftReadCharArray", header->n_dims);
    
    iftCharArray *arr = NULL;
    if (cdtype == IFT_CHAR_TYPE) {
        arr = iftCreateCharArray(header->size);
        iftFree(arr->val);
        arr->val = data;
    }
    else {
        arr = iftConvertToCharArray(data, header->size, cdtype);
        iftFree(data);
    }
    
    iftDestroyNumPyHeader(&header);
    
    return arr;
}


void iftWriteCharArray(const iftCharArray *arr, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);

    long shape[1] = {arr->n};
    iftNumPyHeader *header = iftCreateNumPyHeader(IFT_CHAR_TYPE, shape, 1);
    iftWriteNumPy(header, arr->val, npy_path);
    iftDestroyNumPyHeader(&header);
}


iftCharArray *iftConvertToCharArray(void *data, size_t length, iftCDataType cdtype) {
    if (data == NULL)
        iftError("Input array is NULL", "iftConvertToCharArray");
    if (cdtype != IFT_CHAR_TYPE)
        iftWarning("Array datatype is not CHAR. A conversion will be done and can lost the origin array values.", "iftConvertToCharArray");
    
    iftCharArray *arr = iftCreateCharArray(length);
    
    if (cdtype == IFT_CHAR_TYPE)
        memcpy(arr->val, data, (size_t) length);
    else if (cdtype == IFT_UCHAR_TYPE)
        #pragma omp parallel for
        for (long i = 0; i < length; i++)
            arr->val[i] = (char) ((uchar *) data)[i];
    else if (cdtype == IFT_SHORT_TYPE)
        #pragma omp parallel for
        for (long i = 0; i < length; i++)
            arr->val[i] = (char) ((short *) data)[i];
    else if (cdtype == IFT_USHORT_TYPE)
        #pragma omp parallel for
        for (long i = 0; i < length; i++)
            arr->val[i] = (char) ((ushort *) data)[i];
    else if (cdtype == IFT_INT_TYPE)
        #pragma omp parallel for
        for (long i = 0; i < length; i++)
            arr->val[i] = (char) ((int *) data)[i];
    else if (cdtype == IFT_UINT_TYPE)
        #pragma omp parallel for
        for (long i = 0; i < length; i++)
            arr->val[i] = (char) ((uint *) data)[i];
    else if (cdtype == IFT_LONG_TYPE)
        #pragma omp parallel for
        for (long i = 0; i < length; i++)
            arr->val[i] = (char) ((long *) data)[i];
    else if (cdtype == IFT_ULONG_TYPE)
        #pragma omp parallel for
        for (long i = 0; i < length; i++)
            arr->val[i] = (char) ((ulong *) data)[i];
    else if (cdtype == IFT_FLT_TYPE)
        #pragma omp parallel for
        for (long i = 0; i < length; i++)
            arr->val[i] = (char) ((float *) data)[i];
    else if (cdtype == IFT_DBL_TYPE)
        #pragma omp parallel for
        for (long i = 0; i < length; i++)
            arr->val[i] = (char) ((double *) data)[i];
    else iftError("Invalid datatype during conversion: %s", "iftConvertToCharArray", iftCDataTypeToString(cdtype));
    
    
    return arr;
}
