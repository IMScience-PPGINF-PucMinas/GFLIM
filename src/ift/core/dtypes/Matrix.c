#include "ift/core/dtypes/Matrix.h"

#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "iftMatrix.h"
#include "iftMemory.h"



iftMatrix *iftCreateMatrix(int ncols, int nrows) {
    iftMatrix *M = (iftMatrix *) iftAlloc(1, sizeof(iftMatrix));
    
    M->ncols = ncols;
    M->nrows = nrows;
    M->tbrow = (long*)iftAlloc(nrows,sizeof(long));
    
    for (long r = 0; r < (long)nrows; r++) {
        M->tbrow[r] = r * ncols;
    }
    M->n = (long) ncols * (long) nrows;
    M->allocated = true;
    
    M->val = iftAllocFloatArray(M->n);
    
    return (M);
}

iftMatrix *iftCreateMatrix_omp(int ncols, int nrows) {
    iftMatrix *M = (iftMatrix *) iftAlloc(1, sizeof(iftMatrix));
    
    M->ncols = ncols;
    M->nrows = nrows;
    M->tbrow = (long*)iftAlloc(nrows,sizeof(long));
    
    #pragma omp parallel for shared(M, nrows, ncols)
    for (long r = 0; r < (long)nrows; r++) {
        M->tbrow[r] = r * ncols;
    }
    
    M->n = (long) ncols * (long) nrows;
    M->allocated = true;
    
    M->val = iftAllocFloatArray(M->n);
    
    return (M);
}


void iftDestroyMatrix(iftMatrix **M) {
    if (M != NULL) {
        iftMatrix *aux = *M;
        
        if (aux != NULL) {
            if (aux->allocated && aux->val != NULL)
                iftFree(aux->val);
            iftFree(aux->tbrow);
            iftFree(aux);
        }
        *M = NULL;
    }
}


iftMatrix *iftReadMatrix(const char *format, ...) {
        va_list args;
        char npy_path[IFT_STR_DEFAULT_SIZE];
        
        va_start(args, format);
        vsprintf(npy_path, format, args);
        va_end(args);
        
        iftNumPyHeader *header = NULL;
        void *data = iftReadNumPy(npy_path, &header);
        iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);
        
        if (header->n_dims != 2)
            iftError("Number of dimensions %ld is != 2", "iftReadMatrix", header->n_dims);
        
        long nrows = header->shape[0];
        long ncols = header->shape[1];
    
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
    
        iftMatrix *mat = iftCreateMatrix(ncols, nrows);
        
        // fortran order: column-major - the matrix is stored inverted as (ncols, nrows)
        if (header->fortran_order) {
            iftMatrix *fortran_mat = iftCreateMatrix(nrows, ncols);
            iftFree(fortran_mat->val);
            fortran_mat->val = arr->val;
            arr->val = NULL;
            
            for (int r = 0; r < nrows; r++)
                for (int c = 0; c < ncols; c++)
                    iftMatrixElem(mat, c, r) = iftMatrixElem(fortran_mat, r, c);
            
            iftDestroyMatrix(&fortran_mat);
        }
        else {
            iftFree(mat->val);
            mat->val = arr->val;
            arr->val = NULL;
        }
        
        iftDestroyFloatArray(&arr);
        iftDestroyNumPyHeader(&header);
        
        return mat;
    }


void iftWriteMatrix(const iftMatrix *mat, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    long shape[2] = {mat->nrows, mat->ncols};
    iftNumPyHeader *header = iftCreateNumPyHeader(IFT_FLT_TYPE, shape, 2);
    iftWriteNumPy(header, mat->val, npy_path);
    iftDestroyNumPyHeader(&header);
}


void iftPrintMatrix(const iftMatrix *M) {
    int i, c, r;
    
    i = 0;
    fprintf(stdout, "\n");
    for (r = 0; r < M->nrows; r++) {
        for (c = 0; c < M->ncols; c++) {
            fprintf(stdout, "%lf ", M->val[i]);
            i++;
        }
        fprintf(stdout, "\n");
    }
}

