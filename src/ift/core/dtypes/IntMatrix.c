#include "ift/core/dtypes/IntMatrix.h"

#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftCommon.h"
#include "iftMatrix.h"
#include "iftMemory.h"


iftIntMatrix *iftCreateIntMatrix(int ncols, int nrows) {
    iftIntMatrix *iM = (iftIntMatrix *) iftAlloc(1, sizeof(iftIntMatrix));
    
    iM->ncols = ncols;
    iM->nrows = nrows;
    
    iM->tbrow = iftAllocIntArray(nrows);
    for (int r = 0; r < nrows; r++) {
        iM->tbrow[r] = r * ncols;
    }
    
    iM->n = ncols * nrows;
    iM->val = iftAllocIntArray(iM->n);
    
    return iM;
}


void iftDestroyIntMatrix(iftIntMatrix **iM) {
    iftIntMatrix *aux = *iM;
    
    if (aux != NULL) {
        if (aux->val != NULL)
            iftFree(aux->val);
        iftFree(aux->tbrow);
        iftFree(aux);
        *iM = NULL;
    }
}


iftIntMatrix *iftReadIntMatrix(const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    iftNumPyHeader *header = NULL;
    void *data = iftReadNumPy(npy_path, &header);
    iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);
    
    if (header->n_dims != 2)
        iftError("Number of dimensions %ld is != 2", "iftReadIntMatrix", header->n_dims);
    
    long nrows = header->shape[0];
    long ncols = header->shape[1];
    
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
    
    iftIntMatrix *mat = iftCreateIntMatrix(ncols, nrows);
    
    
    // fortran order: column-major - the matrix is stored inverted as (ncols, nrows)
    if (header->fortran_order) {
        iftIntMatrix *fortran_mat = iftCreateIntMatrix(nrows, ncols);
        iftFree(fortran_mat->val);
        fortran_mat->val = arr->val;
        arr->val = NULL;
        
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                iftMatrixElem(mat, c, r) = iftMatrixElem(fortran_mat, r, c);
        
        iftDestroyIntMatrix(&fortran_mat);
    }
    else {
        iftFree(mat->val);
        mat->val = arr->val;
        arr->val = NULL;
    }
    
    
    iftDestroyIntArray(&arr);
    iftDestroyNumPyHeader(&header);
    
    return mat;
}


void iftWriteIntMatrix(const iftIntMatrix *mat, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    long shape[2] = {mat->nrows, mat->ncols};
    iftNumPyHeader *header = iftCreateNumPyHeader(IFT_INT_TYPE, shape, 2);
    iftWriteNumPy(header, mat->val, npy_path);
    iftDestroyNumPyHeader(&header);
}



iftIntMatrix *iftCopyIntMatrix(int *val, int nrows, int ncols) {
    iftIntMatrix *mat = iftCreateIntMatrix(ncols, nrows);
    
    for (int i = 0; i < mat->n; i++)
        mat->val[i] = val[i];
    
    return mat;
}


void iftResizeIntMatrix(iftIntMatrix **M, int ncols, int nrows) {
    if((ncols == 0) || (nrows ==0))
        iftError("The new quantity of columns and rows must be greater than 0", "iftResizeIntMatrix");
    
    if((*M) == NULL) {
        (*M) = iftCreateIntMatrix(ncols, nrows);
    } else {
        if((ncols != (*M)->ncols) || (nrows != (*M)->nrows)) {
            iftIntMatrix *M1 = iftCopyIntMatrix((*M)->val, (*M)->nrows, (*M)->ncols);
            iftDestroyIntMatrix(M);
            (*M) = iftCreateIntMatrix(ncols, nrows);
            
            for(int r = 0; r < iftMin(M1->nrows, (*M)->nrows); r++) {
                for(int c = 0; c < iftMin(M1->ncols, (*M)->ncols); c++) {
                    iftMatrixElem((*M), c, r) = iftMatrixElem(M1, c, r);
                }
            }
            iftDestroyIntMatrix(&M1);
        }
    }
}


char *iftIntMatrixAsString(const iftIntMatrix *mat) {
    char *str = iftCopyString("[");
    char *aux = NULL;
    
    for (long i = 0; i <= (mat->nrows-2); i++) {
        // 2 chars for square brackets + 10 chars per number (at most) + n-1 commas
        char *str_row = iftAlloc(2 + (10 * mat->ncols) + (mat->ncols-1), sizeof(char));
        
        sprintf(str_row, "[");
        for (long j = 0; j <= (mat->ncols-2); j++) {
            sprintf(str_row, "%s%d, ", str_row, iftMatrixElem(mat, j, i));
        }
        sprintf(str_row, "%s%d]", str_row, iftMatrixElem(mat, mat->ncols-1, i));
        
        aux = str;
        str = iftConcatStrings(3, str, str_row, ", ");
        iftFree(aux);
    }
    
    char *str_row = iftAlloc(2 + (10 * mat->ncols) + (mat->ncols-1), sizeof(char));
    sprintf(str_row, "[");
    for (long j = 0; j <= (mat->ncols-2); j++) {
        sprintf(str_row, "%s%d, ", str_row, iftMatrixElem(mat, j, mat->nrows-1));
    }
    sprintf(str_row, "%s%d]", str_row, iftMatrixElem(mat, mat->ncols-1, mat->nrows-1));
    
    aux = str;
    str = iftConcatStrings(3, str, str_row, "]");
    iftFree(aux);
    
    return str;
}


void iftPrintIntMatrix(const iftIntMatrix* M) {
    for (int r = 0; r < M->nrows; r++) {
        printf("[%d] = [", r);
        
        for (int c = 0; c < (M->ncols - 1); c++) {
            printf(" %d, ", iftMatrixElem(M, c, r));
        }
        printf(" %d ]\n", iftMatrixElem(M, M->ncols - 1, r));
    }
}


iftIntMatrix* iftMergeIntMatrixArrayHorizontal(iftIntMatrix **matArray, int n)
{
    /* compute the number of cols of the new matrix */
    int newColSize = 0, rowSize = 0, lastRowSize = IFT_NIL;
    for(int i = 0; i < n; i++) {
        if(matArray[i] != NULL) {
            newColSize += matArray[i]->ncols;
            rowSize = matArray[i]->nrows;
            if(lastRowSize != IFT_NIL && lastRowSize != matArray[i]->nrows) {
                iftError("Number of rows in all the matrices must be the same.", "iftMergeIntMatrixArrayHorizontal");
            }
            else {
                lastRowSize = matArray[i]->nrows;
            }
        }
    }
    
    /* join the matrices into a new one */
    iftIntMatrix *newMat = iftCreateIntMatrix(newColSize, rowSize);
    int row = 0, col = 0;
    for(int i = 0; i < n; i++) {
        if(matArray[i] != NULL) {
            for(int c = 0; c < matArray[i]->ncols; c++) {
                row = 0;
                for(int r = 0; r < matArray[i]->nrows; r++) {
                    iftMatrixElem(newMat, col, row) = iftMatrixElem(matArray[i], c, r);
                    row++;
                }
                col++;
            }
        }
    }
    
    return newMat;
}

