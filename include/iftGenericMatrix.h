//
// Created by deangeli on 5/15/17.
//

#ifndef _IFTGENERICMATRIX_H
#define _IFTGENERICMATRIX_H

#include "iftCommon.h"
#include "iftGenericVector.h"

typedef struct _iftGenericMatrix{
    iftGenericVector* matrixData;
    size_t numberRows;
    size_t numberColumns;
    size_t numberElements;
}iftGenericMatrix;

iftGenericMatrix* iftCreateGenericMatrix(size_t nrows, size_t ncols, size_t dataSize_bytes);
iftGenericMatrix* iftCreateGenericMatrixGivenVector(iftGenericVector* vector);
iftGenericMatrix* iftStackVerticallyGenericMatrices(iftGenericMatrix* matrix1,iftGenericMatrix* matrix2);
iftGenericMatrix* iftStackVerticallyGenericMatricesGivenVector(iftGenericMatrix* matrix1, iftGenericVector* vector);
void iftReshapeGenericMatrix(iftGenericMatrix* matrix, size_t nrows, size_t ncols);
void iftDestroyGenericMatrix(iftGenericMatrix **pMatrix);
iftGenericMatrix* iftCopyGenericMatrix(iftGenericMatrix* matrix);
void iftSwapGenericMatrixRows(iftGenericMatrix* matrix, size_t index1, size_t index2);
void iftRemoveGenericMatrixRow(iftGenericMatrix* matrix, size_t rowIndex);
iftGenericMatrix* iftGetGenericMatrixRows(iftGenericMatrix* matrix, int startRow,int endRow);
iftGenericMatrix* iftGetGenericMatrixRowsGivenIndices(iftGenericMatrix* matrix, int* indices,int n);
void iftDestroyGenericMatrixVoidPointer(void* matrix);

static inline void iftCopyGenericMatrixRowGivenLines(iftGenericMatrix* matrix1,iftGenericMatrix* matrix2, size_t rowIndex1,size_t rowIndex2){
    size_t nbytes =  matrix1->matrixData->elementSize*matrix1->numberColumns;
    size_t memoryShift1 = matrix1->numberColumns*rowIndex1*matrix1->matrixData->elementSize;
    size_t memoryShift2 = matrix2->numberColumns*rowIndex2*matrix2->matrixData->elementSize;
    memmove( ((unsigned char*)matrix1->matrixData->data)+memoryShift1,
           ((unsigned char*)matrix2->matrixData->data)+memoryShift2,
           nbytes);
}

static inline void iftSetRowValueGivenVector(iftGenericMatrix* matrix, iftGenericVector* vector, size_t rowIndex){
    size_t shift =  matrix->matrixData->elementSize*(matrix->numberColumns*rowIndex);
    size_t nbytes = vector->elementSize*vector->size;
    memcpy((unsigned char*)matrix->matrixData->data+shift,(unsigned char*)vector->data, nbytes);
}

//inline void getSubmetrix(Matrix* matrix,size_t startRow,size_t endRow, size_t startCol, size_t endRow){
//
//}

static inline void iftCopyGenericMatrixRow(iftGenericMatrix* matrix1,iftGenericMatrix* matrix2, size_t rowIndex1){
    size_t rowIndex2 = rowIndex1;
    iftCopyGenericMatrixRowGivenLines(matrix1,matrix2, rowIndex1,rowIndex2);
}



#define MATRIX_GET_ELEMENT_PO_AS(type, matrix, rowIndex,colIndex) iftGenericVectorGetElementAt(type,matrix->matrixData, MATRIX_COMPUTEINDEX(matrix,rowIndex,colIndex))

#define MATRIX_GET_ELEMENT_BI_AS(type, matrix, index) VECTOR_GET_ELEMENT_AS(type,matrix->matrixData, index)

#define MATRIX_COMPUTEINDEX(matrix,rowIndex,colIndex) ((rowIndex*matrix->numberColumns)+colIndex)

#define MATRIX_ADD_ROWS_AS(type,matrix1,matrix2,index1,index2) \
    { \
        type* aux1 = ( (type*) matrix1->matrixData->data) + (index1*matrix1->numberColumns); \
        type* aux2 = ( (type*) matrix2->matrixData->data) + (index2*matrix2->numberColumns); \
        for(size_t index=0; index < matrix1->numberColumns; index++) { \
            aux1[index] += aux2[index]; \
        } \
    } \

#define MATRIX_PRINT_AS(type,symbol,matrix) \
    for(size_t i = 0; i< matrix->numberRows; i++){\
        for(size_t j = 0; j< matrix->numberColumns; j++){\
            printf(symbol, MATRIX_GET_ELEMENT_PO_AS(type,matrix,i,j) ); \
        } \
        printf("\n"); \
    }
#endif //_MATRIX_H
