#include "iftGenericMatrix.h"

#include "ift/core/io/Stream.h"

iftGenericMatrix* iftCreateGenericMatrix(size_t nrows, size_t ncols, size_t dataSize_bytes){
    iftGenericMatrix* matrix = (iftGenericMatrix*)calloc(1,sizeof(iftGenericMatrix));
    size_t nElements = nrows*ncols;
    matrix->numberColumns = ncols;
    matrix->numberRows = nrows;
    matrix->numberElements = nElements;
    matrix->matrixData = iftCreateGenericNullVector(nElements, dataSize_bytes);
    return matrix;
}

iftGenericMatrix* iftCreateGenericMatrixGivenVector(iftGenericVector* vector){
    iftGenericMatrix* matrix = (iftGenericMatrix*)calloc(1,sizeof(iftGenericMatrix));
    size_t nElements = vector->size;
    matrix->numberColumns = vector->size;
    matrix->numberRows = 1;
    matrix->numberElements = nElements;
    matrix->matrixData = iftCreateGenericNullVector(nElements, vector->elementSize);
    memcpy(matrix->matrixData->data, vector->data, matrix->matrixData->elementSize*matrix->numberElements);
    return matrix;
}

iftGenericMatrix* iftCopyGenericMatrix(iftGenericMatrix* matrix){
    iftGenericMatrix* output = iftCreateGenericMatrix(matrix->numberRows,matrix->numberColumns,matrix->matrixData->elementSize);
    output->numberColumns = matrix->numberColumns;
    output->numberRows = matrix->numberRows;
    output->numberElements = matrix->numberElements;
    memmove(output->matrixData->data, matrix->matrixData->data, matrix->matrixData->elementSize*matrix->numberElements);
    return output;
}

void iftReshapeGenericMatrix(iftGenericMatrix* matrix, size_t nrows, size_t ncols){
    if(nrows*ncols != matrix->numberElements){
        iftError("desired shape not match","iftReshapeGenericMatrix");
        return;
    }
    matrix->numberRows = nrows;
    matrix->numberColumns = ncols;
}



iftGenericMatrix* iftStackVerticallyGenericMatrices(iftGenericMatrix* matrix1, iftGenericMatrix* matrix2){

    if(matrix1 == NULL && matrix2 == NULL){
        return NULL;
    }
    if(matrix1 == NULL){
        return iftCopyGenericMatrix(matrix2);
    }

    if(matrix2 == NULL){
        return iftCopyGenericMatrix(matrix1);
    }

    if(matrix1->numberColumns != matrix2->numberColumns){
        iftError("matrices dimension mismatch","iftStackVerticallyGenericMatrices");
        return NULL;
    }

    if(matrix1->matrixData->elementSize != matrix2->matrixData->elementSize){
        iftError("matrices have different data types","iftStackVerticallyGenericMatrices");
        return NULL;
    }
    size_t totalRows = matrix1->numberRows + matrix2->numberRows;
    size_t nCols = matrix1->numberColumns;
    iftGenericMatrix* output = iftCreateGenericMatrix(totalRows,nCols,matrix1->matrixData->elementSize);

    //first matrix
    size_t nbytes1 =  matrix1->matrixData->elementSize*matrix1->numberElements;
    memcpy(output->matrixData->data, matrix1->matrixData->data, nbytes1);
    //second matrix
    size_t nbytes2 =  matrix2->matrixData->elementSize*matrix2->numberElements;
    memcpy((unsigned char*)output->matrixData->data+nbytes1,(unsigned char*)matrix2->matrixData->data, nbytes2);
    return output;
}

iftGenericMatrix* iftStackVerticallyGenericMatricesGivenVector(iftGenericMatrix* matrix1, iftGenericVector* vector){
    if(matrix1 == NULL && vector == NULL){
        return NULL;
    }
    if(vector == NULL){
        return iftCopyGenericMatrix(matrix1);
    }

    if(matrix1 == NULL){
        return iftCreateGenericMatrixGivenVector(vector);
    }



    if(matrix1->numberColumns != vector->size){
        iftError("matrices dimension mismatch","iftStackVerticallyGenericMatricesGivenVector");
        return NULL;
    }

    if(matrix1->matrixData->elementSize != vector->elementSize){
        iftError("matrices have different data types","iftStackVerticallyGenericMatricesGivenVector");
        return NULL;
    }
    size_t totalRows = matrix1->numberRows + 1;
    size_t nCols = matrix1->numberColumns;
    iftGenericMatrix* output = iftCreateGenericMatrix(totalRows,nCols,matrix1->matrixData->elementSize);

    //first matrix
    size_t nbytes1 =  matrix1->matrixData->elementSize*matrix1->numberElements;
    memcpy(output->matrixData->data, matrix1->matrixData->data, nbytes1);
    //second matrix
    size_t nbytes2 =  vector->elementSize*vector->size;
    memcpy((unsigned char*)output->matrixData->data+nbytes1,(unsigned char*)vector->data, nbytes2);
    return output;
}



void iftDestroyGenericMatrix(iftGenericMatrix **pMatrix){
    if ( (*pMatrix) == NULL || pMatrix == NULL){
        return;
    }
    iftDestroyGenericVector(&((*pMatrix)->matrixData));
    free(*pMatrix);
    *pMatrix= NULL;
}

void iftSwapMatrixRows(iftGenericMatrix* matrix, size_t index1, size_t index2){
    unsigned char* auxBuffer = (unsigned char*)calloc(matrix->matrixData->elementSize*matrix->numberRows,sizeof(unsigned char));
    unsigned char* row1 = (unsigned char*)matrix->matrixData->data + (index1 * matrix->matrixData->elementSize * matrix->numberColumns);
    unsigned char* row2 = (unsigned char*)matrix->matrixData->data + (index2 * matrix->matrixData->elementSize * matrix->numberColumns);
    memmove(auxBuffer,row1,matrix->matrixData->elementSize*matrix->numberColumns);
    memmove(row1,row2,matrix->matrixData->elementSize*matrix->numberColumns);
    memmove(row2,auxBuffer,matrix->matrixData->elementSize*matrix->numberColumns);
    free(auxBuffer);
}

void iftRemoveGenericMatrixRow(iftGenericMatrix* matrix, size_t rowIndex){
    if(rowIndex == matrix->numberRows-1){
        matrix->numberRows -= 1;
        matrix->matrixData->size -= matrix->numberColumns;
        return;
    }

    unsigned char* rowAddress = (unsigned char*)matrix->matrixData->data + (rowIndex * matrix->matrixData->elementSize * matrix->numberColumns);
    unsigned char* nextRowAddress = (unsigned char*)matrix->matrixData->data + ( (rowIndex+1) * matrix->matrixData->elementSize * matrix->numberColumns);
    size_t totalBytes2Move = ((matrix->numberRows)-(rowIndex+1))*matrix->matrixData->elementSize * matrix->numberColumns;
    memmove(rowAddress,nextRowAddress,totalBytes2Move);
    matrix->numberRows -= 1;
    matrix->matrixData->size -= matrix->numberColumns;
}

iftGenericMatrix* iftGetMatrixRows(iftGenericMatrix* matrix, int startRow,int endRow){
    int delta = endRow-startRow;
    if(delta < 0){
        iftError("star row is greather or equal than end row","iftGetMatrixRows");
        return NULL;
    }
    size_t shift_bytes = (delta*matrix->matrixData->elementSize*matrix->numberColumns);
    iftGenericMatrix* outputMatrix = iftCreateGenericMatrix(delta,matrix->numberColumns,matrix->matrixData->elementSize);
    unsigned char* source = (unsigned char*)matrix->matrixData->data + (startRow*matrix->matrixData->elementSize*matrix->numberColumns);
    memmove(outputMatrix->matrixData->data,source,shift_bytes);
    return outputMatrix;
}

//iftGenericMatrix* iftGetGenericMatrixRowsGivenIndices  (iftGenericMatrix* matrix, int* indices,int n){
//    if(n <= 0){
//        iftError("vector of indices is empty","iftGetGenericMatrixRowsGivenIndices");
//        return NULL;
//    }
//    iftGenericMatrix* output = iftCreateGenericMatrix(n,matrix->numberColumns,matrix->matrixData->elementSize);
//    int index;
//    for (int i = 0; i < n; ++i) {
//        index = indices[i];
//        iftCopyGenericMatrixRowGivenLines(output,matrix,i,index);
//    }
//    return output;
//}

void iftDestroyMatrixVoidPointer(void* matrix){
    iftGenericMatrix** aux =  (iftGenericMatrix**) matrix;
    iftDestroyGenericMatrix(aux);
    aux = NULL;
}



