#include "iftMatrix.cuh"

#include "iftMemory.cuh"
#include <iftMatrix.h>
#include "ift/core/tools/Dialog.h"

// these includes CANNOT be inside an extern "C" {
#include <cuda.h>
#include <cublas.h>


void iftMultMatricesInPlaceGPU(const iftMatrix *A, const iftMatrix *B, bool transposeA, bool transposeB, iftMatrix **C) {

    double alpha = 1.0, beta = 0.0;

    int colsA, rowsB, rowsA, colsB;

    colsA = transposeA ? A->nrows : A->ncols;
    rowsB = transposeB ? B->ncols : B->nrows;

    rowsA = transposeA ? A->ncols : A->nrows;
    colsB = transposeB ? B->nrows : B->ncols;

    if (colsA != rowsB)
        iftError("Cannot multiply matrices (%d,%d)%s * (%d,%d)%s", "iftMultMatrices", A->nrows, A->ncols, transposeA ? "^T" : "", B->nrows, B->ncols, transposeB ? "^T" : "");

    if ((*C) == NULL) {
        (*C) = iftCreateMatrix(colsB, rowsA);
    }

    if ((*C)->nrows != rowsA || (*C)->ncols != colsB) {
        iftError("Could not perform the matrices multiplication. Dimensions mismatch. (%d,%d)%s * (%d,%d)%s", "iftMultMatrices", A->nrows, A->ncols, transposeA ? "^T" : "", B->nrows, B->ncols, transposeB ? "^T" : "");
        return;
    }

    iftMatrix* M = *C;

    float *dA;

    float *dB;
    float *dM;

    /* Allocation */
    dA = iftAllocFloatArrayGPU(A->n);
    dB = iftAllocFloatArrayGPU(B->n);
    dM = iftAllocFloatArrayGPU(M->n);

    iftCopyFloatArrayToGPU(dA, A->val, A->n);
    iftCopyFloatArrayToGPU(dB, B->val, B->n);

//    /* Kernel */
//    cublasSgemm(transposeA ? 't' : 'n', transposeB ? 'n' : 't', (*C)->nrows, (*C)->ncols, \
//        colsA, alpha, dA, A->ncols, dB, B->ncols, beta, \
//        dM, (*C)->ncols);

    /* Kernel */
    cublasSgemm(transposeB ? 't' : 'n', transposeA ? 't' : 'n', (*C)->ncols, (*C)->nrows, \
        A->ncols, alpha, dB, B->ncols, dA, A->ncols, beta, \
        dM, (*C)->ncols);
//
//    /* Kernel */
//  cublasSgemm(transposeB ? 'n' : 't', transposeA ? 'n' : 't', (*C)->nrows, (*C)->ncols, \
//        rowsA, alpha, dB, B->ncols, dA, A->ncols, beta, \
//        dM, (*C)->ncols);

//    cblas_sgemm(CblasRowMajor, transposeA ? CblasTrans : CblasNoTrans, transposeB ? CblasTrans : CblasNoTrans, (*C)->nrows, (*C)->ncols, \
//        colsA, alpha, A->val, A->ncols, B->val, B->ncols, beta, \
//        (*C)->val, (*C)->ncols);

    cublasStatus_t status = cublasGetError();
    if (status != CUBLAS_STATUS_SUCCESS) {
        iftError("Check Cuda documentation for error %d.\n", "iftCopyFromGPU", status);
    }

    iftCopyFloatArrayFromGPU(M->val, dM, M->n); 


    iftFreeGPU(dA);
    iftFreeGPU(dB);
    iftFreeGPU(dM);
}

void iftCopyTensorToGPU(iftTensor* dst, iftTensor* src, int n) {
    iftCopyIntArrayToGPU(dst->dimension, src->dimension, src->ndimension);
    iftCopyIntArrayToGPU(dst->accDimension, src->accDimension, src->ndimension);
    iftCopyFloatArrayToGPU(dst->val, src->val, src->n);
}

iftTensor* iftAllocTensorGPU(int n) {
    return (iftTensor*) iftAllocGPU(n, sizeof(float));
}


//
//int main() {
//
//    iftMatrix* M = iftRandomMatrix(5,5,0.0,1.0);
//
//    iftMatrix* N = iftRandomMatrix(5,5,0.0,1.0);
//
//    iftPrintMatrix(iftMultMatrices(M, N));
//
//    return 0;
//}

