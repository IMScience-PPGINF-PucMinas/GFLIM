#include "iftMatrix.h"

#include "iftMatrix.cuh"
#include "iftMemory.cuh"
#include "ift/core/io/CSV.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"

//#include "/usr/include/cblas_f77.h" - menotti

/* Decompose A_[m x n] into U_[m x m] S_[m x n] V^t_[n x n], where S
   contains the singular values and the matrices U and V^t are
   orthogonal matrices.
*/

extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );


extern void sgesdd_(char *jobz, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work,
        int *lwork, int *iwork, int *info);

// LU decomoposition of a general matrix
extern int sgetrf_(int *M, int *N, float *A, int *lda, int *IPIV, int *INFO);

// generate inverse of a matrix given its LU decomposition
extern int sgetri_(int *N, float *A, int *lda, int *IPIV, float *WORK, int *lwork, int *INFO);

// / LU decomoposition of a general matrix
extern int dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);

// generate inverse of a matrix given its LU decomposition
extern int dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork, int *INFO);

iftMatrix* iftCreateMatrixPointer(float *val, int ncol, int nrow) {
    iftMatrix* matrix = iftAlloc(1, sizeof(iftMatrix));

    matrix->val = val;
    matrix->nrows = nrow;
    matrix->ncols = ncol;
    matrix->allocated = false;
    matrix->n = nrow*ncol;

    matrix->tbrow = (long*)iftAlloc(nrow,sizeof(long));

    for (int r = 0; r < nrow; ++r) {
        matrix->tbrow[r] = r * ncol;
    }

    return matrix;
}

void iftReshapeMatrix(iftMatrix* matrix, int newCol, int newRow){
    if(newCol*newRow != matrix->n){
        iftError("Dimensions mismatch","iftReshapeMatrix");
        return;
    }
    matrix->nrows = newRow;
    matrix->ncols = newCol;
}

iftMatrix* iftPolynomialExpasion(iftMatrix* matrix, unsigned int degree){
    if(degree != 2){
        iftError("This function does not  supprot this degree","iftPolynomialExpasion");
        return NULL;
    }
    long n_poly0 = 1;
    long n_poly1 =matrix->ncols;
    long n_poly2 = n_poly1*(n_poly1+1)/2;

    long numberElements = n_poly0 + n_poly1 + n_poly2;
    iftMatrix* output = iftCreateMatrix(numberElements,matrix->nrows);
    for (int row = 0; row < matrix->nrows; ++row) {
        iftMatrixElem(output,0,row) = 1;
        int k = 1;
        for (int i = 0; i < n_poly1; ++i) {
            iftMatrixElem(output,k,row) = iftMatrixElem(matrix,i,row);
            k++;
        }

        for (int j = 0; j < n_poly1; ++j) {
            for (int i = j; i < n_poly1; ++i) {
                iftMatrixElem(output,k,row) = iftMatrixElem(matrix,j,row)*iftMatrixElem(matrix,i,row);
                k++;
            }
        }
    }
    return output;
}


void iftDestroyMatrixPointer(iftMatrix** matrix) {
    iftFree((*matrix)->tbrow);
    iftFree(*matrix);
    *matrix = NULL;
}

iftDoubleMatrix* iftCopyDoubleMatrix(const iftDoubleMatrix* m){
    iftDoubleMatrix *M = iftCreateDoubleMatrix(m->ncols, m->nrows);
    iftCopyDoubleArray(M->val, m->val, m->n);
    return (M);
}

iftDoubleMatrix* iftCreateDoubleMatrix(int ncols, int nrows){
    iftDoubleMatrix *M = (iftDoubleMatrix *) iftAlloc(1, sizeof(iftDoubleMatrix));
    M->ncols = ncols;
    M->nrows = nrows;
    M->tbrow = (long*)iftAlloc(nrows,sizeof(long));
    //M->tbrow = iftAllocIntArray(nrows);
    for (long r = 0; r < (long)nrows; r++) {
        M->tbrow[r] = r * ncols;
    }
    M->n = (long) ncols * (long) nrows;
    M->allocated = true;

    M->val = iftAlloc(M->n,sizeof(double));

    return (M);
}

iftDoubleMatrix* iftCreateDoubleMatrixNoAlloc(int ncols, int nrows)
{
    iftDoubleMatrix *M = (iftDoubleMatrix *) iftAlloc(1, sizeof(iftDoubleMatrix));
    M->ncols = ncols;
    M->nrows = nrows;
    M->tbrow = (long*)iftAlloc(nrows,sizeof(long));
    //M->tbrow = iftAllocIntArray(nrows);
    for (long r = 0; r < (long)nrows; r++) {
        M->tbrow[r] = r * ncols;
    }
    M->n = (long) ncols * (long) nrows;
    M->allocated = false;
    M->val = NULL;

    return (M);
}

void iftReallocDoubleMatrix(iftDoubleMatrix *M, int ncols, int nrows)
{
    M->ncols = ncols;
    M->nrows = nrows;
    M->tbrow = (long*)iftRealloc(M->tbrow, nrows * (sizeof *M->tbrow));

    for (long r = 0; r < (long)nrows; r++) {
        M->tbrow[r] = r * ncols;
    }
    M->n = (long) ncols * (long) nrows;
    M->allocated = true;
    M->val = iftRealloc(M->val, M->n * (sizeof *M->val));
}


void iftDestroyDoubleMatrix(iftDoubleMatrix **M){
    if(M == NULL){
        return;
    }
    iftDoubleMatrix *aux = *M;
    if(aux == NULL){
        return;
    }
    if (aux->allocated == true && aux->val != NULL){
        iftFree(aux->val);
    }
    iftFree(aux->tbrow);
    iftFree(aux);
    *M = NULL;
}

iftMatrix *iftReadRawMatrixFILE(FILE *fp) {
    iftMatrix *M = NULL;
    int ncols, nrows;

    if (fread(&ncols, sizeof(int), 1, fp) != 1) iftError("Reading error", "iftReadRawMatrixFILE");
    if (fread(&nrows, sizeof(int), 1, fp) != 1) iftError("Reading error", "iftReadRawMatrixFILE");

    M = iftCreateMatrix(ncols, nrows);

    if (fread((M->val), sizeof(float), M->n, fp) != M->n) iftError("Reading error", "iftReadRawMatrixFILE");

    return M;
}

void iftWriteMatrixCSV(iftMatrix* M, const char* filename) {

    iftCSV* csv = iftCreateCSV(M->nrows, M->ncols);

    for (int r = 0; r < csv->nrows; ++r) {
        for (int c = 0; c < csv->ncols; ++c) {
            sprintf(csv->data[r][c], "%f", iftMatrixElem(M, c, r));
        }
    }

    iftWriteCSV(csv, filename, ',');
    iftDestroyCSV(&csv);
}

iftMatrix *iftReadMatrixCSV(const char* filename) {

    iftCSV* csv = iftReadCSV(filename, ',');
    iftMatrix *M = iftCreateMatrix(csv->ncols, csv->nrows);

    for (int r = 0; r < csv->nrows; ++r) {
        for (int c = 0; c < csv->ncols; ++c) {
            iftMatrixElem(M, c, r) = atof(csv->data[r][c]);
        }
    }

    iftDestroyCSV(&csv);

    return M;
}

void iftWriteRawMatrixFILE(iftMatrix *M, FILE *fp) {
    if (fwrite(&(M->ncols), sizeof(int), 1, fp) != 1) iftError("Writing error", "iftWriteRawMatrixFILE");
    if (fwrite(&(M->nrows), sizeof(int), 1, fp) != 1) iftError("Writing error", "iftWriteRawMatrixFILE");

    if (fwrite((M->val), sizeof(float), M->n, fp) != M->n) iftError("Writing error", "iftWriteRawMatrixFILE");
}

iftMatrix *iftReadRawMatrix(const char *filename) {
    iftDeprecated("iftReadRawMatrix", "iftReadMatrix", "");
    
    FILE *fp = fopen(filename, "rb");

    if(!fp)
        iftError("File %s does not exist", "iftReadRawMatrix", filename);

    iftMatrix *M = NULL;

    M = iftReadRawMatrixFILE(fp);

    fclose(fp);

    return (M);
}

void iftWriteRawMatrix(iftMatrix *M, char *filename) {
    iftDeprecated("iftWriteRawMatrix", "iftWriteMatrix", "");
    
    FILE *fp = fopen(filename, "wb");

    iftWriteRawMatrixFILE(M, fp);

    fclose(fp);
}


void iftWriteBoolMatrix(iftBoolMatrix *M, char *filename)
{
    FILE *fp = fopen(filename, "wb");

    if (fwrite(&(M->ncols), sizeof(int), 1, fp) != 1) iftError("Writing error", "iftWriteBoolMatrix");
    if (fwrite(&(M->nrows), sizeof(int), 1, fp) != 1) iftError("Writing error", "iftWriteBoolMatrix");

    if (fwrite((M->val), sizeof(bool), M->n, fp) != M->n) iftError("Writing error", "iftWriteIntMatrix");

    fclose(fp);
}

iftBoolMatrix *iftReadBoolMatrix(const char *filename)
{
    FILE *fp = fopen(filename, "rb");

    if(!fp)
        iftError("File %s does not exist", "iftReadBoolMatrix", filename);

    iftBoolMatrix *M = NULL;
    int ncols, nrows;

    if (fread(&ncols, sizeof(int), 1, fp) != 1) iftError("Reading error", "iftReadBoolMatrix");
    if (fread(&nrows, sizeof(int), 1, fp) != 1) iftError("Reading error", "iftReadBoolMatrix");

    M = iftCreateBoolMatrix(ncols, nrows);

    if (fread((M->val), sizeof(bool), M->n, fp) != M->n) iftError("Reading error", "iftReadBoolMatrix");

    fclose(fp);

    return (M);
}


void iftLinearCombination(float a, float* A, float b, float* B, float *C, int n) {

	for (int i = 0; i < n; ++i) {
        C[i] = a * A[i] + b * B[i];
    }

//	if(A!=C) {
//		iftCopyFloatArray(C, A, n);
//	}
//	cblas_sscal(n, a, C, 1);
//	cblas_saxpy(n, b, B, 1, C, 1);
}

void iftPrintDoubleMatrix(iftDoubleMatrix *M){
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



iftMatrix *iftMultMatrices(const iftMatrix *A, const iftMatrix *B) {
    iftMatrix *M = NULL; /* M = alpha A*B + beta M */

    if (A->ncols != B->nrows) {
		iftError("Cannot multiply matrices (%d,%d) * (%d,%d)", "iftMultMatrices", A->nrows, A->ncols, B->nrows,
				 B->ncols);
	}

    M = iftCreateMatrix(B->ncols, A->nrows);
    iftMultMatricesInPlace(A, B, 0, 0, &M);
    return (M);
}

void iftMultMatricesInPlace(const iftMatrix *A, const iftMatrix *B, bool transposeA, bool transposeB, iftMatrix **C) {

#ifdef IFT_GPU
    bool run_gpu = FALSE;
    float expectedMemoryAllocMB = sizeof(float)*(A->n + B->n + (B->ncols*A->nrows))/1024.0/1024.0; //A->n + B->n + C->n
    float freeMemoryMB = (float)iftGetFreeMemoryGPU(0)/1024.0/1024.0;
    if (expectedMemoryAllocMB > freeMemoryMB)
        run_gpu = FALSE;
    else
        run_gpu = TRUE;
    if (run_gpu){
        iftMultMatricesInPlaceGPU(A, B, transposeA, transposeB, C);
    } else {
#endif
        double alpha = 1.0, beta = 0.0;

        int colsA, rowsB, rowsA, colsB;

        colsA = transposeA ? A->nrows : A->ncols;
        rowsB = transposeB ? B->ncols : B->nrows;

        rowsA = transposeA ? A->ncols : A->nrows;
        colsB = transposeB ? B->nrows : B->ncols;

        if (colsA != rowsB)
            iftError("Cannot multiply matrices (%d,%d)%s * (%d,%d)%s", "iftMultMatrices", A->nrows, A->ncols,
                     transposeA ? "^T" : "", B->nrows, B->ncols, transposeB ? "^T" : "");

        if ((*C) == NULL) {
            (*C) = iftCreateMatrix(colsB, rowsA);
        }

        if ((*C)->nrows != rowsA || (*C)->ncols != colsB) {
            iftError(
                    "Could not perform the matrices multiplication. Dimensions mismatch. (%d,%d) = (%d,%d)%s * (%d,%d)%s",
                    "iftMultMatrices", (*C)->nrows, (*C)->ncols, A->nrows, A->ncols, transposeA ? "^T" : "", B->nrows,
                    B->ncols, transposeB ? "^T" : "");
            return;
        }

        /* Compute multiplication between matrices */
        cblas_sgemm(CblasRowMajor, transposeA ? CblasTrans : CblasNoTrans, transposeB ? CblasTrans : CblasNoTrans,
                    (*C)->nrows, (*C)->ncols, \
        colsA, alpha, A->val, A->ncols, B->val, B->ncols, beta, \
        (*C)->val, (*C)->ncols);
#ifdef IFT_GPU
    }
#endif
}


iftMatrix *iftMultMatricesChain(int n, ...) {
    int i;
    va_list arguments;

    /* Initializing arguments to store all values after num */
    va_start (arguments, n);

    iftMatrix *cur, *next, *temp;

    cur = iftCopyMatrix(va_arg (arguments, iftMatrix*));

    for (i = 1; i < n; ++i) {
        temp = va_arg (arguments, iftMatrix*);
        next = iftMultMatrices(cur, temp);

        iftDestroyMatrix(&cur);
        cur = next;
    }

    va_end (arguments);                  // Cleans up the list

    return cur;
}


iftMatrix *iftTransposeMatrix(const iftMatrix *A) {
    iftMatrix *B = NULL;
    int c, r;

    B = iftCreateMatrix(A->nrows, A->ncols);
    for (r = 0; r < B->nrows; r++) {
        for (c = 0; c < B->ncols; c++) {
            iftMatrixElem(B, c, r) = iftMatrixElem(A, r, c);
        }
    }
    return (B);
}

void iftTransposeMatrixInPlace(iftMatrix *A) {
    iftMatrix *B = iftTransposeMatrix(A);
    iftCopyFloatArray(A->val, B->val, A->n);

    if(A->nrows < A->ncols) {
        iftFree(A->tbrow);
        A->tbrow = (long*)iftAlloc(A->ncols,sizeof(long));
    }

    iftCopySizeTArray(A->tbrow, B->tbrow, B->nrows);

	iftSwap(A->nrows, A->ncols);
	iftDestroyMatrix(&B);
}

iftMatrix *iftRotationMatrix(char axis, float theta) {
    iftMatrix *A;
    float cos_theta, sin_theta;

    A = iftCreateMatrix(4, 4);
    theta = theta * IFT_PI / 180.0;
    cos_theta = cosf(theta);
    sin_theta = sinf(theta);

    switch (axis) {

        case IFT_AXIS_X:

            A->val[iftGetMatrixIndex(A, 0, 0)] = 1.0;
            A->val[iftGetMatrixIndex(A, 1, 0)] = 0.0;
            A->val[iftGetMatrixIndex(A, 2, 0)] = 0.0;
            A->val[iftGetMatrixIndex(A, 3, 0)] = 0.0;

            A->val[iftGetMatrixIndex(A, 0, 1)] = 0.0;
            A->val[iftGetMatrixIndex(A, 1, 1)] = cos_theta;
            A->val[iftGetMatrixIndex(A, 2, 1)] = -sin_theta;
            A->val[iftGetMatrixIndex(A, 3, 1)] = 0.0;

            A->val[iftGetMatrixIndex(A, 0, 2)] = 0.0;
            A->val[iftGetMatrixIndex(A, 1, 2)] = sin_theta;
            A->val[iftGetMatrixIndex(A, 2, 2)] = cos_theta;
            A->val[iftGetMatrixIndex(A, 3, 2)] = 0.0;

            A->val[iftGetMatrixIndex(A, 0, 3)] = 0.0;
            A->val[iftGetMatrixIndex(A, 1, 3)] = 0.0;
            A->val[iftGetMatrixIndex(A, 2, 3)] = 0.0;
            A->val[iftGetMatrixIndex(A, 3, 3)] = 1.0;

            break;

        case IFT_AXIS_Y:

            A->val[iftGetMatrixIndex(A, 0, 0)] = cos_theta;
            A->val[iftGetMatrixIndex(A, 1, 0)] = 0.0;
            A->val[iftGetMatrixIndex(A, 2, 0)] = sin_theta;
            A->val[iftGetMatrixIndex(A, 3, 0)] = 0.0;

            A->val[iftGetMatrixIndex(A, 0, 1)] = 0.0;
            A->val[iftGetMatrixIndex(A, 1, 1)] = 1.0;
            A->val[iftGetMatrixIndex(A, 2, 1)] = 0.0;
            A->val[iftGetMatrixIndex(A, 3, 1)] = 0.0;

            A->val[iftGetMatrixIndex(A, 0, 2)] = -sin_theta;
            A->val[iftGetMatrixIndex(A, 1, 2)] = 0.0;
            A->val[iftGetMatrixIndex(A, 2, 2)] = cos_theta;
            A->val[iftGetMatrixIndex(A, 3, 2)] = 0.0;

            A->val[iftGetMatrixIndex(A, 0, 3)] = 0.0;
            A->val[iftGetMatrixIndex(A, 1, 3)] = 0.0;
            A->val[iftGetMatrixIndex(A, 2, 3)] = 0.0;
            A->val[iftGetMatrixIndex(A, 3, 3)] = 1.0;

            break;

        case IFT_AXIS_Z:

            A->val[iftGetMatrixIndex(A, 0, 0)] = cos_theta;
            A->val[iftGetMatrixIndex(A, 1, 0)] = -sin_theta;
            A->val[iftGetMatrixIndex(A, 2, 0)] = 0.0;
            A->val[iftGetMatrixIndex(A, 3, 0)] = 0.0;

            A->val[iftGetMatrixIndex(A, 0, 1)] = sin_theta;
            A->val[iftGetMatrixIndex(A, 1, 1)] = cos_theta;
            A->val[iftGetMatrixIndex(A, 2, 1)] = 0.0;
            A->val[iftGetMatrixIndex(A, 3, 1)] = 0.0;

            A->val[iftGetMatrixIndex(A, 0, 2)] = 0.0;
            A->val[iftGetMatrixIndex(A, 1, 2)] = 0.0;
            A->val[iftGetMatrixIndex(A, 2, 2)] = 1.0;
            A->val[iftGetMatrixIndex(A, 3, 2)] = 0.0;

            A->val[iftGetMatrixIndex(A, 0, 3)] = 0.0;
            A->val[iftGetMatrixIndex(A, 1, 3)] = 0.0;
            A->val[iftGetMatrixIndex(A, 2, 3)] = 0.0;
            A->val[iftGetMatrixIndex(A, 3, 3)] = 1.0;

            break;

        default:
            iftError("Invalid option for the axis", "iftRotationMatrix");
    }

    return (A);
}

iftMatrix *iftTranslationMatrix(iftVector T) {
    iftMatrix *A;

    A = iftCreateMatrix(4, 4);

    A->val[iftGetMatrixIndex(A, 0, 0)] = 1.0;
    A->val[iftGetMatrixIndex(A, 0, 1)] = 0.0;
    A->val[iftGetMatrixIndex(A, 0, 2)] = 0.0;
    A->val[iftGetMatrixIndex(A, 0, 3)] = 0.0;

    A->val[iftGetMatrixIndex(A, 1, 0)] = 0.0;
    A->val[iftGetMatrixIndex(A, 1, 1)] = 1.0;
    A->val[iftGetMatrixIndex(A, 1, 2)] = 0.0;
    A->val[iftGetMatrixIndex(A, 1, 3)] = 0.0;

    A->val[iftGetMatrixIndex(A, 2, 0)] = 0.0;
    A->val[iftGetMatrixIndex(A, 2, 1)] = 0.0;
    A->val[iftGetMatrixIndex(A, 2, 2)] = 1.0;
    A->val[iftGetMatrixIndex(A, 2, 3)] = 0.0;

    A->val[iftGetMatrixIndex(A, 3, 0)] = T.x;
    A->val[iftGetMatrixIndex(A, 3, 1)] = T.y;
    A->val[iftGetMatrixIndex(A, 3, 2)] = T.z;
    A->val[iftGetMatrixIndex(A, 3, 3)] = 1.0;

    return (A);
}

iftMatrix *iftScaleMatrix(float sx, float sy, float sz) {
    iftMatrix *A;

    A = iftCreateMatrix(4, 4);

    A->val[iftGetMatrixIndex(A, 0, 0)] = sx;
    A->val[iftGetMatrixIndex(A, 1, 0)] = 0.0;
    A->val[iftGetMatrixIndex(A, 2, 0)] = 0.0;
    A->val[iftGetMatrixIndex(A, 3, 0)] = 0.0;

    A->val[iftGetMatrixIndex(A, 0, 1)] = 0.0;
    A->val[iftGetMatrixIndex(A, 1, 1)] = sy;
    A->val[iftGetMatrixIndex(A, 2, 1)] = 0.0;
    A->val[iftGetMatrixIndex(A, 3, 1)] = 0.0;

    A->val[iftGetMatrixIndex(A, 0, 2)] = 0.0;
    A->val[iftGetMatrixIndex(A, 1, 2)] = 0.0;
    A->val[iftGetMatrixIndex(A, 2, 2)] = sz;
    A->val[iftGetMatrixIndex(A, 3, 2)] = 0.0;

    A->val[iftGetMatrixIndex(A, 0, 3)] = 0.0;
    A->val[iftGetMatrixIndex(A, 1, 3)] = 0.0;
    A->val[iftGetMatrixIndex(A, 2, 3)] = 0.0;
    A->val[iftGetMatrixIndex(A, 3, 3)] = 1.0;

    return (A);
}

iftMatrix *iftCopyMatrix(const iftMatrix *A) {
    iftMatrix *B;
    int i;

    B = iftCreateMatrix(A->ncols, A->nrows);

    for (i = 0; i < A->n; i++)
        B->val[i] = A->val[i];

    return (B);
}


iftVector iftTransformVector(iftMatrix *M, iftVector v) {
    iftVector v1;
    // Since iftMultMatrices uses cblas, it is faster to simply do the transformation
    // directly because iftVector only has 3 coordinates
    v1.x = M->val[iftGetMatrixIndex(M, 0, 0)] * v.x + M->val[iftGetMatrixIndex(M, 1, 0)] * v.y +
           M->val[iftGetMatrixIndex(M, 2, 0)] * v.z;
    v1.y = M->val[iftGetMatrixIndex(M, 0, 1)] * v.x + M->val[iftGetMatrixIndex(M, 1, 1)] * v.y +
           M->val[iftGetMatrixIndex(M, 2, 1)] * v.z;
    v1.z = M->val[iftGetMatrixIndex(M, 0, 2)] * v.x + M->val[iftGetMatrixIndex(M, 1, 2)] * v.y +
           M->val[iftGetMatrixIndex(M, 2, 2)] * v.z;

    return (v1);
}

iftPoint iftTransformPoint(iftMatrix *M, iftPoint v) {
    iftPoint v1;

    // Since iftMultMatrices uses cblas, it is faster to simply do the transformation
    // directly because iftPoint only has 3 coordinates
    v1.x = M->val[iftGetMatrixIndex(M, 0, 0)] * v.x + M->val[iftGetMatrixIndex(M, 1, 0)] * v.y +
           M->val[iftGetMatrixIndex(M, 2, 0)] * v.z + M->val[iftGetMatrixIndex(M, 3, 0)];
    v1.y = M->val[iftGetMatrixIndex(M, 0, 1)] * v.x + M->val[iftGetMatrixIndex(M, 1, 1)] * v.y +
           M->val[iftGetMatrixIndex(M, 2, 1)] * v.z + M->val[iftGetMatrixIndex(M, 3, 1)];
    v1.z = M->val[iftGetMatrixIndex(M, 0, 2)] * v.x + M->val[iftGetMatrixIndex(M, 1, 2)] * v.y +
           M->val[iftGetMatrixIndex(M, 2, 2)] * v.z + M->val[iftGetMatrixIndex(M, 3, 2)];

    return (v1);
}

iftVoxel iftTransformVoxel(const iftMatrix *M, iftVoxel v) {
    iftVoxel vt;

    // Since iftMultMatrices uses cblas, it is faster to simply do the transformation
    // directly because iftPoint only has 3 coordinates
    vt.x = iftRound(M->val[iftGetMatrixIndex(M, 0, 0)] * v.x + M->val[iftGetMatrixIndex(M, 1, 0)] * v.y +
                    M->val[iftGetMatrixIndex(M, 2, 0)] * v.z + M->val[iftGetMatrixIndex(M, 3, 0)]);
    vt.y = iftRound(M->val[iftGetMatrixIndex(M, 0, 1)] * v.x + M->val[iftGetMatrixIndex(M, 1, 1)] * v.y +
                    M->val[iftGetMatrixIndex(M, 2, 1)] * v.z + M->val[iftGetMatrixIndex(M, 3, 1)]);
    vt.z = iftRound(M->val[iftGetMatrixIndex(M, 0, 2)] * v.x + M->val[iftGetMatrixIndex(M, 1, 2)] * v.y +
                    M->val[iftGetMatrixIndex(M, 2, 2)] * v.z + M->val[iftGetMatrixIndex(M, 3, 2)]);

    return vt;
}


double iftMatrixDeterminant(const iftMatrix *A) {
    int i, j, k, l, N = A->ncols;
    double det = 0;
    iftMatrix *B = NULL;

    if (A->ncols != A->nrows)
        iftError("Matrix is not square", "iftMatrixDeterminant");

    if (N == 1) { /* trivial case */
        det = A->val[0];
    } else { /* 2 x 2 matrix */
        if (N == 2) {
            det = (A->val[0] * A->val[3]) - (A->val[2] * A->val[1]);
        } else { /* N x N matrix */
            det = 0;
            for (k = 0; k < N; k++) {
                B = iftCreateMatrix(N - 1, N - 1);
                for (i = 1; i < N; i++) {
                    l = 0;
                    for (j = 0; j < N; j++) {
                        if (j != k) {
                            B->val[iftGetMatrixIndex(B, l, i - 1)] =
                                    A->val[iftGetMatrixIndex(A, j, i)];
                            l++;
                        }
                    }
                }
                det += pow(-1.0, k + 2.0) * A->val[iftGetMatrixIndex(A, k, 0)] *
                       iftMatrixDeterminant(B);
                iftDestroyMatrix(&B);
            }
        }
    }


    return (det);
}

iftMatrix *iftCoFactorMatrix(const iftMatrix *A) {
    int r1, c1, r2, c2, r3, c3, N = A->ncols;
    double det;
    iftMatrix *B, *C;

    if (A->ncols != A->nrows)
        iftError("Matrix is not square", "iftCoFactorMatrix");

    B = iftCreateMatrix(N, N);
    C = iftCreateMatrix(N - 1, N - 1);

    for (c1 = 0; c1 < N; c1++) {
        for (r1 = 0; r1 < N; r1++) {

            r3 = 0;
            for (r2 = 0; r2 < N; r2++) {
                if (r2 != r1) {
                    c3 = 0;
                    for (c2 = 0; c2 < N; c2++) {
                        if (c2 != c1) {
                            C->val[iftGetMatrixIndex(C, c3, r3)] = A->val[iftGetMatrixIndex(A, c2, r2)];
                            c3++;
                        }
                    }
                    r3++;
                }
            }
            det = iftMatrixDeterminant(C);
            B->val[iftGetMatrixIndex(B, c1, r1)] = pow(-1.0, r1 + c1 + 2.0) * det;
        }
    }

    iftDestroyMatrix(&C);
    return (B);
}

iftMatrix *iftInvertMatrix(const iftMatrix *A) {
    iftMatrix *B = NULL;
    double det;

    if (A->ncols != A->nrows)
        iftError("Matrix is not square", "iftInvertMatrix");

    det = iftMatrixDeterminant(A);

    if (fabs(det) < IFT_EPSILON) {
        printf("%f\n", fabs(det));
        iftError("Matrix is not invertible", "iftInvertMatrix");
    } else {
        B = iftPseudoInvertMatrix(A);
    }

    return (B);
}

iftMatrix *iftInvertDiagonalMatrix(const iftMatrix *A) {
    if (A->ncols != A->nrows)
        iftError("Matrix is not square", "iftInvertDiagonalMatrix");

    /* validate the diagonal matrix and compute the inverse */
    iftMatrix *B = iftCreateMatrix(A->ncols, A->nrows);
    for(int c = 0; c < A->ncols; c++) {
        for(int r = 0; r < A->nrows; r++) {
            if(c==r) {
                iftMatrixElem(B, c, r) = 1.0/iftMatrixElem(A, c, r);
            }
            else {
                if(iftMatrixElem(A, c, r) != 0)
                    iftError("Matrix is not a diagonal matrix", "iftInvertDiagonalMatrix");
            }
        }
    }

    return (B);
}

iftMatrix *iftIdentityMatrix(int ncols) {
    iftMatrix *i = iftCreateMatrix(ncols, ncols);
    int idx;

    for (idx = 0; idx < ncols; idx++)
        iftMatrixElem(i, idx, idx) = 1.0f;

    return (i);
}

iftMatrix *iftPseudoInvertMatrix(const iftMatrix *A) {

    if (A->ncols != A->nrows)
        iftError("Matrix is not square", "iftPseudoInvertMatrix");

    iftMatrix *AI = iftCopyMatrix(A);

    int n = AI->nrows;
    int *ipiv = iftAllocIntArray(n + 1);
    int lwork = n * n;
    float *work = iftAllocFloatArray(lwork);
    int info;

    sgetrf_(&n, &n, AI->val, &n, ipiv, &info);
    sgetri_(&n, AI->val, &n, ipiv, work, &lwork, &info);

    iftFree(ipiv);
    iftFree(work);

    if(info>0) {
        iftWarning("Matrix is not invertable", "iftPseudoInvertMatrix");
    }

    return AI;
}


iftDoubleMatrix *iftPseudoInvertDMatrix(const iftDoubleMatrix *A) {

    if (A->ncols != A->nrows)
        iftError("Matrix is not square", "iftPseudoInvertMatrix");

    iftDoubleMatrix *AI = iftCopyDoubleMatrix(A);

    int n = AI->nrows;
    int *ipiv = iftAllocIntArray(n + 1);
    int lwork = n * n;
    double *work = iftAllocDoubleArray(lwork);
    int info;

    dgetrf_(&n, &n, AI->val, &n, ipiv, &info);
    dgetri_(&n, AI->val, &n, ipiv, work, &lwork, &info);

    iftFree(ipiv);
    iftFree(work);

    if(info>0) {
        iftWarning("Matrix is not invertable", "iftPseudoInvertMatrix");
    }

    return AI;
}


/* This function aims at aligning a vector v with the direction of the
   axis Z. First, it rotates v by theta_x around the axis X until it
   falls on the xz plane. The value of theta_x can be obtained from
   the projection v_{yz} of v onto the yz plane as atan(vy/vz). It
   then rotates v (already on the XZ plane) of -theta_y until it
   aligns with the axis Z. The value of theta_y can be obtained from
   the xz plane as atan(vx/vz'), where vz' = vz /
   cos(theta_x). Special cases when vz=0 are treated as well. */

iftMatrix *iftRotationMatrixToAlignVectorWithZ(iftVector v) {
    float theta_x, theta_y, theta_z, m = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    iftMatrix *Rx = NULL, *Ry = NULL, *Rz = NULL, *M = NULL;

    if (m == 0.0)
        iftError("There is no vector", "iftRotationMatrixToAlignVectorWithZ");

    v.x /= m;
    v.y /= m;
    v.z /= m;

    if (iftAlmostZero(v.z)) {
        if (iftAlmostZero(v.x)) {
            if (v.y > 0.0)
                theta_x = 90.0;
            else
                theta_x = -90.0;

            printf("theta_x %f \n", theta_x);

            M = iftRotationMatrix(IFT_AXIS_X, theta_x);
        } else {
            if (iftAlmostZero(v.y)) {
                if (v.x > 0.0)
                    theta_y = -90.0;
                else
                    theta_y = 90.0;

                printf("theta_y %f \n", theta_y);

                M = iftRotationMatrix(IFT_AXIS_Y, theta_y);
            } else { /* v is on the xy plane */
                if (v.y > 0.0)
                    theta_z = -acosf(v.x) * 180.0 / IFT_PI;
                else
                    theta_z = +acosf(v.x) * 180.0 / IFT_PI;

                theta_y = -90.0;
                Rz = iftRotationMatrix(IFT_AXIS_Z, theta_z);
                Ry = iftRotationMatrix(IFT_AXIS_Y, theta_y);
                printf("theta_z %f theta_y %f \n", theta_z, theta_y);
                M = iftMultMatrices(Ry, Rz);
                iftDestroyMatrix(&Ry);
                iftDestroyMatrix(&Rz);
            }
        }
    } else {
        theta_x = atanf(v.y / v.z);
        float vz1 = v.z / cosf(theta_x);
        theta_x = theta_x * 180.0 / IFT_PI;
        theta_y = -atanf(v.x / vz1) * 180.0 / IFT_PI;

        if (v.z < 0) {
            if (!iftAlmostZero(theta_x))
                theta_x = theta_x - 180;
            if (!iftAlmostZero(theta_y))
                theta_y = theta_y - 180;
        }

        printf("theta_x %f theta_y %f \n", theta_x, theta_y);

        Rx = iftRotationMatrix(IFT_AXIS_X, theta_x);
        Ry = iftRotationMatrix(IFT_AXIS_Y, theta_y);
        M = iftMultMatrices(Ry, Rx);
        iftDestroyMatrix(&Ry);
        iftDestroyMatrix(&Rx);
    }

    return (M);
}

/* This function computes the alignment of a vector with the z axis
   and applies the inverse. */

iftMatrix *iftRotationMatrixToAlignZWithVector(iftVector v) {
    iftMatrix *M1 = iftRotationMatrixToAlignVectorWithZ(v);
    iftMatrix *M2 = iftTransposeMatrix(M1);

    iftDestroyMatrix(&M1);

    return (M2);
}


iftMatrix *iftExtendMatrix(const iftMatrix *M, int ncols, int nrows) {
    iftMatrix *E = iftCreateMatrix(ncols, nrows);
    int r, c;

    if ((ncols < M->ncols) && (nrows < M->nrows))
        iftError("Cannot extend the dimensions of the matrix", "iftExtendedMatrix");

    for (r = 0; r < E->nrows; r++)
        for (c = 0; c < E->ncols; c++) {
            if ((r < M->nrows) && (c < M->ncols))
                E->val[iftGetMatrixIndex(E, c, r)] = M->val[iftGetMatrixIndex(M, c, r)];
            else { // set 1 to the main diagonal
                if (c == r)
                    E->val[iftGetMatrixIndex(E, c, r)] = 1.0;
            }
        }

    return (E);

}

void iftResizeMatrix(iftMatrix **M, int ncols, int nrows) {
    if((ncols == 0) || (nrows ==0))
        iftError("The new quantity of columns and rows must be greater than 0", "iftResizeMatrix");
    
    if((*M) == NULL) {
        (*M) = iftCreateMatrix(ncols, nrows);
    } else {
        if((ncols != (*M)->ncols) || (nrows != (*M)->nrows)) {
            iftMatrix *M1 = iftCopyMatrix(*M);
            iftDestroyMatrix(M);
            (*M) = iftCreateMatrix(ncols, nrows);

            for(int r = 0; r < iftMin(M1->nrows, (*M)->nrows); r++) {
                for(int c = 0; c < iftMin(M1->ncols, (*M)->ncols); c++) {
                    iftMatrixElem((*M), c, r) = iftMatrixElem(M1, c, r);
                }
            }
            iftDestroyMatrix(&M1);
        }
    }
}

void iftResizeDoubleMatrix(iftDoubleMatrix **M, int ncols, int nrows) {
    if((ncols == 0) || (nrows ==0))
        iftError("The new quantity of columns and rows must be greater than 0", "iftResizeDoubleMatrix");
    
    if((*M) == NULL) {
        (*M) = iftCreateDoubleMatrix(ncols, nrows);
    } else {
        if((ncols != (*M)->ncols) || (nrows != (*M)->nrows)) {
            iftDoubleMatrix *M1 = iftCopyDoubleMatrix(*M);
            iftDestroyDoubleMatrix(M);
            (*M) = iftCreateDoubleMatrix(ncols, nrows);

            for(int r = 0; r < iftMin(M1->nrows, (*M)->nrows); r++) {
                for(int c = 0; c < iftMin(M1->ncols, (*M)->ncols); c++) {
                    iftMatrixElem((*M), c, r) = iftMatrixElem(M1, c, r);
                }
            }
            iftDestroyDoubleMatrix(&M1);
        }
    }
}


void iftResizeStrMatrix(iftStrMatrix **M, int ncols, int nrows) {
    if((ncols == 0) || (nrows ==0))
        iftError("The new quantity of columns and rows must be greater than 0", "iftResizeStrMatrix");
    
    if((*M) == NULL) {
        (*M) = iftCreateStrMatrix(ncols, nrows);
    } else {
        if((ncols != (*M)->ncols) || (nrows != (*M)->nrows)) {
            iftStrMatrix *M1 = iftCopyStrMatrix((*M)->val, (*M)->nrows, (*M)->ncols);
            iftDestroyStrMatrix(M);
            (*M) = iftCreateStrMatrix(ncols, nrows);

            for(int r = 0; r < iftMin(M1->nrows, (*M)->nrows); r++) {
                for(int c = 0; c < iftMin(M1->ncols, (*M)->ncols); c++) {
                    iftFree(iftMatrixElem((*M), c, r));
                    iftMatrixElem((*M), c, r) = iftCopyString(iftMatrixElem(M1, c, r));
                }
            }
            iftDestroyStrMatrix(&M1);
        }
    }
}

/*
 * Small documentations in:
 * http://dlib.net/dlib/matrix/lapack/gesdd.h.html
 * http://www.netlib.org/lapack/lapack-3.1.1/html/sgesdd.f.html
 *
 * Reminder: In U the eigenvectors are the columns and in Vt the
 * eigenvectors are the rows
 */
void iftSingleValueDecomp(iftMatrix *A, iftMatrix **U, iftMatrix **S, iftMatrix **Vt) {
    int m = A->nrows, n = A->ncols, lda = m, ldu = m;
    int k, D = iftMin(m, n);
    int ldvt = D, info, lwork, *iwork;
    float *work, wkopt;
    float *s, *u, *vt;
    float *a;

    /* Alloc work memory space */

    iwork = iftAllocIntArray(8 * D);
    s = iftAllocFloatArray(D);
    u = iftAllocFloatArray(ldu * D);
    vt = iftAllocFloatArray(ldvt * n);
    a = iftAllocFloatArray(lda * n);

    /* Load input matrix A (Note that Fortran saves array information
     column by column) */

    k = 0;
    for (int j = 0; j < A->ncols; j++)
      for (int i = 0; i < A->nrows; i++)
            a[k++] = A->val[iftGetMatrixIndex(A, j, i)];

    /* Query and allocate the optimal workspace */

    lwork = IFT_NIL;
    sgesdd_("Singular vectors", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt,
            &lwork, iwork, &info);
    lwork = (int) wkopt;

	if(info>0)
		iftError("SVD did not converge.", "iftSingleValueDecomp");


    // minimum space required for work array - check the documentation
    int min_space = 3 * D * D + iftMax(iftMax(m, n), (4 * D * D + 4 * D));
    if (lwork < min_space)
        lwork = min_space;

    work = iftAllocFloatArray(lwork);

    /* Compute SVD */

    sgesdd_("Singular vectors", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work,
            &lwork, iwork, &info);

    iftFree(work);
    iftFree(iwork);
    iftFree(a);


    /* Check for convergence */
    if (info > 0) {
        iftError("The algorithm computing SVD failed to converge", "iftSingleValueDecomp");
    }

    /* Get matrices after decomposition */
    *U = iftCreateMatrix(D, m);
    *S = iftCreateMatrix(D, 1);
    *Vt = iftCreateMatrix(n, D);

    k = 0;
    for (int j = 0; j < (*U)->ncols; j++)
        for (int i = 0; i < (*U)->nrows; i++)
            (*U)->val[iftGetMatrixIndex((*U), i, j)] = u[k++];

    for (k = 0; k < (*S)->n; k++)
        (*S)->val[k] = s[k];

    k = 0;
    for (int j = 0; j < (*Vt)->ncols; j++)
        for (int i = 0; i < (*Vt)->nrows; i++)
            (*Vt)->val[iftGetMatrixIndex((*Vt), j, i)] = vt[k++];

    iftFree(s);
    iftFree(u);
    iftFree(vt);
}

iftMatrix *iftGeneralizedEigenvalues(iftMatrix *A)
{
    iftMatrix *eigenvalues = NULL;

    if (iftIsMatrixSymmetric(A)){

        int n = A->ncols, lda = A->nrows, k, lwork, info; // matrix order
        double wkopt;
        double w[n];
        double *work = NULL;
        double a[lda*n];

        /* Load input matrix A (Note that Fortran saves array information
     column by column) */
        k = 0;
        for (int j = 0; j < A->ncols; j++)
            for (int i = 0; i < A->nrows; i++)
                a[k++] = A->val[iftGetMatrixIndex(A, j, i)];

        /* Query and allocate the optimal workspace */
        lwork = -1;
        dsyev_( "Vectors", "Upper", &n, a, &lda, w, &wkopt, &lwork, &info );
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        /* Solve eigenproblem */
        dsyev_( "Vectors", "Upper", &n, a, &lda, w, work, &lwork, &info );
        /* Check for convergence */

        //dsyev_("N","L",&n,a,&lda,w,work,&lwork,&info);

        if (info < 0) {
            iftError("Parameter %d is invalid.","iftGeneralizedEigenvalues",info);
        }

        if (info > 0){
            iftError("Eigenvalue convergence failed.","iftGeneralizedEigenvalues");
        }
        iftFree(work);

        eigenvalues = iftCreateMatrix(1, n);

        for (k = 0; k < n; k++)
            eigenvalues->val[k] = w[k];

    }

    return eigenvalues;
}

bool iftIsMatrixSymmetric(iftMatrix *A)
{
    iftMatrix *At = iftTransposeMatrix(A);
    if (iftIsMatricesEqual(A,At))
        return TRUE;
    else
        return FALSE;
}


void iftNormalizeMatrix(iftMatrix *M) {
    double norm = iftFrobeniusNorm(M);

    if (norm > 0.0) {
        iftMultMatrixByScalar(M, 1.0/norm);
    }
}

void iftNormalizeMatrixInRange(iftMatrix *M, float minVal, float maxVal)
{
    float curMinVal = iftGetMinimumValueInMatrix(M);
    float curMaxVal = iftGetMaximumValueInMatrix(M);

    for(int i = 0; i < M->n; i++)
        M->val[i] = (maxVal - minVal) * (M->val[i] - curMinVal) / (curMaxVal - curMinVal) + minVal;
}

double iftFrobeniusNorm(iftMatrix *M) {
    int i;
    double norm = 0.0;

    for (i = 0; i < M->n; i++)
        norm += M->val[i] * M->val[i];

    norm = sqrt(norm);
    return (norm);
}

void iftMultMatrixByScalar(iftMatrix *A, double scalar) {
    for (int i = 0; i < A->n; i++)
        A->val[i] *= scalar;
}

/**
@brief Computes the multiplication of N matrices.\n

Dependency: #include <stdarg.h>.\n
This function computes the product of the matrices given as a parameter from left to right.\n
Example: iftComputeTransformation(5, A, B, C, D, E);\n
Will compute the product in this order: ((((B * A) * C) * D) * E)\n


@param  n_args                Integer. Number of matrices.
@param  ...                   iftMatrix. Pointers to matrices.
@return A iftMatrix resultant from the product of all matrices passed as parameters.
*/
iftMatrix *iftComputeTransformation(int n_args, ...) {
    int i;
    iftMatrix *A = NULL, *B = NULL, *T = NULL, *Aux = NULL;

    if (n_args <= 1)
        iftError("Need at least 2 matrices", "iftComputeTransformation");

    va_list ap;
    va_start(ap, n_args);

    A = va_arg(ap, iftMatrix*);
    B = va_arg(ap, iftMatrix*);
    Aux = iftMultMatrices(B, A);
    for (i = 2; i < n_args; i++) {
        A = va_arg(ap, iftMatrix*);
        T = iftMultMatrices(A, Aux);
        iftDestroyMatrix(&Aux);
        Aux = T;
    }
    va_end(ap);
    if (n_args == 2)
        return Aux;
    return T;
}



char *iftMatrixAsString(const iftMatrix *mat) {
    char *str = iftCopyString("[");
    char *aux = NULL;

    for (long i = 0; i <= (mat->nrows-2); i++) {
        // 2 chars for square brackets + 16 chars per number (at most) + n-1 commas
        char *str_row = iftAlloc(2 + (16 * mat->ncols) + (mat->ncols-1), sizeof(char));
	
        sprintf(str_row, "[");
        for (long j = 0; j <= (mat->ncols-2); j++) {
	  sprintf(str_row, "%s%f, ", str_row, iftMatrixElem(mat, j, i));
	  
        }
        sprintf(str_row, "%s%f]", str_row, iftMatrixElem(mat, mat->ncols-1, i));

        aux = str;
        str = iftConcatStrings(3, str, str_row, ", ");
        iftFree(aux);
    }

    char *str_row = iftAlloc(2 + (16 * mat->ncols) + (mat->ncols-1), sizeof(char));
    sprintf(str_row, "[");
    for (long j = 0; j <= (mat->ncols-2); j++) {
        sprintf(str_row, "%s%f, ", str_row, iftMatrixElem(mat, j, mat->nrows-1));
    }
    sprintf(str_row, "%s%f]", str_row, iftMatrixElem(mat, mat->ncols-1, mat->nrows-1));

    aux = str;
    str = iftConcatStrings(3, str, str_row, "]");
    iftFree(aux);
    
    return str;
}


char *iftStrMatrixAsString(const iftStrMatrix *mat) {
    char *str = iftCopyString("[");
    char *aux = NULL;

    for (long i = 0; i <= (mat->nrows-2); i++) {
        // 2 chars for square brackets + 2050 chars per number (at most) + n-1 commas
        char *str_row = iftAlloc(2 + (2050 * mat->ncols) + (mat->ncols-1), sizeof(char));        

        sprintf(str_row, "[");
        for (long j = 0; j <= (mat->ncols-2); j++) {
            sprintf(str_row, "%s\"%s\", ", str_row, iftMatrixElem(mat, j, i));
        }
        sprintf(str_row, "%s\"%s\"]", str_row, iftMatrixElem(mat, mat->ncols-1, i));

        aux = str;
        str = iftConcatStrings(3, str, str_row, ", ");
        iftFree(aux);
    }

    char *str_row = iftAlloc(2 + (2050 * mat->ncols) + (mat->ncols-1), sizeof(char));        
    sprintf(str_row, "[");
    for (long j = 0; j <= (mat->ncols-2); j++) {
        sprintf(str_row, "%s\"%s\", ", str_row, iftMatrixElem(mat, j, mat->nrows-1));
    }
    sprintf(str_row, "%s\"%s\"]", str_row, iftMatrixElem(mat, mat->ncols-1, mat->nrows-1));

    aux = str;
    str = iftConcatStrings(3, str, str_row, "]");
    iftFree(aux);

    return str;
}


iftStrMatrix *iftCreateStrMatrix(int ncols, int nrows) {
    iftStrMatrix *sM = (iftStrMatrix *) iftAlloc(1, sizeof(iftStrMatrix));

    sM->ncols = ncols;
    sM->nrows = nrows;

    sM->tbrow = iftAllocIntArray(nrows);
    for (int r = 0; r < nrows; r++) {
        sM->tbrow[r] = r * ncols;
    }

    sM->n = ncols * nrows;
    sM->val = (char **) iftAlloc(sM->n, sizeof(char *));
    for (int i = 0; i < sM->n; i++)
        sM->val[i] = iftAllocCharArray(1024);

    return sM;
}


void iftDestroyStrMatrix(iftStrMatrix **sM) {
    iftStrMatrix *aux = *sM;

    if (aux != NULL) {
        iftFree(aux->tbrow);
        if (aux->val != NULL) {
            for (int i = 0; i < aux->n; i++)
                iftFree(aux->val[i]);
            iftFree(aux->val);
        }
        iftFree(aux);
        *sM = NULL;
    }
}


iftStrMatrix *iftCopyStrMatrix(char **val, int nrows, int ncols) {
    iftStrMatrix *mat = iftCreateStrMatrix(ncols, nrows);

    for (int i = 0; i < mat->n; i++)
        mat->val[i] = iftCopyString(val[i]);

    return mat;
}


iftBoolMatrix *iftCreateBoolMatrix(int ncols, int nrows) {
    iftBoolMatrix *M = (iftBoolMatrix *) iftAlloc(1, sizeof(iftBoolMatrix));

    M->ncols = ncols;
    M->nrows = nrows;

    M->tbrow = iftAllocIntArray(nrows);
    for (int r = 0; r < nrows; r++) {
        M->tbrow[r] = r * ncols;
    }

    M->n = ncols * nrows;
    M->val = iftAllocBoolArray(M->n);

    return M;
}


void iftDestroyBoolMatrix(iftBoolMatrix **M) {
    iftBoolMatrix *aux = *M;

    if (aux != NULL) {
        if (aux->val != NULL)
            iftFree(aux->val);
        iftFree(aux->tbrow);
        iftFree(aux);
        *M = NULL;
    }
}


iftBoolMatrix *iftCopyBoolMatrix(bool *val, int nrows, int ncols) {
    iftBoolMatrix *mat = iftCreateBoolMatrix(ncols, nrows);

    for (int i = 0; i < mat->n; i++)
        mat->val[i] = val[i];

    return mat;
}


iftVoxel *iftTransformVoxels(iftMatrix *S, iftVoxel *points, int n) {
    iftVoxel *rescaled_points;

    rescaled_points = (iftVoxel *) iftAlloc(n, sizeof(iftVoxel));

    for (int i = 0; i < n; i++) {
        rescaled_points[i] = iftTransformVoxel(S, points[i]);
    }

    return rescaled_points;
}


iftMatrix *iftLeastSquares(const iftMatrix *A, const iftMatrix *B) {

    if (A->nrows != B->nrows) {
        iftError("A and B must have the same number of rows.", "iftLeatSquares");
    }

    iftMatrix *At = iftTransposeMatrix(A);
    iftMatrix *AtA = iftMultMatrices(At, A);
    iftMatrix *AtB = iftMultMatrices(At, B);


    iftMatrix *AtAi = iftPseudoInvertMatrix(AtA);

    iftMatrix *coef = iftMultMatrices(AtAi, AtB);

    iftDestroyMatrix(&At);
    iftDestroyMatrix(&AtA);
    iftDestroyMatrix(&AtAi);
    iftDestroyMatrix(&AtB);

    return coef;
}

double iftGetMinimumValueInMatrix(iftMatrix *matrix) {
    unsigned int i, N;
    double minimumValue = matrix->val[0];
    N = matrix->ncols * matrix->nrows;
    for (i = 0; i < N; i++) {
        if (matrix->val[i] < minimumValue) {
            minimumValue = matrix->val[i];
        }
    }
    return minimumValue;
}

double iftGetMaximumValueInMatrix(iftMatrix *matrix) {
    unsigned int i, N;
    double maximumValue = matrix->val[0];
    N = matrix->ncols * matrix->nrows;
    for (i = 0; i < N; i++) {
        if (matrix->val[i] > maximumValue) {
            maximumValue = matrix->val[i];
        }
    }
    return maximumValue;
}


iftMatrix *iftMatrixMeanColumn(const iftMatrix *matrix) {
    iftMatrix *means = iftCreateMatrix(matrix->ncols, 1);
    int i, j, k;
    k = 0;

    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            means->val[j] += matrix->val[k];
            k++;
        }
    }

    for (j = 0; j < matrix->ncols; j++)
        means->val[j] /= matrix->nrows;

    return means;
}

iftMatrix *iftMatrixSumColumn(const iftMatrix *matrix) {
    iftMatrix *sum = iftCreateMatrix(matrix->ncols, 1);

    for (int i = 0; i < matrix->nrows; i++) {
        for (int j = 0; j < matrix->ncols; j++) {
            sum->val[j] += iftMatrixElem(matrix, j, i);
        }
    }

    return sum;
}

iftMatrix *iftMatrixSumRow(const iftMatrix *matrix) {
    iftMatrix *sum = iftCreateMatrix(1, matrix->nrows);

    for (int i = 0; i < matrix->nrows; i++) {
        for (int j = 0; j < matrix->ncols; j++) {
            sum->val[i] += iftMatrixElem(matrix, j, i);
        }
    }

    return sum;
}

iftMatrix *iftMatrixMaxColumn(const iftMatrix *matrix) {
    iftMatrix *max = iftCreateMatrix(matrix->ncols, 1);

    for (int j = 0; j < matrix->ncols; j++) {
        float maxVal = IFT_INFINITY_FLT_NEG;
        for (int i = 0; i < matrix->nrows; i++) {
            if(iftMatrixElem(matrix, j, i) > maxVal)
                maxVal = iftMatrixElem(matrix, j, i);
        }
        max->val[j] = maxVal;
    }

    return max;
}

iftMatrix *iftMatrixMaxRow(const iftMatrix *matrix) {
    iftMatrix *max = iftCreateMatrix(matrix->nrows, 1);

    for (int i = 0; i < matrix->nrows; i++) {
        float maxVal = IFT_INFINITY_FLT_NEG;
        for (int j = 0; j < matrix->ncols; j++) {
            if(iftMatrixElem(matrix, j, i) > maxVal)
                maxVal = iftMatrixElem(matrix, j, i);
        }
        max->val[i] = maxVal;
    }

    return max;
}

float iftMatrixMaxOneColumn(const iftMatrix *matrix, int c) {

    float maxVal = IFT_INFINITY_FLT_NEG;
    for (int i = 0; i < matrix->nrows; i++) {
        if(iftMatrixElem(matrix, c, i) > maxVal)
            maxVal = iftMatrixElem(matrix, c, i);
    }

    return maxVal;
}

float iftMatrixMaxOneRow(const iftMatrix *matrix, int r) {

    float maxVal = IFT_INFINITY_FLT_NEG;
    for (int j = 0; j < matrix->ncols; j++) {
        if(iftMatrixElem(matrix, j, r) > maxVal)
            maxVal = iftMatrixElem(matrix, j, r);
    }

    return maxVal;
}

void iftComputeOperationBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order, char operationSymbol) {
    unsigned int i, N;
    N = matrix->ncols * matrix->nrows;
    if (order == 'f') {
        switch (operationSymbol) {
            case '+':
                for (i = 0; i < N; i++) {
                    matrix->val[i] = matrix->val[i] + scalar;
                }
                break;
            case '-':
                for (i = 0; i < N; i++) {
                    matrix->val[i] = matrix->val[i] - scalar;
                }
                break;
            case '*':
                for (i = 0; i < N; i++) {
                    matrix->val[i] = matrix->val[i] * scalar;
                }
                break;
            case '/':
                for (i = 0; i < N; i++) {
                    matrix->val[i] = matrix->val[i] / scalar;
                }
                break;
            default:
                printf("Unrecognizable operation symbol %c", operationSymbol);
                break;
        }
    } else if (order == 'b') {
        switch (operationSymbol) {
            case '+':
                for (i = 0; i < N; i++) {
                    matrix->val[i] = scalar + matrix->val[i];
                }
                break;
            case '-':
                for (i = 0; i < N; i++) {
                    matrix->val[i] = scalar - matrix->val[i];
                }
                break;
            case '*':
                for (i = 0; i < N; i++) {
                    matrix->val[i] = scalar * matrix->val[i];
                }
                break;
            case '/':
                for (i = 0; i < N; i++) {
                    matrix->val[i] = scalar / matrix->val[i];
                }
                break;
            default:
                printf("Unrecognizable operation symbol %c", operationSymbol);
                break;
        }
    } else {
        printf("Unrecognizable order %c", order);
    }
}

void iftComputeAdditionBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order) {
    unsigned int i, N;
    N = matrix->ncols * matrix->nrows;
    if (order == 'f') {
        for (i = 0; i < N; i++) {
            matrix->val[i] = matrix->val[i] + scalar;
        }
    } else if (order == 'b') {
        for (i = 0; i < N; i++) {
            matrix->val[i] = scalar + matrix->val[i];
        }
    } else {
        printf("Unrecognizable order %c", order);
    }
}

void iftComputeSubtractionBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order) {
    unsigned int i, N;
    N = matrix->ncols * matrix->nrows;
    if (order == 'f') {
        for (i = 0; i < N; i++) {
            matrix->val[i] = matrix->val[i] - scalar;
        }
    } else if (order == 'b') {
        for (i = 0; i < N; i++) {
            matrix->val[i] = scalar - matrix->val[i];
        }
    } else {
        printf("Unrecognizable order %c", order);
    }
}

void iftComputeMultiplicationBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order) {
    unsigned int i, N;
    N = matrix->ncols * matrix->nrows;
    if (order == 'f') {
        for (i = 0; i < N; i++) {
            matrix->val[i] = matrix->val[i] * scalar;
        }
    } else if (order == 'b') {
        for (i = 0; i < N; i++) {
            matrix->val[i] = scalar * matrix->val[i];
        }
    } else {
        printf("Unrecognizable order %c", order);
    }
}

void iftComputeDivisionBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order) {
    unsigned int i, N;
    N = matrix->ncols * matrix->nrows;
    if (order == 'f') {
        for (i = 0; i < N; i++) {
            matrix->val[i] = matrix->val[i] / scalar;
        }
    } else if (order == 'b') {
        for (i = 0; i < N; i++) {
            matrix->val[i] = scalar / matrix->val[i];
        }
    } else {
        printf("Unrecognizable order %c", order);
    }
}



double iftMatrixSum(const iftMatrix *matrix) {
    long i;
    double sum = 0.0;
    for (i = 0; i < matrix->n; i++) {
        sum += matrix->val[i];
    }
    return sum;
}

iftMatrix *iftCovarianceMatrix(const iftMatrix *X) {
    iftMatrix *cov, *Xt;

    float K = 1.0 / ((float) X->nrows - 1);

    iftMatrix *Xc = iftCopyMatrix(X);
    iftMatrix *mean = iftMatrixMeanColumn(X);

    for (int r = 0; r < X->nrows; ++r)
        for (int c = 0; c < X->ncols; ++c)
            iftMatrixElem(Xc, c, r) -= mean->val[c];

    Xt = iftTransposeMatrix(Xc);
    iftMultMatrixByScalar(Xc, K);
    cov = iftMultMatrices(Xt, Xc);

    iftDestroyMatrix(&Xc);
    iftDestroyMatrix(&Xt);
    iftDestroyMatrix(&mean);

    return (cov);
}

void iftSetDiagonalValue(iftMatrix *matrix, double val) {
    for (int i = 0; i < matrix->nrows; ++i) {
        iftMatrixElem(matrix, i, i) = val;
    }
}

void iftSumMatricesInPlace(const iftMatrix *m1, const iftMatrix *m2, iftMatrix *out) {

    if (m1->nrows != m2->nrows || m1->nrows != out->nrows ||
        m1->ncols != m2->ncols || m1->ncols != out->ncols) {
        iftError("Matrices must have the same dimensions.", "iftSumMatrices");
    }

    for (int i = 0; i < out->n; ++i) {
        out->val[i] = m1->val[i] + m2->val[i];
    }
}

iftMatrix *iftSumMatrices(const iftMatrix *m1, const iftMatrix *m2) {

    if (m1->nrows != m2->nrows ||
        m1->ncols != m2->ncols) {
        iftError("Matrices must have the same dimensions.", "iftSumMatrices");
    }

    iftMatrix *sum = iftCreateMatrix(m1->ncols, m2->nrows);
    iftSumMatricesInPlace(m1, m2, sum);
    return sum;
}

void iftSubtractMatricesInPlace(const iftMatrix *m1, const iftMatrix *m2, iftMatrix *out) {

    if (m1->nrows != m2->nrows || m1->nrows != out->nrows ||
        m1->ncols != m2->ncols || m1->ncols != out->ncols) {
        iftError("Matrices must have the same dimensions.", "iftSumMatrices");
    }

    for (int i = 0; i < out->n; ++i) {
        out->val[i] = m1->val[i] - m2->val[i];
    }
}

iftMatrix *iftSubtractMatrices(const iftMatrix *m1, const iftMatrix *m2) {

    if (m1->nrows != m2->nrows ||
        m1->ncols != m2->ncols) {
        iftError("Matrices must have the same dimensions.", "iftSumMatrices");
    }

    iftMatrix *subtract = iftCreateMatrix(m1->ncols, m2->nrows);
    iftSubtractMatricesInPlace(m1, m2, subtract);
    return subtract;
}

void
iftComputeOperationBetweenMatrixVector(iftMatrix *matrix, iftDblArray *vector, char vectorType, char operationSymbol) {
    unsigned int i, j, k;
    k = 0;
    if (vectorType == 'r') {
        if (matrix->ncols != vector->n) {
            printf("Could not perform the operation between the vector and matrix");
            return;
        }
        switch (operationSymbol) {
            case '+':
                for (i = 0; i < matrix->nrows; i++) {
                    for (j = 0; j < matrix->ncols; j++) {
                        matrix->val[k] += vector->val[j];
                        k++;
                    }
                }
                break;
            case '-':
                for (i = 0; i < matrix->nrows; i++) {
                    for (j = 0; j < matrix->ncols; j++) {
                        matrix->val[k] -= vector->val[j];
                        k++;
                    }
                }
                break;
            case '*':
                for (i = 0; i < matrix->nrows; i++) {
                    for (j = 0; j < matrix->ncols; j++) {
                        matrix->val[k] *= vector->val[j];
                        k++;
                    }
                }
                break;
            case '/':
                for (i = 0; i < matrix->nrows; i++) {
                    for (j = 0; j < matrix->ncols; j++) {
                        matrix->val[k] /= vector->val[j];
                        k++;
                    }
                }
                break;
            default:
                printf("Unrecognizable operation symbol %c", operationSymbol);
                break;
        }
    } else if (vectorType == 'c') {
        if (matrix->nrows != vector->n) {
            printf("Could not perform the operation between the vector and matrix");
            return;
        }
        switch (operationSymbol) {
            case '+':
                for (i = 0; i < matrix->nrows; i++) {
                    for (j = 0; j < matrix->ncols; j++) {
                        matrix->val[k] += vector->val[i];
                        k++;
                    }
                }
                break;
            case '-':
                for (i = 0; i < matrix->nrows; i++) {
                    for (j = 0; j < matrix->ncols; j++) {
                        matrix->val[k] -= vector->val[i];
                        k++;
                    }
                }
                break;
            case '*':
                for (i = 0; i < matrix->nrows; i++) {
                    for (j = 0; j < matrix->ncols; j++) {
                        matrix->val[k] *= vector->val[i];
                        k++;
                    }
                }
                break;
            case '/':
                for (i = 0; i < matrix->nrows; i++) {
                    for (j = 0; j < matrix->ncols; j++) {
                        matrix->val[k] /= vector->val[i];
                        k++;
                    }
                }
                break;
            default:
                printf("Unrecognizable operation symbol %c", operationSymbol);
                break;
        }

    } else {
        printf("Unrecognizable vector type %c", vectorType);
    }
}

void iftAdditionMatrixVectorByColumn(iftMatrix *matrix, iftMatrix *vector) {
    if (matrix->nrows != vector->nrows || vector->ncols > 1) {
        iftError("Could not perform the operation between the vector and matrix. Dimensions must agree",
                 "iftAdditionMatrixVectorByColumn");
    }
    unsigned int i, j, k;
    k = 0;
    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            matrix->val[k] += vector->val[i];
            k++;
        }
    }
}

void iftAdditionMatrixVectorByRow(iftMatrix *matrix, iftMatrix *vector) {
    if (matrix->ncols != vector->ncols || vector->nrows > 1) {
        iftError("Could not perform the operation between the vector and matrix. Dimensions must agree",
                 "iftAdditionMatrixVectorByRow");
    }
    unsigned int i, j, k;
    k = 0;
    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            matrix->val[k] += vector->val[j];
            k++;
        }
    }
}

void iftSubtractionMatrixVectorByColumn(iftMatrix *matrix, iftMatrix *vector) {
    if (matrix->nrows != vector->nrows || vector->ncols > 1) {
        iftError("Could not perform the operation between the vector and matrix. Dimensions must agree",
                 "iftSubtractionMatrixVectorByColumn");
    }
    unsigned int i, j, k;
    k = 0;
    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            matrix->val[k] -= vector->val[i];
            k++;
        }
    }
}

void iftSubtractionMatrixVectorByRow(iftMatrix *matrix, iftMatrix *vector) {
    if (matrix->ncols != vector->ncols || vector->nrows > 1) {
        iftError("Could not perform the operation between the vector and matrix. Dimensions must agree",
                 "iftSubtractionMatrixVectorByRow");
    }
    unsigned int i, j, k;
    k = 0;
    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            matrix->val[k] -= vector->val[j];
            k++;
        }
    }
}

void iftMultiplicationMatrixVectorByColumn(iftMatrix *matrix, iftMatrix *vector) {
    if (matrix->nrows != vector->nrows || vector->ncols > 1) {
        iftError("Could not perform the operation between the vector and matrix. Dimensions must agree",
                 "iftMultiplicationMatrixVectorByColumn");
    }
    unsigned int i, j, k;
    k = 0;
    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            matrix->val[k] *= vector->val[i];
            k++;
        }
    }
}

void iftMultiplicationMatrixVectorByRow(iftMatrix *matrix, iftMatrix *vector) {
    if (matrix->ncols != vector->ncols || vector->nrows > 1) {
        iftError("Could not perform the operation between the vector and matrix. Dimensions must agree",
                 "iftMultiplicationMatrixVectorByRow");
    }
    unsigned int i, j, k;
    k = 0;
    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            matrix->val[k] *= vector->val[j];
            k++;
        }
    }
}

void iftDivisionMatrixVectorByColumn(iftMatrix *matrix, iftMatrix *vector) {
    if (matrix->nrows != vector->nrows || vector->ncols > 1) {
        iftError("Could not perform the operation between the vector and matrix. Dimensions must agree",
                 "iftDivisionMatrixVectorByColumn");
    }
    unsigned int i, j, k;
    k = 0;
    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            matrix->val[k] /= vector->val[i];
            k++;
        }
    }
}

void iftDivisionMatrixVectorByRow(iftMatrix *matrix, iftMatrix *vector) {
    if (matrix->ncols != vector->ncols || vector->nrows > 1) {
        iftError("Could not perform the operation between the vector and matrix. Dimensions must agree",
                 "iftDivisionMatrixVectorByRow");
    }
    unsigned int i, j, k;
    k = 0;
    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            matrix->val[k] /= vector->val[j];
            k++;
        }
    }
}

iftMatrix *iftComputeOperationBetweenMatricesInPointWise(iftMatrix *A, iftMatrix *B, char operationSymbol) {
    unsigned int i, N;
    if (A->ncols != B->ncols || A->nrows != B->nrows) {
        iftError("Could not perform the matrix operation. Matrices have different dimensions",
                 "iftComputeOperationBetweenMatricesInPointWise");
        return NULL;
    }
    iftMatrix *C = iftCreateMatrix(A->ncols, A->nrows);
    N = A->ncols * A->nrows;
    switch (operationSymbol) {
        case '+':
            for (i = 0; i < N; i++) {
                C->val[i] = A->val[i] + B->val[i];
            }
            break;
        case '-':
            for (i = 0; i < N; i++) {
                C->val[i] = A->val[i] - B->val[i];
            }
            break;
        case '*':
            for (i = 0; i < N; i++) {
                C->val[i] = A->val[i] * B->val[i];
            }
            break;
        case '/':
            for (i = 0; i < N; i++) {
                C->val[i] = A->val[i] / B->val[i];
            }
            break;
        default:
            printf("Unrecognizable operation symbol %c", operationSymbol);
            break;
    }
    return C;
}

iftMatrix *iftMatricesAdditionPointWise(iftMatrix *A, iftMatrix *B) {
    if (A->ncols != B->ncols || A->nrows != B->nrows) {
        iftError("Could not perform the matrix pointwise addition. Matrices have different dimensions",
                 "iftMatricesAdditionPointWise");
        return NULL;
    }
    iftMatrix *C = iftCreateMatrix(A->ncols, A->nrows);
    unsigned int i, N;
    N = A->ncols * A->nrows;
    for (i = 0; i < N; i++) {
        C->val[i] = A->val[i] + B->val[i];
    }
    return C;
}

iftMatrix *iftMatricesSubtractionPointWise(iftMatrix *A, iftMatrix *B) {
    if (A->ncols != B->ncols || A->nrows != B->nrows) {
        iftError("Could not perform the matrix pointwise subtraction. Matrices have different dimensions",
                 "iftMatricesSubtractionPointWise");
        return NULL;
    }
    iftMatrix *C = iftCreateMatrix(A->ncols, A->nrows);
    unsigned int i, N;
    N = A->ncols * A->nrows;
    for (i = 0; i < N; i++) {
        C->val[i] = A->val[i] - B->val[i];
    }
    return C;
}


iftMatrix *iftMatricesMultiplicationPointWise(iftMatrix *A, iftMatrix *B) {
    if (A->ncols != B->ncols || A->nrows != B->nrows) {
        iftError("Could not perform the matrix pointwise multiplication. Matrices have different dimensions",
                 "iftMatricesMultiplicationPointWise");
        return NULL;
    }
    iftMatrix *C = iftCreateMatrix(A->ncols, A->nrows);
    unsigned int i, N;
    N = A->ncols * A->nrows;
    for (i = 0; i < N; i++) {
        C->val[i] = A->val[i] * B->val[i];
    }
    return C;
}

iftMatrix *iftMatricesDivisionPointWise(iftMatrix *A, iftMatrix *B) {
    if (A->ncols != B->ncols || A->nrows != B->nrows) {
        iftError("Could not perform the matrix pointwise division. Matrices have different dimensions",
                 "iftMatricesMultiplicationPointWise");
        return NULL;
    }
    iftMatrix *C = iftCreateMatrix(A->ncols, A->nrows);
    unsigned int i, N;
    N = A->ncols * A->nrows;
    for (i = 0; i < N; i++) {
        C->val[i] = A->val[i] / B->val[i];
    }
    return C;
}

void iftMatricesAdditionPointWiseInPlace(iftMatrix *A, iftMatrix *B, iftMatrix **C) {

    if (A->ncols != B->ncols || A->nrows != B->nrows) {
        iftError("Could not perform the matrix pointwise addition. Matrices have different dimensions",
                 "iftMatricesAdditionPointWiseInPlace");
        return;
    }
    if (*C == NULL) {
        (*C) = iftCreateMatrix(A->ncols, A->nrows);
    }
    if ((A->ncols != (*C)->ncols || A->nrows != (*C)->nrows)) {
        iftDestroyMatrix(C);
        (*C) = iftCreateMatrix(A->ncols, A->nrows);
    }

    unsigned int i, N;
    N = A->ncols * A->nrows;
    for (i = 0; i < N; i++) {
        (*C)->val[i] = A->val[i] + B->val[i];
    }

}

void iftMatricesSubtractionPointWiseInPlace(iftMatrix *A, iftMatrix *B, iftMatrix **C) {
    if (A->ncols != B->ncols || A->nrows != B->nrows) {
        iftError("Could not perform the matrix pointwise subtraction. Matrices have different dimensions",
                 "iftMatricesSubtractionPointWiseInPlace");
        return;
    }
    if (*C == NULL) {
        (*C) = iftCreateMatrix(A->ncols, A->nrows);
    }
    if ((A->ncols != (*C)->ncols || A->nrows != (*C)->nrows)) {
        iftDestroyMatrix(C);
        (*C) = iftCreateMatrix(A->ncols, A->nrows);
    }

    unsigned int i, N;
    N = A->ncols * A->nrows;
    for (i = 0; i < N; i++) {
        (*C)->val[i] = A->val[i] - B->val[i];
    }
}

void iftMatricesMultiplicationPointWiseInPlace(iftMatrix *A, iftMatrix *B, iftMatrix **C) {
    if (A->ncols != B->ncols || A->nrows != B->nrows) {
        iftError("Could not perform the matrix pointwise multiplication. Matrices have different dimensions",
                 "iftMatricesMultiplicationPointWiseInPlace");
        return;
    }
    if (*C == NULL) {
        (*C) = iftCreateMatrix(A->ncols, A->nrows);
    }
    if ((A->ncols != (*C)->ncols || A->nrows != (*C)->nrows)) {
        iftDestroyMatrix(C);
        (*C) = iftCreateMatrix(A->ncols, A->nrows);
    }

    unsigned int i, N;
    N = A->ncols * A->nrows;
    for (i = 0; i < N; i++) {
        (*C)->val[i] = A->val[i] * B->val[i];
    }
}

void iftMatricesDivisionPointWiseInPlace(iftMatrix *A, iftMatrix *B, iftMatrix **C) {
    if (A->ncols != B->ncols || A->nrows != B->nrows) {
        iftError("Could not perform the matrix pointwise multiplication. Matrices have different dimensions",
                 "iftMatricesMultiplicationPointWiseInPlace");
        return;
    }
    if (*C == NULL) {
        (*C) = iftCreateMatrix(A->ncols, A->nrows);
    }
    if ((A->ncols != (*C)->ncols || A->nrows != (*C)->nrows)) {
        iftDestroyMatrix(C);
        (*C) = iftCreateMatrix(A->ncols, A->nrows);
    }
    unsigned int i, N;
    N = A->ncols * A->nrows;
    for (i = 0; i < N; i++) {
        (*C)->val[i] = A->val[i] / B->val[i];
    }
}

void iftComputeOperationBetweenMatricesInPointWiseInPlace(iftMatrix *A, iftMatrix *B, iftMatrix *C,
                                                          char operationSymbol) {
    unsigned int i, N;

    N = A->ncols * A->nrows;
    switch (operationSymbol) {
        case '+':
//#pragma omp parallel for private(i) shared(A,B,C,N) schedule (dynamic,1024)
            for (i = 0; i < N; i++) {
                C->val[i] = A->val[i] + B->val[i];
            }
            break;
        case '-':
            for (i = 0; i < N; i++) {
                C->val[i] = A->val[i] - B->val[i];
            }
            break;
        case '*':
            for (i = 0; i < N; i++) {
                C->val[i] = A->val[i] * B->val[i];
            }
            break;
        case '/':
            for (i = 0; i < N; i++) {
                C->val[i] = A->val[i] / B->val[i];
            }
            break;
        default:
            printf("Unrecognizable operation symbol %c", operationSymbol);
            break;
    }
}

void iftComputeOperationBetweenMatricesInPointWiseInPlaceSafe(iftMatrix *A, iftMatrix *B, iftMatrix **C,
                                                              char operationSymbol) {
    unsigned int i, N;

    if (*C == NULL) {
        if ((A->ncols != B->ncols || A->nrows != B->nrows)) {
            printf("Could not perform the matrix addition. Matrices have different dimensions");
            return;
        }
        (*C) = iftCreateMatrix(A->ncols, A->nrows);
    } else if ((A->ncols != B->ncols || A->nrows != B->nrows) || (A->ncols != (*C)->ncols || A->nrows != (*C)->nrows)) {
        printf("Could not perform the matrix addition. Matrices have different dimensions");
        return;
    }
    N = A->ncols * A->nrows;
    switch (operationSymbol) {
        case '+':
//#pragma omp parallel for private(i) shared(A,B,C,N) schedule (dynamic,1024)
            for (i = 0; i < N; i++) {
                (*C)->val[i] = A->val[i] + B->val[i];
            }
            break;
        case '-':
            for (i = 0; i < N; i++) {
                (*C)->val[i] = A->val[i] - B->val[i];
            }
            break;
        case '*':
            for (i = 0; i < N; i++) {
                (*C)->val[i] = A->val[i] * B->val[i];
            }
            break;
        case '/':
            for (i = 0; i < N; i++) {
                (*C)->val[i] = A->val[i] / B->val[i];
            }
            break;
        default:
            printf("Unrecognizable operation symbol %c", operationSymbol);
            break;
    }
}


bool iftIsMatricesEqual(const iftMatrix *A, const iftMatrix *B) {
    if (A->ncols != B->ncols)
        return false;    
    if (A->nrows != B->nrows)
        return false;    

    for (int i = 0; i < A->n; i++) {
        if (A->val[i] != B->val[i])
            return false;
    }

    return true;
}


void iftLogarithmMatrix(iftMatrix *matrix) {
    unsigned int i, N;
    N = matrix->nrows * matrix->ncols;

    for (i = 0; i < N; i++) {
        if (matrix->val[i] < 0) {
            printf("Could not compute the logarithm of the element [%d] because is a negative number\n", i);
        }
        else if (matrix->val[i] < IFT_EPSILON) {
            matrix->val[i] = IFT_INFINITY_DBL_NEG;
        } else {
            matrix->val[i] = log(matrix->val[i]);
        }
    }
}


iftMatrix *iftRandomMatrix(unsigned int ncols, unsigned int nrows, double mean, double variance) {
    iftMatrix *randomMatrix = iftCreateMatrix(ncols, nrows);
    unsigned int i;
    unsigned int N = randomMatrix->n;
    for (i = 0; i < N; i++) {
        //TODO: this is not semi random generation
        randomMatrix->val[i] = ((((double) i) / N) * variance) + mean;
    }
    return randomMatrix;
}

iftMatrix *iftCreateSetMatrix(unsigned int ncols, unsigned int nrows, double value) {
    iftMatrix *matrix = iftCreateMatrix(ncols, nrows);
    iftSetMatrix(matrix, value);
    return matrix;
}

void iftSetMatrix(iftMatrix *matrix, double value) {
    unsigned int k;
    for (k = 0; k < matrix->n; k++) {
        matrix->val[k] = value;
    }
}


void iftPointWiseMatrixPowInPlace(iftMatrix *matrix, double power) {
  
    for (int i = 0; i < matrix->n; i++)
      matrix->val[i] = pow(matrix->val[i], power);

}

iftMatrix *iftPointWiseMatrixPow(iftMatrix *matrix, double power) {
    iftMatrix *pow = iftCopyMatrix(matrix);
    iftPointWiseMatrixPowInPlace(pow, power);
    return pow;
}

iftDblArray *iftComputeSumVector(iftMatrix *matrix, char vectorType) {
    unsigned int i, j, k;
    iftDblArray *outputVector;
    k = 0;
    outputVector = NULL;


    if (vectorType == 'r') {
        outputVector = iftCreateDblArray(matrix->ncols);
        for (i = 0; i < matrix->ncols; i++) {
            outputVector->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                outputVector->val[j] += matrix->val[k];
                k++;
            }
        }
    } else if (vectorType == 'c') {
        outputVector = iftCreateDblArray(matrix->nrows);
        for (i = 0; i < matrix->nrows; i++) {
            outputVector->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                outputVector->val[i] += matrix->val[k];
                k++;
            }
        }

    } else {
        printf("Unrecognizable vector type %c", vectorType);
    }

    return outputVector;
}

iftMatrix *iftComputeMeanVectorBycolumn(iftMatrix *matrix) {
    unsigned int i, j, k;
    iftMatrix *meanVector = iftCreateMatrix(matrix->ncols, 1);
    k = 0;

    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            meanVector->val[j] += matrix->val[k];
            k++;
        }
    }

    for (j = 0; j < matrix->ncols; ++j) {
        meanVector->val[j] /= matrix->nrows;
    }

    return meanVector;
}

iftMatrix *iftComputeSumVectorBycolumn(iftMatrix *matrix) {
    unsigned int i, j, k;
    iftMatrix *sumVector = iftCreateMatrix(matrix->ncols, 1);
    k = 0;

    for (j = 0; j < matrix->ncols; j++) {
        sumVector->val[j] = 0;
    }
    k = 0;

    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            sumVector->val[j] += matrix->val[k];
            k++;
        }
    }
    return sumVector;
}

void iftComputeSumVectorBycolumnInPlace(iftMatrix *matrix, iftMatrix **vector) {
    unsigned int i, j, k;
    k = 0;
    if ((*vector) == NULL) {
        (*vector) = iftCreateMatrix(matrix->ncols, 1);
    }

    for (j = 0; j < matrix->ncols; j++) {
        (*vector)->val[j] = 0;
    }
    k = 0;

    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            (*vector)->val[j] += matrix->val[k];
            k++;
        }
    }
}

iftMatrix *iftComputeSumVectorByRow(iftMatrix *matrix) {
    unsigned int i, j, k;
    iftMatrix *sumVector = iftCreateMatrix(1, matrix->nrows);
    k = 0;

    for (j = 0; j < matrix->nrows; j++) {
        sumVector->val[j] = 0;
    }
    k = 0;

    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            sumVector->val[i] += matrix->val[k];
            k++;
        }
    }
    return sumVector;
}

void iftComputeSumVectorByRowInPlace(iftMatrix *matrix, iftMatrix **vector) {
    unsigned int i, j, k;
    if ((*vector) == NULL) {
        (*vector) = iftCreateMatrix(1, matrix->nrows);
    } else {
        (*vector)->ncols = 1;
        (*vector)->nrows = matrix->nrows;
    }

    k = 0;

    for (j = 0; j < matrix->nrows; j++) {
        (*vector)->val[j] = 0;
    }
    k = 0;

    for (i = 0; i < matrix->nrows; i++) {
        for (j = 0; j < matrix->ncols; j++) {
            (*vector)->val[i] += matrix->val[k];
            k++;
        }
    }
}

void iftTransposeVectorMatrix(iftMatrix *vectorMatrix) {
    iftSwap(vectorMatrix->nrows, vectorMatrix->ncols);
}

void iftComputeSumVectorInplace(iftMatrix *matrix, iftDblArray *array, char vectorType) {
    unsigned int i, j, k;
    k = 0;

    if (vectorType == 'r') {
        if (array->n != matrix->ncols) {
            printf("could not perform operation. The number of columns mismatch the array size");
        }
        for (i = 0; i < matrix->ncols; i++) {
            array->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                array->val[j] += matrix->val[k];
                k++;
            }
        }
    } else if (vectorType == 'c') {
        if (array->n != matrix->nrows) {
            printf("could not perform operation. The number of rows mismatch the array size");
        }
        for (i = 0; i < matrix->nrows; i++) {
            array->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                array->val[i] += matrix->val[k];
                k++;
            }
        }

    } else {
        printf("Unrecognizable vector type %c", vectorType);
    }
}

void iftComputeSumVectorInplaceSafe(iftMatrix *matrix, iftDblArray **array, char vectorType) {
    unsigned int i, j, k;
    k = 0;

    if (vectorType == 'r') {
        if (*array == NULL) {
            *array = iftCreateDblArray(matrix->ncols);
        } else if ((*array)->n != matrix->ncols) {
            printf("could not perform operation. The number of columns mismatch the array size");
        }
        for (i = 0; i < matrix->ncols; i++) {
            (*array)->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                (*array)->val[j] += matrix->val[k];
                k++;
            }
        }
    } else if (vectorType == 'c') {
        if (*array == NULL) {
            *array = iftCreateDblArray(matrix->nrows);
        }
        if ((*array)->n != matrix->nrows) {
            printf("could not perform operation. The number of rows mismatch the array size");
        }
        for (i = 0; i < matrix->nrows; i++) {
            (*array)->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                (*array)->val[i] += matrix->val[k];
                k++;
            }
        }

    } else {
        printf("Unrecognizable vector type %c", vectorType);
    }
}

iftDblArray *iftComputeMeanVector(iftMatrix *matrix, char vectorType) {
    unsigned int i, j, k;
    iftDblArray *outputVector;
    k = 0;
    outputVector = NULL;
    if (vectorType == 'r') {
        outputVector = iftCreateDblArray(matrix->ncols);
        for (i = 0; i < matrix->ncols; i++) {
            outputVector->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                outputVector->val[j] += matrix->val[k];
                k++;
            }
        }
        for (j = 0; j < matrix->ncols; j++) {
            outputVector->val[j] /= matrix->nrows;
        }

    } else if (vectorType == 'c') {
        outputVector = iftCreateDblArray(matrix->nrows);
        for (i = 0; i < matrix->nrows; i++) {
            outputVector->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                outputVector->val[i] += matrix->val[k];
                k++;
            }
        }

        for (i = 0; i < matrix->nrows; i++) {
            outputVector->val[i] /= matrix->ncols;
        }


    } else {
        printf("Unrecognizable vector type %c", vectorType);
    }

    return outputVector;
}

void iftComputeMeanVectorInPlace(iftMatrix *matrix, iftDblArray **array, char vectorType) {
    //TODO: safe check
    unsigned int i, j, k;
    k = 0;
    if (vectorType == 'r') {

        for (i = 0; i < matrix->ncols; i++) {
            (*array)->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                (*array)->val[j] += matrix->val[k];
                k++;
            }
        }
        for (j = 0; j < matrix->ncols; j++) {
            (*array)->val[j] /= matrix->nrows;
        }
    } else if (vectorType == 'c') {
        for (i = 0; i < matrix->nrows; i++) {
            (*array)->val[i] = 0;
        }

        for (i = 0; i < matrix->nrows; i++) {
            for (j = 0; j < matrix->ncols; j++) {
                (*array)->val[i] += matrix->val[k];
                k++;
            }
        }
        for (i = 0; i < matrix->nrows; i++) {
            (*array)->val[i] /= matrix->ncols;
        }
    }

}

iftMatrix *iftCreateDiagonalMatrix(double *vector, long n) {
    iftMatrix *outputMatrix = iftCreateMatrix(n, n);
    unsigned int i, shift;
    for (i = 0, shift = 0; i < outputMatrix->nrows; i++, shift += outputMatrix->ncols) {
        outputMatrix->val[shift + i] = vector[i];
    }
    return outputMatrix;
}

void iftCreateDiagonalMatrixInPlace(double *vector, long n, iftMatrix **matrix) {
    if ((*matrix) == NULL) {
        (*matrix) = iftCreateMatrix(n, n);
    } else if ((*matrix)->nrows != n || (*matrix)->ncols != n) {
        iftDestroyMatrix(matrix);
        (*matrix) = iftCreateMatrix(n, n);
    }
    unsigned int i, shift;
    for (i = 0, shift = 0; i < (*matrix)->nrows; i++, shift += (*matrix)->ncols) {
        (*matrix)->val[shift + i] = vector[i];
    }
}

iftMatrix *iftDiagonalFromMatrix(iftMatrix *matrix, char type)
{
    if(matrix->ncols != matrix->nrows)
        iftError("The diagonal can only be extracted from a square matrix", "iftDiagonalFromMatrix");
        
    iftMatrix *diag = NULL;
    if(type == 'r') {
        diag = iftCreateMatrix(matrix->ncols, 1);
        for(int i = 0; i < matrix->ncols; i++)
            iftMatrixElem(diag, i, 0) = iftMatrixElem(matrix, i, i);
    }
    else if(type == 'c') {
        diag = iftCreateMatrix(1, matrix->ncols);
        for(int i = 0; i < matrix->ncols; i++)
            iftMatrixElem(diag, 0, i) = iftMatrixElem(matrix, i, i);
    }

    return diag;
}

void iftMatrixMinSaturation(iftMatrix* m, float val) {
    #pragma omp parallel for
    for (int i = 0; i < m->n; ++i)
        m->val[i] = iftMax(val, m->val[i]);
}

void iftMatrixElementsSaturation(iftMatrix *matrix, double saturationValue, char saturationType) {
    unsigned int i, N;
    N = matrix->nrows * matrix->ncols;
    //upperBound
    if (saturationType == 'u') {
        for (i = 0; i < N; i++) {
            if (matrix->val[i] > saturationValue) {
                matrix->val[i] = saturationValue;
            }
        }
    } else if (saturationType == 'l') {//lowerbound
        for (i = 0; i < N; i++) {
            if (matrix->val[i] < saturationValue) {
                matrix->val[i] = saturationValue;
            }
        }
    }
}

iftMatrix *iftSubMatrix(iftMatrix *matrix, int startRowIndex, int endRowIndex,
                        int startColumnIndex, int endColumnIndex) {
    int numberRows = endRowIndex - startRowIndex;
    int numberColumns = endColumnIndex - startColumnIndex;
    if (numberRows < 0 || startRowIndex < 0 || endRowIndex > matrix->nrows) {
        printf("The input row indexes are not consistent");
        return NULL;
    }
    if (numberColumns < 0 || startColumnIndex < 0 || endColumnIndex > matrix->ncols) {
        printf("The end Column index is smaller than the start one");
        return NULL;
    }
    iftMatrix *outputMatrix = iftCreateMatrix(numberColumns, numberRows);
    int i, j;
    for (i = startRowIndex; i < endRowIndex; i++) {
        for (j = startColumnIndex; j < endColumnIndex; j++) {
            iftMatrixElem(outputMatrix, j-startColumnIndex, i-startRowIndex) = iftMatrixElem(matrix, j, i);
        }
    }
    return outputMatrix;
}


void iftCopyMatrixInPlace(iftMatrix *A, iftMatrix *B) {
    int N, k;
    N = A->ncols * A->nrows;
    for (k = 0; k < N; k++) {
        B->val[k] = A->val[k];
    }
}

void iftCopyMatrixInPlaceSafe(iftMatrix *A, iftMatrix **B) {
    if (*B == NULL) {
        (*B) = iftCreateMatrix(A->ncols, A->nrows);
    } else if (A->ncols != (*B)->ncols || A->nrows != (*B)->nrows) {
        iftDestroyMatrix(B);
        (*B) = iftCreateMatrix(A->ncols, A->nrows);
    }
    int N, k;
    N = A->ncols * A->nrows;
    for (k = 0; k < N; k++) {
        (*B)->val[k] = A->val[k];
    }
}

iftMatrix* iftMergeMatrixArrayVertical(iftMatrix **matArray, int n)
{
    /* compute the number of rows of the new matrix */
    int newRowSize = 0, colSize = 0, lastColSize = IFT_NIL;
    for(int i = 0; i < n; i++) {
        if(matArray[i] == NULL) {
            iftError("One of the matrices in the array is NULL.", "iftMergeMatrixArrayVertical");
        } else {
            newRowSize += matArray[i]->nrows;
            colSize = matArray[i]->ncols;
            if(lastColSize != IFT_NIL && lastColSize != matArray[i]->ncols) {
                iftError("Number of columns in all the matrices must be the same.", "iftMergeMatrixArrayVertical");
            }
            else {
                lastColSize = matArray[i]->ncols;
            }
        }
    }

    /* join the matrices into a new one */
    iftMatrix *newMat = iftCreateMatrix(colSize, newRowSize);
    int row = 0;
    for(int i = 0; i < n; i++) {
        for(int r = 0; r < matArray[i]->nrows; r++) {
            for(int c = 0; c < matArray[i]->ncols; c++) {
                iftMatrixElem(newMat, c, row) = iftMatrixElem(matArray[i], c, r);
            }
            row++;
        }
    }

    return newMat;
}

iftMatrix* iftMergeMatrixArrayHorizontal(iftMatrix **matArray, int n)
{
    /* compute the number of cols of the new matrix */
    int newColSize = 0, rowSize = 0, lastRowSize = IFT_NIL;
    for(int i = 0; i < n; i++) {
        if(matArray[i]) {
            newColSize += matArray[i]->ncols;
            rowSize = matArray[i]->nrows;
            if(lastRowSize != IFT_NIL && lastRowSize != matArray[i]->nrows) {
                iftError("Number of rows in all the matrices must be the same.", "iftMergeMatrixArrayHorizontal");
            }
            else {
                lastRowSize = matArray[i]->nrows;
            }
        }
    }

    /* join the matrices into a new one */
    iftMatrix *newMat = iftCreateMatrix(newColSize, rowSize);
    int col = 0;
    for(int i = 0; i < n; i++) {
        if(matArray[i]) {
            for(int c = 0; c < matArray[i]->ncols; c++) {
                for(int r = 0; r < matArray[i]->nrows; r++) {
                    iftMatrixElem(newMat, col, r) = iftMatrixElem(matArray[i], c, r);
                }
                col++;
            }
        }
    }

    return newMat;
}


iftMatrix *iftCentralizeMatrix(const iftMatrix *M) {
    iftMatrix *M_central = iftCreateMatrix(M->ncols, M->nrows);

    for (int col = 0; col < M->ncols; col++) {
        float mean = 0.0f;

        #pragma omp parallel for reduction(+:mean)
        for (int row = 0; row < M->nrows; row++)
            mean += iftMatrixElem(M, col, row);
        
        mean /= M->nrows;
        
        // centralize the coefficients
        #pragma omp parallel for
        for (int row = 0; row < M->nrows; row++)
            iftMatrixElem(M_central, col, row) = iftMatrixElem(M, col, row) - mean;
    }

    return M_central;
}


iftMatrix *iftUnitNormMatrix(const iftMatrix *M) {
    iftMatrix *M_norm = iftCreateMatrix(M->ncols, M->nrows);

    // for each column-vector
    for (int col = 0; col < M->ncols; col++) {
        float norm = 0;

        #pragma omp parallel for reduction(+:norm)
        for (int row = 0; row < M->nrows; row++)
            norm += iftMatrixElem(M, col, row) * iftMatrixElem(M, col, row);
        norm = sqrt(norm);
        
        if (norm > IFT_EPSILON) {
            #pragma omp parallel for
            for (int row = 0; row < M->nrows; row++)
                iftMatrixElem(M_norm, col, row) = iftMatrixElem(M, col, row) / norm;
        }
    }

    return M_norm;
}

iftDblArray **iftMatrixToDblArraySetByRow(iftMatrix *M)
{
    iftDblArray **arrays = (iftDblArray**) iftAlloc(M->nrows, sizeof(iftDblArray*));
    for(int r = 0; r < M->nrows; r++) {
        arrays[r] = iftCreateDblArray(M->ncols);
        for(int c = 0; c < M->ncols; c++)
            (arrays[r])->val[c] = iftMatrixElem(M, c, r);
    }

    return arrays;
}

iftFloatArray **iftMatrixToFloatArraySetByRow(iftMatrix *M)
{
    iftFloatArray **arrays = (iftFloatArray**) iftAlloc(M->nrows, sizeof(iftFloatArray*));
    for(int r = 0; r < M->nrows; r++) {
        arrays[r] = iftCreateFloatArray(M->ncols);
        for(int c = 0; c < M->ncols; c++)
            (arrays[r])->val[c] = iftMatrixElem(M, c, r);
    }

    return arrays;
}

iftDblArray **iftMatrixToDblArraySetByColumn(iftMatrix *M)
{
    iftDblArray **arrays = (iftDblArray**) iftAlloc(M->ncols, sizeof(iftDblArray*));
    for(int c = 0; c < M->ncols; c++) {
        arrays[c] = iftCreateDblArray(M->nrows);
        for(int r = 0; r < M->nrows; r++)
            (arrays[c])->val[r] = iftMatrixElem(M, r, c);
    }

    return arrays;
}

iftFloatArray **iftMatrixToFloatArraySetByColumn(iftMatrix *M)
{
    iftFloatArray **arrays = (iftFloatArray**) iftAlloc(M->ncols, sizeof(iftFloatArray*));
    for(int c = 0; c < M->ncols; c++) {
        arrays[c] = iftCreateFloatArray(M->nrows);
        for(int r = 0; r < M->nrows; r++)
            (arrays[c])->val[r] = iftMatrixElem(M, r, c);
    }

    return arrays;
}


double *iftIdentityDouble(int d)
{
    double *eye = iftAlloc(d * d, sizeof *eye);
    for (int i = 0; i < d; i++) {
        eye[i * d + i] = 1;
    }
    return eye;
}

void iftExpandMatrixInPlace(iftMatrix **M, char axis, int n, float value)
{
    iftMatrix *out = NULL;
    if (axis == 0) { // increase number of rows
        out = iftCreateMatrix((*M)->ncols,(*M)->nrows+n);
    } else if (axis == 1) { // increase number of cols
        out = iftCreateMatrix((*M)->ncols+n,(*M)->nrows);
    } else {
        iftError("Axis must be 0 for rows or 1 for columns","iftExpandMatrixInPlace");
    }
    //iftSetMatrix(out,value);
    for (int r = 0; r < out->nrows; r++) {
        for (int c = 0; c < out->ncols; c++) {
            if ((r < (*M)->ncols) && (c < (*M)->ncols)) {
                iftMatrixElem(out, c, r) = iftMatrixElem(*M, c, r);
            } else {
                iftMatrixElem(out, c, r) = value;
            }
        }
    }
    iftDestroyMatrix(M);
    (*M) = out;
}

void iftConcatMatricesInPlace(iftMatrix **M1, iftMatrix *M2, char axis)
{
    if (M2 == NULL)
        return;
    if ((*M1) == NULL){
        (*M1) = iftCopyMatrix(M2);
        return;
    }

    if (axis == 0) { // concatenate cols
        if ((*M1)->nrows != M2->nrows)
            iftError("To concatenate columns, the number of rows must be equal", "iftConcatMatrisesInPlace");
    }else if (axis == 1) { // concatenate rows
        if ((*M1)->ncols != M2->ncols)
            iftError("To concatenate rows, the number of columns must be equal", "iftConcatMatrisesInPlace");
    }else {
        iftError("Axis must be either 0 or 1","iftConcatMatrisesInPlace");
    }

    iftMatrix *M = NULL;
    if (axis == 0)
    {
        M = iftCreateMatrix((*M1)->ncols+M2->ncols,(*M1)->nrows);
        for (int r = 0; r < (*M1)->nrows; r++)
            for (int c = 0; c < (*M1)->ncols; c++)
                iftMatrixElem(M,c,r) = iftMatrixElem((*M1),c,r);
        for (int r = 0; r < (*M1)->nrows; r++)
            for (int c = (*M1)->ncols, k = 0; c < (*M1)->ncols+M2->ncols; c++, k++)
                iftMatrixElem(M,c,r) = iftMatrixElem(M2,k,r);
    } else if (axis == 1){
        M = iftCreateMatrix((*M1)->ncols,(*M1)->nrows+M2->nrows);
        for (int r = 0; r < (*M1)->nrows; r++)
            for (int c = 0; c < (*M1)->ncols; c++)
                iftMatrixElem(M,c,r) = iftMatrixElem((*M1),c,r);
        for (int r = (*M1)->nrows, k = 0; r < (*M1)->nrows+M2->nrows; r++, k++)
            for (int c = 0; c < (*M1)->ncols; c++)
                iftMatrixElem(M,c,r) = iftMatrixElem(M2,k,r);
    }

    iftDestroyMatrix(M1);
    (*M1) = M;
}
