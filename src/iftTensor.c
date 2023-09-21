#include "iftTensor.h"

// #include "iftTensor.cuh"
#include "ift/core/io/CSV.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"

void iftReshapeTensorWithDimArray(iftTensor *T, int ndimension, int *dimension)
{
    int* accDimension = iftAllocIntArray(ndimension);

    long nelems = 1;
    for (int i = 0; i < ndimension; ++i) {
        nelems *= dimension[i];
    }

    accDimension[ndimension-1] = 1;
    for (int i = ndimension-2; i >= 0; --i) {
        accDimension[i] = accDimension[i+1] * dimension[i+1];
    }

    if(nelems != T->n) {
        iftPrintFormattedIntArray(T->dimension, T->ndimension, "From shape = (", ")", ", ");
        iftPrintFormattedIntArray(dimension, ndimension, " to shape = (", ")\n", ", ");
        iftError("Invalid dimension for tensor with %d elements.\n", "iftReshapeTensor", T->n);
        iftFree(dimension);
        iftFree(accDimension);
    }

    T->ndimension = ndimension;

    iftFree(T->dimension);
    T->dimension = iftAllocIntArray(ndimension);
    iftCopyIntArray(T->dimension, dimension, ndimension);

    iftFree(T->accDimension);
    T->accDimension = iftAllocIntArray(ndimension);
    iftCopyIntArray(T->accDimension, accDimension, ndimension);
    iftFree(accDimension);
}

void iftReshapeTensor(iftTensor *T, int ndimension, ...)
{
    /* create the dimension array with the variable arguments */
    va_list arguments;
    int* dimension = iftAllocIntArray(ndimension);
    va_start (arguments, ndimension);
    for (int i = 0; i < ndimension; ++i)
        dimension[i] = va_arg (arguments, int);
    va_end(arguments);

    iftReshapeTensorWithDimArray(T, ndimension, dimension);
    iftFree(dimension);
}

iftTensor* iftCreateTensorWithDimArray(int ndim, int* dim) 
{
  long nelems = 1;
  for (int i = 0; i < ndim; ++i)
    nelems *= dim[i];

  float *val = iftAllocFloatArray(nelems);
  return iftCreateTensorPointerWithDimArray(val, ndim, dim);
}

iftTensor* iftCreateTensorPointerWithDimArray(float *val, int ndim, int *dim)
{
  int* accDimension = iftAllocIntArray(ndim);

  long nelems = 1;
  for (int i = 0; i < ndim; ++i)
    nelems *= dim[i];

  accDimension[ndim-1] = 1;
  for (int i = ndim-2; i >= 0; --i)
    accDimension[i] = accDimension[i+1] * dim[i+1];

  iftTensor* tensor = (iftTensor*) iftAlloc(1, sizeof(iftTensor));

  tensor->ndimension = ndim;
  tensor->dimension = iftAllocIntArray(ndim);
  iftCopyIntArray(tensor->dimension, dim, ndim);
  tensor->accDimension = accDimension;
  tensor->allocated = true;
  tensor->val = val;
  tensor->n = nelems;

  return tensor;
}

iftTensor* iftCreateTensor(int ndimension, ...)
{
    /* create the dimension array with the variable arguments */
    va_list arguments;
    int* dimension = iftAllocIntArray(ndimension);
    va_start (arguments, ndimension);
    for (int i = 0; i < ndimension; ++i)
        dimension[i] = va_arg (arguments, int);
    va_end(arguments);

    iftTensor *tensor = iftCreateTensorWithDimArray(ndimension, dimension);
    iftFree(dimension);

    return tensor;
}

iftTensor* iftCreateTensorPointer(float *val, int ndimension, ...)
{
    /* create the dimension array with the variable arguments */
    va_list arguments;
    int* dimension = iftAllocIntArray(ndimension);
    va_start (arguments, ndimension);
    for (int i = 0; i < ndimension; ++i)
        dimension[i] = va_arg (arguments, int);
    va_end(arguments);

    iftTensor *tensor = iftCreateTensorPointerWithDimArray(val, ndimension, dimension);
    iftFree(dimension);

    return tensor;
}

void iftDestroyTensor(iftTensor** t)
{
    iftTensor* tensor = *t;
    if(tensor!=NULL) {
	    if(tensor->allocated)
            iftFree(tensor->val);
        iftFree(tensor->dimension);
        iftFree(tensor->accDimension);
        iftFree(*t);
        *t = NULL;
    }
}

void iftDestroyTensorPointer(iftTensor** t)
{
  assert(t != NULL);
  if ((*t) != NULL) {
    (*t)->val = NULL;
    iftDestroyTensor(t);
  }
}

iftTensor* iftCreateTensorCopy(iftTensor* src)
{
	iftTensor *dst = (iftTensor*) iftAlloc(1, sizeof(iftTensor));
	
	dst->ndimension = src->ndimension;
    dst->allocated = true;
    dst->n = src->n;

	dst->dimension = iftAllocIntArray(src->ndimension);
    dst->accDimension = iftAllocIntArray(src->ndimension);

	iftCopyIntArray(dst->dimension, src->dimension, src->ndimension);	
	iftCopyIntArray(dst->accDimension, src->accDimension, src->ndimension);

	dst->val = iftAllocFloatArray(src->n);
	iftCopyFloatArray(dst->val, src->val, src->n);

	return dst;
}

int iftGetTensorIndexDebug(iftTensor* tensor, int ndimension, ...)
{
    int dim;
    va_list arguments;
    va_start(arguments, ndimension);

    int* dimension= iftAllocIntArray(ndimension);

    for (int i = 0; i < ndimension; ++i) {
        dimension[i] = va_arg (arguments, int);
    }

    if(ndimension!=tensor->ndimension){
        printf("Position = "); iftPrintIntArray(dimension, ndimension);
        printf("Shape = "); iftPrintIntArray(tensor->dimension, tensor->ndimension);
        iftError("Invalid Index. Dimension mismatch", "iftTensorElem");
    }

    int index = 0;
    for (int i = 0; i < ndimension; ++i) {
        dim = dimension[i];

        if(dim<0 || dim>=tensor->dimension[i]) {
            printf("Position = "); iftPrintIntArray(dimension, ndimension);
            printf("Shape = "); iftPrintIntArray(tensor->dimension, tensor->ndimension);
            iftError("Invalid Index. %dth dimension = %d", "iftTensorElem", i, dim);
        }
        index += dim * tensor->accDimension[i];
    }

	iftFree(dimension);
    va_end(arguments);
    return index;
}

int iftGetTensorIndexWithDimArray(iftTensor* tensor, int *dimension)
{
    int index = 0;
    for (int i = 0; i < tensor->ndimension; ++i)
        index += dimension[i] * tensor->accDimension[i];

    return index;
}

int iftGetTensorIndex(iftTensor* tensor, int ndimension, ...)
{
    va_list arguments;
    int index = 0;
    va_start (arguments, ndimension);
    for (int i = 0; i < tensor->ndimension; ++i)
        index += (va_arg(arguments, int) * tensor->accDimension[i]);
    va_end(arguments);

    return index;
}

iftTensor* iftReadTensorCsv(char* tensorFilename, char* dimension, char separator)
{
	iftCSV* dim = iftReadCSV(dimension, separator);
	iftCSV* tensorCsv = iftReadCSV(tensorFilename, separator);
	iftTensor* tensor = (iftTensor *) iftAlloc(1, sizeof(iftTensor));
	float f;
	int d;

	tensor->ndimension = dim->nrows;
	tensor->n = tensorCsv->nrows;
	tensor->allocated = true;

	tensor->dimension = iftAllocIntArray(dim->nrows);
	tensor->accDimension = iftAllocIntArray(dim->nrows);
	tensor->val = iftAllocFloatArray(tensorCsv->nrows);

	for(int i = 0; i < dim->nrows; i++){
		sscanf(dim->data[i][0], "%d", &d);
		tensor->dimension[i] = d;
	}

	tensor->accDimension[dim->nrows-1] = 1;
    for (int i = tensor->ndimension-2; i >= 0; --i) {
        tensor->accDimension[i] = tensor->accDimension[i+1] * tensor->dimension[i+1];
	}

	for(int i = 0; i < tensorCsv->nrows; i++){
		sscanf(tensorCsv->data[i][0], "%f", &f);
		tensor->val[i] = f;
	}

	iftDestroyCSV(&dim);
	iftDestroyCSV(&tensorCsv);

	return tensor;

}

void iftWriteTensorCsv(iftTensor *tensor, char *tensorFilename, char *dimArrayFilename, char separator)
{
    iftCSV* dimCsv = iftCreateCSV(tensor->ndimension, 1);
    for(int i = 0; i < tensor->ndimension; i++)
		sprintf(dimCsv->data[i][0], "%d", tensor->dimension[i]);
    iftWriteCSV(dimCsv, dimArrayFilename, ';');

	iftCSV* tensorCsv = iftCreateCSV(tensor->n, 1);
	for(int i = 0; i < tensor->n; i++)
		sprintf(tensorCsv->data[i][0], "%f", tensor->val[i]);
    iftWriteCSV(tensorCsv, tensorFilename, separator);

	iftDestroyCSV(&dimCsv);
	iftDestroyCSV(&tensorCsv);
}

/* ------------------------------------------------------------------------------------------------------------------
WARNING: This function needs to be optimized since it performs the batched matrix multiplications one by one.
Several matrix multiplications could be joined in one big matrix to speed up the whole process
--------------------------------------------------------------------------------------------------------------------- */
iftTensor *iftMultTensors(iftTensor *A, iftTensor *B)
{
    /* verify if there are empty tensors */
    if(A->ndimension == 0 || B->ndimension == 0) {
        iftError("One of the tensors is empty (A.shape is %s and B.shape is %s)",
                 "iftMultTensors", iftTensorDimStr(A), iftTensorDimStr(B));
    }
    /* verify that the number of dimensions is correct */
    else if(A->ndimension < B->ndimension) {
        iftError("The number of dimensions of A must be equal or higher than the dimensions of B (A.shape is %s and B.shape is %s)",
                 "iftMultTensors", iftTensorDimStr(A), iftTensorDimStr(B));
    }
    /* verify 1D tensor multiplication */
    else if(A->ndimension == 1 && B->ndimension == 1) {
        if(A->dimension[0] != B->dimension[0])
            iftError("Error in 1D tensor multiplication: dimension 0 in A and B must be equal (A.shape is %s and B.shape is %s)",
                 "iftMultTensors", iftTensorDimStr(A), iftTensorDimStr(B));
    }
    /* verify 2D tensor multiplication */
    else if(A->ndimension == 2 && B->ndimension == 2) {
        if(A->dimension[1] != B->dimension[0])
            iftError("Error in 2D tensor multiplication: dimension 1 in A must be equal to dimension 0 in B (A.shape is %s and B.shape is %s)",
                 "iftMultTensors", iftTensorDimStr(A), iftTensorDimStr(B));
    }
    /* verify nD tensor multiplication */
    else {
        if(A->dimension[A->ndimension-1] != B->dimension[B->ndimension-2])
            iftError("Error in %dD tensor multiplication: dimension %d in A must be equal to dimension %d in B (A.shape is %s and B.shape is %s)",
                 "iftMultTensors", A->ndimension, A->ndimension-1, B->ndimension-2, iftTensorDimStr(A), iftTensorDimStr(B));

        for(int d = B->ndimension-3; d >= 0; d--) {
            int dA = d + (A->ndimension - B->ndimension);
            if(A->dimension[dA] != B->dimension[d] && B->dimension[d] != 1) {
                iftError("Error in %dD tensor multiplication: dimension %d in B must be either equal to dimension %d in A or equal to 1 (A.shape is %s and B.shape is %s)",
                            "iftMultTensors", A->ndimension, d, dA, iftTensorDimStr(A), iftTensorDimStr(B));
            }
        }
    }
    
    /* create output tensor */
    int *outDim = iftAllocIntArray(A->ndimension);
    if(A->ndimension == 1 && B->ndimension == 1) {
        outDim[0] = 1;
    }
    else {
        iftCopyIntArray(outDim, A->dimension, A->ndimension);
        outDim[A->ndimension-1] = B->dimension[B->ndimension-1];
    }
    iftTensor *C = iftCreateTensorWithDimArray(A->ndimension, outDim);

    /* perform matrix multiplication according to the number of dimensions in the tensors */
    iftMatrix *Amat = NULL, *Bmat = NULL, *Cmat = NULL;
    
    /* 1D tensor multiplication */
    if(A->ndimension == 1 && B->ndimension == 1) {
        Amat = iftCreateMatrix(A->dimension[0], 1);
        Bmat = iftCreateMatrix(1, B->dimension[0]);

        for(int i = 0; i < A->dimension[0]; i++) {
            iftMatrixElem(Amat, i, 0) = iftTensorElem(A, i);
            iftMatrixElem(Bmat, 0, i) = iftTensorElem(B, i);
        }

        Cmat = iftMultMatrices(Amat, Bmat);
        iftTensorElem(C, 0) = iftMatrixElem(Cmat, 0, 0);
        
        iftDestroyMatrix(&Amat);
        iftDestroyMatrix(&Bmat);
        iftDestroyMatrix(&Cmat);
    }
    /* 2D tensor multiplication */
    else if(A->ndimension == 2 && B->ndimension == 2) {
        Amat = iftCreateMatrix(A->dimension[1], A->dimension[0]);
        Bmat = iftCreateMatrix(B->dimension[1], B->dimension[0]);

        for(int i = 0; i < A->dimension[0]; i++)
            for(int j = 0; j < A->dimension[1]; j++)
                iftMatrixElem(Amat, j, i) = iftTensorElem(A, i, j);

        for(int i = 0; i < B->dimension[0]; i++)
            for(int j = 0; j < B->dimension[1]; j++)
                iftMatrixElem(Bmat, j, i) = iftTensorElem(B, i, j);

        Cmat = iftMultMatrices(Amat, Bmat);

        for(int i = 0; i < A->dimension[0]; i++)
            for(int j = 0; j < B->dimension[1]; j++)
                iftTensorElem(C, i, j) = iftMatrixElem(Cmat, j, i);

        iftDestroyMatrix(&Amat);
        iftDestroyMatrix(&Bmat);
        iftDestroyMatrix(&Cmat);
    }
    /* nD tensor multiplication (batched matrix multiplication) */
    else {
        /* if tensor B has less dimensions than A, 1 is used to fill the remaining dimensions in B */
        if(A->ndimension != B->ndimension) {
            int *newDim = iftAllocIntArray(A->ndimension);
            iftSetIntArray(newDim, A->ndimension, 1);
            iftCopyIntArray(newDim+(A->ndimension-B->ndimension), B->dimension, B->ndimension);
            iftReshapeTensorWithDimArray(B, A->ndimension, newDim);
        }

        /* batched matrix multiplication */
        int nMatMult = 1;
        for(int d = 0; d < A->ndimension-2; d++)
            nMatMult *= A->dimension[d];

        int *dimIdxCount = iftAllocIntArray(A->ndimension-2);

        for(int m = 0; m < nMatMult; m++) {
            /* slice tensor A */
            int **sliceIdxA = iftAllocIntMatrix(2, A->ndimension);
            for(int d = 0; d < A->ndimension-2; d++) {
                sliceIdxA[d][0] = dimIdxCount[d];
                sliceIdxA[d][1] = dimIdxCount[d]+1;
            }
            sliceIdxA[A->ndimension-2][0] = 0;
            sliceIdxA[A->ndimension-2][1] = A->dimension[A->ndimension-2];
            sliceIdxA[A->ndimension-1][0] = 0;
            sliceIdxA[A->ndimension-1][1] = A->dimension[A->ndimension-1];

            iftTensor *Aprime = iftSliceTensor(A, sliceIdxA);

            /* slice tensor B */
            int **sliceIdxB = iftAllocIntMatrix(2, B->ndimension);
            for(int d = 0; d < B->ndimension-2; d++) {
                sliceIdxB[d][0] = dimIdxCount[d];
                sliceIdxB[d][1] = dimIdxCount[d]+1;
            }
            sliceIdxB[B->ndimension-2][0] = 0;
            sliceIdxB[B->ndimension-2][1] = B->dimension[B->ndimension-2];
            sliceIdxB[B->ndimension-1][0] = 0;
            sliceIdxB[B->ndimension-1][1] = B->dimension[B->ndimension-1];

            iftTensor *Bprime = iftSliceTensor(B, sliceIdxB);

            /* multiply the sliced 2D tensors */
            iftTensor *Cprime = iftMultTensors(Aprime, Bprime);

            int *dimIdxC = iftAllocIntArray(A->ndimension);
            for(int d = 0; d < A->ndimension-2; d++)
                dimIdxC[d] = dimIdxCount[d];

            /* copy the result to the corresponding sliced portion of tensor C */
            for(int i = 0; i < C->dimension[C->ndimension-2]; i++) {
                for(int j = 0; j < C->dimension[C->ndimension-1]; j++) {
                    dimIdxC[C->ndimension-2] = i;
                    dimIdxC[C->ndimension-1] = j;
                    iftTensorElemWithDimArray(C, dimIdxC) = iftTensorElem(Cprime, i, j);
                }
            }

            /* update the dimension index counter */
            dimIdxCount[A->ndimension-3]++;
            for(int d = A->ndimension-4; d >=0 ; d--) {
                if(dimIdxCount[d+1] == A->dimension[d+1]) {
                    dimIdxCount[d]++;
                    dimIdxCount[d+1] = 0;
                }
                else {
                    break;
                }
            }

            iftDestroyTensor(&Aprime);
            iftDestroyTensor(&Bprime);
            iftDestroyTensor(&Cprime);
        }
    }

    return C;
}

iftTensor *iftSliceTensor(iftTensor *A, int **sliceIdx)
{
    /* create the output tensor */
    int *dimArray = iftAllocIntArray(A->ndimension);
    for(int d = 0; d < A->ndimension; d++) {
        if(sliceIdx[d][1] <= sliceIdx[d][0]) {
            iftError("The slice index for dimension %d contains invalid values a, b = (%d, %d). \"a\" must be greater than \"b\".",
                     "iftSliceTensor", d, sliceIdx[d][0], sliceIdx[d][1]);
        }
        dimArray[d] = sliceIdx[d][1] - sliceIdx[d][0];
    }

    iftTensor *Aprime = iftCreateTensorWithDimArray(A->ndimension, dimArray);

    /* compute the number of copy operations that must be performed */
    int nOper = 1;
    for(int d = 0; d < A->ndimension; d++)
        nOper *= Aprime->dimension[d];

    /* initialize the dimension index counters */
    int *dimIdxCountA = iftAllocIntArray(A->ndimension);
    int *dimIdxCountAprime = iftAllocIntArray(A->ndimension);
    for(int d = 0; d < A->ndimension; d++)
        dimIdxCountA[d] = sliceIdx[d][0];

    /* slice the tensor */
    for(int i = 0; i < nOper; i++) {
        iftTensorElemWithDimArray(Aprime, dimIdxCountAprime) = iftTensorElemWithDimArray(A, dimIdxCountA);

        /* update the dimension index counter for tensor A */
        dimIdxCountA[A->ndimension-1]++;
        for(int d = A->ndimension-2; d >=0 ; d--) {
            if(dimIdxCountA[d+1] == sliceIdx[d+1][1]) {
                dimIdxCountA[d]++;
                dimIdxCountA[d+1] = 0;
            }
            else {
                break;
            }
        }

        /* update the dimension index counter for tensor Aprime */
        dimIdxCountAprime[A->ndimension-1]++;
        for(int d = A->ndimension-2; d >=0 ; d--) {
            if(dimIdxCountAprime[d+1] == Aprime->dimension[d+1]) {
                dimIdxCountAprime[d]++;
                dimIdxCountAprime[d+1] = 0;
            }
            else {
                break;
            }
        }
    }

    /* if possible, reduce the dimensions that are equal to 1 */
    int nDimRed = 0;
    for(int d = 0; d < Aprime->ndimension; d++) {
        if(Aprime->dimension[d] == 1)
            nDimRed++;
        else
            break;
    }
    int *newDimArray = iftAllocIntArray(Aprime->ndimension-nDimRed);
    for(int d = Aprime->ndimension-nDimRed; d >=0; d--)
        newDimArray[d] = Aprime->dimension[d+nDimRed];

    iftReshapeTensorWithDimArray(Aprime, Aprime->ndimension-nDimRed, newDimArray);

    return Aprime;
}

void iftPrintTensorDim(iftTensor* T)
{
	if(T==NULL) {
		printf("()");
		return;
	}

	printf("(");
	for (int d = 0; d < T->ndimension; ++d) {
		printf("%d", T->dimension[d]);
		if(d<(T->ndimension-1))
			printf(", ");
	}
	printf(")\n");
}

char *iftTensorDimStr(iftTensor* T)
{
    char *str = NULL, *aux = NULL, dim[256];
	
    if(T==NULL) {
		str = iftCopyString("()");
		return str;
	}

	str = iftCopyString("(");
	for (int d = 0; d < T->ndimension; ++d) {
        aux = iftCopyString(str);
        iftFree(str);
        sprintf(dim, "%d", T->dimension[d]);
        str = iftConcatStrings(2, aux, dim);
        iftFree(aux);
		if(d<(T->ndimension-1)) {
            aux = iftCopyString(str);
            iftFree(str);
            str = iftConcatStrings(2, aux, ", ");
            iftFree(aux);
        }
	}
    aux = iftCopyString(str);
    iftFree(str);
    str = iftConcatStrings(2, aux, ")");
    iftFree(aux);

    return str;
}
