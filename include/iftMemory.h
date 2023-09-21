/**
 * @brief Memory Management functions.
 * @author Peixinho
 * @date Oct 8, 2015
 * @ingroup Memory
 */

#ifndef IFT_IFTMEMORY_H
#define IFT_IFTMEMORY_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if !defined(__APPLE__)
    #include <malloc.h>
#endif


#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Invoke garbage collector.
 * @author Peixinho
 * @ingroup Memory
 * @date Nov, 2016
 */
void iftGC();

/**
 * @brief The amount of memory currently being used by this process, in bytes.
 * @author Peixinho
 * @ingroup Memory
 */
long iftMemoryUsed();

long iftGetPhysicalSystemMemory();

/**
 * @brief Returns the amount of free system memory, in bytes.
 * @author Cesar Castelo
 * @ingroup Memory
 * @date Feb 07, 2019
 */
long iftGetFreePhysicalSystemMemory();

/**
 * @brief Returns the amount of available system memory, in bytes.
 * @author Matheus Cerqueira 
 * @ingroup Memory
 * @date Aug 08, 2022
 */
float iftGetAvailableMemory(void);

/**
 * @brief Shows the current number of allocated objects, for debugging purpose. See also iftAllocObjectsCount().
 * @author Peixinho
 * @ingroup Memory
 */
void iftPrintAllocObjectsCount();

/**
 * @brief Prints a memory size, with an appropriate metric (KB, MB, GB).
 * @param mem Memory size in bytes
 * @author Peixinho
 * @date Jun, 2017
 */
void iftPrintMemoryFormatted(long mem);

/**
 * @brief Count the current number of allocated objects, for debugging purpose.
 * @author Peixinho
 * @ingroup Memory
 */
long iftAllocObjectsCount();

/**
 * @brief Prints the memory amount used by a array.
 * @author Cesar Castelo
 * @param arraySize Size of the array
 * @param elemSize Size of each element
 * @param unit A char indicating the unit of measure (b: bytes, k: kbytes, m: mbytes, g: gbytes)
 * @ingroup Memory
 */
double iftMemoryUsedByArray(long arraySize, int elemSize, char unit);







/**
 * @brief  Verify memory deallocation. 
 * @author Alexandre Falcao
 * @date   Feb 23, 2016
 * @ingroup Memory
 */
bool iftVerifyMemory(long MemDinInicial, long MemDinFinal);

/**
 * @brief Swap the content of two variables.
 * @warning In pointer, this function only swaps the pointers themselves, not their memory content.
 * @author Peixinho
 * @date Mar, 2016
 */
#define iftSwap(x, y) do { __typeof__(x) _IFT_SWAP_ = x; x = y; y = _IFT_SWAP_; } while (0)

/**
 * @brief Swap the content of two strings.
 * @author Peixinho
 * @date Mar, 2016
 */
void iftSwapString(char *a, char *b);


/**
 * @brief Allocates a Char Array (string) with <b>n</b> positions
 * @ingroup Memory
 * @{
 */
bool *iftAllocBoolArray(long n);
char *iftAllocCharArray(long n);
char *iftAllocString(long n);
uchar *iftAllocUCharArray(long n);
short *iftAllocShortArray(long n);
ushort *iftAllocUShortArray(long n);
int *iftAllocIntArray(long n);
uint *iftAllocUIntArray(long n);
long *iftAllocLongIntArray(long n);
#ifndef  __cplusplus
long long *iftAllocLongLongIntArray(long n);
#endif
ulong *iftAllocULongArray(long n);
#ifndef  __cplusplus
ullong *iftAllocULLongArray(long n);
#endif
float *iftAllocFloatArray(long n);
double *iftAllocDoubleArray(long n);
long double *iftAllocLongDoubleArray(long n);
iftComplex *iftAllocComplexArray(long n);

// Aligned memory allocation
uchar *iftAllocAlignedUCharArray(long n, long alignment);
int *iftAllocAlignedIntArray(long n, long alignment);
float *iftAllocAlignedFloatArray(long n, long alignment);
double *iftAllocAlignedDoubleArray(long n, long alignment);
/** @} */


/**
 * @brief Functions to allocate/deallocate a matrix with <b>n</b> positions of different data types.
 * The matrix is allocated as a single array of size (c x r) with the corresponding pointers to the beginning of each row
 * @ingroup Memory
 * @author Cesar Castelo
 * @date Jul 16, 2018
 * @warning These pointers MUST be set free with iftFreeMatrix. A call to iftFree DOES NOT free the memory properly
 * @{
 */
bool **iftAllocBoolMatrix(long c, long r);
char **iftAllocCharMatrix(long c, long r);
uchar **iftAllocUCharMatrix(long c, long r);
short **iftAllocShortMatrix(long c, long r);
ushort **iftAllocUShortMatrix(long c, long r);
int **iftAllocIntMatrix(long c, long r);
uint **iftAllocUIntMatrix(long c, long r);
long **iftAllocLongMatrix(long c, long r);
#ifndef  __cplusplus
long long **iftAllocLongLongMatrix(long c, long r);
#endif
ulong **iftAllocULongMatrix(long c, long r);
#ifndef  __cplusplus
ullong **iftAllocULLongMatrix(long c, long r);
#endif
float **iftAllocFloatMatrix(long c, long r);
double **iftAllocDoubleMatrix(long c, long r);
long double **iftAllocLongDoubleMatrix(long c, long r);
iftComplex **iftAllocComplexMatrix(long c, long r);

#define iftFreeMatrix(M, NROWS) {\
    if(M!=NULL) {\
        for(int i = 0; i < NROWS; i++) {\
            if(M[i]!=NULL) {\
                iftFree(M[i]);\
            }\
        }\
        iftFree(M);\
    }\
}
/** @} */









/**
 * @brief Sets all values in <array_dst> to <val>.
 *
 * @param array_dst The array.
 * @param nelems The number of elements.
 * @param val value to set array
 * @ingroup Memory
 */
void iftSetFloatArray(float *array_dst, int nelems, float val);

/**
 * @brief Copies the float <array_src> to <array_dst>. We expect the number of elements to match.
 *
 * @param array_dst The destination array.
 * @param array_src The source array.
 * @param nelems The number of elements.
 * @ingroup Memory
 */
void iftCopyFloatArray(float *array_dst, float *array_src, int nelems);

/**
 * @brief Copies the float <array_src> to the double <array_dst>. We expect the number of elements to match.
 *
 * @param array_dst The destination double array.
 * @param array_src The source float array.
 * @param nelems The number of elements.
 * @ingroup Memory
 * 
 * @author Samuel Martins
 * @date May 15, 2018
 */
void iftCopyFloatArrayToDblArray(double *array_dst, float *array_src, int nelems);

/**
 * @brief Copies the double <array_src> to the float <array_dst>. We expect the number of elements to match.
 *
 * @param array_dst The destination float array.
 * @param array_src The source double array.
 * @param nelems The number of elements.
 * @ingroup Memory
 * 
 * @author Samuel Martins
 * @date May 15, 2018
 */
void iftCopyDblArrayToFloatArray(float *array_dst, double *array_src, int nelems);

/**
 * @brief Copies the double <array_src> to <array_dst>. We expect the number of elements to match.
 *
 * @param array_dst The destination array.
 * @param array_src The source array.
 * @param nelems The number of elements.
 * @ingroup Memory
 */
void iftCopyDoubleArray(double *array_dst, double *array_src, int nelems);
/**
 * @brief Copies the int <array_src> to <array_dst>. We expect the number of elements to match.
 *
 * @param array_dst The destination array.
 * @param array_src The source array.
 * @param nelems The number of elements.
 * @ingroup Memory
 */
void iftCopyIntArray(int *array_dst, const int *array_src, int nelems);

#ifndef  __cplusplus
/**
 * @brief Copies the long long int <array_src> to <array_dst>. We expect the number of elements to match.
 *
 * @param array_dst The destination array.
 * @param array_src The source array.
 * @param nelems The number of elements.
 * @ingroup Memory
 */
void iftCopyLongLongIntArray(long long *array_dst, const long long *array_src, int nelems);
#endif

/**
 * @brief Concatenates two int arrays.
 *
 * @param array1 The first array.
 * @param n1 The number of elements in the first array.
 * @param array2 The first array.
 * @param n2 The number of elements in the second array.
 * @param nelems Returns the number of elements in the new array.
 * @return The concatenated arrays.
 */
int *iftConcatIntArray(int *array1, int n1, int *array2, int n2, int *nelems);

/**
 * @brief Sets all values in <array_dst> to <val>.
 *
 * @param array_dst The array.
 * @param nelems The number of elements.
 * @param val value to set array
 * @ingroup Memory
 */
void iftSetIntArray(int *array_dst, int nelems, int val);



void iftCopySizeTArray(long *array_dst, const long *array_src, long nelems);





/**
 * @author Cesar Castelo
 * @date Jan 23, 2018
 * @brief Initializes a float array with a given value
 * @param array Float array
 * @param value Value to initialize the array
 * @param nelems Size of the array
 * @return
 */
void iftInitFloatArray(float *array, float value, int nelems);

/**
 * @author Cesar Castelo
 * @date Jan 23, 2018
 * @brief Initializes a single position in a float array with a given value
 * @param array Float array
 * @param pos Position to be initialized
 * @param value Value to initialize the array
 * @param nelems Size of the array
 * @return
 */
void iftInitFloatArrayPos(float *array, int pos, float value, int nelems);

/**
 * @author Cesar Castelo
 * @date Jan 23, 2018
 * @brief Adds two float arrays and save the result in a third array
 * @param array1 Float array to be added
 * @param array2 Float array to be added
 * @param array3 Float array to store the result
 * @param nelems Size of the array
 * @return
 */
void iftSumFloatArrays(float *array1, float *array2, float *array3, int nelems);


/**
 * @brief Sums all elements from the float array <farr>.
 * @author Samuel Martins
 * @date Aug 21, 2018
 */
float iftSumFloatArrayElems(float *farr, int n_elems);


/**
 * @author Cesar Castelo
 * @date Ago 21, 2018
 * @brief Adds two float arrays and save the result in the first array
 * @param array1 Float array to be added and to store the result
 * @param array2 Float array to be added
 * @param nelems Size of the array
 * @return
 */
void iftSumFloatArraysInPlace(float *array1, float *array2, int nelems);

/**
 * @author Cesar Castelo
 * @date Jan 23, 2018
 * @brief Scale (multiply) a float array by a given value
 * @param array Float array to be scaled
 * @param value Value to scale the array
 * @param arrayOut Float array to store the result
 * @param nelems Size of the array
 * @return
 */
void iftScaleFloatArray(float* array, float value, float* arrayOut, int nelems);

/**
 * @author Cesar Castelo
 * @date Oct 4, 2019
 * @brief Substract two float arrays and save the result in a third array
 * @param array1 Float array to substract from
 * @param array2 Float array to be substracted
 * @param array3 Float array to store the result
 * @param nelems Size of the array
 * @return
 */
void iftSubstractFloatArrays(float *array1, float *array2, float *array3, int nelems);

/**
 * @author Cesar Castelo
 * @date Feb 01, 2018
 * @brief Merges two int arrays and stores the result in the first array
 * @param array1 Array to be merged (which will store the result)
 * @param array2 Array to be merged
 * @param n1 Size of array1
 * @param n2 Size of array2
 * @return
 */
void iftMergeIntArraysInPlace(int **array1, int *array2, int n1, int n2);


#ifdef __cplusplus
}
#endif

#endif //IFT_IFTMEMORY_H
