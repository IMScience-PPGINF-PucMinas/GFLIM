#include "iftMemory.h"

#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"

#ifdef __linux__
# include <sys/sysinfo.h>
#include <malloc.h>
#endif

#ifdef __APPLE__
# include <mach/task.h>
# include <mach/mach_init.h>
#endif

#ifdef _WINDOWS
# include <windows.h>
#else
# include <sys/resource.h>
#include <iftCommon.h>

#endif

void iftGC() {
    //GC_gcollect();
}

long iftMemoryUsed()
{

    iftGC();

#if defined(__linux__)
    struct mallinfo info;
    info = mallinfo();
    long mem = info.uordblks;
    return mem;

#elif defined(__APPLE__)
    // Inspired by:
    // http://miknight.blogspot.com/2005/11/resident-set-size-in-mac-os-x.html
    // By default, returns the full virtual arena, but if resident=true, it will report just the resident set in RAM (if supported on that OS).

    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
    long size = t_info.resident_size;
    return size;

#elif defined(_WINDOWS)
    // According to MSDN...
    PROCESS_MEMORY_COUNTERS counters;
    if (GetProcessMemoryInfo (GetCurrentProcess(), &counters, sizeof (counters)))
        return counters.PagefileUsage;
    else return 0;

#else
    // No idea what platform this is
    //are you running this library in a lasagna?
    return 0;   // Punt
#endif
}

long iftGetPhysicalSystemMemory()
{
#if defined(__WIN32) || defined(__WIN64)
    MEMORYSTATUSEX status;
		status.dwLength = sizeof(status);
		GlobalMemoryStatusEx(&status);
		return status.ullTotalPhys;
#elif defined(__linux__)
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
#endif // defined
    return 0;
}

long iftGetFreePhysicalSystemMemory()
{
#if defined(__WIN32) || defined(__WIN64)
    MEMORYSTATUSEX status;
		status.dwLength = sizeof(status);
		GlobalMemoryStatusEx(&status);
		return status.ullAvailPhys;
#elif defined(__linux__)
    //long pages = sysconf(_SC_AVPHYS_PAGES);
    //long page_size = sysconf(_SC_PAGE_SIZE);
    //return pages * page_size;
    return iftGetAvailableMemory();
#endif // defined
    return 0;
}


float iftGetAvailableMemory(){
    char * file_path = "/proc/meminfo"; 

    float AvailMem = 0.0;

    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    fp = fopen(file_path, "r");
    if (fp==NULL){
        iftError("/proc/meminfo Not found",
	       "iftGetAvailableMemory");
    }
    int line_count = 0;
    //get the third line, which corresponds to AvailableMem
    while (((read = getline(&line, &len, fp)) != -1) && line_count<2) {
        line_count++;
    }

    //Now it's time to parse only digits
    int power = 0;
    int tmpi=0;
    for(int i=read-2;i>=0;i--){
        if (line[i] <= 57 && line[i] >= 48){
            tmpi = (int) line[i] - 48;
            AvailMem+= tmpi*pow(10,power);
            power+=1;
        }
    }

    fclose(fp);
    if (line)
        free(line);

    return AvailMem*1024.0;
}


long objectCount = 0;

void iftPrintAllocObjectsCount() {
    printf("You have %lu current allocated objects in memory.\n", iftAllocObjectsCount());
}

long iftAllocObjectsCount() {
    return objectCount;
}

double iftMemoryUsedByArray(long arraySize, int elemSize, char unit)
{
    if(arraySize <= 0)
        iftError("The array must have size equal or greater than 1", "iftMemoryUsedByArray");

    double size = 0.0;
    switch(unit) {
        case 'b':
            size = (double)elemSize*(double)arraySize;
            break;
        
        case 'k':
            size = (double)elemSize*(double)arraySize/1024.0;
            break;

        case 'm':
            size = (double)(elemSize)*(double)arraySize/(1024.0*1024.0);
            break;

        case 'g':
            size = (double)(elemSize)*(double)arraySize/(1024.0*1024.0*1024.0);
            break;
        
        default:
            size = (double)elemSize*(double)arraySize;
    }
    return size;
}







void iftPrintMemoryFormatted(long mem) {

    if(mem<pow(2,20)) {
        printf("%4.2fKB", (float)mem/pow(2.0,10));
    }
    else if(mem<pow(2,30)) {
        printf("%4.2fMB", (float)mem/pow(2.0,20));
    }
    else {
        printf("%4.2fGB", (float)mem/pow(2.0,30));
    }

}

bool iftVerifyMemory(long MemDinInitial, long MemDinFinal)
{  
  if (MemDinInitial!=MemDinFinal){
    printf("\nDynamic memory was not completely deallocated (%ld, %ld). ",
	   MemDinInitial,MemDinFinal);
    iftPrintMemoryFormatted(MemDinFinal-MemDinInitial);
    printf(" not released.\n");
    return(false);
  }

  return(true);
}

void iftSwapString(char *a, char *b) {
    long size = strlen(a)+1;
    char c[size];
    strcpy(c,a);
    strcpy(a,b);
    strcpy(b,c);
}




///////////////////// DATATYPE ALLOCATION FUNCTIONS /////////////////////
bool *iftAllocBoolArray(long n) {
    bool *v = NULL;

    v = (bool *) iftAlloc(n, sizeof(bool));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocBoolArray");
    return (v);
}


char *iftAllocCharArray(long n) {
    char *v = NULL;

    v = (char *) iftAlloc(n, sizeof(char));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocCharArray");
    return (v);
}


char *iftAllocString(long n) {
    return (iftAllocCharArray(n+1));
}


uchar *iftAllocUCharArray(long n) {
    uchar *v = NULL;

    v = (uchar *) iftAlloc(n, sizeof(uchar));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocUCharArray");
    return (v);
}


short *iftAllocShortArray(long n) {
    short *v = NULL;

    v = (short *) iftAlloc(n, sizeof(short));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocShortArray");
    return(v);
}


ushort *iftAllocUShortArray(long n) {
    ushort *v = NULL;

    v = (ushort *) iftAlloc(n, sizeof(ushort));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocUShortArray");
    return(v);
}


int *iftAllocIntArray(long n) {
    int *v = NULL;

    v = (int *) iftAlloc(n, sizeof(int));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocIntArray");
    return(v);
}


uint *iftAllocUIntArray(long n) {
    uint *v = NULL;

    v = (uint *) iftAlloc(n, sizeof(uint));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocUIntArray");
    return(v);
}


long *iftAllocLongIntArray(long n) {
    long *v = NULL;

    v = (long *) iftAlloc(n, sizeof(long));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocLongIntArray");

    return v;
}

#ifndef  __cplusplus
long long *iftAllocLongLongIntArray(long n) {
    long long *v = NULL;

    v = (long long *) iftAlloc(n, sizeof(long long));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocLongLongIntArray");

    return v;
}
#endif

ulong *iftAllocULongArray(long n) {
    ulong *v = NULL;

    v = (ulong *) iftAlloc(n, sizeof(ulong));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocULongArray");

    return (v);
}

#ifndef  __cplusplus

ullong *iftAllocULLongArray(long n) {
    ullong *v = NULL;

    v = (ullong *) iftAlloc(n, sizeof(ullong));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocULLongArray");
    return (v);
}

#endif

float *iftAllocFloatArray(long n) {
    float *v = NULL;
    v = (float *) iftAlloc(n, sizeof(float));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocFloatArray");
    return(v);
}


double *iftAllocDoubleArray(long n) {
    double *v = NULL;

    v = (double *) iftAlloc(n, sizeof(double));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocDoubleArray");
    return (v);
}


long double *iftAllocLongDoubleArray(long n) {
    long double *v = NULL;

    v = (long double *) iftAlloc(n, sizeof(long double));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocLongDoubleArray");
    return(v);
}


iftComplex *iftAllocComplexArray(long n)  {
    iftComplex *v = NULL;

    v = (iftComplex *) iftAlloc(n, sizeof(iftComplex));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocComplexArray");
    return(v);
}


/////////////// ALIGNED MEMORY ALLOCATION /////////////////
uchar *iftAllocAlignedUCharArray(long n, long alignment) {
    uchar *mem = (uchar *) _mm_malloc(n * sizeof(uchar), alignment);
    memset(mem, 0, n * sizeof(uchar));

    return mem;
}


int *iftAllocAlignedIntArray(long n, long alignment) {
    int *mem = (int *) _mm_malloc(n * sizeof(int), alignment);
    memset(mem, 0, n * sizeof(int));

    return mem;
}


float *iftAllocAlignedFloatArray(long n, long alignment) {
    float *mem = (float *) _mm_malloc(n * sizeof(float), alignment);
    memset(mem, 0, n * sizeof(float));

    return mem;
}


double *iftAllocAlignedDoubleArray(long n, long alignment) {
    double *mem = (double *) _mm_malloc(n * sizeof(double), alignment);
    memset(mem, 0, n * sizeof(double));

    return mem;
}


///////////////////// DATATYPE ALLOCATION FUNCTIONS (MATRIX) /////////////////////
bool **iftAllocBoolMatrix(long c, long r) {
    bool **m = NULL;
    m = (bool **) iftAlloc(r, sizeof(bool*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocBoolMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (bool *) iftAlloc(c, sizeof(bool));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocBoolMatrix");
    }

    return (m);
}

char **iftAllocCharMatrix(long c, long r) {
    char **m = NULL;
    m = (char **) iftAlloc(r, sizeof(char*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocCharMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (char *) iftAlloc(c, sizeof(char));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocCharMatrix");
    }

    return (m);
}

uchar **iftAllocUCharMatrix(long c, long r) {
    uchar **m = NULL;
    m = (uchar **) iftAlloc(r, sizeof(uchar*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocUCharMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (uchar *) iftAlloc(c, sizeof(uchar));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocUCharMatrix");
    }

    return (m);
}

short **iftAllocShortMatrix(long c, long r) {
    short **m = NULL;
    m = (short **) iftAlloc(r, sizeof(short*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocShortMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (short *) iftAlloc(c, sizeof(short));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocShortMatrix");
    }

    return (m);
}

ushort **iftAllocUShortMatrix(long c, long r) {
    ushort **m = NULL;
    m = (ushort **) iftAlloc(r, sizeof(ushort*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocUShortMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (ushort *) iftAlloc(c, sizeof(ushort));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocUShortMatrix");
    }

    return (m);
}

int **iftAllocIntMatrix(long c, long r) {
    int **m = NULL;
    m = (int **) iftAlloc(r, sizeof(int*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocIntMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (int *) iftAlloc(c, sizeof(int));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocIntMatrix");
    }

    return (m);
}

uint **iftAllocUIntMatrix(long c, long r) {
    uint **m = NULL;
    m = (uint **) iftAlloc(r, sizeof(uint*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocUIntMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (uint *) iftAlloc(c, sizeof(uint));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocUIntMatrix");
    }

    return (m);
}

long **iftAllocLongMatrix(long c, long r) {
    long **m = NULL;
    m = (long **) iftAlloc(r, sizeof(long*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocLongMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (long *) iftAlloc(c, sizeof(long));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocLongMatrix");
    }

    return (m);
}

#ifndef  __cplusplus
long long **iftAllocLongLongMatrix(long c, long r) {
    long long **m = NULL;
    m = (long long **) iftAlloc(r, sizeof(long long*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocLongLongMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (long long *) iftAlloc(c, sizeof(long long));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocLongLongMatrix");
    }

    return (m);
}
#endif

ulong **iftAllocULongMatrix(long c, long r) {
    ulong **m = NULL;
    m = (ulong **) iftAlloc(r, sizeof(ulong*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocULongMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (ulong *) iftAlloc(c, sizeof(ulong));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocULongMatrix");
    }

    return (m);
}

#ifndef  __cplusplus
ullong **iftAllocULLongMatrix(long c, long r) {
    ullong **m = NULL;
    m = (ullong **) iftAlloc(r, sizeof(ullong*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocULLongMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (ullong *) iftAlloc(c, sizeof(ullong));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocULLongMatrix");
    }

    return (m);
}
#endif

float **iftAllocFloatMatrix(long c, long r) {
    float **m = NULL;
    m = (float **) iftAlloc(r, sizeof(float*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocFloatMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (float *) iftAlloc(c, sizeof(float));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocFloatMatrix");
    }

    return (m);
}

double **iftAllocDoubleMatrix(long c, long r) {
    double **m = NULL;
    m = (double **) iftAlloc(r, sizeof(double*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocDoubleMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (double *) iftAlloc(c, sizeof(double));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocDoubleMatrix");
    }

    return (m);
}

long double **iftAllocLongDoubleMatrix(long c, long r) {
    long double **m = NULL;
    m = (long double **) iftAlloc(r, sizeof(long double*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocLongDoubleMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (long double *) iftAlloc(c, sizeof(long double));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocLongDoubleMatrix");
    }

    return (m);
}

iftComplex **iftAllocComplexMatrix(long c, long r) {
    iftComplex **m = NULL;
    m = (iftComplex **) iftAlloc(r, sizeof(iftComplex*));
    if (m == NULL)
        iftError("Cannot allocate memory space", "iftAllocComplexMatrix");
    
    for(long i = 0; i < r; i++) {
        m[i] = (iftComplex *) iftAlloc(c, sizeof(iftComplex));
        if (m[i] == NULL)
            iftError("Cannot allocate memory space", "iftAllocComplexMatrix");
    }

    return (m);
}


///////////////////// IFT DATATYPE ALLOCATION FUNCTIONS /////////////////////










void iftSetFloatArray(float *array_dst, int nelems, float val) {
    int i;

    for (i = 0; i < nelems; i++) {
        array_dst[i] = val;
    }

}

void iftCopyFloatArray(float *array_dst, float *array_src, int nelems) {
//    int i;
//
//    for (i = 0; i < nelems; i++) {
//      array_dst[i] = array_src[i];
//    }

    memmove(array_dst, array_src, nelems*sizeof(float));
}

void iftCopyFloatArrayToDblArray(double *array_dst, float *array_src, int nelems) {
    #pragma omp parallel for
   for (int i = 0; i < nelems; i++) {
        array_dst[i] = array_src[i];
   }
}

void iftCopyDblArrayToFloatArray(float *array_dst, double *array_src, int nelems) {
    #pragma omp parallel for
   for (int i = 0; i < nelems; i++) {
	    array_dst[i] = array_src[i];
   }
}

void iftCopyDoubleArray(double *array_dst, double *array_src, int nelems) {
    #pragma omp parallel for
    for (int i = 0; i < nelems; i++)
        array_dst[i] = array_src[i];
}

void iftCopyIntArray(int *array_dst, const int *array_src, int nelems) {
    #pragma omp parallel for
    for (int i = 0; i < nelems; i++) {
        array_dst[i] = array_src[i];
    }
}

#ifndef  __cplusplus
void iftCopyLongLongIntArray(long long *array_dst, const long long *array_src, int nelems) {
    #pragma omp parallel for
    for (int i = 0; i < nelems; ++i) {
        array_dst[i] = array_src[i];
    }
}
#endif

void iftCopySizeTArray(long *array_dst, const long *array_src, long nelems) {
    long i;
    for (i = 0; i < nelems; i++) {
        array_dst[i] = array_src[i];
    }

}








int *iftConcatIntArray(int *array1, int n1, int *array2, int n2, int *nelems) {
    int m = n1+n2;
    int *v = iftAllocIntArray(m);

    for (int i = 0; i < n1; i++)
        v[i] = array1[i];

    for (int i = n1; i < m; i++) {
        int idx = i - n1;
        v[i] = array2[idx];
    }

    *nelems = m;

    return v;
}

void iftSetIntArray(int *array_dst, int nelems, int val) {
    int i;

    for (i = 0; i < nelems; i++) {
        array_dst[i] = val;
    }

}









void iftInitFloatArray(float *array, float value, int nelems)
{
    int i;

    for (i=0; i < nelems; i++)
        array[i]=value;
}

void iftInitFloatArrayPos(float *array, int pos, float value, int nelems)
{
    iftInitFloatArray(array,0.0,nelems);
    array[pos]=value;
}

void iftSumFloatArrays(float *array1, float *array2, float *array3, int nelems)
{
    int i;
    for (i=0; i < nelems; i++)
        array3[i] = array1[i] + array2[i];
}

void iftSumFloatArraysInPlace(float *array1, float *array2, int nelems)
{
    int i;
    for (i=0; i < nelems; i++)
        array1[i] = array1[i] + array2[i];
}

float iftSumFloatArrayElems(float *farr, int n_elems) {
    float sum = 0.0f;

    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n_elems; i++)
        sum += farr[i];

    return sum;
}

void iftSubstractFloatArrays(float *array1, float *array2, float *array3, int nelems)
{
    int i;
    for (i=0; i < nelems; i++)
        array3[i] = array1[i] - array2[i];
}

void iftScaleFloatArray(float* array, float value, float* arrayOut, int nelems)
{
    int i=0;
    for (i = 0; i < nelems; ++i)
        arrayOut[i] = array[i]*value;
}

void iftMergeIntArraysInPlace(int **array1, int *array2, int n1, int n2)
{
    int *auxArray = NULL;
    if((*array1) != NULL) {
        auxArray = iftAllocIntArray(n1);
        iftCopyIntArray(auxArray, (*array1), n1);
        iftFree(*array1);
    }
    (*array1) = iftAllocIntArray(n1+n2);

    int c = 0;
    for(int i = 0; i < n1; i++) {
        (*array1)[c] = auxArray[i];
        c++;
    }

    for(int i = 0; i < n2; i++) {
        (*array1)[c] = array2[i];
        c++;
    }

    iftFree(auxArray);
}
