#include "iftCommon.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/io/Stream.h"

#include <iftSort.h>
//#include <MacTypes.h>

#ifdef IFT_GPU
#include <iftMemory.cuh>
#endif

/**
@file
@brief A description is missing here
*/
/************ PRIVATE ************/
//Used in qsort
int compare (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

void iftPlatformInfo() {

#ifdef __unix__
    int cpunumber = sysconf( _SC_NPROCESSORS_ONLN );
    char cpuname[1024];
    FILE* out = popen("lscpu | grep \"Model name\" | sed \"s/Model name://g;s/^[ \\t]*//g\"", "r");
    if(fgets(cpuname, 1024, out)) {
        cpuname[strlen(cpuname)-1]='\0';//remove fgets newline
        printf("CPU: %s. %d cores available. %.fGB memory available\n\n", cpuname, cpunumber, iftGetPhysicalSystemMemory()/pow(1024.0,3));
    }
    else
        printf("CPU: Could not retrive info.\n\n");
#else
    printf("CPU: Could not retrive info. (Platform not supported)\n");
#endif


#ifdef IFT_GPU
    iftStartGPU(0);
#endif
}

/*! \file
 * Author: Alexandre Falcao
 * Date: Aug 22, 2011
 * Last Update: Aug 22, 2011
 */

// based on: https://en.wikipedia.org/wiki/Primality_test
bool iftIsPrime(long n) {
    if (n < 0)
        iftError("Number is Negative (not Natural): %ld... Try Natural Numbers", "iftIsPrime", n);

    if (n <= 1)
        return false;
    else
    if ((n == 2) || (n == 3))
        return true;
    else
    if (((n % 2) == 0) || ((n % 3) == 0))
        return false;

    long sqrt_of_n = (long) sqrt(n); // floor of sqrt(n)

    for (long d = 5; d <= sqrt_of_n; d = d + 6)
        if (((n % d) == 0) || ((n % (d+2)) == 0))
            return false;

    return true;
}


double iftFastPow(double base, unsigned b) {
    double result = 1;
    double t;

    switch (b) {
        case 1:
            return base;
        case 2:
            return base*base;
        case 3:
            return base*base*base;
        case 4:
            t = base*base;
            return t*t;
        case 5:
            t = base*base;
            return t*t*base;
        case 6:
            t = base*base;
            return t*t*t;
        case 7:
            t = base*base*base;
            return t*t*base;
        case 8:
            t = base*base;
            t = t*t;
            return t*t;
        case 9:
            t = base*base;
            t = t*t;
            return t*t*base;
        case 10:
            t = base*base;
            t = t*t;
            t = t*t;
            return t*base*base;
        default:
            return pow(base, b);
    }

    return result;
}

void iftPrintBoolArray(bool *v, long n)
{
    for(long i = 0; i < n; i++)
        fprintf(stdout,"%d ", v[i]);
    fprintf(stdout,"\n");
}

void iftPrintFormattedBoolArray(bool *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%d", v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s",suffix);
}

void iftPrintCharArray(char *v, long n)
{
    for(long i = 0; i < n; i++)
        fprintf(stdout,"%c ", v[i]);
    fprintf(stdout,"\n");
}

void iftPrintFormattedCharArray(char *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%c", v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s",suffix);
}

void iftPrintUCharArray(uchar *v, long n)
{
    for(long i = 0; i < n; i++)
        printf("%X ", v[i]); // prints in hexadecimal format
    printf("\n");
}

void iftPrintFormattedUCharArray(uchar *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%X",v[i]); // prints in hexadecimal format
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s", suffix);
}

void iftPrintShortArray(short *v, long n)
{
    for(long i = 0; i < n; i++)
        fprintf(stdout,"%d ", v[i]);
    fprintf(stdout,"\n");
}

void iftPrintFormattedShortArray(short *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%d", v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s",suffix);
}

void iftPrintUShortArray(ushort *v, long n)
{
    for(long i = 0; i < n; i++)
        fprintf(stdout,"%d ", v[i]);
    fprintf(stdout,"\n");
}

void iftPrintFormattedUShortArray(ushort *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%d", v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s",suffix);
}

void iftPrintIntArray(int *v, long n)
{
    for(long i = 0; i < n; i++)
        fprintf(stdout,"%d ", v[i]);
    fprintf(stdout,"\n");
}

void iftPrintFormattedIntArray(int *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%d", v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s",suffix);
}

void iftPrintUIntArray(uint *v, long n)
{
    for(long i = 0; i < n; i++)
        fprintf(stdout,"%d ", v[i]);
    fprintf(stdout,"\n");
}

void iftPrintFormattedUIntArray(uint *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%d", v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s",suffix);
}

void iftPrintLongArray(long *v, long n)
{
    for(long i = 0; i < n; i++)
        printf("%ld ", v[i]);
    printf("\n");
}

void iftPrintFormattedLongArray(long *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%ld",v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s", suffix);
}

#ifndef  __cplusplus
void iftPrintLongLongIntArray(long long *v, long n)
{
    for(long i = 0; i < n; i++)
        printf("%lld ", v[i]);
    printf("\n");
}

void iftPrintFormattedLongLongIntArray(long long *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%lld",v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s", suffix);
}
#endif

void iftPrintULongArray(ulong *v, long n)
{
    for(long i = 0; i < n; i++)
        printf("%lu ", v[i]);
    printf("\n");
}

void iftPrintFormattedULongArray(ulong *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%lu",v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s", suffix);
}

#ifndef  __cplusplus
void iftPrintULLongArray(ullong *v, long n)
{
    for(long i = 0; i < n; i++)
        printf("%llu ", v[i]);
    printf("\n");
}

void iftPrintFormattedULLongArray(ullong *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%llu",v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s", suffix);
}
#endif

void iftPrintFloatArray(float *v, long n)
{
    for(long i = 0; i < n; i++)
        fprintf(stdout,"%f ", v[i]);
    fprintf(stdout,"\n");
}

void iftPrintFormattedFloatArray(float *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%.4f",v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s", suffix);
}

void iftPrintDoubleArray(double *v, long n)
{
    for(long i = 0; i < n; i++)
        fprintf(stdout,"%lf ", v[i]);
    fprintf(stdout,"\n");
}

void iftPrintFormattedDoubleArray(double *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%.4lf",v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s", suffix);
}

void iftPrintLongDoubleArray(long double *v, long n)
{
    for(long i = 0; i < n; i++)
        fprintf(stdout,"%Lf ", v[i]);
    fprintf(stdout,"\n");
}

void iftPrintFormattedLongDoubleArray(long double *v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        fprintf(stdout,"%Lf",v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s", suffix);
}

void iftPrintStrArray(char **v, long n)
{
    for(long i = 0; i < n; i++)
        printf("%s ", v[i]);
    printf("\n");
}

void iftPrintFormattedStrArray(char **v, long n, const char *prefix, const char *suffix, const char *separator)
{
    printf("%s",prefix);
    for(long i = 0; i < n; i++) {
        printf("%s", v[i]);
        if(i < n-1)
            printf("%s", separator);
    }
    printf("%s", suffix);
}

float iftRandomUniform(float low, float high) {
    double d;
    d = ((double) rand()) / ((double) RAND_MAX);
    return low + d * (high - low);
}


/*! \brief
 * Randomly selects an integer from low to high
 */

int iftRandomInteger (int low, int high){
    int k;
    double d;

    d = (double) rand () / ((double) RAND_MAX + 0.5);
    k =  iftMin((int)(d * (high - low + 1.0) ) + low, high);

    return k;
}

/*
 * Randomly selects nelems of the set [low, high]
 */
int *iftRandomIntegers(int low, int high, int nelems) {
    char msg[512];

    if (low > high) {
        sprintf(msg, "Low is greater than High: (%d, %d)", low, high);
        iftError(msg, "iftRandomIntegers");
    }

    int total_of_elems = (high - low + 1);

    if (nelems > total_of_elems) {
        sprintf(msg, "Nelems = %d is greater than the total of integer number in the range: [%d, %d]",
                nelems, low, high);
        iftError(msg, "iftRandomIntegers");
    }

    int *selected = iftAllocIntArray(nelems);
    int *values   = iftAllocIntArray(total_of_elems);
    int *count    = iftAllocIntArray(total_of_elems);

    int t = 0;
    for (int i = low; i <= high; i++) {
        values[t] = i;
        count[t]  = 100;
        t++;
    }

    if (nelems == total_of_elems) {
        iftFree(count);
        iftFree(selected);

        return values;
    }

    // Randomly select samples
    t = 0;
    int roof = total_of_elems - 1;
    while (t < nelems) {
        int i = iftRandomInteger(0, roof);
        int v = values[i];

        if (count[i] == 0) {
            selected[t] = v;
            iftSwap(values[i], values[roof]);
            iftSwap(count[i], count[roof]);
            t++;
            roof--;
        } else {
            count[i]--;
        }
    }
    iftFree(values);
    iftFree(count);

    return selected;
}

float iftRandomNormalFloatAux(void) {
    float v1,v2,s;

    do {
        v1 = 2.0 * ((float) rand()/RAND_MAX) - 1;
        v2 = 2.0 * ((float) rand()/RAND_MAX) - 1;

        s = v1*v1 + v2*v2;
    } while ( s >= 1.0 );

    if (s == 0.0)
        return 0.0;
    else
        return (v1*sqrt(-2.0 * log(s) / s));
}

/*! \brief
 * Randomly selects a normal distributed (N(0,1)) float number
 */

float iftRandomNormalFloat(float mean, float variance) {
    return variance*iftRandomNormalFloatAux() + mean;
}


double iftRandomNormalDouble(){
    double d;
    d = (double) rand () / ((double) RAND_MAX + 1);
    return d;
}


double iftLog(double val, double base) {
    return (log(val) / log(base));
}


/*! \brief
 * Returns the initial time
 */

timer *iftTic(){
    timer *tic=NULL;
    tic = (timer *)iftAlloc(1, sizeof(timer));
    gettimeofday(tic,NULL);
    return(tic);
}

/*! \brief
 * Returns the final time
 */

timer *iftToc(){
    timer *toc=NULL;
    toc = (timer *)iftAlloc(1, sizeof(timer));
    gettimeofday(toc,NULL);
    return(toc);
}

/*! \brief
 * Prints the computational time from tstart to tend in ms.
 */


void iftPrintCompTime(timer *tstart, timer *tend, const char *msg, ...)
{
    va_list args;
    char    final_msg[IFT_STR_DEFAULT_SIZE];
    float   t=iftCompTime(tstart,tend);

    va_start(args,msg);
    vsprintf(final_msg, msg, args);
    va_end(args);

    printf("\n %s: %f ms\n",final_msg,t);
}


/*! \brief
 * Writes a timer to a given file, using the following format %s: %f, where %s is the
 * given information corresponding to the time and %lf is the current time in milliseconds.
 * @param tic is freed by the function.
 */
void iftWriteTimerToFile(const char *filename, const char *information, timer *tic)
{
    FILE *f = fopen(filename, "a+");
    fprintf(f,"%s: %lf\n", information, (tic->tv_sec)*1000.0 + (tic->tv_usec)*0.001);
    fclose(f);
    iftFree(tic);
}

/*! \brief
 * Returns the difference from the initial to the final times
 */

float iftCompTime(timer *tic, timer *toc)
{
    float t=0.0;
    if ((tic!=NULL)&&(toc!=NULL)){
        t = (toc->tv_sec-tic->tv_sec)*1000.0 +
            (toc->tv_usec-tic->tv_usec)*0.001;
        iftFree(tic);iftFree(toc);
    }
    return(t);
}


char *iftFormattedTime(float runtime) {
    char *formatted = iftAllocCharArray(128);

    if (runtime < 1.0)
        sprintf(formatted, "%f ms", runtime);
    else {
        char msg[512];
        int days, hours, mins, secs;
        int ms = runtime; // ignores decimal miliseconds

        days  = ms/86400000;
        ms    = ms%86400000;

        hours = ms/3600000;
        ms    = ms%3600000;

        mins  = ms/60000;
        ms    = ms%60000;

        secs  = ms/1000;
        ms    = ms%1000;

        // builds the formatted time
        if (days) {
            sprintf(msg, "%d day(s) ", days);
            strcat(formatted, msg);
        }
        if (hours) {
            sprintf(msg, "%d hour(s) ", hours);
            strcat(formatted, msg);
        }
        if (mins) {
            sprintf(msg, "%d min(s) ", mins);
            strcat(formatted, msg);
        }
        if (secs) {
            sprintf(msg, "%d sec(s) ", secs);
            strcat(formatted, msg);
        }
        if (ms) {
            sprintf(msg, "%d ms", ms);
            strcat(formatted, msg);
        }
    }

    return formatted;
}


void iftPrintFormattedTime(float runtime) {
    char *form_time = iftFormattedTime(runtime);
    fprintf(stdout, "-> Time elapsed: %s\n", form_time);
    iftFree(form_time);
}


/*! \brief
 * Generates a new seed for rand(), used in iftRandomInteger()
 */

void iftRandomSeed(unsigned int seed)
{

    srand(seed);

}

/*! \brief
 * Returns the factorial of a number or NIL in case of overflow
 */

long double iftFactorial(int n)
{
    long double fact=1;
    int i;


    for (i=1; i <= n; i++) {
        if (fact > (IFT_INFINITY_LDBL / i)){ // overflow
            return(IFT_NIL);
        }
        fact=fact*i;
    }
    return(fact);
}

/*! \brief
 * Returns the limit to avoid overflow in the factorial computation
 */

int iftFactorialLimit()
{
    long double fact=1.0;
    int i;

    i = 1;
    while (fact < (IFT_INFINITY_LDBL / i)){
        fact=fact*i;
        i++;
    }
    return(i-1);
}

iftVector iftVectorScalarDiv(iftVector vec, double s) {
    return ((!iftAlmostZero(s))? (iftVector) {(vec).x / (s), (vec).y / (s), (vec).z / (s)} : vec);
}


bool iftVoxelsAreEqual(iftVoxel u1, iftVoxel u2)
{
    return ((u1.x==u2.x)&&(u1.y==u2.y)&&(u1.z==u2.z));
}

iftVoxel iftMeanVoxel(iftVoxel u1, iftVoxel u2) {
    iftVoxel v;

    v.x = (u1.x + u2.x) / 2;
    v.y = (u1.y + u2.y) / 2;
    v.z = (u1.z + u2.z) / 2;

    return v;
}

iftVector iftNormalizeVector(iftVector v)
{
    iftVector u;
    float     m = iftVectorMagnitude(v);

    u.x = v.x; u.y = v.y; u.z = v.z;

    if (!iftAlmostZero(m)){
        u.x = v.x/m;
        u.y = v.y/m;
        u.z = v.z/m;
    }

    return(u);
}


iftVector iftProjectVector(iftVector U, iftVector V)
{
    float V_norm = iftVectorMagnitude(V);
    iftVector V_n = iftNormalizeVector(V);

    if(iftAlmostZero(V_norm))
        iftError("Impossible to project 0-length vector U onto V",
                 "iftProjectVector");

    return (iftVector)iftVectorScalarProd(V, iftVectorInnerProduct(U, V_n) / V_norm);
}


iftVector iftOrthogonalizeVector(iftVector u, iftVector v) {
    float u_mag = iftVectorMagnitude(u);
    float v_mag = iftVectorMagnitude(v);

    if (iftAlmostZero(u_mag) || iftAlmostZero(v_mag))
        iftWarning("Vector u or Vector v is a zero-vector. Resulting orthogonal vector is a zero-vecto", "iftOrthogonalizeVector");

    iftVector proj_u_on_v = iftProjectVector(u, v);
    iftVector u_orth = iftVectorSub(u, proj_u_on_v);

    float u_orth_mag = iftVectorMagnitude(u_orth);
    if (iftAlmostZero(u_orth_mag))
        iftWarning("Vector u and v are parallel. Resulting orthogonal vector is a zero-vector", "iftOrthogonalizeVector");

    float inner_prod = iftVectorInnerProduct(u_orth, v);

    double epsilon = 1e-4;
    // if the inner product is not almost zero (considering a larger epsilon than that used on function iftAlmostZero)
    if (!((inner_prod >= -epsilon) && (inner_prod <= epsilon)))
        iftError("Error when orthogonalizing vector u onto v. Resulting orthogonal vector is not actually orthogonal",
            "iftOrthogonalizeVector");


    return u_orth;
}


double iftVoxelLineDist2D(iftVoxel P0, iftVoxel P1, iftVoxel P2, double P1P2)
{
    return(abs((P2.x - P1.x)*(P1.y - P0.y) - (P1.x - P0.x)*(P2.y - P1.y))/P1P2);
}

iftVoxel iftClosestVoxelOnLine2D(iftVoxel P0, iftVoxel P1, iftVoxel P2)
{
    iftVector Vproj;
    iftVector V12 = iftVectorSub(P2, P1);
    iftVector V10 = iftVectorSub(P0, P1);

    // We project the vector P1->P0 on P1->P2
    Vproj = iftProjectVector(V10, V12);

    // And then convert the result into the closest voxel
    // on the line segment P1P2 voxel
    return (iftVoxel)iftVectorRound((iftVector) iftVectorSum(Vproj, P1));
}


int iftVoxelLinePosition2D (iftVoxel P0, iftVoxel P1, iftVoxel P2)
{
    return (P2.x - P1.x) * (P0.y - P1.y) - (P2.y - P1.y) * (P0.x - P1.x);
}


void iftRemoveCarriageReturn(char *line)
{
    int i;

    for (i=0; i < strlen(line); i++){
        if ((line[i] == '\n') || (line[i] == '\r')){
            line[i]='\0';
        }
    }
}

void iftSwitchVoxels(iftVoxel *u, iftVoxel *v)
{
    iftVoxel w;

    w  = *v;
    *v = *u;
    *u = w;
}

float iftIncrementalMean(float mean, int n, float x) {
    float diff = (float)(n)/(n+1);
    return mean*diff + x/(n+1);
}

float iftDecrementalMean(float mean, int n, float x) {
    double diff = (double)(n)/(n-1);
    return mean*diff - x/(n-1);
}

float iftIncrementalVariance(float mean, float var, int n, float x) {
    double e2x = iftFastPow(mean, 2);//square of x expectation
    double ex2 = var + e2x;//expectation of squared x

    e2x = iftFastPow(iftIncrementalMean(mean, n, x), 2);
    ex2 = iftIncrementalMean(ex2, n, iftFastPow(x,2));

    return (float)(ex2 - e2x);

}

float iftDecrementalVariance(float mean, float var, int n, float x) {
    double e2x = iftFastPow(mean, 2);//square of x expectation
    double ex2 = var + e2x;//expectation of squared x

    e2x = iftFastPow(iftDecrementalMean(mean, n, x), 2);
    ex2 = iftDecrementalMean(ex2, n, iftFastPow(x,2));

    return (float)(ex2 - e2x);
}

float iftSigmoidalValue(float value, float alfa)
{
    float result;

    result = 1 / (1 + exp(-alfa * value));

    return(result);
}

void iftVerifyToken(FILE *fp, char *token, char *function)
{
    char label[30];

    if ((fscanf(fp,"%s",label)!=1)||(strcmp(label,token)!=0)){
        char msg[100];
        sprintf(msg,"Could not find token %s",token);
        iftError(msg, function);
    }

}

void iftReadIntValue(FILE *fp, int *value, char *token, char *function)
{

    iftVerifyToken(fp, token, function);
    if (fscanf(fp,"%d\n",value)!=1){
        char msg[100];
        sprintf(msg,"Could not read %s value",token);
        iftError(msg, function);
    }
}

void iftReadIntValues(FILE *fp, int **value, int nvalues, char *token, char *function)
{

    iftVerifyToken(fp, token, function);
    for (int i = 0; i < nvalues; i++) {
        if (fscanf(fp,"%d",&(*value)[i])!=1){
            char msg[100];
            sprintf(msg,"Could not read %s value",token);
            iftError(msg, function);
        }
    }
    if (fscanf(fp,"\n")!=0){
        iftError("Could not read \\n in file", "iftReadIntValues");
    }
}

void iftWriteIntValue(FILE *fp, int value, char *token)
{
    fprintf(fp,"%s %d\n",token,value);
}

void iftWriteIntValues(FILE *fp, int *value, int nvalues, char *token)
{
    fprintf(fp,"%s",token);
    for (int i=0; i < nvalues; i++) {
        fprintf(fp," %d",value[i]);
    }
    fprintf(fp,"\n");
}

void iftReadFloatValue(FILE *fp, float *value, char *token, char *function)
{
    iftVerifyToken(fp, token, function);
    if (fscanf(fp,"%f\n",value)!=1){
        char msg[100];
        sprintf(msg,"Could not read %s value",token);
        iftError(msg, function);
    }
}

void iftReadFloatValues(FILE *fp, float **value, int nvalues, char *token, char *function)
{

    iftVerifyToken(fp, token, function);
    for (int i = 0; i < nvalues; i++) {
        if (fscanf(fp,"%f",&(*value)[i])!=1){
            char msg[100];
            sprintf(msg,"Could not read %s value",token);
            iftError(msg, function);
        }
    }
    if (fscanf(fp,"\n")!=0){
        iftError("Could not read \\n in file", "iftReadFloatValues");
    }
}

void iftWriteFloatValue(FILE *fp, float value, char *token)
{
    fprintf(fp,"%s %f\n",token,value);
}

void iftWriteFloatValues(FILE *fp, float *value, int nvalues, char *token)
{
    fprintf(fp,"%s",token);
    for (int i=0; i < nvalues; i++) {
        fprintf(fp," %f",value[i]);
    }
    fprintf(fp,"\n");
}

void iftReadDoubleValue(FILE *fp, double *value, char *token, char *function)
{

    iftVerifyToken(fp, token, function);
    if (fscanf(fp,"%lf\n",value)!=1){
        char msg[100];
        sprintf(msg,"Could not read %s value",token);
        iftError(msg, function);
    }
}

void iftReadDoubleValues(FILE *fp, double **value, int nvalues, char *token, char *function)
{

    iftVerifyToken(fp, token, function);
    for (int i = 0; i < nvalues; i++) {
        if (fscanf(fp,"%lf",&(*value)[i])!=1){
            char msg[100];
            sprintf(msg,"Could not read %s value",token);
            iftError(msg, function);
        }
    }
    if (fscanf(fp,"\n")!=0){
        iftError("Could not read \\n in file", "iftReadDoubleValues");
    }
}

void iftWriteDoubleValue(FILE *fp, double value, char *token)
{
    fprintf(fp,"%s %lf\n",token,value);
}

void iftWriteDoubleValues(FILE *fp, double *value, int nvalues, char *token)
{
    fprintf(fp,"%s",token);
    for (int i=0; i < nvalues; i++) {
        fprintf(fp," %lf",value[i]);
    }
    fprintf(fp,"\n");
}

void iftSkipComments(FILE *fp)
{
    //skip for comments
    while (fgetc(fp) == '#') {
        while (fgetc(fp) != '\n');
    }
    fseek(fp,-1,SEEK_CUR);

}


void iftUnitNorm(float *feats, int nelems){
    int i;
    float sum = 0.0;

    for (i = 0; i < nelems; ++i) {
        sum+= feats[i]*feats[i];
    }
    sum = sqrtf(sum);
    if (sum > 0.001)
      for (i = 0; i < nelems; ++i) {
        feats[i]/=sum;
      }
}

void iftNormalizeFeatures(float *feats, int nelems)
{
    int i;
    float maximum = IFT_INFINITY_FLT_NEG;
    for (i = 0; i < nelems; ++i)
    {
        maximum = iftMax(maximum, feats[i]);
    }

    for (i = 0; i < nelems; ++i)
    {
        feats[i]/=maximum;
    }
}

void iftStandardizeFeatures(float *feats, int nelems)
{
    int i;
    float mean = 0;
    for (i = 0; i < nelems; ++i)
    {
        mean += feats[i];
    }
    mean /= (float)nelems;

    float stdev = 0;
    for (i = 0; i < nelems; ++i)
    {
        stdev += powf(feats[i] - mean, 2.0);
    }
    stdev /= (float)nelems;
    stdev = sqrtf(stdev);

    for (i = 0; i < nelems; ++i)
    {
        feats[i] = (feats[i] - mean) / stdev;
    }
}

iftIntArray *iftIntArrayUnique(const int *array, int n) {
    if (array == NULL)
        iftError("Array is NULL", "iftIntArrayUnique");
    if (n <= 0)
        iftError("Number of Elements %d is <= 0", "iftIntArrayUnique", n);

    iftIntArray *v = iftCreateIntArray(n);
    iftCopyIntArray(v->val, array, n);

    qsort(v->val, n, sizeof(int), compare);

    int m = 1;

    for (int i = 1; i < n; i++) {
        if (v->val[i] > v->val[m - 1]) {
            v->val[m++] = v->val[i];
        }
    }

    v->val = (int*) iftRealloc(v->val, m * sizeof(int));
    v->n = m;

    return v;
}


int iftCountUniqueIntElems(const int *array, int n) {
    int *copy = iftAllocIntArray(n);
    iftCopyIntArray(copy, array, n);

    iftIntArray *unique = iftIntArrayUnique(copy, n);
    int count           = unique->n;

    iftDestroyIntArray(&unique);
    iftFree(copy);

    return count;
}

int iftCountElemsPerUniqueIntElem(int *array, int n, int **elems_out, int **quant_out)
{
	int maxVal = iftMaxValIntArray(array, n, NULL);
	int *count = iftAllocIntArray(maxVal + 1);

    /* count the number of repetitions for each element */
	for(int i = 0; i < n; i++) {
        count[array[i]]++;
	}

    /* count the ones that are actually present */
    int numbElems = 0;
    for(int i = 0; i < maxVal+1; i++)
        if(count[i] > 0)
            numbElems++;

    /* create the return arrays */
    (*elems_out) = iftAllocIntArray(numbElems);
    (*quant_out) = iftAllocIntArray(numbElems);
    int s = 0;
    for (int i = 0; i < maxVal+1; i++) {
        if (count[i] > 0) {
            (*elems_out)[s] = i;
            (*quant_out)[s] = count[i];
            s++;
        }
    }

	return numbElems;
}


bool iftIntArrayContainValue(const int *array, int n, int value) {
    for (size_t i = 0; i < n; i++) {
        if (array[i] == value)
            return true;
    }

    return false;
}


void iftWriteRawIntArray(char *filename, int *array, int n){
    FILE *fp  = fopen(filename,"wb");
    fwrite(array,sizeof(int),n,fp);
    fclose(fp);
}

int* iftReadRawIntArray(char *filename, int n){
    FILE *fp = fopen(filename, "rb");
    int *array = iftAllocIntArray(n);
    if(fread(array, sizeof(int),n,fp) != n)
        iftError("Reading error", "iftReadRawIntArray");
    fclose(fp);
    return array;
}


float iftIntMean(const int *x, int n) {
    float x_mean = 0.0;
    
    for (int i = 0; i < n; i++)
        x_mean += x[i];
    x_mean /= (float) n;
    
    return x_mean;
}


/**
 * @brief Compute the Mean value of the float vector x.
 *
 * @param x Float vector.
 * @param n Size of the vector.
 * @return Mean value.
*/
float iftMean(float *x, int n) {
    float x_mean = 0.0;

    for (int i = 0; i < n; i++)
        x_mean += x[i];
    x_mean /= (float) n;

    return x_mean;
}

/**
 * @brief Compute the Variance value of the float vector x.
 *
 * @param x Float vector.
 * @param n Size of the vector.
 * @return Variance value.
*/
float iftVar(float *x, int n) {
    float x_mean = iftMean(x, n);
    float var = 0.0;

    for (int i = 0; i < n; i++)
        var += (x[i]-x_mean)*(x[i]-x_mean);
    var /= (float) n;

    return var;
}


float iftStd(float *x, int n) {
    return sqrtf(iftVar(x, n));
}



/**
 * @brief Compute the Covariance of the float vectors x and y.
 *
 * @param x Float vector.
 * @param y Float vector.
 * @param n Size of the vectors.
 * @return Covariance value.
*/
float iftCov(float *x, float *y, int n) {
    float x_mean = iftMean(x, n);
    float y_mean = iftMean(y, n);
    float cov = 0.0;

    for (int i = 0; i < n; i++)
        cov += (x[i]-x_mean)*(y[i]-y_mean);
    cov /= (float) n;

    return cov;
}

int iftAlmostZero(double x){
    return (x >= -IFT_EPSILON) && (x <= IFT_EPSILON);
}

inline int iftSafeMod(int a, int n)
{
    int r = a % n;

    return (r >= 0) ? r : n+r;
}

float iftSquaredFeatDistance(float *A, float *B, int n)
{
    float dist=0.0;
    int    i;
    for (i=0; i < n; i++)
        dist += (A[i]-B[i])*(A[i]-B[i]);

    return(dist);
}

float iftFeatDistance(float *A, float *B, int n)
{
    float dist=0.0;
    int    i;
    for (i=0; i < n; i++)
        dist += (A[i]-B[i])*(A[i]-B[i]);

    return(sqrt(dist));
}


float iftSignedManhattanDistance(float *A, float *B, int n) {
    float dist = 0;
    for (int i = 0; i < n; i++)
        dist += (B[i] - A[i]);

    return dist;
}


long iftNormalizationValue(long maxval) {
    long norm_val = 1;
    
    if (maxval < 0)
        iftError("Input value %ld < 0", "iftNormalizationValue", maxval);
    else if (maxval <= 1)
        norm_val = 1;
    else if (maxval <= 255)
        norm_val = 255;
    else if (maxval <= 4095)
        norm_val = 4095;
    else if (maxval <= 65535)
        norm_val = 65535;
    else if (maxval <= 4294967295)
        norm_val = 4294967295;
    else iftError("Invalid maxval number %ld with number of bits > 32. It only supports values within [0, 2ˆn_bits -1], " \
                  "where n_bits in {1, 8, 12, 16, 32}", "iftNormalizationValue", maxval);
    
    return norm_val;
}


int iftArgMin(const int *x, long n)
{
    int i, best_i = IFT_NIL;
    int x_min = IFT_INFINITY_INT;

    for(i = 0; i < n; i++)
    {
        if(x[i] < x_min)
        {
            best_i = i;
            x_min = x[i];
        }
    }

    return best_i;
}

int iftArgMax(const int *x, long n)
{
    int i, best_i = IFT_NIL;
    int x_max = IFT_INFINITY_INT_NEG;

    for(i = 0; i < n; i++)
    {
        if(x[i] > x_max)
        {
            best_i = i;
            x_max = x[i];
        }
    }

    return best_i;
}

int iftDArgmin(const double *x, int n)
{
    int i, best_i = IFT_NIL;
    double x_min = IFT_INFINITY_DBL;

    for(i = 0; i < n; i++)
    {
        if(x[i] < x_min)
        {
            best_i = i;
            x_min = x[i];
        }
    }

    return best_i;
}

int iftDArgmax(const double *x, int n)
{
    int i, best_i = IFT_NIL;
    double x_max = IFT_INFINITY_DBL_NEG;

    for(i = 0; i < n; i++)
    {
        if(x[i] > x_max)
        {
            best_i = i;
            x_max = x[i];
        }
    }

    return best_i;
}

int iftFArgmin(const float *x, int n)
{
    int i, best_i = IFT_NIL;
    float x_min = IFT_INFINITY_FLT;

    for(i = 0; i < n; i++)
    {
        if(x[i] < x_min)
        {
            best_i = i;
            x_min = x[i];
        }
    }

    return best_i;
}

int iftFArgmax(const float *x, int n)
{
    int i, best_i = IFT_NIL;
    float x_max = IFT_INFINITY_FLT_NEG;

    for(i = 0; i < n; i++)
    {
        if(x[i] > x_max)
        {
            best_i = i;
            x_max = x[i];
        }
    }

    return best_i;
}

int iftFindIntArrayElementIndex(int *x, int n, int elem) {
    for (int idx = 0; idx < n; ++idx) {
        if(x[idx]==elem)
            return idx;
    }
    return IFT_NIL;
}



double *iftGaussPDF(double *vals, size_t n_vals, double mean, double stdev) {
    double variance = stdev * stdev;
    double denom    = (sqrtf(2.0 * IFT_PI * variance));
    double var_x_2  = 2.0 * variance;

    double *pdf = iftAllocDoubleArray(n_vals);

#pragma omp parallel for
    for (size_t i = 0; i < n_vals; i++)
        pdf[i] = exp(-( powf(vals[i]-mean, 2) / var_x_2)) / denom;

    return pdf;
}


double iftGaussCDF(double x, double mean, double stdev) {
    double area = 0.5 * (1 + erf((x - mean) / (stdev * sqrtf(2.0))));

    return area;
}


int iftMedianIntVal(const iftIntArray *arr) {
    iftIntArray *copy = iftCreateIntArray(arr->n);
    iftCopyIntArray(copy->val, arr->val, arr->n);
    iftIntArray *idx = iftIntRange(0, copy->n-1, 1);

    iftQuickSort(copy->val, idx->val, 0, copy->n-1, IFT_INCREASING);

    int median = copy->val[copy->n / 2];

    iftDestroyIntArray(&copy);
    iftDestroyIntArray(&idx);

    return median;
}


double iftBetaFunction(double x, double y) {
    //gamma function deprecated in MacOS, changed to tgamma
    return (tgamma(x)*tgamma(y))/tgamma(x+y);
}

double iftFPDF(double x, double d1, double d2) {

    return sqrt((pow(d1*x, d1)*pow(d2,d2)) / pow((d1*x + d2), d1+d2)) / (x*iftBetaFunction(d1/2, d2/2));
}

double iftPearson4PDF(double x, double m, double v, double a,  double l, double k) {
    double xx = (x-l)/a;
    return k*pow(1.0+iftFastPow(xx, 2), -m)*exp(-v*atan(xx));
}

double iftComplexGamma(double x,double y) {
/* returns abs(gamma(x+iy)/gamma(x))^2 */
    const double y2=y*y, xmin = (2*y2>10.0) ? 2*y2 : 10.0;
    double r=1, s=1, p=1, f=0;
    while(x<xmin) {
        const double t = y/x++;
        r *= 1 + t*t;
    }
    while (p > s*DBL_EPSILON) {
        p *= y2 + f*f;
        p /= x++ * ++f;
        s += p;
    }
    return 1.0/(r*s);
}

double iftPearson4Constant(double m,double v,double a) {
/* returns k */
    //iftAssert(m>0.5 && a>0, "Invalid value.", "iftPearson4Constant");

    double f0 = 0.5*M_2_SQRTPI;
    double f1 = iftComplexGamma(m,0.5*v);
    double f2 = exp(lgamma(m)-lgamma(m-0.5))/a;
    return f0*f1*f2;
}

double iftDiracDeltaFunction(double x) {
    return iftAlmostZero(x)? 1.0 : 0.0;
}

iftRandomLinearCongruentialGenerator iftCreateRandomLinearCongruentialGenerator(){
    iftRandomLinearCongruentialGenerator randomLCG;
    randomLCG.X = 1;
    randomLCG.a  = 1103515245;
    randomLCG.b = 1;
    randomLCG.c = 12345;
    randomLCG.m = 2147483648;
    randomLCG.randMax = randomLCG.m -1;
    return randomLCG;
}

int iftMaxValIntArray(const int *array, int n, const int *mask)
{
    int maxVal = IFT_INFINITY_INT_NEG;
    
    if (mask) {
        for(int i = 0; i < n; i++) {
            if(mask[i] && array[i] > maxVal) {
                maxVal = array[i];
            }
        }
    }
    else {
        for(int i = 0; i < n; i++) {
            if(array[i] > maxVal) {
                maxVal = array[i];
            }
        }
    }
    
    return maxVal;
}


float iftMaxValFloatArray(float *array, int n) {
    float maxVal = IFT_INFINITY_FLT_NEG;

    for (int i = 0; i < n; i++)
        if (array[i] > maxVal)
            maxVal = array[i];

    return maxVal;
}

int iftMinValIntArray(const int *array, int n, const int *mask)
{
    int minVal = IFT_INFINITY_INT;
    
    if (mask) {
        for(int i = 0; i < n; i++) {
            if(mask[i] && array[i] < minVal) {
                minVal = array[i];
            }
        }
    }
    else {
        for(int i = 0; i < n; i++) {
            if(array[i] < minVal) {
                minVal = array[i];
            }
        }
    }
    
    return minVal;
}


float iftMinValFloatArray(float *array, int n) {
    float minVal = IFT_INFINITY_FLT;

    for (int i = 0; i < n; i++)
        if (array[i] < minVal)
            minVal = array[i];

    return minVal;
}
