/**
 * @file iftCommon.h
 * @brief Common definitions and functions to remaining modules.
 */

#ifndef IFT_COMMON_H
#define IFT_COMMON_H

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <cblas.h>
#include <ctype.h>
#include <dirent.h>
#include <float.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <mm_malloc.h>
#include <regex.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#ifdef __linux
#include <dirent.h>
#include <omp.h>
#endif


#if !defined(__APPLE__)
    #include <malloc.h>
#endif
#if defined(__linux)
#endif

#include "ift/core/dtypes/BasicDataTypes.h"

#include "iftMemory.h"
#include "iftDataTransference.h"

/**
 * Common definitions
 */

#ifndef PI
#define PI IFT_PI
#endif

typedef enum ift_axis_order {
  XYZ, XZY, YXZ, YZX, ZXY, ZYX
} iftAxisOrder;

typedef void (*FreeFunction)(void *);

//! swig()
typedef enum {
    IFT_INTERIOR=0, IFT_EXTERIOR=1, IFT_BOTH=2
} iftSide;

#define IFT_RANDOM_SEED (unsigned int) 213344
#define IFT_MAXWEIGHT     4095.0
#define IFT_PI      3.14159265358979323846
#define IFT_WHITE       0
#define IFT_GRAY        1
#define IFT_BLACK       2
#define IFT_NIL        -1
#define IFT_INCREASING  1
#define IFT_DECREASING  0
#define IFT_EPSILON     1E-07

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

//! swig()
typedef enum {
    IFT_AXIS_X, IFT_AXIS_Y, IFT_AXIS_Z
} iftAxis;


/**
 * Common operations
 */

#define IFT_PRINT_BOOL(bool_val) ((bool_val) ? "true" : "false")

/**
 * @brief Converts a bool name (str) into its corresponding bool value
 * @author Cesar Castelo
 * @date Mar 19, 2019
 * @ingroup Common
 */
#define iftBoolStrToBool(bool_str) (!strcmp(bool_str, "true") ? true : false);

#ifndef iftMax
#define iftMax(x,y) (((x) > (y))?(x):(y))
#endif

#ifndef iftMin
#define iftMin(x,y) (((x) < (y))?(x):(y))
#endif

#define iftRound(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))
#define iftRoundHalfDown(x) ((x < 0) ? (int) (x - (0.5 + IFT_EPSILON)) : (int) (x + (0.5 - IFT_EPSILON)))
#define iftSign(x) ((x >= 0) ? 1 : -1)

#define iftPowerOfTwo(x) ((x) * (x))

/**
 * @brief Fixes a float by assigning to 0.0 if it is NaN or infinite.
 * @author Samuel Martins
 * @date Sep 13, 2018
 */
#define iftFixedFloat(x) (isfinite(x) ? (x) : 0)


/**
 * @brief Check if the number x belongs to the interval [a, b]
 * @author Samuka Martins
 * @date Jun 28, 2018
 */
#define iftInClosedInterval(x, a, b) (( (x) >= (a) && (x) <= (b) ))

/** @brief Euclidean Distance between voxels or points */
#define iftPointDistance(u,v) (sqrtf((u.x-v.x)*(u.x-v.x)+(u.y-v.y)*(u.y-v.y)+(u.z-v.z)*(u.z-v.z)))

/** @brief Euclidean Distance between voxels or porints */
#define iftVoxelDistance(u,v) iftPointDistance(u,v)

/** @brief Computes the Squared Distance between two voxels or points */
#define iftSquaredVoxelDistance(u,v) ((u.x-v.x)*(u.x-v.x)+(u.y-v.y)*(u.y-v.y)+(u.z-v.z)*(u.z-v.z))

/* @brief Computes the Smooth Euclidean Distances from a pre-computed Squared Euclidean distance among 26-neighbors
 * @param squared_dist The pre-computed Squared Euclidean Distance among 26-neighbors. The details are in N. Kiryati and G. Sze kely, "Estimating shortest paths and minimal distances on digitized three-dimensional surfaces," Pattern Recognition, vol. 26, pp. 1623-1637, 1993.
 * @return The resulting Smooth Euclidean Distance
 */
#define iftSmoothEuclideanDistance(squared_dist) ((fabs(squared_dist-1)<1.0e-6)? 0.9016 : (fabs(squared_dist-2)<1.0e-6)? 1.289 : 1.615 )

/**
 * @brief Converts a Point to a Voxel
 * @author Samuel
 * @date May 24, 2016
 */
#define iftPointToVoxel(point) ((iftVoxel) {(int) (point).x, (int) (point).y, (int) (point).z})

/**
 * @brief Converts a Voxel to a Point
 * @author Samuel
 * @date May 24, 2016
 */
#define iftVoxelToPoint(v) ((iftPoint) {(v).x, (v).y, (v).z})

#define iftImageCenter(img) ((iftVoxel){((img)->xsize)/2, ((img)->ysize)/2, ((img)->zsize)/2})


/** @brief Constant used to align memory */
#define IFT_MEMORY_ALIGNMENT 16

/**
 * @brief Checks if the long positive number is prime or not.
 * @author Samuel Martins
 * @date Jan 17th, 2016
 * @note Based on https://en.wikipedia.org/wiki/Primality_test
 */
bool iftIsPrime(long n);

/**
 * @brief Prints information about the running platform.
 * @author Peixinho
 * @date Jun, 2017
 * @warning Still not available for all platforms
 */
void iftPlatformInfo();

/**
 * @brief Computes a fast power for small integer exponents.
 * @param base
 * @param b
 * @return Returns base^b
 * @author Peixinho
 * @date May, 2016
 * @warning The funcion is faster upto b=10, after that the regular function pow is called.
 */
double iftFastPow(double base, unsigned b);

/**
 * @brief Prints an array of lentgh <n>.
 * @param v Array of elements
 * @param n Array size
 * @author Cesar Castelo
 * @date Jul 18, 2018
 * @{
 */
void iftPrintBoolArray(bool *v, long n);
void iftPrintCharArray(char *v, long n);
void iftPrintUCharArray(uchar *v, long n);
void iftPrintShortArray(short *v, long n);
void iftPrintUShortArray(ushort *v, long n);
void iftPrintIntArray(int *v, long n);
void iftPrintUIntArray(uint *v, long n);
void iftPrintLongArray(long *v, long n);
#ifndef  __cplusplus
void iftPrintLongLongIntArray(long long *v, long n);
#endif
void iftPrintULongArray(ulong *v, long n);
#ifndef  __cplusplus
void iftPrintULLongArray(ullong *v, long n);
#endif
void iftPrintFloatArray(float *v, long n);
void iftPrintDoubleArray(double *v, long n);
void iftPrintLongDoubleArray(long double *v, long n);
void iftPrintStrArray(char **v, long n);
/** @} */

/**
 * Prints an array of lentgh <n>, following an specific format.
 * @param v Array of elements
 * @param n Array size
 * @param prefix Printed before array
 * @param suffix Printed after array
 * @param delimiter Printed to separate the elements
 * @author Cesar Castelo
 * @date Jul 18, 2018
 * @{
 */
void iftPrintFormattedBoolArray(bool *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedCharArray(char *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedUCharArray(uchar *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedShortArray(short *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedUShortArray(ushort *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedIntArray(int *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedUIntArray(uint *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedLongArray(long *v, long n, const char *prefix, const char *suffix, const char *separator);
#ifndef  __cplusplus
void iftPrintFormattedLongLongIntArray(long long *v, long n, const char *prefix, const char *suffix, const char *separator);
#endif
void iftPrintFormattedULongArray(ulong *v, long n, const char *prefix, const char *suffix, const char *separator);
#ifndef  __cplusplus
void iftPrintFormattedULLongArray(ullong *v, long n, const char *prefix, const char *suffix, const char *separator);
#endif
void iftPrintFormattedFloatArray(float *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedDoubleArray(double *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedLongDoubleArray(long double *v, long n, const char *prefix, const char *suffix, const char *separator);
void iftPrintFormattedStrArray(char **v, long n, const char *prefix, const char *suffix, const char *separator);
/** @} */


/**
 * @brief Returns a random uniform value in the range [low, high]
 * @author Peixinho
 */
float iftRandomUniform(float low, float high);

/**
 * @brief Returns a random integer number between low and high.
 */
int iftRandomInteger (int low, int high);

/**
 * @brief Randomly selects nelems of the set [low, high]
 */
int *iftRandomIntegers (int low, int high, int nelems);

/**
 * @brief Randomly selects a normal distributed (N(0,1)) float number
 */
float iftRandomNormalFloat(float mean, float var);

double iftRandomNormalDouble();

/**
 * @brief Returns initial time
 */
timer *iftTic(void);

/**
 * @brief Returns final time
 */
timer *iftToc(void);

/**
 * @brief Prints the computational time from tstart to tend, as obtained by iftTic() and iftToc() functions.
 */
void iftPrintCompTime(timer *tstart, timer *tend, const char *msg, ...);

/**
 * @brief Computes the difference in ms from the initial time to the final time.
 *
 * @param tic The initial time.
 * @param toc The final time.
 *
 * @return The computed time difference in milliseconds.
 *
 * @warning The memory allocated for <tic> and <toc> is freed inside this function.
 */
float iftCompTime(timer *tic, timer *toc);


/**
 * @brief Returns a string with the formatted time: days hours mins secs ms.
 *
 * @author Samuel
 *
 * Given a time in ms, this function returns a formatted time: days hours mins secs ms.
 * The runtime in MILISECONDS can be obtained with the function iftCompTime(timer *tic, timer *toc).
 *
 * @param runtime Time in ms.
 * @return The formatted time.
 */
//! swig()
char *iftFormattedTime(float runtime);


/**
 * @brief Prints the runtime in the following format: days hours mins secs ms.
 *
 * @author Samuel Martins
 * @date October 23, 2015
 * @sa iftFormattedTime()
 *
 * Given a time in ms, this function prints the runtime following the format: days hours mins secs ms.
 * The runtime in MILISECONDS can be obtained with the function iftCompTime(timer *tic, timer *toc).
 *
 * @param runtime Time in ms.
 */
void iftPrintFormattedTime(float runtime);


/**
 * @brief Writes a timer to a given file, using the following format %s: %f, where %s is the
 * given information corresponding to the time and %f is the current time in milliseconds.
 *
 * @param tic is freed by the function.
 */
void iftWriteTimerToFile(const char *filename, const char *information, timer *tic);

/**
 * @brief Generates seed for rand(), used in iftRandomInteger.
 */
//! swig()
void iftRandomSeed(unsigned int);

/**
 * @brief Returns the factorial of a number or NIL in case of overflow
 */
long double iftFactorial(int n);

/**
 * @brief Returns the limit to avoid overflow in factorial computation
 */
int iftFactorialLimit(void);


void iftUnitNorm(float *feats, int nelems);
void iftNormalizeFeatures(float *feats, int nelems);
void iftStandardizeFeatures(float *feats, int nelems);

float iftSquaredFeatDistance(float *A, float *B, int n);
float iftFeatDistance(float *A, float *B, int n);

/**
 * @brief Computes the Manhattan Distance between two arrays of float considering the sinal, ie.,
 * wihout the module between the feature difference: dist = (b[0]-a[0]) + (b[1]-a[1]) + ...
 */
float iftSignedManhattanDistance(float *A, float *B, int n);

/**
 * @brief Computes the log(val) in the specified base.
 * @author Peixinho
 */
double iftLog(double val, double base);


/**
 * @brief Vector operations defined as macros so that iftVoxels, iftVectors, and iftPoints may be used interchangeably.
 */
/**
 * @brief Subtracts point vec1 from vec2 (it works for iftVectors/iftVoxels/iftPoints).
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 *
 * @warning The result must be cast as an iftVoxel, iftVector, or iftPoint.
 */
#define iftVectorSub(vec1, vec2) {(vec1).x - (vec2).x, (vec1).y - (vec2).y, (vec1).z - (vec2).z, (vec1).t - (vec2).t}
#define iftVoxelSub(u, v) {(u).x - (v).x, (u).y - (v).y, (u).z - (v).z, (u).t - (v).t}

/**
 * @brief Adds points vec1 and vec2 (it works for iftVectors/iftVoxels/iftPoints).
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 *
 * @warning The result must be cast as an iftVoxel, iftVector, or iftPoint.
 */
#define iftVectorSum(vec1, vec2) {(vec1).x + (vec2).x, (vec1).y + (vec2).y, (vec1).z + (vec2).z, (vec1).t + (vec2).t}
#define iftVoxelSum(u, v) {(u).x + (v).x, (u).y + (v).y, (u).z + (v).z, (u).t + (v).t}

/**
 * @brief Computes the cross product between two voxels, points, or vectors.
 *
 * @author Thiago Vallin Spina (changed the macro function's name)
 * @date Mar 04, 2016
 *
 * @warning The result must be cast as an iftVoxel, iftVector, or iftPoint.
 */
#define iftVectorCrossProd(a, b) {(a).y*(b).z - (a).z*(b).y, (a).z*(b).x - (a).x*(b).z, (a).x*(b).y - (a).y*(b).x}
/**
 * @brief Computes the inner product between two voxels, points, or vectors.
 *
 * @author Thiago Vallin Spina (changed the macro function's name)
 * @date Mar 04, 2016
 *
 * @warning The result must be cast as an iftVoxel, iftVector, or
 * iftPoint. It ignores time.
 */
#define iftVectorInnerProduct(a, b) (((a).x*(b).x + (a).y*(b).y + (a).z*(b).z))

/**
 * @brief Computes the Vector Magnitude from a voxel or point
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 *
 */
#define iftVectorMagnitude(v) (sqrtf(iftVectorInnerProduct((v),(v))))


/**
 * @brief Returns the XY angle of a 2D vector
 * @author Peixinho
 * @date April, 2016
 * @ingroup Geometry
 */
#define iftVectorXYAngle(v) ( (atan2(v.y, v.x)+IFT_PI) * (180.0/IFT_PI) )

/**
 * @brief Rounds a vector to integer coordinates.
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 *
 * @warning The result must be cast as an iftVoxel, iftVector, or iftPoint.
 */
#define iftVectorRound(vec1) {iftRound((vec1).x), iftRound((vec1).y), iftRound((vec1).z), iftRound((vec1).t)}
/**
 * @brief Multiplies a vector by an scalar.
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 *
 * @warning The result must be cast as an iftVoxel, iftVector, or iftPoint.
 */
#define iftVectorScalarProd(vec, s) {(vec).x*(s), (vec).y*(s), (vec).z*(s), (vec).t*(s)}
/**
 * @brief Adds an scalar to a vector.
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 *
 * @warning The result must be cast as an iftVoxel, iftVector, or iftPoint.
 */
#define iftVectorScalarSum(vec, s) {(vec).x+(s), (vec).y+(s), (vec).z+(s), (vec).t+(s)}
/**
 * @brief Divides a vector by an scalar.
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 *
 * @return The vector divided by the scalar, or the vector itself if the scalar is close to 0.0
 *
 * @warning The result must be cast as an iftVoxel, iftVector, or iftPoint.
 */
iftVector iftVectorScalarDiv(iftVector vec, double s);
/**
 * @brief Verifies if a vector has coordinates (0, 0, 0).
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 */
#define iftIsNullVector(vec) (iftAlmostZero((vec).x) && iftAlmostZero((vec).y) && iftAlmostZero((vec).z) && iftAlmostZero((vec).t))
/**
 * @brief Verifies if three points, voxels, or vectors are collinear.
 *
 * This is achieved by testing if the cross product between the vectors P1->P2 and P1->P3 is Null (0,0,0).
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 */
#define iftCollinearPoints(P1, P2, P3) (iftIsNullVector((iftVector)iftVectorCrossProd((iftVector)iftVectorSub((P2), (P1)), (iftVector)iftVectorSub((P3), (P1)))))

/**
 * @brief Tests if two iftVoxels are equal
 *
 * @author Thiago Vallin Spina
 * @date Apr 19, 2016
 *
 * @param u1 The first voxel
 * @param u2 The second voxel
 * @return True if their coordinates are the same
 */
bool iftVoxelsAreEqual(iftVoxel u1, iftVoxel u2);

/**
 * @brief Computes the mean voxel
 *
 * @author Thiago Vallin Spina
 * @date May 1, 2016
 *
 * @param u1 The first voxel
 * @param u2 The second voxel
 * @return The mean voxel (u1 + u2) / 2
 */
iftVoxel iftMeanVoxel(iftVoxel u1, iftVoxel u2);

/**
 * @brief Tests if two iftPoints/iftVectors are equal
 *
 * @author Thiago Vallin Spina
 * @date Apr 19, 2016
 *
 * @param u1 The first point/vector
 * @param u2 The second point/vector
 * @return True if their coordinates are (almost) the same
 */
#define iftVectorsAreEqual(u1, u2) (iftAlmostZero((u1).x-(u2).x) && iftAlmostZero((u1).y-(u2).y) && \
                                      iftAlmostZero((u1).z-(u2).z))

/**
 * @breif Check if the vector u is a zero-vector
 * @author Samuel Martins
 * @date Apr 22, 2018
 */
#define iftIsZeroVector(u) ((iftAlmostZero(u.x) && iftAlmostZero(u.y) && iftAlmostZero(u.z)))


/**
 * @brief Normalizes a vector if its norm is greater than 0.0, otherwise returns the vector itself.
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 */
iftVector iftNormalizeVector(iftVector v);

/**
 * @brief Projects vector U onto vector V.
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 *
 * @warning If vector V has norm close to 0.0 the function issues and iftError.
 */
iftVector iftProjectVector(iftVector U, iftVector V);


/**
 * @brief Orthogonalize the vector u onto vector v based on Gram–Schmidt process.
 *
 * https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
 *
 * @warning If any vector is a zero-vector or if they are parallel, a zero-vector is returned.
 *
 * @author Samuel Martins
 * @date Apr 22, 2018
 * @ingroup Geometry
 */
iftVector iftOrthogonalizeVector(iftVector u, iftVector v);


/**
 * @brief Returns the distance between P0 and the line from P1 to P2, whose size is P1P2
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 */
double iftVoxelLineDist2D(iftVoxel P0, iftVoxel P1, iftVoxel P2, double P1P2);

/**
 * @brief Given the vector from P1 to P0, returns the voxel projected on the line between P1 and P2.
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 */
iftVoxel iftClosestVoxelOnLine2D(iftVoxel P0, iftVoxel P1, iftVoxel P2);

/**
 * @brief Returns the position of P0 with respect to the line from P1 to P2.
 *
 * Negative values indicate left side, 0 indicates on the line,
 * and positive values indicate right side.
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 * @ingroup Geometry
 */
int iftVoxelLinePosition2D(iftVoxel P0, iftVoxel P1, iftVoxel P2);


/**
 * @brief Computes the circle's area given a radius.
 *
 * @author Thiago Vallin Spina
 * @date Mar 04, 2016
 *
 * @ingroup Geometry
 */
static inline double iftCircleArea(double r) {
    return (r*r)*IFT_PI;
}

/* @brief Gets rid of the carriage return and the line feed characteres introduced by DOS systems when reading
 * strings from ASCII files.
 *
 * @author Alexandre Xavier Falcao
 *
 * @param line The string in which all carriage return characters are removed.
 */
void iftRemoveCarriageReturn(char *line);

/**
 * Updates the mean value after removing one element from the set
 * @param mean average value of the original set
 * @param n size of the original set
 * @param x removed value
 * @return the updadated average value
 * @author Peixinho
 * @date Mar, 2017
 */
float iftDecrementalMean(float mean, int n, float x);

/**
 * Print memory amount with the specific quantifier (KB,MB,GB)
 * @author Peixinho
 * @param mem
 */
void iftPrintMemoryFormatted(long mem);

/**
 * Updates the variance value after removing one element from the set
 * @param mean average value of the original set
 * @param var variance value of the original set
 * @param n size of the original set
 * @param x removed value
 * @return
 */
float iftDecrementalVariance(float mean, float var, int n, float x);

/**
 * @brief Evaluates the sigmoid function, with x = value.
 *
 * @param alfa Controls the decay of the function.
 * @author Renzo Phellan
 */
float iftSigmoidalValue(float value, float alfa);

void iftVerifyToken(FILE *fp, char *token, char *function);
void iftReadIntValue(FILE *fp, int *value, char *token, char *function);
void iftReadIntValues(FILE *fp, int **value, int nvalues, char *token, char *function);
void iftWriteIntValue(FILE *fp, int value, char *token);
void iftWriteIntValues(FILE *fp, int *value, int nvalues, char *token);
void iftReadFloatValue(FILE *fp, float *value, char *token, char *function);
void iftReadFloatValues(FILE *fp, float **value, int nvalues, char *token, char *function);
void iftWriteFloatValue(FILE *fp, float value, char *token);
void iftWriteFloatValues(FILE *fp, float *value, int nvalues, char *token);
void iftReadDoubleValue(FILE *fp, double *value, char *token, char *function);
void iftReadDoubleValues(FILE *fp, double **value, int nvalues, char *token, char *function);
void iftWriteDoubleValue(FILE *fp, double value, char *token);
void iftWriteDoubleValues(FILE *fp, double *value, int nvalues, char *token);
void iftSkipComments(FILE *fp);


/* These functions are currently used to communicate with numpy */
void iftWriteRawIntArray(char *filename, int *array, int n);
int* iftReadRawIntArray(char *filename, int n);


/**
 * @brief Computes the mean value of an integer array.
 * @author Samuel Martins
 * @date Nov 21, 2018
 */
float iftIntMean(const int *x, int n);


/**
 * @brief Computes the mean value of a float array.
 */
float iftMean(float *x, int n);
/**
 * @brief Computes the variance of a float array.
 */
float iftVar(float *x, int n);
/**
 * @brief Computes the standard deviation of a float array.
 */
float iftStd(float *x, int n);

/**
 * @brief Computes the covariance of float arrays.
 */
float iftCov(float *x, float *y, int n);

/**
 * @brief Tests if a real (float/double) number is zero. It should be used for comparisons if two floats are equal.
 *
 * @author Thiago Vallin Spina
 */
int iftAlmostZero(double x);
/**
 * @brief Computes the mod operation of <a> by <n>. If the mod value is negative, it is placed within [0,n-1].
 */
int iftSafeMod(int a, int n);


/**
 * @brief Gets an Integer Array with the Unique Elements from an Input Array. The resulting Unique array is SORTED in ASCENDING order;
 *
 * Ex: Given the array = [4, 2, 3, 1, 4, 3, 2, 1], with size n = 8, the resulting array of unique elements will be:
 * [1, 2, 3, 4]
 *
 * @param  array Input Int Array to be scanned.
 * @param  n     Number of Elements from the Input Int Array.
 * @return       A Sorted Integer Array in Ascending Order with the Unique Elements from <b>array</b>.
 */
iftIntArray *iftIntArrayUnique(const int *array, int n);



/**
 * @brief Counts the Number of Unique Elements from an array.
 *
 * @param array Integer Array whose Number of Unique Elements will be counted.
 * @param n     Array Size.
 * @return      The number of Unique Elements.
 */
int iftCountUniqueIntElems(const int *array, int n);

/**
 * @brief Counts the number of elements per each unique element in an array.
 * It returns the number of elements and uses out pointers to return the elements and the number of repetitions for each element
 * @param array Input array
 * @param n Input array Size.
 * @param elems_out Array with the elements that are present in the result
 * @param quant_out NUmber of repetitions of each element
 * @return Number of elements that are present in the result
 * @author Cesar Castelo
 * @date Abr 24, 2018
 */
int iftCountElemsPerUniqueIntElem(int *array, int n, int **elems_out, int **quant_out);

/**
 * @brief Checks if the Integer Array contains the target value.
 *
 * @param  array Interget Array.
 * @param  n     Array Size
 * @param  value Target value
 * @return       true/false
 *
 * @author Samuka
 * @date Dec 19, 2016
 */
bool iftIntArrayContainValue(const int *array, int n, int value);


/**
 * @brief This function finds the upper bound value (normalize value) of the range of a positive value.
 *
 * It finds the the most suitable normalization value to be used during image color/feature conversion/normalization from a given positive value.
 *
 * This number is 2ˆn_bits - 1, where the supported number of bits <n_bits> per voxel (image depth) is in {1, 8, 12, 16, 32}.
 * Therefore, the possible upper bound ranges are respectively: 1, 255, 4095, 65535, 4294967295.
 *
 * @warning The input value must be positive.
 *
 * @param maxval The maximum value of an image, for instance. It must be positive.
 * @return The normalization value.
 *
 * @author Alexandre Xavier Falcao, Samuel Martins
 * @date Dec 5, 2018
 */
//! swig()
long iftNormalizationValue(long maxval);


/**
 * @brief Converts a number in radians to degrees.
 *
 * @author Thiago Vallin Spina
 */
static inline double iftDegrees(double rad) {
  return rad * 180.0 / IFT_PI;
}

/**
 * @brief Converts a number in degrees to radians.
 *
 * @author Thiago Vallin Spina
 */
static inline double iftRadians(double deg) {
  return deg * IFT_PI / 180.0;
}

/**
 * @brief Converts a negative radian number to positive.
 *
 * @author Thiago Vallin Spina
 */
static inline double iftNegativeToPositiveRadians(double rad) {
  return (rad < 0) ? 2 * IFT_PI + rad : rad;
}

/**
 * @brief Converts a negative degree number to positive.
 *
 * @author Thiago Vallin Spina
 */
static inline double iftNegativeToPositiveDegrees(double deg) {
  return (deg < 0) ? 360.0+deg : deg;
}

/**
 * @brief Compute the index of the element in vector x with minimum value
 *
 * @param x Float vector.
 * @param n Size of the vectors.
 * @return index.
*/
int iftArgMin(const int *x, long n);

/**
 * @brief Compute the index of the element in vector x with maximum value
 *
 * @param x Float vector.
 * @param n Size of the vectors.
 * @return index.
*/
int iftArgMax(const int *x, long n);

/**
 * @brief Compute the index of the element in vector x with minimum value
 *
 * @param x Float vector.
 * @param n Size of the vectors.
 * @return index.
*/
int iftFArgmin(const float *x, int n);

/**
 * @brief Compute the index of the element in vector x with maximum value
 *
 * @param x Float vector.
 * @param n Size of the vectors.
 * @return index.
*/
int iftFArgmax(const float *x, int n);

/**
 * @brief Compute the index of the element in vector x with minimum value
 *
 * @param x Double vector.
 * @param n Size of the vectors.
 * @return index.
*/
int iftDArgmin(const double *x, int n);

/**
 * @brief Compute the index of the element in vector x with maximum value
 *
 * @param x Double vector.
 * @param n Size of the vectors.
 * @return index.
*/
int iftDArgmax(const double *x, int n);


/**
 * @brief Finds the index of <elem> in array <x> and returns it.
 *
 * @author Peixinho
 *
 * @param x Unsorted array of elements.
 * @param n The size of the array.
 * @param elem The elem to be found.
 *
 * @return The element's index in <x> or NIL if not found
 */
int iftFindIntArrayElementIndex(int *x, int n, int elem);


/**
 * @brief Compute the Gaussian Probability Density Function for a set of values.
 * @author Samuel Martins
 * @date May 12, 2016
 *
 * @param  vals   Array of Values to compute the Gaussian PDF.
 * @param  n_vals Number of values.
 * @param  mean   Mean used to compute the Gaussian PDF.
 * @param  stdev  Standard Deviation used to compute the Gaussian PDF.
 * @return        Array of probabilities, one for each input value.
 */
double *iftGaussPDF(double *vals, size_t n_vals, double mean, double stdev);


/**
 * @brief Compute the Cumulative Density of the (Gaussian) PDF (Normal Distribution) in a Random Variable with value x.
 * @author Samuel Martins
 * @date May 12, 2016
 *
 * It gives the Area under the Probability Density Function (Normal Distribution) from minus infinity to x. \n
 * Then, it also describes the probability that a random variable X takes on a value less than or
 * equal to a number x. \n
 *
 * @param  x      Value from the Random Variable to compute the CDF.
 * @param  mean   Mean used to compute the Gaussian PDF.
 * @param  stdev  Standard Deviation used to compute the Gaussian PDF.
 * @return        The area under the PDF from minus to x.
 */
double iftGaussCDF(double x, double mean, double stdev);


/**
 * @brief Finds the Median Value from an integer array.
 * @author Samuka Martins
 * @date Dec 5, 2017
 */
int iftMedianIntVal(const iftIntArray *arr);


/**
 * @brief Computes the beta function B(x,y)
 * @param x
 * @param y
 * @author Peixinho
 * @date Apr, 2017
 */
double iftBetaFunction(double x, double y);

/**
 * Fisher–Snedecor F-distribution.
 * @param x
 * @param d1
 * @param d2
 * @return the probability of x
 * @author Peixinho
 * @date Apr, 2017
 */
double iftFPDF(double x, double d1, double d2);

/**
 * Pearson type IV distribution.
 * @param x
 * @param m kurtosis shape parameter
 * @param v skew shape parameter
 * @param a deviation parameter
 * @param l location parameter
 * @return the probability of x
 * @author Peixinho
 * @date Apr, 2017
 */
double iftPearson4PDF(double x, double m,double v,double a, double l, double k);


/**
 *
 * @param m
 * @param v
 * @param a
 * @return
 */
double iftPearson4Constant(double m,double v,double a);

/**
 * @brief Computes the dirac delta function of x.
 * @author Thiago Vallin Spina
 * @date Oct 27, 2016
 *
 * @param x a real-valued number.
 *
 * @return if x is (almost) zero, returns 1.0, otherwise, returns 0.0.
 */
double iftDiracDeltaFunction(double x);

/**
 * @brief Finds the maximum value in an integer array
 * @author Cesar Castelo
 * @date Mar 18, 2019
 * @param array Unsorted array of elements.
 * @param n The size of the array.
 *
 * @return The maximum value
 */
int iftMaxValIntArray(const int *array, int n, const int *mask);


/**
 * @brief Finds the maximum value in a float array
 * @author Cesar Castelo
 * @date Mar 18, 2019
 * @param array Unsorted array of elements.
 * @param n The size of the array.
 *
 * @return The maximum value
 */
float iftMaxValFloatArray(float *array, int n);

/**
 * @brief Finds the minimum value in an integer array
 * @author Cesar Castelo
 * @date Mar 18, 2019
 * @param array Unsorted array of elements.
 * @param n The size of the array.
 *
 * @return The minimum value
 */
int iftMinValIntArray(const int *array, int n, const int *mask);


/**
 * @brief Finds the minimum value in a float array
 * @author Cesar Castelo
 * @date Mar 18, 2019
 * @param array Unsorted array of elements.
 * @param n The size of the array.
 *
 * @return The minimum value
 */
float iftMinValFloatArray(float *array, int n);


//ref: https://en.wikipedia.org/wiki/Linear_congruential_generator
typedef struct _ift_randomLinearCongruentialGenerator {
    unsigned long int X;
    unsigned long int a;
    unsigned long int b;//parameter used in rand.c
    unsigned long int c;
    unsigned long int m;
    unsigned long int randMax;
}iftRandomLinearCongruentialGenerator;

iftRandomLinearCongruentialGenerator iftCreateRandomLinearCongruentialGenerator();
inline unsigned long iftRandomLCG(iftRandomLinearCongruentialGenerator* randomLCG) // RAND_MAX assumed to be 32767
{
    randomLCG->X = (randomLCG->X * randomLCG->a + randomLCG->c) % randomLCG->m;
    return randomLCG->X;
}

#define iftUnusedParameter(paramater) (void)(paramater)

#ifdef __cplusplus
}
#endif

#endif
