#ifndef _IFTMATRIX_H_
#define _IFTMATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ift/core/dtypes/BasicDataTypes.h"
#include "ift/core/dtypes/DblArray.h"
#include "ift/core/dtypes/Matrix.h"
#include "ift/core/dtypes/FloatArray.h"
#include "iftCommon.h"


// Its type definition in iftBasicDataType.h
// iftMatrix



iftMatrix* iftCreateMatrixPointer(float *val, int ncol, int nrow);
/**
 * @brief Reshape the matrix
 * @author Deangeli
 * @date may, 2017
 * @param matrix
 * @param newCol
 * @param newRow
 */
void iftReshapeMatrix(iftMatrix* matrix, int newCol, int newRow);
iftMatrix* iftPolynomialExpasion(iftMatrix* matrix, unsigned int degree);
void iftDestroyMatrixPointer(iftMatrix**);


/**
 * @brief Matrix of double values
 */
//! swig(extend = iftDoubleMatrixExt.i, destroyer = iftDestroyDoubleMatrix)
typedef struct ift_double_matrix {
    double *val;                  /* matrix's values */
    int ncols, nrows; // number of columns, number of rows,
    long  *tbrow; //look-up table to speed up index access
    long n;                        /* ncols * nrows */
    bool allocated;
} iftDoubleMatrix;


/**
 * @brief Matrix of strings of 1024 characters.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 * 
 * @note Its type definition in iftBasicDataType.h
 */
// iftStrMatrix
struct ift_str_matrix {
    /** Matrix of elements. */
    char **val;
    /** Number of Columns, Rows, and look-up table to speed up index access. */
    int ncols, nrows, *tbrow;
    /** Total number of elements = ncols*nrows. */
    int n;
};

typedef struct ift_bool_matrix {
    bool *val;          // matrix's values
    int ncols, nrows;  // number of columns, number of rows,
    int  *tbrow;       //look-up table to speed up index access
    int n;             // ncols * nrows
} iftBoolMatrix;


#define iftGetMatrixCol(m,i) ((i) % (m)->ncols)
#define iftGetMatrixRow(m,i) ((i) / (m)->ncols)
#define iftGetMatrixIndex(m,c,r) ((c)+(m)->tbrow[(r)])
#define iftMatrixRowPointer(m, r) ((m)->val + iftGetMatrixIndex((m), 0, (r)))


/**
 * Macro to access the element m[c,r] in matrix m.
 * @author Peixinho
 * @date Jun, 2016
 */
#define iftMatrixElem(m, c, r) (m)->val[iftGetMatrixIndex((m), (c), (r))]

iftDoubleMatrix* iftCopyDoubleMatrix(const iftDoubleMatrix* m);

//! swig(newobject)
iftDoubleMatrix *iftCreateDoubleMatrix(int ncols, int nrows);
iftDoubleMatrix *iftCreateDoubleMatrixNoAlloc(int ncols, int nrows);

void iftReallocDoubleMatrix(iftDoubleMatrix *M, int ncols, int nrows);

iftMatrix *iftCopyMatrix(const iftMatrix *A);
void       iftDestroyDoubleMatrix(iftDoubleMatrix **M);

//! swig(newobject)
iftMatrix *iftReadRawMatrix(const char *filename);

//! swig()
void        iftWriteRawMatrix(iftMatrix *M, char *filename);

void        iftPrintDoubleMatrix(iftDoubleMatrix *M);

iftBoolMatrix *iftReadBoolMatrix(const char *filename);
void        iftWriteBoolMatrix(iftBoolMatrix *M, char *filename);

/**
 * @brief Solves the least squares problem. Finding the best x, that solves Ax=b.
 */
iftMatrix *iftLeastSquares(const iftMatrix* A, const iftMatrix* B);
iftMatrix *iftMultMatricesChain(int n, ...);
iftMatrix *iftMultMatrices(const iftMatrix *A, const iftMatrix *B);
iftMatrix *iftComputeTransformation(int n_args, ...);
void       iftMultMatrixByScalar(iftMatrix *A, double scalar );
iftMatrix *iftInvertMatrix(const iftMatrix *A);
iftMatrix *iftInvertDiagonalMatrix(const iftMatrix *A);
void iftLinearCombination(float a, float* A, float b, float* B, float *C, int n);

/**
 * @brief Computes the inversion of Matrix A in place, without check for invertibility.
 */
iftMatrix * iftPseudoInvertMatrix(const iftMatrix *A);
iftDoubleMatrix *iftPseudoInvertDMatrix(const iftDoubleMatrix *A);
iftMatrix *iftTransposeMatrix(const iftMatrix *A);
double     iftMatrixDeterminant(const iftMatrix *A);
iftMatrix *iftCoFactorMatrix(const iftMatrix *A);
iftMatrix *iftRotationMatrix(char axis, float theta);
iftMatrix *iftTranslationMatrix(iftVector T);
iftMatrix *iftScaleMatrix(float sx, float sy, float sz);
iftVector  iftTransformVector(iftMatrix *M, iftVector u);
iftPoint   iftTransformPoint(iftMatrix *M, iftPoint  u);


/**
 * @brief Transform a voxel from a 4x4 transformation matrix on homogeneous coordinate.
 * 
 * @param  M 4X4 Transformation Matrix on homogeneous coordinates.
 * @param  v Voxel to be transformed 
 * @return   Transformed voxel.
 *
 * @author Samuel Martins
 * @date Apr 22, 2018
 */
iftVoxel iftTransformVoxel(const iftMatrix *M, iftVoxel v);


iftMatrix *iftIdentityMatrix(int ncols);
iftMatrix *iftRotationMatrixToAlignVectorWithZ(iftVector v);
iftMatrix *iftRotationMatrixToAlignZWithVector(iftVector v);
iftMatrix *iftExtendMatrix(const iftMatrix *M, int ncols, int nrows);
void       iftSingleValueDecomp(iftMatrix *A, iftMatrix **U, iftMatrix **S, iftMatrix **Vt);
void       iftNormalizeMatrix(iftMatrix *M);
double     iftFrobeniusNorm(iftMatrix *M);
iftMatrix  *iftReadRawMatrixFILE(FILE *fp);
void       iftWriteRawMatrixFILE(iftMatrix *M, FILE *fp);


/**
 * @brief Computes the generalized Eigenvalues
 * @param A: Input Matrix
 * @author Azael Sousa
 * @date May 06, 2021
 */
iftMatrix *iftGeneralizedEigenvalues(iftMatrix *A);


/**
 * @brief Normalize the matrix values into the range specified by minVal and maxVal
 * @param M Matrix
 * @param minValue Minimum value
 * @param maxValue Maximum value
 * @author Cesar Castelo
 * @date Set 04, 2018
 */
void iftNormalizeMatrixInRange(iftMatrix *M, float minVal, float maxVal);

/**
 * @brief Reallocates memory for a matrix and copies the original data. The new size could be higher or lower
 * @warning If the original matrix is larger then some data will be lost, i.e. only ncols columns and nrows rows will be copied
 * @param M Matrix pointer
 * @param ncols New number of columns
 * @param nrows New number of rows
 * @author Cesar Castelo
 * @date Jul 20, 2018
 */
void iftResizeMatrix(iftMatrix **M, int ncols, int nrows);
void iftResizeDoubleMatrix(iftDoubleMatrix **M, int ncols, int nrows);
void iftResizeStrMatrix(iftStrMatrix **M, int ncols, int nrows);

/**
 * Saves a iftMatrix as CSV
 * @param M
 * @param filename
 * @author Peixinho
 * @date Jun, 2017
 */
void iftWriteMatrixCSV(iftMatrix* M, const char* filename);

/**
 * Reads a iftMatrix as CSV
 * @param filename
 * @author Cesar Castelo
 * @date Mar 12, 2019
 */
iftMatrix *iftReadMatrixCSV(const char* filename);






/**
 * @brief Rearrange a matrix to a string for printing.
 * @author Samuka
 * @date Sep 16, 2017
 */
char *iftMatrixAsString(const iftMatrix *mat);

/**
 * @brief Rearrange a matrix to a string for printing.
 * @author Samuka
 * @date Sep 16, 2017
 */
char *iftStrMatrixAsString(const iftStrMatrix *mat);


/**
 * @brief Creates a Matrix of Strings with 1024 Characters.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 *
 * @param ncols Number of Columns.
 * @param nrows Number of Rows.
 * @return The created String Matrix.
 */
iftStrMatrix *iftCreateStrMatrix(int ncols, int nrows);


/**
 * @brief Destroys a Matrix of Strings.
 *
 * @author Samuel Martins
 * @date October 20, 2015
 */
void iftDestroyStrMatrix(iftStrMatrix **sM);


/**
 * @brief Copies a string matrix
 * @author Samuka Martins
 * @date Sep 17, 2017
 */
iftStrMatrix *iftCopyStrMatrix(char **val, int nrows, int ncols);

/**
 * @brief Creates a Matrix of Bool values.
 *
 * @author Cesar Castelo
 * @date Set 04, 2018
 *
 * @param ncols Number of Columns.
 * @param nrows Number of Rows.
 * @return The created Bool Matrix.
 */
iftBoolMatrix *iftCreateBoolMatrix(int ncols, int nrows);


/**
 * @brief Destroys a Matrix of Bool values.
 *
 * @author Cesar Castelo
 * @date Set 04, 2018
 */
void iftDestroyBoolMatrix(iftBoolMatrix **iM);


/**
 * @brief Copies a Bool Matrix
 * @author Cesar Castelo
 * @date Set 04, 2018
 */
iftBoolMatrix *iftCopyBoolMatrix(bool *val, int nrows, int ncols);

/**
 * @brief Transforms an array of voxels given.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015

 * @param S      Transformation matrix.
 * @param points Voxel array.
 * @param n      Number of voxels in the array.
 * @return The transformed voxels.
 */
iftVoxel *iftTransformVoxels(iftMatrix *S, iftVoxel *points, int n);

/**
 * @brief Gets the minimum value among all elements in the matrix
 *
 * @author Deangeli Gomes Neves
 * @date july 4, 2016
 *
 * @param matrix input matrix.
 * @return minimum value in matrix.
 */
double iftGetMinimumValueInMatrix(iftMatrix *matrix);

/**
 * @brief Gets the maximum value among all elements in the matrix
 *
 * @author Deangeli Gomes Neves
 * @date july 4, 2016
 *
 * @param matrix input matrix.
 * @return maximumvalue in matrix.
 */
double iftGetMaximumValueInMatrix(iftMatrix *matrix);


/**
 * @brief Computes the mean value for each matrix column
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 4, 2016
 *
 * @param matrix input matrix.
 * @return vector that contains the mean value of each column
 */
iftMatrix* iftMatrixMeanColumn(const iftMatrix *matrix);

/**
 * @brief Computes the sum for each matrix column
 * @author Peixinho
 * @date Jul, 2016
 * @param matrix input matrix.
 * @return Row matrix that contains the sum of each column
 */
iftMatrix* iftMatrixSumColumn(const iftMatrix *matrix);


/**
 * @brief Computes the sum for each matrix row
 * @author Peixinho
 * @date Jul, 2016
 * @param matrix input matrix.
 * @return Column matrix that contains the sum of each row
 */
iftMatrix* iftMatrixSumRow(const iftMatrix *matrix);

/**
 * @brief Obtains the max value for each matrix column
 * @author Cesar Castelo
 * @date Nov 28, 2018
 * @param matrix input matrix.
 * @return Row matrix that contains the max value of each column
 */
iftMatrix* iftMatrixMaxColumn(const iftMatrix *matrix);


/**
 * @brief Obtains the max value for each matrix row
 * @author Cesar Castelo
 * @date Nov 28, 2018
 * @param matrix input matrix.
 * @return Column matrix that contains the max value of each row
 */
iftMatrix* iftMatrixMaxRow(const iftMatrix *matrix);

/**
 * @brief Obtains the max value for one matrix column
 * @author Cesar Castelo
 * @date Nov 28, 2018
 * @param matrix input matrix.
 * @param c column.
 * @return Max value of the chosen column
 */
float iftMatrixMaxOneColumn(const iftMatrix *matrix, int c);


/**
 * @brief Obtains the max value for one matrix row
 * @author Cesar Castelo
 * @date Nov 28, 2018
 * @param matrix input matrix.
 * @param r row
 * @return Max calue of the chosen row
 */
float iftMatrixMaxOneRow(const iftMatrix *matrix, int r);

/**
 * @brief Performs a mathematical operation between a matrix and a scalar
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 5, 2016
 *
 *
 *
 * @param matrix input matrix.
 * @param scalar scalar value
 * @param order 'f' to matrix[i] (operator) by scalar, 'b' to scalar (operator) by operator
 * @param operationSymbol '+' to addition, '-' to subtraction, '*' to multiplication, '/' to division.
 *
 */
void iftComputeOperationBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order, char operationSymbol);

/**
 * @brief Performs a addition operation between a scalar and  all elements of matrix.
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 *
 *
 * @param matrix input matrix.
 * @param scalar scalar value
 * @param order 'f' to matrix[i][j] = matrix[i][j] +  scalar, 'b' to matrix[i][j] = scalar + matrix[i][j]
 *
 */
void iftComputeAdditionBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order);

/**
 * @brief Performs a subtraction operation between a scalar and  all elements of matrix.
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * @param matrix input matrix.
 * @param scalar scalar value
 * @param order 'f' to matrix[i][j] = matrix[i][j] - scalar, 'b' to matrix[i][j] = scalar - matrix[i][j]
 *
 */
void iftComputeSubtractionBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order);

/**
 * @brief Performs a multiplication operation between a scalar and  all elements of matrix.
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * @param matrix input matrix.
 * @param scalar scalar value
 * @param order 'f' to matrix[i][j] = matrix[i][j] * scalar, 'b' to matrix[i][j] = scalar * matrix[i][j]
 *
 */
void iftComputeMultiplicationBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order);

/**
 * @brief Performs a division operation between a scalar and  all elements of matrix.
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * @param matrix input matrix.
 * @param scalar scalar value
 * @param order 'f' to matrix[i][j] = matrix[i][j] / scalar, 'b' to matrix[i][j] = scalar / matrix[i][j]
 *
 */
void iftComputeDivisionBetweenMatrixScalar(iftMatrix *matrix, double scalar, char order);



/**
 * @brief Performs a mathematical operation between a matrix and a vector
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 5, 2016
 *
 * if the matrix dimensons are MxN, and the input vector is a column vector, its dimensions must be
 * Mx1, otherwise 1xN.
 *
 *
 * @param matrix input matrix.
 * @param Vector row or column vector
 * @param vectorType 'r' to row vector, 'c' to column vector
 * @param operationSymbol '+' to addition, '-' to subtraction, '*' to multiplication, '/' to division.
 *
 */
void iftComputeOperationBetweenMatrixVector(iftMatrix *matrix, iftDblArray* vector, char vectorType,
                                            char operationSymbol);

void iftAdditionMatrixVectorByColumn(iftMatrix *matrix, iftMatrix* vector);
void iftAdditionMatrixVectorByRow(iftMatrix *matrix, iftMatrix* vector);
void iftSubtractionMatrixVectorByColumn(iftMatrix *matrix, iftMatrix* vector);
void iftSubtractionMatrixVectorByRow(iftMatrix *matrix, iftMatrix* vector);
void iftMultiplicationMatrixVectorByColumn(iftMatrix *matrix, iftMatrix* vector);
void iftMultiplicationMatrixVectorByRow(iftMatrix *matrix, iftMatrix* vector);
void iftDivisionMatrixVectorByColumn(iftMatrix *matrix, iftMatrix* vector);
void iftDivisionMatrixVectorByRow(iftMatrix *matrix, iftMatrix* vector);






/**
 * @brief Performs a mathematical pointwise operation between two matrices
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 5, 2016
 *
 * The matrices must have the same dimensions. The operation is A "operator" by B
 *
 *
 * @param A input matrix A.
 * @param B input matrix B
 * @param operationSymbol '+' to addition, '-' to subtraction, '*' to multiplication, '/' to division.
 *
 * @return a matrix C that is result of A (pointwise operator) by B
 */
iftMatrix* iftComputeOperationBetweenMatricesInPointWise(iftMatrix *A, iftMatrix *B, char operationSymbol);

/**
 * @brief Performs pointwise addition between two matrices
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * The matrices must have the same dimensions. The operation is A + B
 *
 *
 * @param A input matrix A.
 * @param B input matrix B
 *
 * @return a matrix C that is result of A + B
 */
iftMatrix* iftMatricesAdditionPointWise(iftMatrix *A, iftMatrix *B);

/**
 * @brief Performs pointwise subtraction between two matrices
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * The matrices must have the same dimensions. The operation is A - B
 *
 *
 * @param A input matrix A.
 * @param B input matrix B
 *
 * @return a matrix C that is result of A - B
 */
iftMatrix* iftMatricesSubtractionPointWise(iftMatrix *A, iftMatrix *B);

/**
 * @brief Performs pointwise multiplication between two matrices
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * The matrices must have the same dimensions. The operation is A * B
 *
 *
 * @param A input matrix A.
 * @param B input matrix B
 *
 * @return a matrix C that is result of A * B
 */
iftMatrix* iftMatricesMultiplicationPointWise(iftMatrix *A, iftMatrix *B);

/**
 * @brief Performs pointwise division between two matrices
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * The matrices must have the same dimensions. The operation is "A divided by B"
 *
 *
 * @param A input matrix A.
 * @param B input matrix B
 *
 * @return a matrix C that is result of "A divided by B"
 */
iftMatrix* iftMatricesDivisionPointWise(iftMatrix *A, iftMatrix *B);


/**
 * @brief Performs pointwise addition between the matrices <b>A</b> and <b>B</b>, and stores the result in <b>C</b>
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * The matrices <b>A</b> and <b>B</b> must have the same dimensions. If the matrix <b>C</b> is NULL, the function
 * allocates memory to C. If C has not the same dimensions of A/B, the functions destroy the C pointer, and then
 * creates a new <b>iftMatrix</b> pointer that matches the dimensions requirements. The operation is C = A + B
 *
 *
 * @param A input matrix A.
 * @param B input matrix B
 * @param C pointer to pointer where the result will be stored
 */
void iftMatricesAdditionPointWiseInPlace(iftMatrix *A, iftMatrix *B, iftMatrix **C);

/**
 * @brief Performs pointwise subtraction between the matrices <b>A</b> and <b>B</b>, and stores the result in <b>C</b>
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * The matrices <b>A</b> and <b>B</b> must have the same dimensions. If the matrix <b>C</b> is NULL, the function
 * allocates memory to C. If C has not the same dimensions of A/B, the functions destroy the C pointer, and then
 * creates a new <b>iftMatrix</b> pointer that matches the dimensions requirements. The operation is C = A - B
 *
 *
 * @param A input matrix A.
 * @param B input matrix B
 * @param C pointer to pointer where the result will be stored
 */
void iftMatricesSubtractionPointWiseInPlace(iftMatrix *A, iftMatrix *B, iftMatrix **C);

/**
 * @brief Performs pointwise multiplication between the matrices <b>A</b> and <b>B</b>, and stores the result in <b>C</b>
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * The matrices <b>A</b> and <b>B</b> must have the same dimensions. If the matrix <b>C</b> is NULL, the function
 * allocates memory to C. If C has not the same dimensions of A/B, the functions destroy the C pointer, and then
 * creates a new <b>iftMatrix</b> pointer that matches the dimensions requirements. The operation is C = A * B
 *
 *
 * @param A input matrix A.
 * @param B input matrix B
 * @param C pointer to pointer where the result will be stored
 */
void iftMatricesMultiplicationPointWiseInPlace(iftMatrix *A, iftMatrix *B, iftMatrix **C);

/**
 * @brief Performs pointwise division between the matrices <b>A</b> and <b>B</b>, and stores the result in <b>C</b>
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 12, 2016
 *
 * The matrices <b>A</b> and <b>B</b> must have the same dimensions. If the matrix <b>C</b> is NULL, the function
 * allocates memory to C. If C has not the same dimensions of A/B, the functions destroy the C pointer, and then
 * creates a new <b>iftMatrix</b> pointer that matches the dimensions requirements. The operation is C = A / B
 *
 *
 * @param A input matrix A.
 * @param B input matrix B
 * @param C pointer to pointer where the result will be stored
 */
void iftMatricesDivisionPointWiseInPlace(iftMatrix *A, iftMatrix *B, iftMatrix **C);



/**
 * @brief Performs a mathematical pointwise operation between two matrices <b>A</b> and <b>B</b>, and stores the result
 * in <b>C</b>
 *
 * @author Deangeli Gomes Neves
 * @date july 5, 2016
 *
 * The matrices must have the same dimensions. The operation is C = A "operator" by B
 *
 * @param A input matrix A.
 * @param B input matrix B
 * @param C the matrix where the results will be stored
 * @param operationSymbol '+' to addition, '-' to subtraction, '*' to multiplication, '/' to division.
 *
 */
void iftComputeOperationBetweenMatricesInPointWiseInPlace(iftMatrix *A, iftMatrix *B,iftMatrix *C,
                                                          char operationSymbol);



void iftComputeOperationBetweenMatricesInPointWiseInPlaceSafe(iftMatrix *A, iftMatrix *B,iftMatrix **C,
                                                              char operationSymbol);


/**
 * @brief Sum of matrix elements
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 4, 2016
 *
 *
 * @param matrix input matrix.
 *
 * @return the sum of matrix elements
 */
double  iftMatrixSum(const iftMatrix *matrix);

/**
 * @brief Computes the covariance matrix from an input matrix
 * @author Cesar Castelo
 * @date Oct 8, 2019
 *
 * @param X input matrix
 * @return The covariance matrix
 */
iftMatrix *iftCovarianceMatrix(const iftMatrix *X);

/**
 * @brief Sets all diagonal elements of a matrix, given a value
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 4, 2016
 *
 * The matrix must be a square matrix
 *
 * @param matrix input matrix.
 * @param diagonalValue value to be set in diagonal elements of input matrix
 *
 */
void  iftSetDiagonalValue(iftMatrix *matrix, double diagonalValue);

/**
 *
 * @param m1
 * @param m2
 * @param sum
 * @return
 */
void iftSumMatricesInPlace(const iftMatrix* m1, const iftMatrix* m2, iftMatrix* sum);

iftMatrix* iftSubtractMatrices(const iftMatrix* m1, const iftMatrix* m2);
void iftSubtractMatricesInPlace(const iftMatrix* m1, const iftMatrix* m2, iftMatrix* sub);


/**
 *
 * @param m1
 * @param m2
 * @return
 */
iftMatrix* iftSumMatrices(const iftMatrix* m1, const iftMatrix* m2);

/**
 * @brief Computes the logarithm values of all elementos of matrix
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 4, 2016
 *
 *
 * @param matrix input matrix.
 *
 */
void  iftLogarithmMatrix(iftMatrix *matrix);


/**
 * @brief Checks if the matrices are equal.
 * @author Samuka
 * @date Dec 15, 2016
 */
bool iftIsMatricesEqual(const iftMatrix *A, const iftMatrix *B);

/**
 * @brief Checks if the matrix is symmetric
 * @author Azael Sousa
 * @date May 6, 2021
 */
bool iftIsMatrixSymmetric(iftMatrix *A);



/**
 * @brief Creates a matrix that its elements are generated randomly
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 5, 2016
 *
 * @param numberRows number of rows in matrix.
 * @param numberColumns number of columns in matrix.
 * @param mean mean value amonng the matrix elements
 * @param variance variance
 *
 * @return a matrix that contains random elements values
 */
iftMatrix *iftRandomMatrix(unsigned int numberRows, unsigned int numberColumns, double mean, double variance);

/**
 * @brief Creates a matrix that its elements are all set as <b>value</b>
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 5, 2016
 *
 * @param numberRows number of rows in matrix.
 * @param numberColumns number of columns in matrix.
 * @param value value to be set in all elements
 *
 * @return a matrix that contains a specified value in elements
 */
iftMatrix *iftCreateSetMatrix(unsigned int numberRows,unsigned int numberColumns, double value);


//TODO
void iftSetMatrix(iftMatrix* matrix, double value);


/**
 * @brief Raises by power all elements of matrix
 * @author Deangeli Gomes Neves
 * @date july 5, 2016
 *
 * @param matrix input matrix
 * @param power Power coefficient
 */
void iftPointWiseMatrixPowInPlace(iftMatrix *matrix, double power);

/**
 * @brief Raises by power all elements of matrix
 * @author Peixinho
 * @date Jul, 2016
 *
 * @param matrix input matrix
 * @param power Power coefficient
 * @return The powered matrix
 */
iftMatrix* iftPointWiseMatrixPow(iftMatrix *matrix, double power);

/**
 * @brief Computes the sum of all elements for each row vector or column vector in the input matrix
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 6, 2016
 *
 * @param matrix input matrix
 * @param vectorType 'r' to row vector, 'c' to column vector
 * @return a vector the contains the sum of elements for each row vector or column vector
 */
iftDblArray *iftComputeSumVector(iftMatrix *matrix, char vectorType);
iftMatrix* iftComputeSumVectorBycolumn(iftMatrix* matrix);
void iftComputeSumVectorBycolumnInPlace(iftMatrix* matrix,  iftMatrix **vector);
iftMatrix* iftComputeSumVectorByRow(iftMatrix* matrix);
void iftComputeSumVectorByRowInPlace(iftMatrix* matrix, iftMatrix **vector);
void iftTransposeVectorMatrix(iftMatrix* vectorMatrix);
void iftComputeSumVectorInplace(iftMatrix *matrix,iftDblArray *array , char vectorType);
void iftComputeSumVectorInplaceSafe(iftMatrix *matrix,iftDblArray **array , char vectorType);


/**
 * @brief Computes the mean of all row vectors or column vectors in the input matrix
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 6, 2016
 *
 * @param matrix input matrix
 * @param vectorType 'r' to row vector, 'c' to column vector
 * @return a vector the contains the sum of the row vectors or column vectors
 */
iftDblArray *iftComputeMeanVector(iftMatrix *matrix, char vectorType);
void iftComputeMeanVectorInPlace(iftMatrix *matrix,iftDblArray** array,char vectorType);

/**
 * @brief Create a square matrix where the diagonal elements are the vector entries
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 6, 2016
 *
 * @param vector vector
 * @return a vector the contains the sum of the row vectors or column vectors
 */
iftMatrix *iftCreateDiagonalMatrix(double *vector, long n);

void iftCreateDiagonalMatrixInPlace(double *vector, long n,iftMatrix **matrix);

/**
 * @brief Create a vector matrix containing the elements of the diagonal of the input matrix
 * @param matrix Matrix to extract the diagonal from
 * @param type Type of the resulting vector matrix ('r'=row or 'c'=column)
 * @return A vector matrix containing the diagonal from the input matrix
 *
 * @author Cesar Castelo
 * @date Oct 7, 2019
 */
iftMatrix *iftDiagonalFromMatrix(iftMatrix *matrix, char type);

/**
 * @brief Saturate the elements of matrix. The saturation may be upper bound or lower bound
 *
 *
 * @author Deangeli Gomes Neves
 * @date july 6, 2016
 *
 * @param matrix input matrix
 * @param saturationValue limit that defines the upper or lower bound
 * @param saturationType 'u' to upper bound, 'l' to lower bound.
 * @return a vector the contains the sum of the row vectors or column vectors
 */
void iftMatrixElementsSaturation(iftMatrix *matrix, double saturationValue, char saturationType);

/**
 * @brief Sets the minimum value for a given matrix.
 * @param m Input matrix
 * @param val Minimum value
 * @author Peixinho
 * @date Dec, 2016
 */
void iftMatrixMinSaturation(iftMatrix* m, float val);

iftMatrix* iftSubMatrix(iftMatrix *matrix, int startRow, int endRow, int startColumn, int endColumn);
void iftCopyMatrixInPlace(iftMatrix *A,iftMatrix *B);
void iftCopyMatrixInPlaceSafe(iftMatrix *A,iftMatrix **B);
void iftTransposeMatrixInPlace(iftMatrix *A);
void iftMultMatricesInPlace(const iftMatrix *A, const iftMatrix *B, bool transposeA, bool transposeB, iftMatrix **C);

/**
 * @brief Merges an array of iftMatrix into a new iftMatrix (vertically)
 * @param matArray Array of matrices
 * @param n Size of the array
 * @author Cesar Castelo
 * @date Jan 24, 2018
 */
iftMatrix* iftMergeMatrixArrayVertical(iftMatrix **matArray, int n);

/**
 * @brief Merges an array of iftMatrix into a new iftMatrix (horizontally)
 * @param matArray Array of matrices
 * @param n Size of the array
 * @author Cesar Castelo
 * @date Oct 17, 2019
 */
iftMatrix* iftMergeMatrixArrayHorizontal(iftMatrix **matArray, int n);

/**
 * @brief Centralize a Matrix by computing the mean value for each column (one mean per column)
 * and subtracting the values of each column for its corresponding mean.
 * @param  M Input Matrix.
 * @return   Centralized Matrix.
 */
iftMatrix *iftCentralizeMatrix(const iftMatrix *M);


/**
 * @brief Normalize a Matrix by transforming each column-vector in a unit-vector.
 * @param  M Input Matrix.
 * @return   Normalized Matrix.
 */
iftMatrix *iftUnitNormMatrix(const iftMatrix *M);

/**
 * @brief Functions to convert a matrix into a set of arrays
 * @param  M Input Matrix.
 * @return Array of arrays containing the matrix's data.
 */
iftDblArray **iftMatrixToDblArraySetByRow(iftMatrix *M);
iftFloatArray **iftMatrixToFloatArraySetByRow(iftMatrix *M);
iftDblArray **iftMatrixToDblArraySetByColumn(iftMatrix *M);
iftFloatArray **iftMatrixToFloatArraySetByColumn(iftMatrix *M);

double *iftIdentityDouble(int d);

/**
 * @brief Expand the number of rows (axis = 0) or columns (axis = 1) of a
 * given matrix.
 * @author Azael M. Sousa
 * @param  M Input Matrix.
 * @param  axis Axis to expand matrix.
 * @param  n Number of rows or columns to be added.
 * @param  value Value of the new dimension.
 */
void iftExpandMatrixInPlace(iftMatrix **M, char axis, int n, float value);

/**
 * @brief Concatenate 2 matrices in place. The output Matrix returns in M1.
 * @author Azael M. Sousa
 * @param  M1 first Matrix.
 * @param  M2 second Matrix.
 * @param  axis Axis to expand matrix.
 */
void iftConcatMatricesInPlace(iftMatrix **M1, iftMatrix *M2, char axis);

#ifdef __cplusplus
}
#endif

#endif
