
#ifndef IFT_METRIC_LEARN_CORE
#define IFT_METRIC_LEARN_CORE

#ifdef __cplusplus
extern "C" {
#endif

#include "iftMatrix.h"

double *iftPCATransform(const double *data, int n, int d, int out_d);

double *iftSpaceTransform(const double *data, const double *L, int n, int d, int d_out);

float *iftSpaceTransformFloat(const float *data, const double *L, int n, int d, int d_out);

iftMatrix *iftSpaceTransformMatrix(const iftMatrix *X, const double *L, int d, int d_out);

double *iftCopyTransform(const double *L, int n_row, int n_col);

static inline double iftFeatureDistance(const double *x_1, const double *x_2, int d)
{
    double dist = 0.0;
    for (int i = 0; i < d; i++) {
        double diff = x_1[i] - x_2[i];
        dist += diff * diff;
    }
    return dist;
}


/***** Begin Kernel Functions *****/

typedef double (iftKernelFunction)(const double *x1, const double *x2, int d, double param1, double param2);

/**
 * @details K(x1, x2) = x1^T * x2
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param not_used1
 * @param not_used2
 * @return          K(x1, x2)
 */
double iftKernelLinear(const double *x1, const double *x2, int d, double not_used1, double not_used2);


/**
 * @details K(x1, x2) = (x1^T * x2 + r)^n
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param r         r > 0
 * @param n         n > 0
 * @return          K(x1, x2)
 */
double iftKernelPolynomial(const double *x1, const double *x2, int d, double r, double n);


/**
 * @details K(x1, x2) = exp(-||x1 - x2||^2/(2 * sigma^2))
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param sigma     sigma (standard deviation)
 * @param not_used1
 * @return          K(x1, x2)
 */
double iftKernelGaussian(const double *x1, const double *x2, int d, double sigma, double not_used1);


/**
 * @details K(x1, x2) = exp(-alpha||x1 - x2||)
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param alpha     alpha > 0
 * @param not_used1
 * @return          K(x1, x2)
 */
double iftKernelLaplacian(const double *x1, const double *x2, int d, double alpha, double not_used1);


/**
 * @details K(x1, x2) = tanh(k * x1^T * x2 + b)
 * @brief Similar to the sigmoid function, may result in overfitting, ranges from -1 to 1;
 * @attention This kernel is only conditionally positive definite
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param k         k, usually equal to 1/N, where N is the data dimension
 * @param b         b < 0
 * @return          K(x1, x2)
 */
double iftKernelSigmoidal(const double *x1, const double *x2, int d, double k, double b);


/**
 * @details K(x1, x2) = x1^T * x2 / (||x1|| * ||x2||)
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param not_used1
 * @param not_used2
 * @return          K(x1, x2)
 */
double iftKernelCosine(const double *x1, const double *x2, int d, double not_used1, double not_used2);


/**
 * @details K(x1, x2) = (theta / ||x1 - x2||) * sin(||x1 - x2|| / theta)
 * @attention This kernel is only positive definite on R^3
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param theta     theta
 * @param not_used1
 * @return          K(x1, x2)
 */
double iftKernelWave(const double *x1, const double *x2, int d, double theta, double not_used1);


/**
 * @details K(x1, x2) = - log(||x1 - x2||^n) + 1
 * @attention This kernel is only conditionally positive definite
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param n         n
 * @param not_used1
 * @return          K(x1, x2)
 */
double iftKernelLog(const double *x1, const double *x2, int d, double n, double not_used1);


/**
 * @details K(x1, x2) = 1 / (1 + ||x1 - x2||^n)
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param n         n
 * @param not_used1
 * @return          K(x1, x2)
 */
double iftKerneltStudent(const double *x1, const double *x2, int d, double n, double not_used1);


/**
 * @details K(x1, x2) = sum_i^d (x1i - x2i)^2 / (x1i + x2i)
 * @cite http://www.robots.ox.ac.uk/~vedaldi/assets/pubs/vedaldi11efficient.pdf
 * @param x1        d-dimensional array
 * @param x2        d-dimensional array
 * @param d         d
 * @param not_used1
 * @param not_used2
 * @return          K(x1, x2)
 */
double iftKernelChiSqr(const double *x1, const double *x2, int d, double not_used1, double not_used2);


/***** End Kernel Functions *****/

//! swig()
typedef enum ift_ml_kernel {
    KERNEL_LINEAR,
    KERNEL_POLYNOMIAL,
    KERNEL_GAUSSIAN,
    KERNEL_LAPLACIAN,
    KERNEL_SIGMOIDAL,
    KERNEL_COSINE,
    KERNEL_WAVE,
    KERNEL_LOG,
    KERNEL_TSTUDENT,
    KERNEL_CHISQUARED
} iftMetricKernel;

//! swig()
iftKernelFunction *iftMetricLearnKernel(iftMetricKernel type);


/**
 * @param data          Data (n, d)
 * @param n             Data length
 * @param d             Feature vector length
 * @param K             Kernel Function
 * @param param1        Kernel function first parameter
 * @param param2        Kernel function second parameter
 * @return              Gram matrix (n, n)
 */
double *iftGramianMatrix(const double *data, int n, int d, iftKernelFunction *K, double param1, double param2);

//! swig(newobject)
iftDoubleMatrix* iftGramMatrix(const iftDoubleMatrix *data, iftKernelFunction *K, double param1, double param2);

/**
 * @param train_data    Data used to train kernel method (n_train, d)
 * @param new_data      Data to be project into kernel space (n_new, d)
 * @param n_train       Train data length
 * @param n_new         New data length
 * @param d             Feature vector length
 * @param K             Kernel Function
 * @param param1        Kernel function first parameter
 * @param param2        Kernel function second parameter
 * @return              Feature matrix (n_new, n_train)
 */
double *iftKernelFeatures(const double *train_data, const double *new_data, int n_train, int n_new, int d,
                          iftKernelFunction *K, double param1, double param2);

//! swig(newobject)
iftDoubleMatrix *iftGramTransform(const iftDoubleMatrix *gram, const iftDoubleMatrix *omega);

#ifdef __cplusplus
}
#endif

#endif
