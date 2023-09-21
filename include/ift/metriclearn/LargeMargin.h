#ifndef IFT_LARGE_MARGIN
#define IFT_LARGE_MARGIN

#ifdef __cplusplus
extern "C" {
#endif

#include "iftMatrix.h"
#include "iftMImage.h"

typedef struct ift_impostor {
    int example;
    int target;
    int impostor;
    struct ift_impostor *next;
} iftImpostorSet;

void iftInsertImpostorSet(iftImpostorSet **S, int example, int target, int impostor);

void iftDestroyImpList(iftImpostorSet ***S, int n);

void iftCheckNumberOfTargets(const int *label, int n,  int *k_targets);

double *iftUpdateOmega(const double *omega, const double *grad, int n, int dim_out, double learn_rate);


static inline void iftUpdateMatrixE(double *E, const double *diff, int n, int i, int j)
{
    int ii = i * n;
    int jj = j * n;
    for (int k = 0; k < n; k++) {
        E[ii + k] += diff[k];
        E[jj + k] -= diff[k];
    }
}

static inline void iftDiffKiKj(const double *gram, double *diff, int n, int i, int j, double c)
{
    int ii = i * n;
    int jj = j * n;
    for (int k = 0; k < n; k++) {
        diff[k] = c * (gram[ii + k] - gram[jj + k]);
    }
}


/**
 * @brief Implementation of Distance Metric Learning for Large Margin Nearest Neighbor Classification
 * http://www.jmlr.org/papers/v10/weinberger09a.html
 * @author Jordão Bragantini
 * @param data          [n x d] matrix of features
 * @param label         [n x 1] vector of labels
 * @param L             [d_out x d] starting distance transform
 * @param n             number of samples
 * @param d             dimension of original space
 * @param d_out         dimension of project space
 * @param k_targets     number of targets
 * @param learn_rate    learning rate for gradient descent
 * @param iterations    maximum number of iterations
 * @param verbose       print information or not
 * @return optimal distance transform matrix
 */

//! swig(newobject)
double *iftLMCompAnalysis(const double *data, const int *label, const double *L_in, int n, int d, int d_out,
                          int k_targets, double learn_rate, int iterations, bool verbose);

/**
 * @param gram          Gram Matrix of kernel of choice, (n, n)
 * @param label         Label array (n)
 * @param n             n
 * @param dim_out       dim_out
 * @param k_targets     k_targets
 * @param c             scale of impostor penalization
 * @param learn_rate    learning rate
 * @param iterations    iterations
 * @param verbose       verbose
 * @return              Matrix of weights omega (dim_out, n), equivalent to L
 */
//! swig(newobject)
double *iftKernelLMCA(const double *gram, const int *label, int n, int dim_out,
                      int k_targets, double c, double learn_rate, int iterations, bool verbose);

//! swig(newobject)
iftDoubleMatrix *iftKLMCA(const iftDoubleMatrix *gram, const iftIntArray *label, int dim_out,
                          int k_targets, double c, double learn_rate, int iterations, bool verbose);


#ifdef __cplusplus
}
#endif

#endif //IFT_LARGE_MARGIN
