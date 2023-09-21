#include "ift/metriclearn/MetricLearnCore.h"
#include "ift/core/io/Stream.h"

extern void dgesdd_(char *jobz, int *m, int *n, double* a, int *lda, double* s, double* u, int *ldu, double* vt,
                    int* ldvt, double *work, int *lwork, int *iwork, int *info);

double *iftPCATransform(const double *data, int n, int d, int out_d) {

    if (out_d > d) {
        iftWarning("Invalid Out Dimension Input", "iftPCATransform");
        return iftAlloc(d * out_d, sizeof(double));
    }

    double *mean = iftAlloc(d, sizeof *mean);
    double *centralized = iftAlloc(n * d, sizeof *centralized);

    for (int i = 0; i < n; i++) {
        int row = i * d;
        for (int j = 0; j < d; j++) {
            mean[j] += data[row + j];
        }
    }

    for (int i = 0; i < d; i++) {
        mean[i] /= n;
    }

    for (int i = 0; i < n; i++) {
        int row = i * d;
        for (int j = 0; j < d; j++) {
            centralized[row + j] = data[row + j] - mean[j];
        }
    }

    iftFree(mean);

    double *cov_matrix = iftAlloc(d * d, sizeof *cov_matrix);
    double scale = 1.0 / (n - 1);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, d, d, n,
                scale, centralized, d, centralized, d, 0.0, cov_matrix, d);

    iftFree(centralized);

    char jobz = 'A';
    double *S = iftAlloc(d, sizeof *S);
    double *U = iftAlloc(d * d, sizeof *U);
    double *Vt = iftAlloc(d * d, sizeof *Vt);
    int lwork = -1;
    double work_param;
    int *iwork = iftAlloc(d * 8, sizeof *iwork);
    int info;

    dgesdd_(&jobz, &d, &d, cov_matrix, &d, S, U, &d, Vt, &d, &work_param, &lwork, iwork, &info);

    lwork = (int) work_param;
    double *work = iftAlloc(lwork, sizeof *work);

    dgesdd_(&jobz, &d, &d, cov_matrix, &d, S, U, &d, Vt, &d, work, &lwork, iwork, &info);

    if (info != 0)
        iftError("SVD was not successful", "iftPCAReduTransf");

    iftFree(cov_matrix);
    iftFree(S);
    iftFree(Vt);
    iftFree(work);
    iftFree(iwork);

    double *out = iftAlloc(d * out_d, sizeof *out);
    for (int i = 0; i < out_d; i++) {
        int row = i * d;
        for (int j = 0; j < d; j++) {
            out[row + j] = U[row + j];
        }
    }

    iftFree(U);

    return out;
}

double *iftSpaceTransform(const double *data, const double *L, int n, int d, int d_out)
{
    double *new_data = iftAlloc(n * d_out, sizeof *new_data);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, d_out, d, 1.0, data, d, L, d, 0.0, new_data, d_out);

    return new_data;
}


float *iftSpaceTransformFloat(const float *data, const double *L, int n, int d, int d_out)
{
    float *new_data = iftAlloc(n * d_out, sizeof *new_data);
    float *Lf = iftAlloc(d * d_out, sizeof *Lf);
    for (int i = 0; i < d * d_out; i++)
    {
        Lf[i] = (float) L[i];
    }

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, d_out, d, 1.0, data, d, Lf, d, 0.0, new_data, d_out);

    iftFree(Lf);

    return new_data;
}

iftMatrix *iftSpaceTransformMatrix(const iftMatrix *X, const double *L, int d, int d_out)
{
    if (X->ncols != d) {
        iftError("Number of Matrix Columns does not match with transformation", "iftSpaceTransformMatrix");
    }
    iftMatrix *out = iftCreateMatrix(d_out, X->nrows);
    float *aux = iftSpaceTransformFloat(X->val, L, X->nrows, X->ncols, d_out);
    iftFree(out->val);
    out->val = aux;

    return out;
}

double *iftCopyTransform(const double *L, int n_row, int n_col)
{
    double *new_L = iftAlloc(n_row * n_col, sizeof *new_L);
    for (int i = 0; i < n_row * n_col; i++)
    {
        new_L[i] = L[i];
    }

    return new_L;
}


/***** Kernel Functions *****/

double iftKernelLinear(const double *x1, const double *x2, int d, double not_used1, double not_used2)
{
    double dot = 0;
    for (int i = 0; i < d; i++) {
        dot += x1[i] * x2[i];
    }
    return dot;
}


double iftKernelPolynomial(const double *x1, const double *x2, int d, double r, double n)
{
    double dot = 0;
    for (int i = 0; i < d; i++) {
        dot += x1[i] * x2[i];
    }
    return pow((dot + r), n);
}


double iftKernelGaussian(const double *x1, const double *x2, int d, double sigma, double not_used1)
{
    double norm = 0;
    for (int i = 0; i < d; i++) {
        double diff = x1[i] * x2[i];
        norm += diff * diff;
    }
    return exp(- norm/(2 * sigma * sigma));
}


double iftKernelLaplacian(const double *x1, const double *x2, int d, double alpha, double not_used1)
{
    double norm = 0;
    for (int i = 0; i < d; i++) {
        double diff = x1[i] * x2[i];
        norm += diff * diff;
    }
    norm = sqrt(norm);
    return exp(- alpha * norm);
}


double iftKernelSigmoidal(const double *x1, const double *x2, int d, double k, double b)
{
    double dot = 0;
    for (int i = 0; i < d; i++) {
        dot += x1[i] * x2[i];
    }
    return tanh(k * dot + b);
}


double iftKernelCosine(const double *x1, const double *x2, int d, double not_used1, double not_used2)
{
    double dot = 0;
    double norm1 = 0;
    double norm2 = 0;
    for (int i = 0; i < d; i++) {
        dot += x1[i] * x2[i];
        norm1 += x1[i] * x1[i];
        norm2 += x2[i] * x2[i];
    }

    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);
    return dot / (norm1 * norm2);
}


double iftKernelWave(const double *x1, const double *x2, int d, double theta, double not_used1)
{
    double norm = 0;
    for (int i = 0; i < d; i++) {
        double diff = x1[i] * x2[i];
        norm += diff * diff;
    }
    norm = sqrt(norm);
    return theta / norm * sin(norm / theta);
}


double iftKernelLog(const double *x1, const double *x2, int d, double n, double not_used1)
{
    double norm = 0;
    for (int i = 0; i < d; i++) {
        double diff = x1[i] * x2[i];
        norm += diff * diff;
    }
    norm = sqrt(norm);
    norm = pow(norm, n);
    return -log(norm) + 1;
}


double iftKerneltStudent(const double *x1, const double *x2, int d, double n, double not_used1)
{
    double norm = 0;
    for (int i = 0; i < d; i++) {
        double diff = x1[i] * x2[i];
        norm += diff * diff;
    }
    norm = sqrt(norm);
    norm = pow(norm, n);
    return 1 / (1 + norm);
}


double iftKernelChiSqr(const double *x1, const double *x2, int d, double not_used1, double not_used2)
{
    double sum = 0;
    for (int i = 0; i < d; i++) {
        double diff = x1[i] - x2[i];
        sum += diff * diff / (x1[i] + x2[i]);
    }
    return sum;
}



iftKernelFunction *iftMetricLearnKernel(iftMetricKernel type)
{
    switch (type)
    {
        case KERNEL_LINEAR:
            return &iftKernelLinear;
        case KERNEL_POLYNOMIAL:
            return &iftKernelPolynomial;
        case KERNEL_GAUSSIAN:
            return &iftKernelGaussian;
        case KERNEL_LAPLACIAN:
            return &iftKernelLaplacian;
        case KERNEL_SIGMOIDAL:
            return &iftKernelSigmoidal;
        case KERNEL_COSINE:
            return &iftKernelCosine;
        case KERNEL_WAVE:
            return &iftKernelWave;
        case KERNEL_LOG:
            return &iftKernelLog;
        case KERNEL_TSTUDENT:
            return &iftKerneltStudent;
        case KERNEL_CHISQUARED:
            return &iftKernelChiSqr;
        default:
            iftError("Kernel type not found", "iftMetricLearnKernel");
    }
    return NULL;
}

double *iftGramianMatrix(const double *data, int n, int d, iftKernelFunction *K, double param1, double param2)
{
    double *gram = iftAlloc(n * n, sizeof *gram);
    if (!gram)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftGramMatrix");

    // currently the gram matrix is symmetric this might not be true in the future
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        int row = i * n;
        int ii = i * d;
        for (int j = i; j < n; j++) {
            int jj = j * d;
            gram[row + j] = (*K)(&data[ii], &data[jj], d, param1, param2);
            gram[j * n + i] = gram[row + j];
        }
    }

    return gram;
}

iftDoubleMatrix* iftGramMatrix(const iftDoubleMatrix *data, iftKernelFunction *K, double param1, double param2)
{
    iftDoubleMatrix *gram = iftCreateDoubleMatrixNoAlloc(data->nrows, data->nrows);
    gram->val = iftGramianMatrix(data->val, data->nrows, data->ncols, K, param1, param2);
    gram->allocated = true;
    return gram;
}


double *iftKernelFeatures(const double *train_data, const double *new_data, int n_train, int n_new, int d,
                          iftKernelFunction *K, double param1, double param2)
{
    double *feat = iftAlloc(n_train * n_new, sizeof *feat);
    if (!feat)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftKernelFeatures");

    #pragma omp parallel for
    for (int i = 0; i < n_new; i++) {
        int row = i * n_train;
        int ii = i * d;
        for (int j = 0; j < n_train; j++) {
            int jj = j * d;
            feat[row + j] = (*K)(&new_data[ii], &train_data[jj], d, param1, param2);
        }
    }

    return feat;
}


iftDoubleMatrix *iftGramTransform(const iftDoubleMatrix *gram, const iftDoubleMatrix *omega)
{
    iftDoubleMatrix *out = iftCreateDoubleMatrixNoAlloc(omega->nrows, gram->nrows);
    out->val = iftSpaceTransform(gram->val, omega->val, gram->nrows, omega->ncols, omega->nrows);
    out->allocated = true;
    return out;
}
