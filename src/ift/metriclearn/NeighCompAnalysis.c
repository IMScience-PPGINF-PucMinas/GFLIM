#include "ift/metriclearn/MetricLearnCore.h"
#include "ift/metriclearn/NeighCompAnalysis.h"
#include "ift/core/io/Stream.h"

#include "iftMatrix.h"


static inline double iftTransformedL2Norm(const double *L, const double *x, int d, int d_out, double *pre_allocated_Lx)
{
    double norm = 0;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, d_out, d, 1.0, x, d, L, d, 0.0, pre_allocated_Lx, d_out);

    for (int i = 0; i < d_out; i++) {
        norm += pre_allocated_Lx[i] * pre_allocated_Lx[i];
    }

    return norm;
}

double *iftScaleDoubleMatrix(const double *data, int n, int d)
{
    double *min = iftAlloc(d, sizeof *min);
    double *max = iftAlloc(d, sizeof *max);

    for (int i = 0; i < d; i++) {
        min[i] = IFT_INFINITY_DBL;
        max[i] = IFT_INFINITY_DBL_NEG;
    }

    for (int i = 0; i < n; i++) {
        int row = i * d;
        for (int j = 0; j < d; j++) {
            if (data[row + j] < min[j]) {
                min[j] = data[row + j];
            }
            if (data[row + j] > max[j]) {
                max[j] = data[row + j];
            }
        }
    }

    double *L = iftAlloc(d * d, sizeof *L);
    for (int i = 0; i < d; i++) {
        L[i * d + i] = 1.0 / (max[i] - min[i] + IFT_EPSILON);
    }

    iftFree(min);
    iftFree(max);

    return L;
}

void iftShuffleArray(int *array, int n)
{
    /* not very random */
    for (int i = 0; i < n; i++) {
        int j = (int) (i + random() / (RAND_MAX / (n - i) + 1));
        int t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

double *iftNeighCompAnalysis(const double *data, const int *label, const double *L_in,
                             int n, int d, int d_out, int iterations, double learn_rate)
{
    double *L = NULL;
    if (L_in == NULL)
    {
        if (d == d_out)
            L = iftScaleDoubleMatrix(data, n, d);
        else
            L = iftPCATransform(data, n, d, d_out);
    } else {
        L = iftCopyTransform(L, d_out, d);
    }

    int *idx = iftAlloc(n, sizeof *idx);
    for (int i = 0; i < n; i++) idx[i] = i;
    iftShuffleArray(idx, n);

    double *diff = iftAlloc(d, sizeof *diff);
    double *Lx = iftAlloc(d_out, sizeof *Lx);
    double *softmax = iftAlloc(n, sizeof *softmax);

    double *first = iftAlloc(d * d, sizeof *first);
    double *second = iftAlloc(d * d, sizeof *second);

    double *grad = iftAlloc(d * d_out, sizeof *grad);

    for (int ite = 0; ite < iterations; ite++)
    {
        int i = idx[ite % n];
        int i_row = i * d;

        double softmax_norm = IFT_EPSILON;
        for (int k = 0; k < n; k++) {
            if (k == i) {
                softmax[k] = 0.0;
                continue;
            }

            int k_row = k * d;
            for (int j = 0; j < d; j++) {
                diff[j] = data[i_row + j] - data[k_row + j];
            }

            softmax[k] = exp(-iftTransformedL2Norm(L, diff, d, d_out, Lx));
            softmax_norm += softmax[k];
        }

        for (int k = 0; k < n; k++) {
            softmax[k] /= softmax_norm;
        }

        double p_ik = 0.0;
        for (int k = 0; k < n; k++) {
            if (label[i] == label[k]) p_ik += softmax[k];
        }

        for (int k = 0; k < d * d; k++) {
            first[k] = 0.0;
            second[k] = 0.0;
        }

        for (int di = 0; di < d; di++) {
            int di_row = di * d;
            for (int dj = 0; dj < d; dj++) {
                for (int k = 0; k < n; k++) {
                    int k_row = k * d;
                    double xii = data[i_row + di] * data[i_row + dj];
                    double xkk = data[k_row + di] * data[k_row + dj];
                    double xik = -data[i_row + di] * data[k_row + dj];
                    double xki = -data[k_row + di] * data[i_row + dj];

                    first[di_row + dj] += softmax[k] * (xii + xkk + xik + xki);

                    if (label[k] == label[i]) {
                        second[di_row + dj] += softmax[k] * (xii + xkk + xik + xki);
                    }
                }
            }
        }

        for (int j = 0; j < d * d; j++) {
            first[j] = first[j] * p_ik - second[j];
        }

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, d_out, d, d,
                    2 * learn_rate, L, d, first, d, 0.0, grad, d);

        for (int j = 0; j < d * d_out; j++) {
            L[j] += grad[j];
        }
    }

    iftFree(idx);
    iftFree(diff);
    iftFree(Lx);
    iftFree(softmax);
    iftFree(first);
    iftFree(second);
    iftFree(grad);

    return L;
}
