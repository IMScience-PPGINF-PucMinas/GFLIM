#include "ift/metriclearn/MetricLearnCore.h"
#include "ift/metriclearn/LargeMargin.h"

#include "ift/core/io/Stream.h"
#include <omp.h>
#include <ift/core/dtypes/IntArray.h>

void iftInsertImpostorSet(iftImpostorSet **S, int example, int target, int impostor)
{
    iftImpostorSet *i = iftAlloc(1, sizeof *i);
    i->example = example;
    i->target = target;
    i->impostor = impostor;
    i->next = *S;
    *S = i;
}


iftImpostorSet *iftCopyImpostorSet(const iftImpostorSet *S)
{
    iftImpostorSet *C = NULL;
    for (const iftImpostorSet *s = S; s != NULL; s = s->next) {
        iftInsertImpostorSet(&C, s->example, s->target, s->impostor);
    }

    return C;
}


int iftImpostorSetSize(const iftImpostorSet *S)
{
    int i = 0;
    for (const iftImpostorSet *s = S; s != NULL; s = s->next) {
        i++;
    }
    return i;
}


void iftDestroyImpostorSet(iftImpostorSet **S)
{
    iftImpostorSet *i = NULL;
    while (*S != NULL)
    {
        i = *S;
        *S = i->next;
        iftFree(i);
    }
    *S = NULL;
}


void iftDestroyImpList(iftImpostorSet ***S, int n)
{
    iftImpostorSet **aux = *S;
    for (int i = 0; i < n; i++)
        iftDestroyImpostorSet(&aux[i]);
    iftFree(aux);
    aux = NULL;
}

int iftImpSetListSize(iftImpostorSet **S, int n)
{
    int count = 0;
    for (int i = 0; i < n; i++) {
        for (iftImpostorSet* s = S[i]; s != NULL; s = s->next)
            count++;
    }

    return count;
}


void iftCheckNumberOfTargets(const int *label, int n,  int *k_targets)
{
    int n_labels = iftMaxValIntArray(label, n, NULL) + 1;
    int *n_per_class = iftAlloc(n_labels, sizeof *n_per_class);
    for (int i = 0; i < n; i++)
    {
        int lb = label[i];
        n_per_class[lb]++;
    }

    int min_k_targets = IFT_INFINITY_INT;
    for (int i = 0; i < n_labels; i++)
    {
        if (min_k_targets > n_per_class[i])
            min_k_targets = n_per_class[i];
    }

    iftFree(n_per_class);

    if ((min_k_targets-1) < *k_targets)
    {
        *k_targets = (min_k_targets - 1);
        char buf[IFT_STR_DEFAULT_SIZE];
        sprintf(buf, "Number of target neighbors reduced to %d, number of samples with same class not sufficient.", *k_targets);
        iftWarning(buf, "iftCheckNumberOfTargets");
    }
}

int *iftSelectTargetNeighbors(const double *data, const int *label, int n, int d, int k_targets)
{
    int *targets = iftAlloc(k_targets * n, sizeof *targets);
    double *t_dist = iftAlloc(k_targets * n, sizeof *t_dist);
    double *max_dist = iftAlloc(n, sizeof *max_dist);

    for (int i = 0; i < k_targets * n; i++) {
        t_dist[i] = IFT_INFINITY_DBL;
        targets[i] = -1;
    }

    for (int i = 0; i < n; i++) {
        max_dist[i] = IFT_INFINITY_DBL;
    }

    double *dist_tab = NULL;
    // checking if its able to allocate memory
    if ((((long) n) * ((long) n)) == (n * n)) {
        dist_tab = iftAlloc(n * n, sizeof *dist_tab);
    }


    if (dist_tab != NULL) {
        #pragma omp parallel for
        for (int i = 0; i < (n - 1); i++) {
            int row_idx = i * n;
            for (int j = i + 1; j < n; j++) {
                if (label[i] == label[j]) {
                    double dist = iftFeatureDistance(&data[i * d], &data[j * d], d);
                    dist_tab[row_idx + j] = dist;
                    dist_tab[j * n + i] = dist;
                }
            }
        }
    }

    # pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            if (label[i] == label[j] && i != j)
            {
                double dist = ((dist_tab != NULL) ? dist_tab[i * n + j] :
                               iftFeatureDistance(&data[i * d], &data[j * d], d));

                if (dist < max_dist[i]) {
                    int first = k_targets * i;
                    int last = k_targets * (i + 1) - 1;
                    int idx = last;
                    while (idx > first && dist < t_dist[idx - 1]) {
                        t_dist[idx] = t_dist[idx - 1];
                        targets[idx] = targets[idx - 1];
                        idx--;
                    }

                    t_dist[idx] = dist;
                    targets[idx] = j;
                    max_dist[i] = t_dist[last];
                }
            }
        }
    }

    iftFree(t_dist);
    iftFree(max_dist);
    iftFree(dist_tab);

    return targets;
}


iftImpostorSet *iftFindImpostors(const double *Ldata, const int *label, const int *targets, int n, int d, int k_targets)
{
    double *t_dist = iftAlloc(n * k_targets, sizeof *t_dist);

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        int t_row = i * k_targets;
        for (int j = 0; j < k_targets; j++)
        {
            int tk = targets[t_row + j];
            t_dist[t_row + j] = iftFeatureDistance(&Ldata[i * d], &Ldata[tk * d], d);
        }
    }

    const int n_threads = 8;
    iftImpostorSet **S_ary = iftAlloc(n_threads, sizeof *S_ary);
    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < n; i++) {
        int t_id = omp_get_thread_num();
        int t_row = i * k_targets;
        for (int ik = 0; ik < n; ik++)
        {
            if (label[i] != label[ik])
            {
                double i_dist = iftFeatureDistance(&Ldata[i * d], &Ldata[ik * d], d);
                for (int j = 0; j < k_targets; j++)
                {
                    double diff = t_dist[t_row + j] + 1 - i_dist;
                    if (diff > 0) {
                        int tk = targets[t_row + j];
                        iftInsertImpostorSet(&S_ary[t_id], i, tk, ik);
                    }
                }
            }
        }
    }

    iftFree(t_dist);

    iftImpostorSet *S = S_ary[0];
    for (int i = 1; i < n_threads; i++) {
        iftImpostorSet *last = NULL;
        for (iftImpostorSet *s = S_ary[i]; s != NULL; s = s->next) {
            last = s;
        }
        if (last != NULL) {
            last->next = S;
            S = S_ary[i];
        }
    }

    iftFree(S_ary);

    return S;
}


iftImpostorSet *iftFindKImpostors(const double *Ldata, const int *label, const int *targets,
                                  int n, int d, int k_targets, int k_impostors)
{
    double *t_dist = iftAlloc(n * k_targets, sizeof *t_dist);

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        int t_row = i * k_targets;
        for (int j = 0; j < k_targets; j++)
        {
            int tk = targets[t_row + j];
            t_dist[t_row + j] = iftFeatureDistance(&Ldata[i * d], &Ldata[tk * d], d);
        }
    }

    const int n_threads = 8;
    iftImpostorSet **S_ary = iftAlloc(n_threads, sizeof *S_ary);
    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < n; i++)
    {
        int count = 0;
        int t_id = omp_get_thread_num();
        int t_row = i * k_targets;
        for (int ik = 0; ik < n; ik++)
        {
            if (count >= k_impostors) break;

            if (label[i] != label[ik])
            {
                double i_dist = iftFeatureDistance(&Ldata[i * d], &Ldata[ik * d], d);
                bool is_imp = false;
                for (int j = 0; j < k_targets; j++)
                {
                    double diff = t_dist[t_row + j] + 1 - i_dist;
                    if (diff > 0) {
                        int tk = targets[t_row + j];
                        iftInsertImpostorSet(&S_ary[t_id], i, tk, ik);
                        is_imp = true;
                    }
                }
                if (is_imp == true) count++;
            }
        }
    }

    iftFree(t_dist);

    iftImpostorSet *S = S_ary[0];
    for (int i = 1; i < n_threads; i++) {
        iftImpostorSet *last = NULL;
        for (iftImpostorSet *s = S_ary[i]; s != NULL; s = s->next) {
            last = s;
        }
        if (last != NULL) {
            last->next = S;
            S = S_ary[i];
        }
    }

    iftFree(S_ary);

    return S;
}


iftImpostorSet *iftFindPartialImpostors(const iftImpostorSet *old_imp, const double *Ldata, int d, iftImpostorSet **extra)
{
    iftImpostorSet *S = NULL;
    for (const iftImpostorSet *s = old_imp; s != NULL; s = s->next)
    {
        int i = s->example;
        int tk = s->target;
        int ik = s->impostor;
        double t_dist = iftFeatureDistance(&Ldata[i * d], &Ldata[tk * d], d);
        double i_dist = iftFeatureDistance(&Ldata[i * d], &Ldata[ik * d], d);
        double diff = t_dist + 1 - i_dist;
        if (diff > 0) {
            iftInsertImpostorSet(&S, i, tk, ik);
        } else {
            iftInsertImpostorSet(extra, i, tk, ik);
        }
    }

    return S;
}


double *iftUpdateTransform(const double *grad_imp, const double *grad_target, const double *L, double learn_rate, int d, int d_out)
{
    double *new_L = iftAlloc(d * d_out, sizeof *new_L);
    double *sum_grad = iftAlloc(d * d, sizeof *sum_grad);
    double *prod_grad = iftAlloc(d * d_out, sizeof *prod_grad);

    for (int i = 0; i < d * d; i++) {
        sum_grad[i] = grad_imp[i] + grad_target[i];
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, d_out, d, d, 2 * learn_rate,
                L, d, sum_grad, d, 0.0, prod_grad, d);

    for (int i = 0; i < d * d_out; i++) {
        new_L[i] = L[i] - prod_grad[i];
    }

    iftFree(sum_grad);
    iftFree(prod_grad);

    return new_L;
}


double *iftTargetGradient(const double *data, const int *targets, int n, int d, int k_targets)
{
    double *grad = iftAlloc(d * d, sizeof *grad);

    #pragma omp parallel for
    for (int di = 0; di < d; di++)
    {
        int row_idx = di * d;
        for (int dj = 0; dj < d; dj++)
        {
            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < k_targets; k++)
                {
                    int tk = targets[i * k_targets + k];

                    grad[row_idx + dj] +=  + data[i * d + di] * data[i * d + dj]
                                           + data[tk * d + di] * data[tk * d + dj]
                                           - data[i * d + di] * data[tk * d + dj]
                                           - data[tk * d + di] * data[i * d + dj];
                }
            }
        }
    }

    return grad;
}


double *iftImpostorGradient(const double *data, const iftImpostorSet *S, int d)
{
    double *grad = iftAlloc(d * d, sizeof *grad);

    #pragma omp parallel for
    for (int di = 0; di < d; di++)
    {
        int row_idx = di * d;
        for (int dj = 0; dj < d; dj++)
        {
            for (const iftImpostorSet *s = S; s != NULL; s = s->next)
            {
                int i = s->example;
                int tk = s->target;
                int ik = s->impostor;

                /* impostor */
                grad[row_idx + dj] += - data[i * d + di] * data[i * d + dj]
                                      - data[ik * d + di] * data[ik * d + dj]
                                      + data[i * d + di] * data[ik * d + dj]
                                      + data[ik * d + di] * data[i * d + dj];

                /* target */
                grad[row_idx + dj] += + data[i * d + di] * data[i * d + dj]
                                      + data[tk * d + di] * data[tk * d + dj]
                                      - data[i * d + di] * data[tk * d + dj]
                                      - data[tk * d + di] * data[i * d + dj];
            }
        }
    }

    return grad;
}


double *iftPartialImpGradient(const double *old_grad, const double *data, const iftImpostorSet *missing,
                              const iftImpostorSet *extra, int d)
{
    double *grad = iftAlloc(d * d, sizeof *grad);

    for (int i = 0; i < d * d; i++) {
        grad[i] = old_grad[i];
    }

    #pragma omp parallel for
    for (int di = 0; di < d; di++)
    {
        int row_idx = di * d;
        for (int dj = 0; dj < d; dj++)
        {
            for (const iftImpostorSet *s = missing; s != NULL; s = s->next)
            {
                int i = s->example;
                int tk = s->target;
                int ik = s->impostor;

                /* impostor */
                grad[row_idx + dj] += - data[i * d + di] * data[i * d + dj]
                                      - data[ik * d + di] * data[ik * d + dj]
                                      + data[i * d + di] * data[ik * d + dj]
                                      + data[ik * d + di] * data[i * d + dj];

                /* target */
                grad[row_idx + dj] += + data[i * d + di] * data[i * d + dj]
                                      + data[tk * d + di] * data[tk * d + dj]
                                      - data[i * d + di] * data[tk * d + dj]
                                      - data[tk * d + di] * data[i * d + dj];
            }

            for (const iftImpostorSet *s = extra; s != NULL; s = s->next)
            {
                int i = s->example;
                int tk = s->target;
                int ik = s->impostor;

                /* impostor */
                grad[row_idx + dj] += + data[i * d + di] * data[i * d + dj]
                                      + data[ik * d + di] * data[ik * d + dj]
                                      - data[i * d + di] * data[ik * d + dj]
                                      - data[ik * d + di] * data[i * d + dj];

                /* target */
                grad[row_idx + dj] += - data[i * d + di] * data[i * d + dj]
                                      - data[tk * d + di] * data[tk * d + dj]
                                      + data[i * d + di] * data[tk * d + dj]
                                      + data[tk * d + di] * data[i * d + dj];
            }
        }
    }

    return grad;
}


// L_row = d_out, L_col = d
double iftGradientLoss(const double *grad, const double *L, int L_nrow, int L_ncol)
{
    double *M = iftAlloc(L_ncol * L_ncol, sizeof *M);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, L_ncol, L_ncol, L_nrow, 1.0, L,
                L_ncol, L, L_ncol, 0.0, M, L_ncol);

    double loss = 0;
    for (int i = 0; i < L_ncol * L_ncol; i++) {
        loss += M[i] * grad[i];
    }

    iftFree(M);

    return loss;
}


int iftCompareImpostors(const iftImpostorSet *a, const iftImpostorSet *b)
{
    if (a->example > b->example) return 1;
    if (a->example < b->example) return -1;

    if (a->impostor > b->impostor) return 1;
    if (a->impostor < b->impostor) return -1;

    if (a->target > b->target) return 1;
    if (a->target < b->target) return -1;

    return 0;
}


void iftImpostorsDiff(const iftImpostorSet *imp, const iftImpostorSet *next_imp,
                      iftImpostorSet **missing, iftImpostorSet **extra)
{
    const iftImpostorSet *old = imp;
    const iftImpostorSet *new = next_imp;
    /* sorted decreasing order */
    while (old != NULL && new != NULL) {
        int cmp = iftCompareImpostors(old, new);
        if (cmp == 0) {
            old = old->next;
            new = new->next;
        } else if (cmp > 0) {
            iftInsertImpostorSet(extra, old->example, old->target, old->impostor);
            old = old->next;
        } else {
            iftInsertImpostorSet(missing, new->example, new->target, new->impostor);
            new = new->next;
        }
    }

    for ( ; old != NULL; old = old->next) {
        iftInsertImpostorSet(extra, old->example, old->target, old->impostor);
    }

    for ( ; new != NULL; new = new->next) {
        iftInsertImpostorSet(missing, new->example, new->target, new->impostor);
    }
}


double *iftLMCompAnalysis(const double *data, const int *label, const double *L_in, int n, int d, int d_out,
                          int k_targets, double learn_rate, int iterations, bool verbose)
{
    double *L = NULL;
    if (L_in != NULL) {
        L = iftCopyTransform(L_in, d_out, d);
    } else {
        L = iftPCATransform(data, n, d, d_out);
	}

    int *targets = iftSelectTargetNeighbors(data, label, n, d, k_targets);
    double *grad_target = iftTargetGradient(data, targets, n, d, k_targets);

    if (verbose) {
        printf("Targets neighbors found and static gradient computed\n");
    }

    double *Ldata = iftSpaceTransform(data, L, n, d, d_out);

    iftImpostorSet *imp = iftFindImpostors(Ldata, label, targets, n, d_out, k_targets);
    iftFree(Ldata);

    if (imp == NULL) {
        iftWarning("No impostor found", "iftLMCompAnalysis");
        goto no_impostor;
    }

    double *grad_imp = iftImpostorGradient(data, imp, d);

    int count_imp = iftImpostorSetSize(imp);
    double loss = iftGradientLoss(grad_target, L, d_out, d) + iftGradientLoss(grad_imp, L, d_out, d) + count_imp;

    int ite = 0;
    double delta;
    do {
        ite++;
        do {
            double *next_L = iftUpdateTransform(grad_imp, grad_target, L, learn_rate, d, d_out);
            double *next_Ldata = iftSpaceTransform(data, next_L, n, d, d_out);

            iftImpostorSet *next_imp = NULL;
            iftImpostorSet *missing = NULL;
            iftImpostorSet *extra = NULL;

            /* great optimization but produces approximate results and was not converging
             * i think it only works when computing the gradient with the M matrix and not L, given M = L^t * L
             */
//            if (ite % 10 == 0) {
//                next_imp = iftFindImpostors(next_Ldata, label, targets, n, d_out, k_targets);
//                iftImpostorsDiff(imp, next_imp, &missing, &extra);
//            } else {
//                next_imp = iftFindPartialImpostors(imp, next_Ldata, d_out, &extra);
//            }

            next_imp = iftFindImpostors(next_Ldata, label, targets, n, d_out, k_targets);

            iftImpostorsDiff(imp, next_imp, &missing, &extra);
            double *next_grad_imp = iftPartialImpGradient(grad_imp, data, missing, extra, d);
            iftDestroyImpostorSet(&missing);
            iftDestroyImpostorSet(&extra);

            count_imp = iftImpostorSetSize(next_imp);

            double next_loss = iftGradientLoss(next_grad_imp, next_L, d_out, d) +
                               iftGradientLoss(grad_target, next_L, d_out, d) + count_imp;
            delta = next_loss - loss;

            if (delta > 0) {
                learn_rate /= 2;
                iftFree(next_Ldata);
                iftFree(next_L);
                iftDestroyImpostorSet(&next_imp);
                iftFree(next_grad_imp);
            } else {
                loss = next_loss;
                learn_rate *= 1.01;
                iftFree(next_Ldata);
                iftFree(L);
                iftDestroyImpostorSet(&imp);
                iftFree(grad_imp);
                L = next_L;
                imp = next_imp;
                grad_imp = next_grad_imp;
            }
        } while (learn_rate > 1e-22 && delta > 0);

        if (verbose) {
            printf("---------------------\n"
                   "Iteration: %d\nLoss: %lf\nDelta: %lf\nActive Impostors: %d\nLearning Rate: %e\n",
                   ite, loss, delta, count_imp, learn_rate);
        }

    } while (ite < iterations && count_imp > 0 && fabs(delta) > IFT_EPSILON && learn_rate > 1e-22);

    iftFree(grad_imp);
    no_impostor: iftDestroyImpostorSet(&imp);
    iftFree(grad_target);
    iftFree(targets);

    return L;
}


int *iftSelectTargetKLMCA(const double *gram, const int *label, int n, int k_targets)
{
    int *targets = iftAlloc(k_targets * n, sizeof *targets);
    double *t_dist = iftAlloc(k_targets * n, sizeof *t_dist);
    double *max_dist = iftAlloc(n, sizeof *max_dist);

    for (int i = 0; i < k_targets * n; i++) {
        t_dist[i] = IFT_INFINITY_DBL;
        targets[i] = -1;
    }

    for (int i = 0; i < n; i++) {
        max_dist[i] = IFT_INFINITY_DBL;
    }

    # pragma omp parallel for
    for (int i = 0; i < n; i++) {
        int row = i * n;
        for (int j = 0; j < n; j++)
        {
            if (label[i] == label[j] && i != j)
            {
                double dist = gram[row + j]; // K(xi, xj)
                if (dist < max_dist[i]) {
                    int first = k_targets * i;
                    int last = k_targets * (i + 1) - 1;
                    int idx = last;
                    while (idx > first && dist < t_dist[idx - 1]) {
                        t_dist[idx] = t_dist[idx - 1];
                        targets[idx] = targets[idx - 1];
                        idx--;
                    }

                    t_dist[idx] = dist;
                    targets[idx] = j;
                    max_dist[i] = t_dist[last];
                }
            }
        }
    }

    iftFree(t_dist);
    iftFree(max_dist);

    return targets;
}


double *iftGradientKLMCA(const double *omega, const double *gram, iftImpostorSet **imp, const int *targets,
                         int n, int dim_out, int k_targets, double c)
{
    /*
     * Note: the matrices E are transpose because it's more cache efficient
     */
    double *E = iftAlloc(n * n, sizeof *E);
    if (!E)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftTargetGradKLMCA");

    double *diff = iftAlloc(n, sizeof *diff);

    for (int i = 0; i < n; i++)
    {
        int count = 0;
        for (iftImpostorSet *s = imp[i]; s != NULL; s = s->next)
        {
            count++;
            iftDiffKiKj(gram, diff, n, i, s->impostor, - c); // c is negative because this is the opposite of the target
            iftUpdateMatrixE(E, diff, n, i, s->impostor);
        }

        int row = i * k_targets;
        for (int j = 0; j < k_targets; j++) {
            int tk = targets[row + j];
            iftDiffKiKj(gram, diff, n, i, tk, count * c + 1);
            iftUpdateMatrixE(E, diff, n, i, tk);
        }
    }

    iftFree(diff);
    double *grad = iftAlloc(n * dim_out, sizeof *grad);
    if (!grad)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftTargetGradKLMCA");

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                dim_out, n, n, 2.0, omega, n, E, n, 0.0, grad, n);

    iftFree(E);

    return grad;
}


double *iftUpdateOmega(const double *omega, const double *grad, int n, int dim_out, double learn_rate)
{
    double *next_omega = iftAlloc(n * dim_out, sizeof *next_omega);

    for (int i = 0; i < n * dim_out; i++)
        next_omega[i] = omega[i] - learn_rate * grad[i];

    return next_omega;
}

iftImpostorSet **iftFindImpKLMCA(const double *Ldata, const int *label, const int *targets,
                                 int n, int dim_out, int k_targets, double c, double *out_loss, int *imp_count)
{

    double *t_dist = iftAlloc(n * k_targets, sizeof *t_dist);

    double loss = 0;
    #pragma omp parallel for reduction(+:loss)
    for (int i = 0; i < n; i++)
    {
        int i_targ = i * k_targets;
        int ii = i * dim_out;
        for (int j = 0; j < k_targets; j++)
        {
            int tk = targets[i * k_targets + j];
            int jj = tk * dim_out;
            t_dist[i_targ + j] = iftFeatureDistance(&Ldata[ii], &Ldata[jj], dim_out);
            loss += t_dist[i_targ + j];
        }
    }


    int count = 0;
    iftImpostorSet **imp = iftAlloc(n, sizeof *imp);
    #pragma omp parallel for reduction(+:loss)
    for (int i = 0; i < n; i++)
    {
        int d_row = i * dim_out;
        for (int ik = 0; ik < n; ik++)
        {
            if (label[i] != label[ik])
            {
                double i_dist = iftFeatureDistance(&Ldata[d_row], &Ldata[ik * dim_out], dim_out);
                for (int j = 0; j < k_targets; j++)
                {
                    int tt = i * k_targets + j;
                    double diff = t_dist[tt] + 1 - i_dist;
                    if (diff > 0) {
                        int tk = targets[tt];
                        iftInsertImpostorSet(&imp[i], i, tk, ik);
                        loss += c * diff;
                        #pragma omp atomic
                        count++;
                    }
                }
            }
        }
    }

    *out_loss = loss;
    *imp_count = count;

    iftFree(t_dist);

    return imp;
}


double *iftKernelLMCA(const double *gram, const int *label, int n, int dim_out,
                      int k_targets, double c, double learn_rate, int iterations, bool verbose)
{
    iftCheckNumberOfTargets(label, n, &k_targets);

    double *omega = iftPCATransform(gram, n, n, dim_out);
	
	if (k_targets <= 0)
		return omega;

    int *targets = iftSelectTargetKLMCA(gram, label, n, k_targets);

    if (verbose) {
        printf("Targets neighbors found and static gradient computed\n");
    }

    double *Ldata = iftSpaceTransform(gram, omega, n, n, dim_out);

    double loss = 0;
    int imp_count = 0;
    iftImpostorSet **imp = iftFindImpKLMCA(Ldata, label, targets, n, dim_out, k_targets, c, &loss, &imp_count);
    iftFree(Ldata);

    if (!imp_count) {
        iftWarning("No impostor found", "iftLMCompAnalysis");
        goto no_impostor;
    }

    double *grad = iftGradientKLMCA(omega, gram, imp, targets, n, dim_out, k_targets, c);

    int ite = 0;
    double delta;
    do {
        ite++;
        do {
            double *next_omega = iftUpdateOmega(omega, grad, n, dim_out, learn_rate);
            double *next_Ldata = iftSpaceTransform(gram, next_omega, n, n, dim_out);

            iftDestroyImpList(&imp, n);
            double next_loss = 0;
            imp = iftFindImpKLMCA(next_Ldata, label, targets, n, dim_out, k_targets, c, &next_loss, &imp_count);
            iftFree(next_Ldata);

            if (!imp_count) {
                if (verbose)
                    printf("No impostor found, done!\n");
                goto done;
            }

            delta = next_loss - loss;
            if (delta > 0) {
                learn_rate /= 2;
                iftFree(next_omega);
            } else {
                loss = next_loss;
                learn_rate *= 1.01;
                iftFree(omega);
                iftFree(grad);
                omega = next_omega;
                grad = iftGradientKLMCA(omega, gram, imp, targets, n, dim_out, k_targets, c);
            }
        } while (learn_rate > 1e-22 && delta > 0);

        if (verbose) {
            printf("---------------------\n"
                   "Iteration: %d\nLoss: %lf\nDelta: %lf\nActive Impostors: %d\nLearning Rate: %e\n",
                   ite, loss, delta, imp_count, learn_rate);
        }

    } while (ite < iterations && imp && fabs(delta) > IFT_EPSILON && learn_rate > 1e-22);

    done: iftFree(grad);
    no_impostor: iftDestroyImpList(&imp, n);
    iftFree(targets);

	return omega;
}


iftDoubleMatrix *iftKLMCA(const iftDoubleMatrix *gram, const iftIntArray *label, int dim_out,
                          int k_targets, double c, double learn_rate, int iterations, bool verbose)
{
    if (label->n != gram->nrows) {
        iftError("Gram matrix and labels must have the same length\n", "iftKLMCA");
    }

    iftDoubleMatrix *omega = iftCreateDoubleMatrixNoAlloc(gram->nrows, dim_out);
    omega->val = iftKernelLMCA(gram->val, label->val, gram->nrows, dim_out, k_targets, c, learn_rate, iterations, verbose);
    omega->allocated = true;
    return omega;
}
