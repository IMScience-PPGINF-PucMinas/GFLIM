
#ifndef IFT_NEIGH_COMP_ANALYSIS
#define IFT_NEIGH_COMP_ANALYSIS

#ifdef __cplusplus
extern "C" {
#endif

double *iftNeighCompAnalysis(const double *data, const int *label, const double *L_in,
                             int n, int d, int d_out, int iterations, double learn_rate);

#ifdef __cplusplus
}
#endif

#endif
