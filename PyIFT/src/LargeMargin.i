%include "iftMatrix.i"
%include "iftMImage.i"

%newobject iftLMCompAnalysis;
%feature("autodoc", "2");
double *iftLMCompAnalysis(const double *data, const int *label, const double *L_in, int n, int d, int d_out,
                          int k_targets, double learn_rate, int iterations, bool verbose);

%newobject iftKernelLMCA;
%feature("autodoc", "2");
double *iftKernelLMCA(const double *gram, const int *label, int n, int dim_out,
                      int k_targets, double c, double learn_rate, int iterations, bool verbose);

%newobject iftKLMCA;
%feature("autodoc", "2");
iftDoubleMatrix *iftKLMCA(const iftDoubleMatrix *gram, const iftIntArray *label, int dim_out,
                          int k_targets, double c, double learn_rate, int iterations, bool verbose);

