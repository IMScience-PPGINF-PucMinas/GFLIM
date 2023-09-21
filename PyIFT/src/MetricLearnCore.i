%include "iftMatrix.i"



typedef enum {
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

%feature("autodoc", "2");
iftKernelFunction *iftMetricLearnKernel(iftMetricKernel type);

%newobject iftGramMatrix;
%feature("autodoc", "2");
iftDoubleMatrix* iftGramMatrix(const iftDoubleMatrix *data, iftKernelFunction *K, double param1, double param2);

%newobject iftGramTransform;
%feature("autodoc", "2");
iftDoubleMatrix *iftGramTransform(const iftDoubleMatrix *gram, const iftDoubleMatrix *omega);

