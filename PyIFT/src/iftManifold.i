%include "iftDataSet.i"
%include "iftGenericMatrix.i"
%include "iftClassification.i"

%newobject iftDimReductionByTSNE;
%feature("autodoc", "2");
iftDataSet* iftDimReductionByTSNE(iftDataSet* Z, int ndim, double perplexity, size_t max_iter);

