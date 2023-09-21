%include "iftCommon.i"
%include "iftImage.i"
%include "iftMatrix.i"
%include "iftMetrics.i"
%include "iftSeeds.i"

%feature("autodoc", "2");
double iftDiceSimilarity(const iftImage *bin_source, const iftImage *bin_target);

%newobject iftASSD;
%feature("autodoc", "2");
double iftASSD(const iftImage *bin_source, const iftImage *bin_target);

%newobject iftGeneralBalancedCoeff;
%feature("autodoc", "2");
double iftGeneralBalancedCoeff(const iftImage *bin_source, const iftImage *bin_target);

%newobject iftShannonEntropy;
%feature("autodoc", "2");
double          iftShannonEntropy(const iftImage *image);

