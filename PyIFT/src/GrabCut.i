%include "iftImage.i"
%include "iftMImage.i"
%include "iftDataSet.i"
%include "iftKmeans.i"
%include "iftMaxflow.i"

%newobject iftUncertMapToLabel;
%feature("autodoc", "2");
iftImage *iftUncertMapToLabel(const iftImage *map);

%newobject iftLabelToUncertMap;
%feature("autodoc", "2");
iftImage *iftLabelToUncertMap(const iftImage *label, const iftAdjRel *A);

%newobject iftMaybeForeground;
%feature("autodoc", "2");
iftImage *iftMaybeForeground(const iftImage *label);

%newobject iftGrabCut;
%feature("autodoc", "2");
iftImage *iftGrabCut(const iftMImage *mimg, const iftImage *regions, double beta, int n_iters);

%newobject iftGMMDataSetDist;
%feature("autodoc", "2");
iftFImage *iftGMMDataSetDist(iftDataSet *train, const iftMImage *mimg, int n_comps);

