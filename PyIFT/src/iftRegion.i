%include "iftImage.i"
%include "iftCommon.i"
%include "iftImage.i"
%include "iftAdjacency.i"
%include "iftDataSet.i"
%include "iftMathMorph.i"

%newobject iftRelabelRegions;
%feature("autodoc", "2");
iftImage* iftRelabelRegions(iftImage* labelled, iftAdjRel* adj_rel);

%newobject iftGeodesicCenters;
%feature("autodoc", "2");
iftLabeledSet* iftGeodesicCenters(const iftImage* label_image);

%newobject iftJoinProbabilityRegions;
%feature("autodoc", "2");
iftImage *iftJoinProbabilityRegions(iftMImage *prob, iftImage *label, int norm_val, bool decrement);

