%include "iftMemory.i"
%include "iftFImage.i"
%include "iftMImage.i"
%include "iftRadiometric.i"

%newobject iftGraphCut;
%feature("autodoc", "2");
iftImage* iftGraphCut(const iftImage *gradient, iftLabeledSet *seeds, int beta);

%newobject iftGraphCutFromMImage;
%feature("autodoc", "2");
iftImage* iftGraphCutFromMImage(const iftMImage * mimg, iftLabeledSet* seeds, int beta);

%newobject iftGraphCutWithObjmap;
%feature("autodoc", "2");
iftImage* iftGraphCutWithObjmap(const iftMImage * mimg, const iftImage *objmap, iftLabeledSet* seeds, float alpha, int beta);

%newobject iftOptimumPathGraphCut;
%feature("autodoc", "2");
iftImage* iftOptimumPathGraphCut(const iftImage *gradient, iftLabeledSet* seeds, int beta);

%newobject iftGraphCutBeta;
%feature("autodoc", "2");
double iftGraphCutBeta(iftMImage *mimg);

