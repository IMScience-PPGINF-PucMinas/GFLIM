%include "iftAdjacency.i"
%include "iftCommon.i"
%include "iftCurve.i"
%include "iftDataSet.i"
%include "iftFImage.i"
%include "iftImage.i"
%include "iftImageMath.i"
%include "iftInterpolation.i"
%include "iftMathMorph.i"
%include "iftMImage.i"
%include "iftPlane.i"
%include "iftRepresentation.i"
%include "iftSeeds.i"

%feature("autodoc", "2");
void iftDrawLabels(iftImage *img, const iftImage *label_img, const iftColorTable *cmap,
                   const iftAdjRel *A, bool fill, float alpha);

%feature("autodoc", "2");
void       iftDrawObject(iftImage *img, iftImage *bin, iftColor YCbCr, iftAdjRel *B);

%newobject iftDraw2DFeatureSpace;
%feature("autodoc", "2");
iftImage *iftDraw2DFeatureSpace(iftDataSet *Z, iftFeatureLabel opt, iftSampleStatus status);

%newobject iftColorizeComp;
%feature("autodoc", "2");
iftImage *iftColorizeComp(iftImage *label);

%newobject iftColorizeImageByLabels;
%feature("autodoc", "2");
iftImage *iftColorizeImageByLabels(iftImage *orig, iftImage *label, iftColorTable *ctb);

%newobject iftColorizeCompOverImage;
%feature("autodoc", "2");
iftImage *iftColorizeCompOverImage(iftImage *orig, iftImage *label);

%feature("autodoc", "2");
void iftDrawBoundingBoxBordersInPlace(iftImage *img, iftBoundingBox bb, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint);

