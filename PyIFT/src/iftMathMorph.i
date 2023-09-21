%include "iftCommon.i"
%include "iftImage.i"
%include "iftImageMath.i"
%include "iftSeeds.i"
%include "iftAdjacency.i"
%include "iftFiltering.i"
%include "iftCompTree.i"

%newobject iftDilate;
%feature("autodoc", "2");
  iftImage *iftDilate(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftErode;
%feature("autodoc", "2");
  iftImage *iftErode(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftClose;
%feature("autodoc", "2");
  iftImage *iftClose(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftOpen;
%feature("autodoc", "2");
  iftImage *iftOpen(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftAsfCO;
%feature("autodoc", "2");
  iftImage *iftAsfCO(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftAsfOC;
%feature("autodoc", "2");
  iftImage *iftAsfOC(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftDilateBin;
%feature("autodoc", "2");
  iftImage *iftDilateBin(const iftImage *bin, iftSet **seed, float radius);

%newobject iftErodeBin;
%feature("autodoc", "2");
iftImage *iftErodeBin(const iftImage *bin, iftSet **seed, float radius);

%newobject iftErodeLabelImage;
%feature("autodoc", "2");
iftImage *iftErodeLabelImage(const iftImage *label_img, float radius);

%newobject iftCloseBin;
%feature("autodoc", "2");
iftImage *iftCloseBin(const iftImage *bin, float radius);

%newobject iftOpenBin;
%feature("autodoc", "2");
iftImage *iftOpenBin(const iftImage *bin, float radius);

%newobject iftAsfCOBin;
%feature("autodoc", "2");
  iftImage *iftAsfCOBin(const iftImage *bin, float radius);

%newobject iftAsfOCBin;
%feature("autodoc", "2");
  iftImage *iftAsfOCBin(const iftImage *bin, float radius);

%newobject iftCloseRecBin;
%feature("autodoc", "2");
  iftImage *iftCloseRecBin(const iftImage *bin, float radius);

%newobject iftOpenRecBin;
%feature("autodoc", "2");
  iftImage *iftOpenRecBin(const iftImage *bin, float radius);

%newobject iftMorphGrad;
%feature("autodoc", "2");
  iftImage *iftMorphGrad(const iftImage *img, const iftAdjRel *A);

%newobject iftSuperiorRec;
%feature("autodoc", "2");
  iftImage *iftSuperiorRec(const iftImage *img, const iftImage *marker, const iftImage *mask);

%newobject iftInferiorRec;
%feature("autodoc", "2");
  iftImage *iftInferiorRec(const iftImage *img, const iftImage *marker, const iftImage *mask);

%newobject iftOpenRec;
%feature("autodoc", "2");
  iftImage *iftOpenRec(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftCloseRec;
%feature("autodoc", "2");
  iftImage *iftCloseRec(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftAsfCORec;
%feature("autodoc", "2");
  iftImage *iftAsfCORec(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftAsfOCRec;
%feature("autodoc", "2");
  iftImage *iftAsfOCRec(const iftImage *img, const iftAdjRel *A, const iftImage *mask);

%newobject iftCloseBasins;
%feature("autodoc", "2");
  iftImage *iftCloseBasins(const iftImage *img, iftSet *seed, const iftImage *mask);

%newobject iftOpenDomes;
%feature("autodoc", "2");
  iftImage *iftOpenDomes(const iftImage *img, iftSet *seed, const iftImage *mask);

%newobject iftHClose;
%feature("autodoc", "2");
  iftImage *iftHClose(const iftImage *img, int H, const iftImage *mask);

%newobject iftHOpen;
%feature("autodoc", "2");
  iftImage *iftHOpen(const iftImage *img, int H, const iftImage *mask);

%newobject iftFastAreaClose;
%feature("autodoc", "2");
  iftImage *iftFastAreaClose(const iftImage *img, int thres);

%newobject iftFastAreaOpen;
%feature("autodoc", "2");
iftImage *iftFastAreaOpen(const iftImage *img, int thres);

%newobject iftAreaClose;
%feature("autodoc", "2");
  iftImage *iftAreaClose(const iftImage *img, int area_thres, iftCompTree *ctree);

%newobject iftVolumeClose;
%feature("autodoc", "2");
  iftImage *iftVolumeClose(const iftImage *img, int volume_thres, iftCompTree *ctree);

%newobject iftAreaOpen;
%feature("autodoc", "2");
  iftImage *iftAreaOpen(const iftImage *img, int area_thres, iftCompTree *ctree);

%newobject iftVolumeOpen;
%feature("autodoc", "2");
  iftImage *iftVolumeOpen(const iftImage *img, int volume_thres, iftCompTree *ctree);

