%include "iftAdjacency.i"
%include "iftClassification.i"
%include "iftClustering.i"
%include "iftCommon.i"
%include "iftDataSet.i"
%include "iftFiltering.i"
%include "iftFImage.i"
%include "iftImage.i"
%include "iftImageForest.i"
%include "iftImageMath.i"
%include "iftMathMorph.i"
%include "iftMetrics.i"
%include "iftMImage.i"
%include "iftRadiometric.i"
%include "iftRepresentation.i"
%include "iftSeeds.i"

%newobject iftEnhanceEdges;
%feature("autodoc", "2");
iftImage      *iftEnhanceEdges(iftImage *img, iftAdjRel *A, iftLabeledSet *seed, float alpha);

%newobject iftWatershed;
%feature("autodoc", "2");
iftImage *iftWatershed(const iftImage *basins, iftAdjRel *Ain, iftLabeledSet *seeds, iftSet *forbidden);

%newobject iftWaterCut;
%feature("autodoc", "2");
iftImage *iftWaterCut(iftMImage *mimg, iftAdjRel *Ain, iftLabeledSet *seeds, iftSet *forbidden);

%newobject iftWaterGray;
%feature("autodoc", "2");
iftImage  *iftWaterGray(iftImage *basins, iftImage *marker, iftAdjRel *A);

%feature("autodoc", "2");
void iftWaterGrayForest(iftImageForest *fst, iftImage *marker);

%newobject iftOrientedWatershed;
%feature("autodoc", "2");
iftImage *iftOrientedWatershed(const iftImage *img, const iftImage *grad, iftAdjRel *Ain,
                               iftFloatArray *alpha, iftLabeledSet *seeds, iftSet *forbidden);

%feature("autodoc", "2");
int iftOtsu(const iftImage *img);

%newobject iftThreshold;
%feature("autodoc", "2");
iftImage *iftThreshold(const iftImage *img, int lowest, int highest, int value);

%newobject iftAboveAdaptiveThreshold;
%feature("autodoc", "2");
  iftImage  *iftAboveAdaptiveThreshold(iftImage *img, iftImage *mask, iftAdjRel *A, float perc, int niters, int value);

%newobject iftBelowAdaptiveThreshold;
%feature("autodoc", "2");
  iftImage  *iftBelowAdaptiveThreshold(iftImage *img, iftImage *mask, iftAdjRel *A, float perc, int niters, int value);

%newobject iftBorderImage;
%feature("autodoc", "2");
iftImage *iftBorderImage(const iftImage *label, bool get_margins);

%newobject iftSelectAndPropagateRegionsAboveAreaByColor;
%feature("autodoc", "2");
iftImage *iftSelectAndPropagateRegionsAboveAreaByColor(iftImage *img, iftImage *label, int area);

%newobject iftSmoothRegionsByDiffusion;
%feature("autodoc", "2");
iftImage *iftSmoothRegionsByDiffusion(const iftImage *label_img, const iftImage *orig_img,
                                      float smooth_factor, int n_iters);

%newobject iftOrientedWaterCut;
%feature("autodoc", "2");
iftImage *iftOrientedWaterCut(const iftImage *img, iftAdjRel *Ain, iftFloatArray *alpha,
                                          iftLabeledSet *seeds, iftSet *forbidden);

%newobject iftOrientedColorWaterCut;
%feature("autodoc", "2");
iftImage *iftOrientedColorWaterCut(iftMImage *mimg, iftImage *orient, iftAdjRel *Ain, float beta, iftLabeledSet *seeds);

%newobject iftEnhancedWaterCut;
%feature("autodoc", "2");
iftImage *iftEnhancedWaterCut(iftMImage *mimg, iftImage *objmap, iftAdjRel *Ain, iftLabeledSet *seeds, float alpha);

%newobject iftEnhancedWatershed;
%feature("autodoc", "2");
iftImage *iftEnhancedWatershed(iftImage *basins, iftImage *objmap, iftAdjRel *Ain, iftLabeledSet *seeds, float alpha);

%newobject iftSuperPixelMajorityVote;
%feature("autodoc", "2");
iftImage *iftSuperPixelMajorityVote(iftImage *comp, iftImage *objmap, float threshold);

%newobject iftBoundingBoxArrayToLabel;
%feature("autodoc", "2");
iftImage *iftBoundingBoxArrayToLabel(const iftImage *img, const iftBoundingBoxArray *bb_ary);

%newobject iftFindOptPath;
%feature("autodoc", "2");
iftSet *iftFindOptPath(const iftMImage *mimg, int src, int dst, float sigma, iftFImage *pathval, iftImage *pred);

%newobject iftFindOptPathWithProb;
%feature("autodoc", "2");
iftSet *iftFindOptPathWithProb(const iftMImage *mimg, const iftFImage *obj, const iftFImage *bkg, int src, int dst,
                               float sigma, float gamma, iftFImage *pathval, iftImage *pred);

%newobject iftSinglePathContour;
%feature("autodoc", "2");
iftImage *iftSinglePathContour(const iftImage *label);

%newobject iftFillContour;
%feature("autodoc", "2");
iftImage *iftFillContour(const iftImage *contour);

%newobject iftDISF;
%feature("autodoc", "2");
iftImage *iftDISF(iftMImage *mimg, iftAdjRel *A, int num_init_seeds, int num_superpixels, iftImage *mask);

