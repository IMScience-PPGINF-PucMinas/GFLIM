%include "iftCommon.i"
%include "iftImage.i"
%include "iftFImage.i"

%newobject iftLinearStretch;
%feature("autodoc", "2");
iftImage *iftLinearStretch(iftImage *img, double f1, double f2, double g1, double g2);

%newobject iftNormalize;
%feature("autodoc", "2");
iftImage *iftNormalize(const iftImage *img, double minval, double maxval);

%newobject iftNormalizeWithNoOutliersInRegion;
%feature("autodoc", "2");
iftImage *iftNormalizeWithNoOutliersInRegion(const iftImage *img, const iftImage *region,
                                            int minval, int maxval, float perc);

%newobject iftNormalizeWithNoOutliers;
%feature("autodoc", "2");
iftImage *iftNormalizeWithNoOutliers(const iftImage *img, int minval, int maxval, float perc);

%newobject iftMatchHistogram;
%feature("autodoc", "2");
iftImage *iftMatchHistogram(const iftImage *img, const iftImage *img_mask,
                            const iftImage *ref, const iftImage *ref_mask);

