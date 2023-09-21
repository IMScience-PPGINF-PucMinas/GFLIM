%include "iftCommon.i"
%include "iftImage.i"
%include "iftMatrix.i"
%include "iftInterpolation.i"
%include "iftDataSet.i"

%newobject iftRotateImage2D;
%feature("autodoc", "2");
iftImage *iftRotateImage2D(iftImage *img, float theta);

%newobject iftAlignObject;
%feature("autodoc", "2");
iftImage *iftAlignObject(iftImage *bin);

%newobject iftRotateImage;
%feature("autodoc", "2");
iftImage *iftRotateImage(const iftImage *img, float theta_x, float theta_y);

%newobject iftScaleImage;
%feature("autodoc", "2");
iftImage *iftScaleImage(const iftImage *img, float sx, float sy, float sz);

%newobject iftFlipImage;
%feature("autodoc", "2");
iftImage *iftFlipImage(const iftImage *img, iftAxis axis);

