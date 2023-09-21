%include "iftImage.i"
%include "iftInterpolation.i"
%include "iftPlane.i"

%newobject iftAlignBrainByMSP;
%feature("autodoc", "2");
iftImage *iftAlignBrainByMSP(const iftImage *img, iftPlane **msp_out);

%newobject iftRotateImageToMSP;
%feature("autodoc", "2");
iftImage *iftRotateImageToMSP(const iftImage *img, const iftPlane *msp, iftInterpolationType interp_type);

%newobject iftMSPAsImage;
%feature("autodoc", "2");
iftImage *iftMSPAsImage(const iftPlane *msp, const iftImage *ref_img);

