%include "iftCommon.i"
%include "iftImage.i"
%include "iftFImage.i"

%newobject iftAdd;
%feature("autodoc", "2");
iftImage *iftAdd(const iftImage *img1, const iftImage *img2);

%newobject iftSub;
%feature("autodoc", "2");
iftImage *iftSub(const iftImage *img1, const iftImage *img2);

%newobject iftSubReLU;
%feature("autodoc", "2");
iftImage *iftSubReLU(const iftImage *img1, const iftImage *img2);

%newobject iftAbsSub;
%feature("autodoc", "2");
iftImage *iftAbsSub(const iftImage *img1, const iftImage *img2);

%newobject iftIntersec;
%feature("autodoc", "2");
iftImage *iftIntersec(const iftImage *label_img1, const iftImage *label_img2);

%newobject iftAnd;
%feature("autodoc", "2");
iftImage *iftAnd(const iftImage *img1, const iftImage *img2);

%newobject iftOr;
%feature("autodoc", "2");
iftImage *iftOr(const iftImage *img1, const iftImage *img2);

%newobject iftMult;
%feature("autodoc", "2");
iftImage *iftMult(const iftImage *img1, const iftImage *img2);

%newobject iftMultByIntScalar;
%feature("autodoc", "2");
iftImage *iftMultByIntScalar(const iftImage *img, int scalar);

%newobject iftAbs;
%feature("autodoc", "2");
iftImage *iftAbs(const iftImage *img);

%newobject iftXor;
%feature("autodoc", "2");
iftImage *iftXor(const iftImage *img1, const iftImage *img2);

%newobject iftComplement;
%feature("autodoc", "2");
iftImage *iftComplement(const iftImage *img);

%newobject iftReLUImage;
%feature("autodoc", "2");
iftImage *iftReLUImage(const iftImage *img);

%newobject iftMask;
%feature("autodoc", "2");
iftImage *iftMask(const iftImage *img, const iftImage *mask);

%newobject iftBinarize;
%feature("autodoc", "2");
iftImage *iftBinarize(const iftImage *label_img);

%newobject iftAddValue;
%feature("autodoc", "2");
iftImage *iftAddValue(const iftImage *img, int val);

%newobject iftBinaryFrom1To255;
%feature("autodoc", "2");
iftImage *iftBinaryFrom1To255(const iftImage *bin);

%feature("autodoc", "2");
float iftMeanValue(const iftImage *img, iftImage *mask);

%newobject iftRoundFImage;
%feature("autodoc", "2");
iftImage *iftRoundFImage(const iftFImage *fimg);

