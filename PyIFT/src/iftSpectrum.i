%include "iftImage.i"
%include "iftMImage.i"

%newobject iftGaborFilter2D;
%feature("autodoc", "2");

iftMImage *iftGaborFilter2D(iftMImage *spec, float u0, float v0, float P, float x0, float y0, float a, float b, float theta);

