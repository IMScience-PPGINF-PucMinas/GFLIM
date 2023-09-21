%include "iftCommon.i"
%include "iftImage.i"
%include "iftFImage.i"
%include "iftMImage.i"
%include "iftDataSet.i"
%include "iftAdjacency.i"
%include "iftPlane.i"
%include "iftMatrix.i"
%include "iftSegmentation.i"
%include "iftRepresentation.i"
%include "iftMSPS.i"



typedef enum {
  IFT_NEAREST_NEIGHBOR_INTERPOLATION, IFT_LINEAR_INTERPOLATION
} iftInterpolationType;

%newobject iftInterpByNearestNeighbor;
%feature("autodoc", "2");
iftImage *iftInterpByNearestNeighbor(const iftImage *img, float sx, float sy, float sz);

%newobject iftInterp2DByNearestNeighbor;
%feature("autodoc", "2");
iftImage *iftInterp2DByNearestNeighbor(const iftImage *img, float sx, float sy);

%newobject iftInterp;
%feature("autodoc", "2");
iftImage  *iftInterp(const iftImage *img, float sx, float sy, float sz);

%newobject iftInterp2D;
%feature("autodoc", "2");
iftImage *iftInterp2D(const iftImage *img, float sx, float sy);

%newobject iftResizeImage;
%feature("autodoc", "2");
iftImage* iftResizeImage(const iftImage *img, int xsize, int ysize, int zsize);

