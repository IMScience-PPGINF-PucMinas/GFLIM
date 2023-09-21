%include "iftCommon.i"
%include "iftAdjacency.i"
%include "iftImage.i"
%include "iftRadiometric.i"
%include "iftFImage.i"
%include "iftSeeds.i"
%include "iftSegmentation.i"

%newobject iftGeodesicDistTrans;
%feature("autodoc", "2");
iftFImage *iftGeodesicDistTrans(const iftSet *S, const iftImage *mask, const iftAdjRel *A);

%newobject iftBorderDistTrans;
%feature("autodoc", "2");
iftImage  *iftBorderDistTrans(const iftImage *label, iftAdjRel *A);

%newobject iftSetDistTrans;
%feature("autodoc", "2");
  iftImage  *iftSetDistTrans(iftSet *S, int xsize, int ysize, int zsize);

%newobject iftMSSkel;
%feature("autodoc", "2");
iftFImage      *iftMSSkel(iftImage *bin);

%newobject iftMSSkel2D;
%feature("autodoc", "2");
iftFImage *iftMSSkel2D(iftImage *label_img, iftAdjRel *A, iftSide side, iftImage **dist_out,
                       iftImage **relabel_img_out);

%newobject iftEuclDistTrans;
%feature("autodoc", "2");
iftImage *iftEuclDistTrans(const iftImage *label_img, const iftAdjRel *Ain, iftSide side, iftImage **root_out,
                           iftImage **edt_label_out, iftImage **pred_out);

%newobject iftGeodesicDistTransFromSet;
%feature("autodoc", "2");
iftImage *iftGeodesicDistTransFromSet(const iftImage *mask, const iftSet *set, iftImage **root_out);

%newobject iftMedialAxisTrans2D;
%feature("autodoc", "2");
iftImage       *iftMedialAxisTrans2D(iftImage *bin, float scale_thres, iftSide side);

%newobject iftShapeReconstruction;
%feature("autodoc", "2");
iftImage       *iftShapeReconstruction(iftImage *medial_axis, int value);

%newobject iftTerminalPoints2D;
%feature("autodoc", "2");
iftImage       *iftTerminalPoints2D(iftImage *skel);

%newobject iftBranchPoints2D;
%feature("autodoc", "2");
iftImage *iftBranchPoints2D(iftImage *skel);

%newobject iftTerminalPointSet2D;
%feature("autodoc", "2");
iftSet         *iftTerminalPointSet2D(iftImage *skel);

%newobject iftContourToArray;
%feature("autodoc", "2");
iftVoxelArray *iftContourToArray(const iftImage *contour);

%newobject iftApproxContour;
%feature("autodoc", "2");
iftVoxelArray *iftApproxContour(const iftImage *contour, double epsilon);

%newobject iftApproxVoxelArray;
%feature("autodoc", "2");
iftVoxelArray *iftApproxVoxelArray(const iftVoxelArray *array, double epsilon);

%newobject iftVoxelArrayToSet;
%feature("autodoc", "2");
iftSet *iftVoxelArrayToSet(const iftImage *image, const iftVoxelArray *array);

%newobject iftNearestInContour;
%feature("autodoc", "2");
iftSet *iftNearestInContour(const iftImage *contour, const iftSet *set, const iftImage *mask);

