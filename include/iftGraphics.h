#ifndef IFT_GRAPHICS_H_
#define IFT_GRAPHICS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftAdjacency.h"
#include "iftCommon.h"
#include "iftCurve.h"
#include "iftDataSet.h"
#include "iftFImage.h"
#include "iftImage.h"
#include "iftImageMath.h"
#include "iftInterpolation.h"
#include "iftMathMorph.h"
#include "iftMImage.h"
#include "iftPlane.h"
#include "iftRepresentation.h"
#include "iftSeeds.h"

#include "ift/core/dtypes/Color.h"
#include "ift/imgproc/dtypes/BoundingBox.h"
#include "ift/imgproc/dtypes/Roi.h"


#define  NUM_OF_NORMALS    65161
#define  SCENE_NORMAL      0
#define  OBJECT_NORMAL     1
#define  RAYCASTING        0
#define  SPLATTING         1

/* Graphical context with the visualization attributes for ray casting
   and voxel splatting, including current viewing direction
   parameters.

       Octants of the Scene
 z
  /
 /
 ---- x
 |
y|
       4 ------- 5
      /        /|
     /        / |
    /        /  |
  0 -------- 1  |
   |  6     |  / 7
   |        | /
   |        |/
  2 -------- 3


The origin of the image to be rendered on the viewing plane is
initially at (0,0,-diag/2). The observer is initially at (diag/2,
diag/2, -infinity) with a single white light source at the same
position. The transformation T is applied to the scene for voxel
splatting or the inverse of T is applied to the vieweing plane for ray
casting. The front-to-back (FTB) voxel access order for splatting
depends on the pair (closest octant, principal axis). The FTB
possibilities (dx,dy,dz,xo,yo,zo,xf,yf,zf) are stored in a 2D array
whose column is the closest octant and line is the principal axis. Any
change in viewing direction requires to update all viewing parameters below.

*/

/* Extra surface rendering buffers, especially required for voxel
   splatting. They all have the same dimensions of the resulting
   rendered image (diag x diag pixels) to avoid clipping. */

typedef struct ift_surface_rendering_buffers {
  float  depth;         /* distance to the viewing plane, whose center
			   is initially positioned at
			   (diag/2,diag/2,-diag/2) */
  float  opacity;       /* accumulated opacity of the voxels projected
			   onto the viewing plane */
  int    voxel;         /* voxels in the scene that have been projected
			   onto the viewing plane */
  int    object;        /* object in the scene whose voxel has been
			   projected onto the viewing plane */
} iftSRBuffers;

/* Object opacity, visibility, and color attributes in the scene */

typedef struct ift_object_attributs {
  float opacity;          /* opacity value in [0,1] */
  float red, green, blue; /* proportions in [0,1] of red, green and
			     blue components in its color */
  char visibility;        /* visibility flag (0/1 for invisible/visible) */
} iftObjectAttributes;

/* Parameters of the Phong's illumination model */

typedef struct ift_phong_model {
  float      ka;     /* ambient reflection constant  (0.1 by default) */
  float      kd;     /* diffuse reflection constant  (0.7 by default) */
  float      ks;     /* specular refrection constant (0.2 by default) */
  float      ns;     /* exponent of the specular reflection component (5 by default) */
  iftVector *normal; /* normal table used to speed up shading */
  float     *Idist;  /* distance table used to convert voxel's distance to the viewing plane into intensity in [0,1] (also known as depth shading) */
  int        ndists; /* number of distance values, which is the size of the scene's diagonal. */
} iftPhongModel;

/* Viewing direction and FTB voxel access order */

typedef struct ift_ftb_access {
  int dx, dy, dz; /* -1 or 1. Default is 1,1,1 */
  int xo, yo, zo; /* the x,y,z original coordinates of the closest  octant: the initial FTB voxel access coordinates */
  int xf, yf, zf; /* the x,y,z original coordinates of the farthest octant: the final FTB voxel access coordinates */
} iftFTBaccess;

typedef struct ift_viewing_direction {
  iftMatrix           *T;            /* scene's transformation matrix: default is spin=tilt=0.0 */
  iftMatrix           *Tinv;         /* inverse of the scene's transformation matrix */
  iftMatrix           *R;            /* rotation matrix, only for vectors */
  iftMatrix           *Rinv;         /* inverse of the rotation matrix, only for vectors */
  char                 paxis;        /* principal axis, 0 for 'x', 1 for 'y', or 2 for 'z', which is the last one to be visited in the voxel splatting loop. This is the most orthogonal axis to the viewing plane. */
  int                  octant;       /* closest octant, 0, 1, ..., 7, to the vieweing plane, which together with the principal axis indicate the front-to-back voxel access order for splatting. */
  iftFTBaccess        *ftb;          /* gives the front-to-back access order from the closest octant to the farthest octant. Used for voxel splatting only. */
  iftVector            V;            /* visualization direction --- vector that points to the observer. */
} iftViewDir;


/* Graphical context of the scene for 3D visualization */

typedef struct ift_graphic_context {
  iftObjectAttributes *object;       /* list of attributes per object
					in the scene */
  int                  nobjs;        /* number of objects in the scene */
  float                overall_opac; /* it indicates the overall opacity of the objects in the scene */
  char                 proj_mode;    /* voxel projection mode: ray casting is 0 and voxel splatting is 1 */
  iftPhongModel       *phong;        /* rendering parameters of the Phong's illumination model */
  iftViewDir          *viewdir;      /* viewing direction parameters */
  iftSRBuffers        *surf_render;  /* extra surface rendering buffers, especially required for voxel splatting */
  iftFImage           *scene;        /* the 3D image property for visualization (intensity, BMD, etc.) */
  iftImage            *label;        /* the 3D label image with objects for visualization, used only for surface rendering */
  iftImage            *border;        /* the 3D labeled borders of the objects in label, used only for surface rendering */
  iftPlane             face[6];      /* planes of the 6 faces of the scene  */
  iftImage            *normal;       /* normal indices of the scene's voxels */
  iftFImage           *opacity;      /* opacity scene used for volume rendering only */
} iftGraphicalContext;

/*
   Graphical Context for 3D visualization of the image properties in
   scene. The object labels i > 0 (0 is background) are used only for
   surface rendering (i.e., label==NULL for volume rendering). The
   normal vectors may be estimated as the gradient of the scene
   (normal_type=SCENE_NORMAL), as the gradient of the EDT of the
   objects (normal_type=OBJECT_NORMAL) or on-the-fly from the index
   buffer (normal_type=NIL).

*/

void                 iftSetSceneNormal(iftGraphicalContext *gc);
void                 iftSetProjectionMode(iftGraphicalContext *gc, char proj_mode);
void                 iftSetObjectNormal(iftGraphicalContext *gc, char normal_type);
void                 iftSetObjectColor(iftGraphicalContext *gc, int object, float red, float green, float blue);
void                 iftSetSceneOpacity(iftGraphicalContext *gc, float min_val, float max_val, iftImage *grad, int grad_thres, float max_opac);
void                 iftSetObjectOpacity(iftGraphicalContext *gc, int object, float opacity);
void                 iftSetObjectVisibility(iftGraphicalContext *gc, int object, char visibility);
void                 iftSetViewDir(iftGraphicalContext *gc, float tilt, float spin);
void                 iftSetAmbientReflection(iftGraphicalContext *gc, float ka);
void                 iftSetDiffuseReflection(iftGraphicalContext *gc, float kd);
void                 iftSetSpecularReflection(iftGraphicalContext *gc, float ks);
char                 iftIntersectionPoints(iftGraphicalContext *gc, iftPoint P0, iftVector n, iftPoint *P1, iftPoint *Pn);
iftImage            *iftSurfaceRender(iftGraphicalContext *gc);
iftImage            *iftVolumeRender(iftGraphicalContext *gc);
iftGraphicalContext *iftCreateGraphicalContext(iftFImage *scene, iftImage *label);
void                 iftDestroyGraphicalContext(iftGraphicalContext *gc);

iftImage  *iftGraphicImage(iftImage *img); /* creates image for overlaying drawings */
iftImage  *iftDrawPath(iftImage *img, iftImage *pred, int last, iftFColor normRGB, iftAdjRel *B); /* Draw a path in the predecessor map with terminus at point last */
void       iftDrawPoint(iftImage *img, iftVoxel u, iftColor YCbCr, iftAdjRel *B, int rangeValue);
void       iftDrawPointAlpha(iftImage *img, iftVoxel u, iftColor YCbCr, iftAdjRel *B, int rangeValue, float alpha);
void       iftDrawPoints(iftImage *img, iftSet *S, iftColor YCbCr, iftAdjRel *B);
void       iftDrawCurve(iftImage *img, iftCurve *curve, iftColor YCbCr, iftAdjRel *B);
void       iftDrawLine(iftImage *img, iftVoxel u1, iftVoxel u2, iftColor YCbCr, iftAdjRel *B);
void       iftDrawBorders(iftImage *img, iftImage *label, iftAdjRel *A, iftColor YCbCr, iftAdjRel *B);


/**
 * @brief Draws the objects of a label image (borders and interiors), with the colors defined in cmap, in an input image.
 * 
 * If fill is true, the interior of the object is drawn in <img>.
 * An input adjacency is required to define and draw the objects' borders in the image. If it is NULL,
 * the border is not drawn.
 * The object border has a higher color than the one from its interior.
 * 
 * @note For a label image with n objects (value in the range 0 (bg) to n), a color map of n + 1 colors
 * is required, so that the first valid color must be at index [1].
 * 
 * @param img Input image where the objects are drawn (in-place).
 * @param label_img Label image with n objects.
 * @param cmap Color map with n + 1 colors, one for each object and bg.
 * @param A Adjacency relation to define and draw the object borders. If NULL, no borders are drawn.
 * @param fill Fill the objects' interiors.
 * @param alpha Alpha factor to draw the colors on the image.
 * 
 * @author Samuel Martins (only a few changes)
 * @date Oct 15, 2018
 */
//! swig(stable)
void iftDrawLabels(iftImage *img, const iftImage *label_img, const iftColorTable *cmap,
                   const iftAdjRel *A, bool fill, float alpha);


//! swig(stable)
void       iftDrawObject(iftImage *img, iftImage *bin, iftColor YCbCr, iftAdjRel *B);
void       iftDrawRoots(iftImage *img, iftImage *root, iftColor YCbCr, iftAdjRel *B);

//! swig(newobject)
iftImage *iftDraw2DFeatureSpace(iftDataSet *Z, iftFeatureLabel opt, iftSampleStatus status);

iftImage *iftDrawVoxelSamples(iftDataSet *Z, iftFeatureLabel opt, char *ref_data_type);
//! swig(newobject)
iftImage *iftColorizeComp(iftImage *label);
//! swig(newobject)
iftImage *iftColorizeImageByLabels(iftImage *orig, iftImage *label, iftColorTable *ctb);
//! swig(newobject)
iftImage *iftColorizeCompOverImage(iftImage *orig, iftImage *label);
iftImage  *iftProjectMaxValue(iftImage *img,  iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMinValue(iftImage *img,  iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMeanValue(iftImage *img, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMaxObjectValue(iftImage *img,  iftImage *obj, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMinObjectValue(iftImage *img,  iftImage *obj, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftProjectMeanObjectValue(iftImage *img, iftImage *obj, iftPlane *cutplane, int xviewsize, int yviewsize, int slabsize);
iftImage  *iftCurvilinearSplatting(iftImage *img, iftFImage *dist, iftPlane *cutplane, int viewsize, float depth);
void       iftTiltSpinFromPrincipalAxis(iftImage *bin, float *tilt, float *spin);

/**
 * @brief Draws a line from px1 to pxn recursively using color YCbCr
 * for a given adjacency relation. The main advantage of this function
 * over iftDrawLine is that the line will be "perfectly" rasterized,
 * guaranteeing that endpoints px1 and pxn are part of it and that
 * there will be no "holes" on the line.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param img colored Image.
 * @param px1 start voxel.
 * @param pxn end voxel.
 * @param YCbCr YCbCr color.
 * @param B adjacency used to draw the line.
*/
void iftDrawLineRecursive(iftImage *img, iftVoxel px1, iftVoxel pxn, iftColor YCbCr, iftAdjRel *B);

/**
 * @brief Draws the lines of a 2D polygon given by a set of points and then fills it to form a binary mask.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param points The set of vertices.
 * @param n      The number of points.
 * @param xsize  The X image dimensions.
 * @param ysize  The Y image dimensions.
 * @return The X tile coordinate.
 */
iftImage *iftMaskImageFromPolygon2D(iftVoxel *points, int n, int xsize, int ysize);

/**
 * @brief Draws the borders of a given Bounding Box in a pre-allocated Image.
 * The border will have value/color <b>YCbCr</b>.
 * @author Samuel Martins
 * @date May 30, 16
 * 
 * @param img               Pre-allocated Image where the borders will be drawn.
 * @param bb                Bounding box.
 * @param YCbCr             Value/Color used to draw the borders.
 * @param A                 Adjacency to determine the width of the bounding box
 * @param drawCentralPoint  Whether or not to draw the central point of the bounding box.
 */
//! swig(stable)  
void iftDrawBoundingBoxBordersInPlace(iftImage *img, iftBoundingBox bb, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint);

/**
 * @brief Draws the borders of a given set of Bounding Boxes in a pre-allocated Image.
 * The border will have value/color <b>YCbCr</b>.
 * @author Cesar Castelo
 * @date Jul 10, 2018
 * 
 * @param img               Pre-allocated Image where the borders will be drawn.
 * @param bbs               Bounding box array.
 * @param YCbCr             Value/Color used to draw the borders.
 * @param drawCentralPoint  Whether or not to draw the central point of the bounding box.
 */
void iftDrawBoundingBoxArrayBordersInPlace(iftImage *img, iftBoundingBoxArray *bbs, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint);

/**
 * @brief Draws the borders of a given ROI in a pre-allocated Image.
 * The border will have value/color <b>YCbCr</b>.
 * @author Cesar Castelo
 * @date Mar 22, 2019
 * 
 * @param img               Pre-allocated Image where the borders will be drawn.
 * @param roi               Region of interest.
 * @param YCbCr             Value/Color used to draw the borders.
 * @param A                 Adjacency to determine the width of the border
 * @param drawCentralPoint  Whether or not to draw the central point of the ROI.
 */
void iftDrawRoiBordersInPlace(iftImage *img, iftRoi *roi, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint);

/**
 * @brief Draws the borders of a given set of ROIs in a pre-allocated Image.
 * The border will have value/color <b>YCbCr</b>.
 * @author Cesar Castelo
 * @date Mar 22, 2019
 * 
 * @param img               Pre-allocated Image where the borders will be drawn.
 * @param roiArray          Array of ROIs.
 * @param YCbCr             Value/Color used to draw the borders.
 * @param drawCentralPoint  Whether or not to draw the central point of each ROI.
 */
void iftDrawRoiArrayBordersInPlace(iftImage *img, iftRoiArray *roiArray, iftColor YCbCr, iftAdjRel *A, bool drawCentralPoint);

/**
 * @brief Draws binary seeds on pre-allocated image
 * @author Jordao
 * @date Oct 23, 17
 *
 * @param img   Pre-allocated Image where the seeds will be drawn.
 * @param seeds Binary seeds
 * @param YCbCr Value/Color used to draw the borders.
 * @param A     Radius of the drawing
 */
void iftDrawBinaryLabeledSeeds(iftImage *img, iftLabeledSet *seeds, iftColor YCbCr, iftAdjRel *A);

/**
 * @brief Draws a legend in an image containing a box with small boxes inside
 * (one per each color in the color table from firstColor to lastColor)
 * @author Cesar Castelo
 * @date Set 18, 2018
 *
 * @param img           Pre-allocated Image where the legend will be drawn.
 * @param ctb           Color table containing the colors to be drawn
 * @param firstColor    Index of the first color to be used from the color table
 * @param lastColor     Index of the last color to be used from the color table
 * @param location      Location for the legend: top-left, top-right, bottom-left, bottom-right
 * @param colorBoxSize  Size of each color box (in pixels)
 */
void iftDrawColorLegend(iftImage *img, iftColorTable *ctb, int firstColor, int lastColor, char *location, int colorBoxSize);

/**
 * @brief Draws a fancier legend in an image containing a white box with small boxes inside
 * (one per each color in the color table from firstColor to lastColor)
 * @author Cesar Castelo
 * @date Nov 16, 2018
 *
 * @param img           Pre-allocated Image where the legend will be drawn.
 * @param ctb           Color table containing the colors to be drawn
 * @param firstColor    Index of the first color to be used from the color table
 * @param lastColor     Index of the last color to be used from the color table
 * @param location      Location for the legend: top-left, top-right, bottom-left, bottom-right
 * @param colorBoxSize  Size of each color box (in pixels)
 */
void iftDrawFancyColorLegend(iftImage *img, iftColorTable *ctb, int firstColor, int lastColor, char *location, int colorBoxSize);


#ifdef __cplusplus
}
#endif

#endif

