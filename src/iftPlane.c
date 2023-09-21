#include "iftPlane.h"

#include "ift/core/dtypes/DblArray.h"
#include "ift/core/dtypes/Dict.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Json.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"


iftPlane *iftCreatePlane()
{
  iftPlane *pl=(iftPlane *)iftAlloc(1,sizeof(iftPlane));

  iftSetPlanePos(pl,0,0,0);
  iftSetPlaneOrient(pl,0.0,0.0,1.0);

  return(pl);
}

void iftDestroyPlane(iftPlane **pl) {
    if (pl) {
        iftPlane *aux = *pl;

        if (aux) {
            iftFree(aux);
            *pl = NULL;
        }
    }
}


iftPlane *iftReadPlane(const char *format, ...) {
    va_list args;
    char path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(path, format, args);
    va_end(args);
    
    if (!iftEndsWith(path, ".json"))
        iftError("Pathname %s is not a json file", "iftReadPlane", path);
    
    iftPlane *pl = iftCreatePlane();
    
    iftDict *json = iftReadJson(path);
    iftDblArray *normal = iftGetDblArrayFromDict("normal", json);
    pl->normal.x = normal->val[0];
    pl->normal.y = normal->val[1];
    pl->normal.z = normal->val[2];
    
    iftDblArray *ref_point = iftGetDblArrayFromDict("reference-point", json);
    pl->pos.x = ref_point->val[0];
    pl->pos.y = ref_point->val[1];
    pl->pos.z = ref_point->val[2];
    
    iftDestroyDict(&json);
    
    return pl;
}

void iftWritePlane(const iftPlane *pl, const char *format, ...) {
    if (pl) {
        va_list args;
        char path[IFT_STR_DEFAULT_SIZE];
    
        va_start(args, format);
        vsprintf(path, format, args);
        va_end(args);
        
        if (!iftEndsWith(path, ".json"))
            iftError("Pathname %s is not a json file", "iftWritePlane", path);
        
        char *parent_dir = iftParentDir(path);
        if (!iftDirExists(parent_dir))
            iftMakeDir(parent_dir);
        iftFree(parent_dir);
        
        iftDict *json = iftCreateDict();
        
        iftDblArray *normal = iftCreateDblArray(3);
        normal->val[0] = pl->normal.x;
        normal->val[1] = pl->normal.y;
        normal->val[2] = pl->normal.z;
        iftInsertIntoDict("normal", normal, json);
        
        iftDblArray *ref_point = iftCreateDblArray(3);
        ref_point->val[0] = pl->pos.x;
        ref_point->val[1] = pl->pos.y;
        ref_point->val[2] = pl->pos.z;
        iftInsertIntoDict("reference-point", ref_point, json);
        
        iftWriteJson(json, path);
        
        iftDestroyDict(&json); // it destroys the normal and ref_points
    }
}


void iftTranslatePlane(iftPlane *pl, float Tx, float Ty, float Tz){
  pl->pos.x += Tx;
  pl->pos.y += Ty;
  pl->pos.z += Tz;
}

void iftRotatePlane(iftPlane *pl, char axis, float theta)
{
  iftMatrix *R=iftRotationMatrix(axis,theta);
  iftVector N;
  
  N = iftTransformVector(R,pl->normal);
  iftDestroyMatrix(&R);

  iftSetPlaneOrient(pl,N.x,N.y,N.z);
}

void iftSetPlanePos(iftPlane *pl, float x, float y, float z)
{
  pl->pos.x = x;
  pl->pos.y = y;
  pl->pos.z = z;
}

void iftSetPlaneOrient(iftPlane *pl, float Nx, float Ny, float Nz)
{
  float m = sqrtf(Nx*Nx + Ny*Ny + Nz*Nz); 

  if (m==0.0)
      iftError("There is no orientation", "iftSetPlaneOrient");

  pl->normal.x = Nx/m;
  pl->normal.y = Ny/m;
  pl->normal.z = Nz/m;
}


iftPlane *iftDefinePlaneFromPoints(iftPoint P1, iftPoint P2, iftPoint P3)
{
  iftPlane *pl=NULL;
  iftVector a, b;

  if (!iftCollinearPoints(P1,P2,P3)){
    pl = iftCreatePlane();
    pl->pos.x = (P1.x+P2.x+P3.x)/3.0;
    pl->pos.y = (P1.y+P2.y+P3.y)/3.0;
    pl->pos.z = (P1.z+P2.z+P3.z)/3.0;
    a.x       = P2.x - P1.x;
    a.y       = P2.y - P1.y;
    a.z       = P2.z - P1.z;
    b.x       = P3.x - P1.x;
    b.y       = P3.y - P1.y;
    b.z       = P3.z - P1.z;
    pl->normal = (iftVector)iftVectorCrossProd(a,b);
  }
  else
  {
     iftWarning("Collinear points were found.", "iftDefinePlaneFromPoints");
  }
  return(pl);
}

iftPlane *iftDefinePlaneFromVoxels(iftVoxel P1, iftVoxel P2, iftVoxel P3)
{
  iftPlane *pl=NULL;
  iftVector a, b;

  if (!iftCollinearPoints(P1,P2,P3)){
    pl = iftCreatePlane();    
    pl->pos.x = (P1.x+P2.x+P3.x)/3.0;
    pl->pos.y = (P1.y+P2.y+P3.y)/3.0;
    pl->pos.z = (P1.z+P2.z+P3.z)/3.0;
    a.x       = P2.x - P1.x;
    a.y       = P2.y - P1.y;
    a.z       = P2.z - P1.z;
    b.x       = P3.x - P1.x;
    b.y       = P3.y - P1.y;
    b.z       = P3.z - P1.z;
    pl->normal = (iftVector)iftVectorCrossProd(a,b);
  }

  return(pl);
}

char iftPlaneSide(iftPlane *pl, iftPoint P)
{
  iftVector a;
  float iprod;

  a.x = (P.x - pl->pos.x); 
  a.y = (P.y - pl->pos.y); 
  a.z = (P.z - pl->pos.z); 

  iprod = iftVectorInnerProduct(a, pl->normal);

  if (fabs(iprod) < IFT_EPSILON) 
    return(0); // point belongs to the plane
  else{
    if (iprod > 0.0)
      return(1); // point is in the positive side
    else
      return(-1); // point is in the negative side
  }
}

float iftDistPlanePoint(iftPlane *pl, iftPoint P)
{
  iftVector a;
  float     dist;

  a.x = (P.x - pl->pos.x); 
  a.y = (P.y - pl->pos.y); 
  a.z = (P.z - pl->pos.z); 

  dist = iftVectorInnerProduct(a, pl->normal);

  return(dist);
}

char iftPlaneSideVoxel(iftPlane *pl, iftVoxel P)
{
  iftVector a;
  float iprod;

  a.x = (P.x - pl->pos.x); 
  a.y = (P.y - pl->pos.y); 
  a.z = (P.z - pl->pos.z); 

  iprod = iftVectorInnerProduct(a, pl->normal);

  if (fabs(iprod) < IFT_EPSILON) 
    return(0); // point belongs to the plane
  else{
    if (iprod > 0.0)
      return(1); // point is in the positive side
    else
      return(-1); // point is in the negative side
  }
}

float iftDistPlaneVoxel(iftPlane *pl, iftVoxel P)
{
  iftVector a;
  float     dist;

  a.x = (P.x - pl->pos.x); 
  a.y = (P.y - pl->pos.y); 
  a.z = (P.z - pl->pos.z); 

  dist = iftVectorInnerProduct(a, pl->normal);

  return(dist);
}



iftPlane *iftGetPlaneFrom3Points(iftPoint A, iftPoint B, iftPoint C, bool normalize) {
    /* A = (x0, y0, z0), B = (x1, y1, z1), C = (x2, y2, z2)
     * AB = (x1 - x0, y1 - y0, z1 - z0) = (AB.x, AB.y, AB.z)
     * AC = (x2 - x0, y2 - y0, z2 - z0) = (AC.x, AC.y, AC.z)
     *
     * norm to the plane = (a, b, c) = AB x AC (cross product)
     *        |  a     b     c   |  a     b   |
     * norm = | AB.x  AB.y  AB.z | AB.x  AB.y | = 0
     *        | AC.x  AC.y  AC.z | AC.x  AC.y |
     *
     * a = (AB.y * AC.z) - (AB.z * AC.y)
     * b = (AB.z * AC.x) - (AB.x * AC.z)
     * c = (AB.x * AC.y) - (AB.y * AC.x)
     *
     * Given a point P = (x, y, z) on the plane
     * plane equation: norm x AP (dot product) = norm * (x-x0, y-y0, z-z0)
     */
    
    float x0 = A.x, y0 = A.y, z0 = A.z;
    float x1 = B.x, y1 = B.y, z1 = B.z;
    float x2 = C.x, y2 = C.y, z2 = C.z;
    
    iftVector AB = {x1 - x0, y1 - y0, z1 - z0};
    iftVector AC = {x2 - x0, y2 - y0, z2 - z0};

    iftPlane *pl = iftCreatePlane();
    pl->normal.x = (AB.y * AC.z) - (AB.z * AC.y);
    pl->normal.y = (AB.z * AC.x) - (AB.x * AC.z);
    pl->normal.z = (AB.x * AC.y) - (AB.y * AC.x);

    if (normalize)
        pl->normal = iftNormalizeVector(pl->normal);
    
    pl->pos.x = x0;
    pl->pos.y = y0;
    pl->pos.z = z0;
    
    return pl;
}


float iftGetBiasFromPlane(const iftPlane *pl) {
    /* P0 = pl.pos = (x0, y0, z0)
     * norm to the plane = pl->normal = (a, b, c)
     * 
     * general plane equation
     * ax + by + cz + d = 0
     * 
     * bias d = - (ax0 + by0 + cz0)
     */
    return (- ((pl->normal.x * pl->pos.x) + (pl->normal.y * pl->pos.y) + (pl->normal.z * pl->pos.z)));
}


iftPoint iftGetPointFromYZPlaneProjection(const iftPlane *pl, float y, float z) {
    /* 
     * norm to the plane = pl->normal = (a, b, c)
     * 
     * general plane equation
     * ax + by + cz + d = 0
     * 
     * Find x from P = (x, y, z), given y and z
     * x = -(by + cz + d) / a
     */
    float d = iftGetBiasFromPlane(pl);
    float x = - ((pl->normal.y * y) + (pl->normal.z * z) + d) / pl->normal.x;

    return (iftPoint) {.x = x, .y = y, .z = z};
}
