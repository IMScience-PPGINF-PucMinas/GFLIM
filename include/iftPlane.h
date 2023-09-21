#ifndef IFT_PLANE_H_
#define IFT_PLANE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"
#include "iftMatrix.h"


/**
 * @brief Struct that defines a Plane.
 * @author Unknown, Samuel Martins (last update)
 * @date Nov 21, 2018 
 */
typedef struct ift_plane {
    /** reference point from the plane */
    iftPoint pos;
    /** normal vector which gives its orientation */
    iftVector normal;
} iftPlane;

/**
 * @brief Creates/allocates a plane.
 * @author Samuel Martins (last update)
 * @date Nov 21, 2018 
 */
iftPlane *iftCreatePlane(void);
void iftDestroyPlane(iftPlane **pl);


/**
 * @brief Reads a plane.
 *
 * @param pathname Json pathname to of the plane.
 * @param ... Optional var args to build the pathname dinamically (like printf).
 *
 * @author Samuel Martins (last update)
 * @date Nov 22, 2018
 */
iftPlane *iftReadPlane(const char *pathname, ...);


/**
 * @brief Writes a plane to a json file <pathname>.
 *
 * @param pl Plane to be written.
 * @param pathname Json pathname to write the plane.
 * @param ... Optional var args to build the pathname dinamically (like printf).
 *
 * @author Samuel Martins (last update)
 * @date Nov 22, 2018
 */
void iftWritePlane(const iftPlane *pl, const char *pathname, ...);


/**
 * @brief Copies a plane.
 * @author Samuel Martins
 * @date Nov 25, 2018
 */
inline iftPlane *iftCopyPlane(const iftPlane *pl) {
    if (pl) {
        iftPlane *copy = iftCreatePlane();
        copy->normal = pl->normal;
        copy->pos = pl->pos;
        return copy;
    }
    else return NULL;
}


void iftTranslatePlane(iftPlane *pl, float Tx, float Ty, float Tz);
void iftRotatePlane(iftPlane *pl, char axis, float theta);

/**
 * @brief Scales a given plane.
 * @param pl Plane to be scaled.
 * @param scaling_factor Scaling factor. 0.5 means scaling of 50%.
 * @return Scaled plane.
 *
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
inline iftPlane *iftScalePlane(const iftPlane *pl, float scaling_factor) {
    iftPlane *pl_scaled = iftCreatePlane();
    pl_scaled->normal = pl->normal;
    pl_scaled->pos.x = pl->pos.x * scaling_factor;
    pl_scaled->pos.y = pl->pos.y * scaling_factor;
    pl_scaled->pos.z = pl->pos.z * scaling_factor;
    
    return pl_scaled;
}

void iftSetPlanePos(iftPlane *pl, float x, float y, float z);
void iftSetPlaneOrient(iftPlane *pl, float Nx, float Ny, float Nz);
void iftSetViewupOrient(iftPlane *pl, float vx, float vy, float vz);
iftPlane *iftDefinePlaneFromPoints(iftPoint P1, iftPoint P2, iftPoint P3);
iftPlane *iftDefinePlaneFromVoxels(iftVoxel P1, iftVoxel P2, iftVoxel P3);
char iftPlaneSide(iftPlane *pl, iftPoint P);
char iftPlaneSideVoxel(iftPlane *pl, iftVoxel P);
float iftDistPlanePoint(iftPlane *pl, iftPoint P);
float iftDistPlaneVoxel(iftPlane *pl, iftVoxel P);


/**
 * @brief Gets a plane from 3 non-collinear points (points not on a single line).
 *
 * @attention There is no checking for the collinearity of the points.
 *
 * @param A First point.
 * @param B Second point.
 * @param C Third point.
 * @param normalize Consider the unit vector of the norm.
 * @return Resulting Plane.
 *
 * @author Samuel Martins
 * @date Nov 21, 2018
 */
iftPlane *iftGetPlaneFrom3Points(iftPoint A, iftPoint B, iftPoint C, bool normalize);


/**
 * @brief Computes the bias from the plane based on its normal vector and reference point.
 * @author Samuel Martins
 * @date Nov 21, 2018
 */
float iftGetBiasFromPlane(const iftPlane *pl);


/**
 * @brief Figures out the x-coordinate of the point on the plane with given y- and z-coordinates.
 * It returns the points (x, y, z).
 * 
 * @param pl Plane.
 * @param y Y-coordinate of the point on the plane.
 * @param z Z-coordinate of the point on the plane.
 * @return Resulting point on the plane.
 * 
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
iftPoint iftGetPointFromYZPlaneProjection(const iftPlane *pl, float y, float z);


#ifdef __cplusplus
}
#endif

#endif
