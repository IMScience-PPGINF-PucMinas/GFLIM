#ifndef IFT_ADJACENCY_H_
#define IFT_ADJACENCY_H_

/**
 * @file iftAdjacency.h
 * @brief Structs and function prototypes for adjacency relation.
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * 
 */

/** @addtogroup Image
 * @{ */

#ifdef __cplusplus
extern "C" {
#endif

#include "iftCommon.h"


/**
 * @brief Defines a neighborhood around a voxel.
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 *
 */
//! swig(extend = iftAdjRelExt.i, destroyer = iftDestroyAdjRel)
typedef struct ift_adjrel {
    int *dx, *dy, *dz, *dt;
    /* displacements to achieve the n adjacent voxels. */
    int n; /* number of adjacent voxels. */
} iftAdjRel;

/**
 * @brief Defines a neighborhood around a voxel, with a faster Iteration, 
 * Disconsidering the image's border
 * 
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @warning Can only be used for images with same size. See also ::iftAdjRel.
 *
 */
typedef struct ift_fastadjrel {
    int n;
    /* number of adjacent voxels */
    int *dq;
    /* displacements to reach adjacent voxels for a given image */
    int bx, by, bz; /* sizes of the image's border to be disconsidered */
} iftFastAdjRel;

/**
 * @brief Allocates memory for a 3D adjacency relation
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 *
 */
//! swig(newobject)
iftAdjRel *iftCreateAdjRel(int n); /*  */


/**
 * @brief Deallocates memory for a adjacency relation object.
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 *
 */
void iftDestroyAdjRel(iftAdjRel **A);

/**
 * @brief Creates a 3D ball of radius <r> as adjacency relation.
 *
 * The displacement vectors are sorted by their magnitude, similarly to the function iftSpheric.
 * However, there is no sense in the order, for example, clockwise, etc.
 * 
 * @param  r Sphere radius
 * @return   Spheric Adjacency relation of radius <r>
 *
 * @author Samuel Martins
 * @date May 2, 2018
 * @ingroup Adjacency
 */
//! swig(newobject)
iftAdjRel *iftSpheric(float r);

/**
 * @brief Creates a 4D ball of radius <r> as adjacency relation.
 *
 * The displacement vectors are sorted by their magnitude, similarly to the function iftSpheric.
 * However, there is no sense in the order, for example, clockwise, etc.
 * 
 * @param  r Sphere radius
 * @return   Spheric Adjacency relation of radius <r>
 *
 * @author Ilan Silva
 * @date Feb 11, 2022
 * @ingroup Adjacency
 */
//! swig(newobject)
iftAdjRel *iftHyperSpheric(float r);

/**
 * @brief Creates a 4D semisphere in the t-axis of radius <r> as adjacency relation.
 *
 * The displacement vectors are sorted by their magnitude, similarly to the function iftSpheric.
 * However, there is no sense in the order, for example, clockwise, etc.
 *
 * @param  r Sphere radius
 * @param  positive signal of the displacements in the t-axis
 * @return   SemiHyperSpheric Adjacency relation of radius <r>
 *
 * @author Ilan Silva
 * @date Mar 23, 2022
 * @ingroup Adjacency
 */
//! swig(newobject)
iftAdjRel *iftSemiHyperSpheric(float r, bool positive);

/**
 * @brief Creates a 3D half-ball of radius @a r as adjacency relation, in the corresponding axis and direction.
 * This adjacency is useful for segmenting a volume in a single direction, e.g., a video-volume where z axis is time.
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @param r Sphere radius.
 * @param axis Axis to be considered.
 * @param direction The direction to follow in @a axis. (-1, 1)
 * @return The hemispheric adjacency.
 * */
//! swig(newobject)
iftAdjRel *iftHemispheric(float r, char axis, int direction);


/**
 * @brief Creates a 3D ball surface of radius r as adjacency relation.
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @warning This function is used for 3D images.
 */
//! swig(newobject)
iftAdjRel *iftSphericEdges(float r);

/**
 * @brief Creates a 2D ball of radius <r> on the xy plane as adjacency relation.
 *
 * The displacement vectors are sorted by their magnitude, similarly to the function iftSpheric.
 * However, there is no sense in the order, for example, clockwise, etc.
 * 
 * @param  r Sphere radius.
 * @return   Spheric Adjacency relation of radius <r>.
 *
 * @author Samuel Martins
 * @date May 2, 2018
 * @ingroup Adjacency
 */
//! swig(newobject, stable)
iftAdjRel *iftCircular(float r);


/**
 * @brief Creates a circunference of radius r on the xy plane as adjacency relation
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @warning This function is used for 2D images.
 * @param r Circle radius.
 * */
//! swig(newobject)
iftAdjRel *iftCircularEdges(float r);

/**
 * @brief Creates a 2D ball of radius @a r on the xy plane as adjacency relation for contour pixel labeling.
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @param r Circle radius.
 **/
//! swig(newobject)
iftAdjRel *iftClockCircular(float r);

/**
 * @brief Creates a 2D adjacency relation with only 
 * the right-hand side of the input adjacency relation
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @param r Circle radius.
 **/
//! swig(newobject)
iftAdjRel *iftRightSide(iftAdjRel *A, float r);

/**
 * @brief Creates a 2D adjacency relation with only 
 * the left-hand side of the input adjacency relation
 *
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @param r Circle radius.
 **/
//! swig(newobject)
iftAdjRel *iftLeftSide(iftAdjRel *A, float r);

/**
 * @author Alexandre Falcao
 * @date Jun 10, 2016
 * @ingroup Adjacency
 * @brief Creates a rectangle adjacency with specified dimensions.
 * @warning This function is used for 2D images.
 * @param xsize Rectangle width.
 * @param ysize Rectangle height.
 * @return The rectangular adjacency.
 */
//! swig(newobject)
iftAdjRel *iftRectangular(int xsize, int ysize);

/**
 * @author Alexandre Falcao
 * @date Nov 18, 2020
 * @ingroup Adjacency
 * @brief Creates a rectangle adjacency with specified dimensions and dilation rates.
 * @warning This function is used for 2D images.
 * @param xsize: Rectangle width.
 * @param ysize: Rectangle height.
 * @param sx > 1: dilation rate along x.
 * @param sy > 1: dilation rate along y.
 * @return Rectangular adjacency.
 */
//! swig(newobject)
iftAdjRel *iftRectangularWithDilation(int xsize, int ysize, int sx, int sy);

  
/**
 * @brief Creates a cross adjacency with specified dimensions.
 * @warning This function is used for 2D images.
 * @param xsize Cross width.
 * @param ysize Cross height.
 * @return The cross adjacency.
 */
//! swig(newobject)
iftAdjRel *iftCross(int xsize, int ysize);
  
/**
 * @brief Creates a Cuboid adjacency with specified dimensions.
 *
 * @warning This function is used for 3D images.
 *
 * @param xsize Rectangle width.
 * @param ysize Rectangle height.
 * @param zsize Rectangle depth.
 * @return The cuboid adjacency.
 */
//! swig(newobject)
iftAdjRel *iftCuboid(int xsize, int ysize, int zsize);

/**
 * @brief Creates a HyperCuboid adjacency with specified dimensions.
 *
 * @warning This function is used for 4D images.
 *
 * @param xsize Rectangle width.
 * @param ysize Rectangle height.
 * @param zsize Rectangle depth.
 * @param tsize Rectangle depth.
 * @return The hyper cuboid adjacency.
 */
//! swig(newobject)
iftAdjRel *iftHyperCuboid(int xsize, int ysize, int zsize, int tsize);

/**
 * @author Azael Sousa
 * @date Jan 09, 2021
 * @ingroup Adjacency
 * @brief Creates a rectangle adjacency with specified dimensions and dilation rates.
 * @warning This function is used for 2D images.
 * @param xsize: Rectangle width.
 * @param ysize: Rectangle height.
 * @param sx > 1: dilation rate along x.
 * @param sy > 1: dilation rate along y.
 * @return Rectangular adjacency.
 */
//! swig(newobject)
iftAdjRel *iftCuboidWithDilation(int xsize, int ysize, int zsize, int sx, int sy, int sz);


/**
 * @brief Creates a new copy of the given adjacency.
 * @param A Adjacency to be copied.
 * @return Copy of adjacency @a A.
 */
//! swig(newobject)
iftAdjRel *iftCopyAdjacency(const iftAdjRel *A);

iftFastAdjRel *iftCreateFastAdjRel(iftAdjRel *A, int *tby, int *tbz);

/* create an adjacency relation to speed up implementations for a given image by computing the displacements to the adjaceny voxels based on the look-up tables tby and tbz of the image. The fast implementation must disconsider the image's border */

void iftDestroyFastAdjRel(iftFastAdjRel **F);

//! swig()
void iftMaxAdjShifts(const iftAdjRel *A, int *dx, int *dy, int *dz);

void iftWriteAdjRel(iftAdjRel *A, char *filename);
void iftPrintAdjRel(iftAdjRel *A);

iftAdjRel *iftReadAdjRel(char *filename);

/** Read an Adjacency Relation from a Binary File */
iftAdjRel *iftReadAdjRelBinFile(const char *path);

/** Writes an Adjacency Relation from a Binary File */
void iftWriteAdjRelBinFile(const iftAdjRel *A, const char *path);

//! swig(newobject)
iftVoxel iftGetAdjacentVoxel(const iftAdjRel *A, iftVoxel u, int adj);

/**
 * @brief Check if the adjacency is suitable for 3D images.
 * @param A The checked adjacency.
 * @return 1 if it is a 3D adjacency, 0 otherwise.
 */
int iftIsAdjRel3D(iftAdjRel *A);


/**
 * @brief Given an adjacency relation <A>, preferably created with iftCircular or iftSpheric, this function extracts
 * all displacements whose magnitude is within the closed interval [radius_min, radius_max] into a
 * new adjacency relation. Hence, this function may be used to create a shell (ring) adjacency relation. NOTE THAT: this
 * function generalizes the concept of ring adjacency to any kind of adjacency relations other than iftCircular and
 * iftSpheric. We could replace it by a simpler functiond called iftCircularRing (2D) and iftSphericRing (3D), but we
 * implemented this one instead for the sake of generalization. Moreover, by filtering an existing adjacency relation
 * we ensure that <A> may be broken into several non-overlapping rings, if we call iftExtractAdjacencyShell with
 * different non-overlapping intervals.
 *
 * @author Thiago Vallin Spina
 * @date Feb 10, 2016
 *
 * @param A The input adjacency relation (preferably a circular or spherical one).
 * @param radius_min The minimum length that the displacement must have to be part of the new adjacency relation.
 * @param radius_max The maximum length that the displacement must have to be part of the new adjacency relation.
 * @param include_center This parameter forces that displacement (0,0,0) to be added in the new adjacency relation, since
 * it may not be contained within [radius_min, radius_max]
 *
 * @return The filtered adjacency relation with displacement (0,0,0) in the first position, if required and/or part of
 * the interval.
 */
iftAdjRel *iftExtractAdjacencyShell(iftAdjRel *A, double radius_min, double radius_max, uchar include_center);


/**
 * @brief Find the Boundaries (extreme displacements) of an Adjacency Relation <A>, return them into
 * a new Adjacency Relation.
 *
 * To find such boundaries (extreme displacements), we create an "image" and to fill <A> inside it,
 * with value 1 for the voxels representing the displacements of A.
 * We then use a second adjacency relation B to find the borders of this "object". If it is NULL,
 * 4-neighborhood (2D) or 6-neighborhood (3D) is used. 
 * There could be some cases in that other adjacencies, such that 8-neighborhood and 26-neighborhood,
 * are better.
 * 
 * @param  A Adjacency whose boundaries (extremities) are found.
 * @param  B Adjacency relation used to 
 * @return   Resulting adjacency relation only with the extreme displacements of A.
 * 
 * @author Samuel Martins
 * @date Apr 27, 2018
 */
iftAdjRel *iftAdjacencyBoundaries(const iftAdjRel *A, const iftAdjRel *B);


/**
 * @brief This was created to define an order invariant to direction for contour extraction,
 *        the idea is similar to always start looking from the same side when going through an
 *        labyrinth to avoid cycles.
 * @details
 * 6-Neighbor ClockCircular sorting:
 *
 *      2 | 3 | 4
 *      1 | 0 | 5
 *      8 | 7 | 6
 *
 * For the current pixel "q" and predecessor "p":
 *
 *        |   |
 *        | p |
 *        |   | q
 *
 * The right side of "p" to "q" is the index 1.
 *
 * @param  pred     Voxel of predecessor of "cur"
 * @param  cur      Current voxel
 * @return  Starting index of ClockCircular with radii sqrtf(2.0) of right side pixel
 *          of vector point from "pred" to "cur"
 *
 * @author Jordão Bragantini
 * @date Aug 12, 2018
 */
int ift6NeighborsClockFromLeft(iftVoxel *pred, iftVoxel *cur);

void iftSidePixels(const iftVoxel *u, const iftVoxel *v, float radii, iftVoxel *left, iftVoxel *right);

void iftSidePixelsWithJitter(const iftVoxel *u, const iftVoxel *v, float radii,
                             int mod, iftVoxel *left, iftVoxel *right);


#ifdef __cplusplus
}
#endif

/** @} */

#endif
