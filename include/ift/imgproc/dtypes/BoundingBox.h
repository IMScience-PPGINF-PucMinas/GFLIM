//
// Created by samuel on 17/01/19.
//

#ifndef IFT_BOUNDING_BOX_H
#define IFT_BOUNDING_BOX_H

#ifdef __cplusplus
extern "C" {
#endif


#include "ift/core/dtypes/BasicDataTypes.h"


/**
 * @brief Image Bounding Box.
 * @author Samuel Martins
 * @date Mar 2, 2016
 * @ingroup DataTypes
 */
//! swig()
typedef struct ift_bounding_box {
    /** The leftmost and topmost point of the Bounding Box */
    iftVoxel begin;
    /** The rightmost and bottommost point of the Bounding Box */
    iftVoxel end;
} iftBoundingBox;


/**
 * @brief Image Bounding Box Array (with no label).
 * @author Cesar Castelo
 * @date Jan 8, 2018
 * @ingroup DataTypes
 */
//! swig(destroyer = iftDestroyBoundingBoxArray)
typedef struct ift_bounding_box_array {
    iftBoundingBox *val;
    long n;
} iftBoundingBoxArray;


/**
 * @brief Label Image Bounding Box.
 * @author Samuka Martins
 * @date Sep 15, 2017
 * @ingroup DataTypes
 */
typedef struct ift_label_bounding_box {
    /** The leftmost and topmost point of the Bounding Box */
    iftVoxel begin;
    /** The rightmost and bottommost point of the Bounding Box */
    iftVoxel end;
    /** Label from the Bounding Box */
    int label;
} iftLabelBoundingBox;


/**
 * @brief Creates an array of Bounding Boxes.
 * @author Cesar Castelo
 * @date Jan 8, 2018
 */
iftBoundingBoxArray *iftCreateBoundingBoxArray(long n);


/**
 * @brief Destroys an array of Bounding Boxes.
 * @author Cesar Castelo
 * @date Jan 8, 2018
 */
void iftDestroyBoundingBoxArray(iftBoundingBoxArray **bbs);

/**
 * @brief Read a CSV file of bounding boxes.
 *
 * Where each row in the CSV file contains a bounding box (6 coordinates) following:
 * begin.x,begin.y,begin.z,end.x,end.y,end.z
 *
 * @author Samuel Martins
 * @date Apr 27, 2018
 */
iftBoundingBoxArray *iftReadBoundingBoxArray(const char *csv_path, ...);

/**
 * @brief Write an array of bounding boxes as a CSV file.
 *
 * Each bounding box (6 coordinates) is store as a row in csv file, following:
 * begin.x,begin.y,begin.z,end.x,end.y,end.z
 *
 * @author Samuel Martins
 * @date Apr 27, 2018
 */
void iftWriteBoundingBoxArray(const iftBoundingBoxArray *bbs, const char *csv_path, ...);


/**
 * @brief Gets the central voxel from a bounding box.
 * @author Samuka
 * @date Jan 5, 2017
 */
iftVoxel iftBoundingBoxCenterVoxel(iftBoundingBox bb);

/**
 * @brief Scale a Bounding Box by a factor of alpha. To enlage the BB, use alpha > 1.0.
 * @param  bb    Bounding Box to be scaled.
 * @param  alpha Scale factor.
 * @return       Scaled Bounding Box.
 *
 * @author Samuka Martins
 * @date Aug 31, 2017
 */
iftBoundingBox iftScaleBoundingBox(iftBoundingBox bb, float alpha);

/**
 * @brief Count the number of voxels inside a bounding box
 * @author Samuel Martins
 * @date Aug 31, 2017
 */
//! swig()
long iftBoundingBoxBoxVolume(iftBoundingBox bb);


/**
 * @brief Print the coordinates from a Bounding Box
 * @author Samuka Martins
 * @date Sep 18. 2017
 */
//! swig()
void iftPrintBoundingBox(iftBoundingBox bb);


/**
 * @brief Centralize a given bounding box to the voxel <new_center>
 * @author Samuka Martins
 * @date Sep 18, 2017
 */
//! swig(newobject)
iftBoundingBox iftCentralizeBoundingBox(iftBoundingBox bb, iftVoxel new_center);


/**
 * @brief Return the size of a bounding box.
 * @author Samuka Martins
 * @date Sep 18, 2017
 */
iftImageDomain iftBoundingBoxSize(iftBoundingBox bb);


/**
 * @brief Find out the minimum bounding box that fits a pair of input bounding boxes.
 * @param  bb1 Bounding Box 1
 * @param  bb2 Bounding Box 2
 * @return     Min Bounding Box that fits the input bounding boxes
 *
 * @author Samuel Martins
 * @date Nov 14, 2017
 */
iftBoundingBox iftMinBoundingBoxOfBoundingBoxes(iftBoundingBox bb1, iftBoundingBox bb2);


/**
 * @brief Gets the min. bounding box that fits all voxels from an array.
 * @param  voxels Array of voxels.
 * @param  n      Size of the Array.
 * @return        Resulting Min. Bouning Box.
 *
 * @author Samuel Martins
 * @date Nov 14, 2017
 */
iftBoundingBox iftMinBoundingBoxOfVoxels(const iftVoxel *voxels, size_t n);


/**
 * @brief Read the Coordinates of a bounding box from a csv file: x0,y0,z0,x1,y1,z1
 * @author Samuka Martins
 * @date Jan 16, 2018
 */
iftBoundingBox iftReadBoundingBox(const char *format, ...);


/**
 * @brief Write the Coordinates of a bounding box into a csv file: x0,y0,z0,x1,y1,z1
 * @author Samuka Martins
 * @date Jan 16, 2018
 */
void iftWriteBoundingBox(iftBoundingBox bb, const char *filename, ...);


/**
 * @brief Fit a bounding box to a given image domain. It is already on such domain, it returns a
 * copy of the input bounding box.
 *
 * @param  bb     Input Bounding Box.
 * @param  domain Image domain where the bounding box is fit.
 * @return        Fit bounding box. If the input bounding box was already on the image domina, a
 *                copy is returned.
 *
 * @author Samuka Martins
 * @date Apr 25, 2018.
 */
iftBoundingBox iftFitBoundingBoxOnImageDomain(iftBoundingBox bb, iftImageDomain domain);



#ifdef __cplusplus
}
#endif

#endif //IFT_BOUNDING_BOX_H
