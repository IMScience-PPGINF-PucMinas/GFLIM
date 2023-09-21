/**
 * @file iftDicom.h
 * @brief Structs and function prototypes for dicom manipulation.
 * 
 * An example of dicom manipulation is found in demo/Miscelaneuous/iftDicom2Scene.c
 * 
 * @author Samuel Martins (sbmmartins)
 */
#ifndef _IFT_DICOM_H_
#define _IFT_DICOM_H_


#ifdef __cplusplus
extern "C" {
#endif


#include "iftCommon.h"
#include "iftImage.h"
#include "iftMatrix.h"
#include "iftRadiometric.h"

typedef struct _ift_dicom_slice {
    char filename[512];
    float z_position; /* Position from the slice in the z-axis */
} iftDicomSlice;


typedef struct _ift_dicom_image {
    int xsize, ysize, zsize; /* all dicom slices have the same xsize and ysize.
                                zsize = number of slices. */
    float dx, dy, dz; /* all dicom slices have the same dx and dy */
    int bits_per_voxel; /* all dicom slices have the same bits_per_voxel (image depth) */
    char *orientation; /* all dicom slices have the same orientation */
    iftDicomSlice **slices; /* vector of iftDicomSlices* */
    float orientation_vectors[2][3];    /* all dicom slices have the same orientation vectors
				       (this variable is a matrix in which the first line 
				       is the first vector and the second line is the 
                                       second vector, both with 3 elements). 
				       These vectors indicate the direction of the x 
				       and y axes as read from the image, which may be 
				       different from the usual (1, 0, 0), (0, 1, 0).*/
} iftDicomImage;


iftDicomSlice *iftCreateDicomSlice();
void iftDestroyDicomSlice(iftDicomSlice **dicom_slice);

iftDicomImage *iftCreateDicomImage(int nslices);
iftDicomImage *iftReadDicomImage(char *dir_pathname);

void iftDestroyDicomImage(iftDicomImage **dicom);

void iftPrintDicomInfo(iftDicomImage *dicom);

iftImage *iftConvertDicomCoordinates(iftImage *img, iftDicomImage *dicom);

iftImage *iftConvertDicom2Scene(iftDicomImage *dicom);


#ifdef __cplusplus
}
#endif

#endif
