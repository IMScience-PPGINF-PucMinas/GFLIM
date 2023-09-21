//
// Created by Samuel Martins on 06/12/18.
//

#ifndef IFT_ELASTIX_H
#define IFT_ELASTIX_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftImage.h"
#include "iftDataSet.h"



/**
 * @brief Register an (moving) image on to a fixed image by Elastix.
 *
 * If more than one elastix parameter files are passed, multiple registration are done, following the same order of
 * the files.
 *
 * @note If a mask for the moving image is passed, it defines the region on the moving image that will be registered.
 * @note Similarly, if a mask for the fixed image is passed, only its region is considered during registration.
 *
 * @param moving_img Moving Image to be registered.
 * @param fixed_img Fixed Image (Standard Space).
 * @param moving_mask Optional Mask to define the region on the moving image that will be registerd.
 * @param fixed_mask Optional Mask to define the region on the fixed image considered during registration.
 * @param elastix_files Elastix Registration Parameter Files.
 * @param def_fields_basename Optional string to define the basename from the resulting deformation fields.
 * @param def_fields Reference to return the resulting deformation fields from the registration. If NULL, the def. fields are deleted.
 * @return Registered image.
 *
 * @author Samuka Martins
 * @date Mar 13, 2018
 */
iftImage *iftRegisterImageByElastix(const iftImage *moving_img, const iftImage *fixed_img, const iftImage *moving_mask,
                                    const iftImage *fixed_mask, const iftFileSet *elastix_files,
                                    const char *def_fields_basename, iftFileSet **def_fields, const char output_file[]);


/**
 * @brief Finds and register a test image into the best normal (not abnormal) space.
 *
 * If more than one elastix parameter files are passed, multiple registrations are done, following the same order of
 * the files.
 *
 * @note If a mask for the moving image is passed, it defines the region on the moving image that will be registered.
 * @note Similarly, if the path for the fixed images masks is passed, only its region is considered during registration.
 *
 * @param moving_img Moving Image to be registered.
 * @param input_mask Moving Mask.
 * @param Z_normal_spaces Clusterized DataSet containing all normal images.
 * @param normal_mask_dir Dir containing.
 * @param transf_files Elastix Registration Parameter Files.
 * @param def_fields_basename Optional string to define the basename from the resulting deformation fields.
 * @param best_space_id Return the id of the sample related to the best normal reference space.
 * @param best_space_name File of the normal reference space.
 * @param registration_errors Return each registration error, according to the MSE error metric.
 * @return Registered image.
 *
 * @author Azael Sousa
 * @date Apr 1, 2021
 */
iftImage *iftRegisterImageByElastixIntoBestNormalSpace(iftImage *moving_img, iftImage *input_mask, iftDataSet *Z_normal_spaces,
                                                       const char *normal_mask_dir, iftFileSet *transf_files,
                                                       const char *def_fields_basename, int *best_space_id, char **best_space_name,
                                                       iftFloatArray **registration_errors);

/**
 * @brief Run the program iftRegisterSetImageByElastix.
 * @author Samuka
 * @date Jan 6, 2017
 * @ingroup Registration
 */
void iftRunProgRegisterImageSetByElastix(const char *moving_img_entry, const char *fixed_img_path, int img_depth,
                                         const char *affine_params_path, const char *bspline_params_path,
                                         const char *out_dir);

/**
 * @brief Transforms/Maps an Image by Transformix Program.
 * @author Samuel Martins
 * @date Dec 9th, 2015
 *
 * @param img Image to be Transformed/Mapped.
 * @param def_fields_path Pathname from the Tranformix Deformartion Fields File.
 * @return The Transformated/Mapped Image.
 */
iftImage *iftTransformImageByTransformix(const iftImage *img, const char *def_fields_path);


/**
 * @brief Run the program iftTransformImageSetByTransformix.
 * @author Samuka
 * @date Jan 9, 2017
 * @ingroup Registration
 */
void iftRunProgTransformImageSetByTransformix(const char *img_entry, const char *affine_def_fields_entry,
                                              const char *bspline_def_fields_entry, const char *out_dir);


/**
 * @brief Moves the Deformation Fields of <def_fields> to the new basename <out_basename>, fixing the new pathnames
 * inside the files.
 *
 * @warning The old deformation fields' files are removed from the disk.
 *
 * @param def_fields Deformation Fields.
 * @param out_basename New output basename.
 * @return File set with the deformation fields' pathnames in the new directory.
 *
 * @author Samuel Martins
 * @date Jan 3, 2019
 */
iftFileSet *iftMoveElastixDefFields(const iftFileSet *def_fields, const char *out_basename);

#ifdef __cplusplus
}
#endif

#endif
