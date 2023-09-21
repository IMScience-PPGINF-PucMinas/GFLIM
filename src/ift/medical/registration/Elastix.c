#include "ift/medical/registration/Elastix.h"

#include "ift/core/dtypes/IntArray.h"
#include "ift/core/io/Dir.h"
#include "ift/core/tools/OS.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"
#include "iftImageMath.h"
#include "iftRegistration.h"



/********************** PRIVATE FUNCTIONS *************************/
iftFileSet *_iftAdjustElastixParameterFiles(const iftFileSet *elastix_files, const iftImage *moving_img,
                                            const iftImage *fixed_img, const char *tmp_dir) {
    iftFileSet *elastix_files_final = iftCreateFileSet(elastix_files->n);
    
    for (int i = 0; i < elastix_files->n; i++) {
        char *filename = iftFilename(elastix_files->files[i]->path, NULL);
        char *out_transf_path = iftJoinPathnames(2, tmp_dir, filename);
        iftCopyFromDisk(elastix_files->files[i]->path, out_transf_path);
        
        elastix_files_final->files[i] = iftCreateFile(out_transf_path);
        
        // changing image dimensions from 3D to 2D
        if (!iftIs3DImage(moving_img))
            iftSed(out_transf_path, "(MovingImageDimension 3)", "(MovingImageDimension 2)");
        if (!iftIs3DImage(fixed_img))
            iftSed(out_transf_path, "(FixedImageDimension 3)", "(FixedImageDimension 2)");
        
        iftFree(filename);
        iftFree(out_transf_path);
    }
    
    return elastix_files_final;
}


/**
 * @brief Runs the Elastix program.
 * @author Samuel
 * @November 17, 2015
 * @ingroup Elastix
 */
iftImage *_iftRunElastix(const char *moving_path, const char *fixed_path, const char *moving_mask_path,
                        const char *fixed_mask_path, const iftFileSet *transf_files, const char *out_dir, const char output_file[]) {
    char args[4096];

    // Registration: Assuming elastix is already in your OS PATH
    // just turn on the BSpline Registration    
    sprintf(args, "-m %s -f %s", moving_path, fixed_path);
    for (size_t i = 0; i < transf_files->n; i++)
        sprintf(args, "%s -p %s", args, transf_files->files[i]->path);
    sprintf(args, "%s -out %s", args, out_dir);

    if (moving_mask_path != NULL)
        sprintf(args, "%s -mMask %s", args, moving_mask_path);

    if (fixed_mask_path != NULL)
        sprintf(args, "%s -fMask %s", args, fixed_mask_path);

    if (output_file != NULL){
        sprintf(args, "%s >> %s", args, output_file);
    }

    printf("\n * args = %s\n", args);

    iftRunProgram("elastix", args);

    char *reg_path = NULL;
    for (int i = 0; i < transf_files->n; i++) {
        char *path = NULL;

        char filename[128];
        sprintf(filename, "result.%d.nii.gz", i);
        path = iftJoinPathnames(2, out_dir, filename);
        printf("path = %s\n", path);
        
        if (iftFileExists(path))
            reg_path = path;
        else iftFree(path);
    }

    if (reg_path == NULL)
        iftError("Resulting Images (result.?.nii.gz) not found", "_iftRunElastix");

    iftImage *reg_img = iftReadImageByExt(reg_path);
    iftFree(reg_path);

    return reg_img;
}

/**
 * @brief Runs the Elastix program.
 * @author Samuel
 * @November 17, 2015
 * @ingroup Elastix
 */
iftImage *_iftRunElastixUnstop(const char *moving_path, const char *fixed_path, const char *moving_mask_path,
                         const char *fixed_mask_path, const iftFileSet *transf_files, const char *out_dir, const char output_file[]) {
    char args[4096];

    // Registration: Assuming elastix is already in your OS PATH
    // just turn on the BSpline Registration
    sprintf(args, "-m %s -f %s", moving_path, fixed_path);
    for (size_t i = 0; i < transf_files->n; i++)
        sprintf(args, "%s -p %s", args, transf_files->files[i]->path);
    sprintf(args, "%s -out %s", args, out_dir);

    if (moving_mask_path != NULL)
        sprintf(args, "%s -mMask %s", args, moving_mask_path);

    if (fixed_mask_path != NULL)
        sprintf(args, "%s -fMask %s", args, fixed_mask_path);

    if (output_file != NULL){
        sprintf(args, "%s >> %s", args, output_file);
    }

    char *cmd = iftConcatStrings(2, "elastix ", args);
    if (system(cmd) != 0){
        iftFree(cmd);
        return NULL;
    }
    iftFree(cmd);

    char *reg_path = NULL;
    for (int i = 0; i < transf_files->n; i++) {
        char *path = NULL;

        char filename[128];
        sprintf(filename, "result.%d.nii.gz", i);
        path = iftJoinPathnames(2, out_dir, filename);

        if (iftFileExists(path))
            reg_path = path;
        else iftFree(path);
    }

    if (reg_path == NULL)
        iftError("Resulting Images (result.?.nii.gz) not found", "_iftRunElastix");

    iftImage *reg_img = iftReadImageByExt(reg_path);
    iftFree(reg_path);

    return reg_img;
}


iftFileSet *_iftGetElatixDefFields(const char *elastix_dir, const char *new_basename) {
    iftFileSet *def_fields = iftLoadFileSetFromDirByRegex(elastix_dir, "^TransformParameters\\.[0-9]+\\.txt$", true);
    
    if (new_basename) {
        char *parent_dir = iftParentDir(new_basename);
        iftMakeDir(parent_dir);
        iftFree(parent_dir);
    }
    
    iftFileSet *def_fields_out = iftCreateFileSet(def_fields->n);
    for (long i = 0; i < def_fields->n; i++) {
        char *path = NULL;
        
        if (new_basename) {
            char *filename = iftCopyString(".%ld.txt", i);
            path = iftConcatStrings(2, new_basename, filename);
            iftFree(filename);
        }
        else {
            char suffix[512];
            sprintf(suffix, ".%ld.txt", i);
            path = iftMakeTempFile("DefFields.", suffix, NULL);
        }
        
        def_fields_out->files[i] = iftCreateFile(path);
        iftCopyFromDisk(def_fields->files[i]->path, def_fields_out->files[i]->path);
    
        // changes the previous transformation pathname inside the current transformation field by the new pathname
        if ((i - 1) >= 0)
            iftSed(def_fields_out->files[i]->path, def_fields->files[i - 1]->path, def_fields_out->files[i - 1]->path);
            
        // Avoids holes when transforming a binary mask
        // http://elastix.isi.uu.nl/FAQ.php#Q_TransformBinarySegmentation
        iftSed(def_fields_out->files[i]->path, "(FinalBSplineInterpolationOrder 3)", "(FinalBSplineInterpolationOrder 0)");
        
        iftFree(path);
    }
    iftDestroyFileSet(&def_fields);
    
    return def_fields_out;
}

/********************** PUBLIC FUNCTIONS *************************/
iftImage *iftRegisterImageByElastix(const iftImage *moving_img, const iftImage *fixed_img, const iftImage *moving_mask,
                                    const iftImage *fixed_mask, const iftFileSet *elastix_files,
                                    const char *def_fields_basename, iftFileSet **def_fields, const char output_file[]) {

    if (iftIsColorImage(moving_img))
        iftError("Moving Image is colored: Our current Elastix's wrap only accepts gray images", "_iftRegisterImageByElastix");
    if (iftIsColorImage(fixed_img))
        iftError("Fixed Image is colored: Our current Elastix's wrap only accepts gray images", "_iftRegisterImageByElastix");
    
    // Saving Images
    char *tmp_dir     = iftMakeTempDir("tmpdir_", NULL, NULL);
    char *moving_path = iftJoinPathnames(2, tmp_dir, "moving_img.nii.gz");
    char *fixed_path  = iftJoinPathnames(2, tmp_dir, "fixed_img.nii.gz");
    iftWriteImageByExt(moving_img, moving_path);
    iftWriteImageByExt(fixed_img, fixed_path);
    
    char *moving_mask_path = NULL;
    if (moving_mask) {
        moving_mask_path = iftJoinPathnames(2, tmp_dir, "moving_mask.nii.gz");
        iftWriteImageByExt(moving_mask, moving_mask_path);
    }
    
    char *fixed_mask_path = NULL;
    if (fixed_mask) {
        fixed_mask_path = iftJoinPathnames(2, tmp_dir, "fixed_mask.nii.gz");
        iftWriteImageByExt(fixed_mask, fixed_mask_path);
    }
    
    iftFileSet *elastix_files_final = _iftAdjustElastixParameterFiles(elastix_files, moving_img, fixed_img, tmp_dir);
    iftImage *reg_img = _iftRunElastix(moving_path, fixed_path, moving_mask_path, fixed_mask_path, elastix_files_final, tmp_dir, output_file);
    
    // elastix can output a shifted image
    int min, max;
    iftMinMaxValues(reg_img, &min, &max);
    
    if (min < 0) {
        printf("- Shifting the Image: [%d, %d] to ", min, max);
        iftAddScalarInPlace(reg_img, -1*min); // right shift
        max = max + (-1 * min);
        min = 0;
        printf("[%d, %d]\n", min, max);
    }
    
    if (def_fields)
        *def_fields = _iftGetElatixDefFields(tmp_dir, def_fields_basename);
    
    // deleting temporary elastix output directory
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
    
    iftFree(moving_path);
    iftFree(fixed_path);
    iftFree(moving_mask_path);
    iftFree(fixed_mask_path);
    iftDestroyFileSet(&elastix_files_final);
    
    return reg_img;
}

iftImage *iftRegisterImageByElastixUnstop(const iftImage *moving_img, const iftImage *fixed_img, const iftImage *moving_mask,
                                    const iftImage *fixed_mask, const iftFileSet *elastix_files,
                                    const char *def_fields_basename, iftFileSet **def_fields, const char output_file[]) {

    if (iftIsColorImage(moving_img))
        iftError("Moving Image is colored: Our current Elastix's wrap only accepts gray images", "_iftRegisterImageByElastix");
    if (iftIsColorImage(fixed_img))
        iftError("Fixed Image is colored: Our current Elastix's wrap only accepts gray images", "_iftRegisterImageByElastix");

    // Saving Images
    char *tmp_dir     = iftMakeTempDir("tmpdir_", NULL, NULL);
    char *moving_path = iftJoinPathnames(2, tmp_dir, "moving_img.nii.gz");
    char *fixed_path  = iftJoinPathnames(2, tmp_dir, "fixed_img.nii.gz");
    iftWriteImageByExt(moving_img, moving_path);
    iftWriteImageByExt(fixed_img, fixed_path);

    char *moving_mask_path = NULL;
    if (moving_mask) {
        moving_mask_path = iftJoinPathnames(2, tmp_dir, "moving_mask.nii.gz");
        iftWriteImageByExt(moving_mask, moving_mask_path);
    }

    char *fixed_mask_path = NULL;
    if (fixed_mask) {
        fixed_mask_path = iftJoinPathnames(2, tmp_dir, "fixed_mask.nii.gz");
        iftWriteImageByExt(fixed_mask, fixed_mask_path);
    }

    iftFileSet *elastix_files_final = _iftAdjustElastixParameterFiles(elastix_files, moving_img, fixed_img, tmp_dir);
    iftImage *reg_img = _iftRunElastixUnstop(moving_path, fixed_path, moving_mask_path, fixed_mask_path, elastix_files_final, tmp_dir, output_file);

    // elastix can output a shifted image
    if (reg_img != NULL) {
        int min, max;
        iftMinMaxValues(reg_img, &min, &max);

        if (min < 0) {
            iftAddScalarInPlace(reg_img, -1 * min); // right shift
            max = max + (-1 * min);
            min = 0;
        }

        if (def_fields){	
            *def_fields = _iftGetElatixDefFields(tmp_dir, def_fields_basename);
	}
    }

    // deleting temporary elastix output directory
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);

    iftFree(moving_path);
    iftFree(fixed_path);
    iftFree(moving_mask_path);
    iftFree(fixed_mask_path);
    iftDestroyFileSet(&elastix_files_final);

    return reg_img;
}

iftImage *iftRegisterImageByElastixIntoBestNormalSpace(iftImage *moving_img, iftImage *input_mask, iftDataSet *Z_normal_spaces,
                                                       const char *normal_mask_dir, iftFileSet *transf_files,
                                                       const char *def_fields_basename, int *best_space_id, char **best_space_name,
                                                       iftFloatArray **registration_errors)
{
    char filename[4096];

    /* Transforming input image into a sample */
    iftSample *s   = (iftSample*)calloc(1,sizeof(iftSample));
    s->feat        = (float*)calloc(Z_normal_spaces->nfeats,sizeof(float));
    s->id          = Z_normal_spaces->nsamples;

    for (int f = 0; f < Z_normal_spaces->nfeats; f++)
        s->feat[f] = moving_img->val[f];

    /* Computing nearest sample on the dataset Z and define the group of the new sample */
    iftIntArray *group_idx = NULL;
    iftFloatArray *dists = iftSampleDistanceToGroups(Z_normal_spaces,s,&group_idx,1,TRUE);
    iftDestroyFloatArray(&dists);
    iftFree(s->feat);
    iftFree(s);

    /* Registering input image into images of the same group */
    iftFileSet *fs     = NULL;
    if (Z_normal_spaces->ref_data_type == IFT_REF_DATA_FILESET)
        fs = Z_normal_spaces->ref_data;
    else
        iftError("DataSet Reference Data must be a fileset","main");


    iftIntArray *samples_per_group       = iftCreateIntArray(Z_normal_spaces->ngroups);
    iftFree(samples_per_group->val);
    samples_per_group->val               = iftCountSamplesPerGroupDataSet(Z_normal_spaces);
    iftImage *best_reg                   = NULL;
    float min_error                      = IFT_INFINITY_DBL;
    int pos                              = 0;
    iftFileSet *best_df                  = NULL;
    char *tmp_dir                        = iftMakeTempDir("temp_dir",NULL,"./");
    int g = 1;
    while ((best_reg == NULL) && (g <= Z_normal_spaces->ngroups)) {
        if (*registration_errors != NULL){
            iftDestroyFloatArray(registration_errors);
        }
        *registration_errors = iftCreateFloatArray(samples_per_group->val[group_idx->val[g]]);
	printf("Running for group %d\n",group_idx->val[g]);
        for (int i = 0; i < Z_normal_spaces->nsamples; i++) {

            if (Z_normal_spaces->sample[i].group != group_idx->val[g])
                continue;

            char output_file[200] = "./temp.txt";

            iftFileSet *def_fields = NULL;
            iftImage *fixed_img = iftReadImageByExt(fs->files[Z_normal_spaces->sample[i].id]->path);
            iftImage *fixed_mask = NULL;

            if (normal_mask_dir != NULL) {
                char *f = iftFilename(fs->files[Z_normal_spaces->sample[i].id]->path, NULL);
                sprintf(filename, "%s/%s", normal_mask_dir, f);
                fixed_mask = iftReadImageByExt(filename);
                iftFree(f);
            }

            iftImage *reg_img = iftRegisterImageByElastixUnstop(moving_img, fixed_img, input_mask, fixed_mask,
                                                                transf_files, def_fields_basename, &def_fields,
                                                                output_file);

            if (reg_img == NULL){
                (*registration_errors)->val[pos] = IFT_INFINITY_DBL;
	    }
            else {
                (*registration_errors)->val[pos] = iftRegistrationRMSE(reg_img, fixed_img);

                if ((*registration_errors)->val[pos] < min_error) {
                    min_error = (*registration_errors)->val[pos];
                    iftDestroyImage(&best_reg);
                    best_reg = iftCopyImage(reg_img);
                    *best_space_name = iftCopyString(fs->files[i]->path);
                    *best_space_id = pos;
                    iftDestroyFileSet(&best_df);
		    puts(tmp_dir);
		    for (int a = 0; a < def_fields->n; a++)
			    puts(def_fields->files[a]->path);
                    best_df = iftMoveElastixDefFields(def_fields, tmp_dir);
                }
            }

            pos++;
            iftDestroyImage(&reg_img);
            iftDestroyImage(&fixed_img);
            iftDestroyImage(&fixed_mask);
            iftDestroyFileSet(&def_fields);
        }
        g++;
    }
    iftDestroyIntArray(&samples_per_group);
    iftDestroyIntArray(&group_idx);

    puts(def_fields_basename);
    if (best_df == NULL)
	    puts("best_df is NULL");
    iftFileSet *aux = iftMoveElastixDefFields(best_df,def_fields_basename);
    iftDestroyFileSet(&aux);
    iftRemoveDir(tmp_dir);
    free(tmp_dir);
    iftDestroyFileSet(&best_df);

    return best_reg;
}


void iftRunProgRegisterImageSetByElastix(const char *moving_img_entry, const char *fixed_img_path, int img_depth,
                                         const char *affine_params_path, const char *bspline_params_path,
                                         const char *out_dir) {
    char opts[4096];
    sprintf(opts, "--moving-image-entry %s --fixed-image %s --image-depth %d --output-dir %s",
            moving_img_entry, fixed_img_path, img_depth, out_dir);
    if (affine_params_path != NULL)
        sprintf(opts, "%s --affine-param-file %s", opts, affine_params_path);
    if (bspline_params_path != NULL)
        sprintf(opts, "%s --bspline-param-file %s", opts, bspline_params_path);

    iftRunProgram("iftRegisterImageSetByElastix", opts);
}


iftImage *iftTransformImageByTransformix(const iftImage *img, const char *def_fields_path) {
    char args[4096];
    char *tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);

    char *img_path = iftJoinPathnames(2, tmp_dir, "image.nii.gz");
    iftWriteImageByExt(img, img_path);

    sprintf(args, "-in %s -out %s -tp %s", img_path, tmp_dir, def_fields_path);
    iftRunProgram("transformix", args);

    char *mapped_path = iftJoinPathnames(2, tmp_dir, "result.nii.gz");
    iftImage *reg_img = iftReadImageByExt(mapped_path);
    iftFree(mapped_path);

    // DESTROYERS
    iftFree(img_path);
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);

    return reg_img;
}


void iftRunProgTransformImageSetByTransformix(const char *img_entry, const char *affine_def_fields_entry,
                                              const char *bspline_def_fields_entry, const char *out_dir) {
    char opts[4096];
    sprintf(opts, "--img-entry %s --affine-def-fields-entry %s --output-dir %s", img_entry, affine_def_fields_entry,
            out_dir);
    if (bspline_def_fields_entry != NULL)
        sprintf(opts, "%s --bspline-def-fields-entry %s", opts, bspline_def_fields_entry);

    iftRunProgram("iftTransformImageSetByTransformix", opts);
}


iftFileSet *iftMoveElastixDefFields(const iftFileSet *def_fields, const char *out_basename) {
    if (def_fields == NULL)
        iftError("Def. Fields is NULL", "iftMoveElastixDefFields");
    if (out_basename == NULL)
        iftError("Output basename is NULL", "iftMoveElastixDefFields");
    
    char *parent_dir = iftParentDir(out_basename);
    if (!iftDirExists(parent_dir))
        iftMakeDir(parent_dir);
    iftFree(parent_dir);
    
    iftFileSet *new_def_fields = iftCreateFileSet(def_fields->n);
    
    for (long i = 0; i < def_fields->n; i++) {
        char *out_def_field = iftCopyString("%s.%ld.txt", out_basename, i);
        new_def_fields->files[i] = iftCreateFile(out_def_field);
        iftCopyFromDisk(def_fields->files[i]->path, new_def_fields->files[i]->path);
        iftRemoveFile(def_fields->files[i]->path);
        
        if (i > 0)
            iftSed(new_def_fields->files[i]->path, def_fields->files[i - 1]->path, new_def_fields->files[i - 1]->path);
        iftFree(out_def_field);
    }
    
    //iftDestroyFileSet(&new_def_fields);
    
    return new_def_fields;
}


