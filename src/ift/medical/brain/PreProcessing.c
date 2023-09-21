#include "ift/medical/brain/PreProcessing.h"

#include "ift/core/tools/OS.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Stream.h"
#include "ift/medical/brain/MSP.h"
#include "iftFiltering.h"
#include "iftInterpolation.h"
#include "iftRadiometric.h"


/********************** PRIVATE FUNCTIONS *************************/



/********************** PUBLIC FUNCTIONS *************************/
void iftSetANTsEnvironment(void) {
    iftOS os = iftGetOS();
    char env[IFT_STR_DEFAULT_SIZE];
    sprintf(env, "%s:%s/externals/ANTs/", getenv("PATH"), getenv("NEWIFT_DIR"));

    if (os == IFT_LINUX)
        strcat(env, "Linux");
    else if (os == IFT_MACOS)
        strcat(env, "MacOS");

    setenv("PATH", env, true);
}


iftImage *iftN4BiasFieldCorrection(const iftImage *img, const iftImage *mask, int shrink_factor,
                                   iftImage **out_bias) {
    iftSetANTsEnvironment();
    
    char *tmp_dir = iftMakeTempDir("tmpdir_", NULL, NULL);
    char *input_img_path = iftJoinPathnames(2, tmp_dir, "input_img.nii.gz");
    char *output_img_path = iftJoinPathnames(2, tmp_dir, "output_img.nii.gz");
    
    iftWriteImageByExt(img, input_img_path);
    
    char args[IFT_STR_DEFAULT_SIZE];
    sprintf(args, "--input-image %s", input_img_path);
    
    
    // If a mask was passed
    if (mask != NULL) {
        if (!iftIsDomainEqual(img, mask)) {
            iftError("Input Image and Mask with different domains\n(%d, %d, %d)\n(%d, %d, %d)",
                     "iftN4BiasFieldCorrection", img->dx, img->dy, img->dz, mask->dx, mask->dy, mask->dz);
        }
        
        char *mask_path = iftJoinPathnames(2, tmp_dir, "mask.nii.gz");
        iftWriteImageByExt(mask, mask_path);
        sprintf(args, "%s --mask-image %s", args, mask_path);
        iftFree(mask_path);
    }
    
    // Check the image dimension
    if (iftIs3DImage(img)) {
        sprintf(args, "%s --image-dimensionality 3", args);
    }
    else sprintf(args, "%s --image-dimensionality 2", args);
    
    sprintf(args, "%s --shrink-factor %d", args, shrink_factor);
    sprintf(args, "%s --convergence [50x50x50x50, 0.0001]", args);
    
    char *output_bias_path = NULL;
    // don't save output bias image
    if (out_bias == NULL) {
        sprintf(args, "%s  --output %s", args, output_img_path);
    }
    else {
        output_bias_path = iftJoinPathnames(2, tmp_dir, "output_bias.nii.gz");
        sprintf(args, "%s  --output [%s, %s]", args, output_img_path, output_bias_path);
    }
    
    iftRunProgram("N4BiasFieldCorrection", args);
    
    iftImage *out_img = iftReadImageByExt(output_img_path);
    
    // Read the resulting output bias
    if (output_bias_path != NULL) {
        *out_bias = iftReadImageByExt(output_bias_path);
    }
    
    iftRemoveDir(tmp_dir);
    iftFree(tmp_dir);
    iftFree(input_img_path);
    iftFree(output_img_path);
    iftFree(output_bias_path);
    
    return out_img;
}


iftImage *iftEffectiveN4BiasFieldCorrection(const iftImage *img, const iftImage *mask, iftImage **out_bias) {
    iftImage *corr = iftN4BiasFieldCorrection(img, mask, 4, NULL);
    iftImage *corr_final = iftN4BiasFieldCorrection(corr, mask, 6, out_bias);
    iftDestroyImage(&corr);
    
    return corr_final;
}


iftImage *iftBrainMRIPreProcessing(const iftImage *mri, int nbits, const iftPlane *msp_in, const iftImage *mask,
                                   const iftImage *ref_mri, const iftImage *ref_mask, bool skip_n4,
                                   bool skip_median_filter, bool skip_msp_alignment, bool skip_hist_matching,
                                   iftPlane **msp_out) {
    iftImage *out = iftCopyImage(mri);
    iftImage *aux = NULL;
    
    // N4 Bias Field Correction
    if (!skip_n4) {
        aux = out;
        out = iftEffectiveN4BiasFieldCorrection(out, NULL, NULL);
        iftDestroyImage(&aux);
    }
    
    // Median Filtering
    if (!skip_median_filter) {
        aux = out;
        out = iftMedianFilter(out, NULL);
        iftDestroyImage(&aux);
    }
    
    // Intensity Normalization
    if (nbits > 0) {
        int min, max;
        iftMinMaxValues(out, &min, &max);
        int new_max = powf(2, nbits) - 1;
        
        if (min != 0 || max != new_max) {
            aux = out;
            out = iftLinearStretch(out, min, max, 0, new_max);
            iftDestroyImage(&aux);
        }
    }
    
    // MSP alignment
    if (!skip_msp_alignment) {
        aux = out;
        
        iftPlane *msp = iftCopyPlane(msp_in);
        
        if (msp_in)
            out = iftRotateImageToMSP(out, msp_in, IFT_LINEAR_INTERPOLATION);
        else
            out = iftAlignBrainByMSP(out, &msp);
        
        if (msp_out)
            *msp_out = msp;
        else iftDestroyPlane(&msp);
        
        iftDestroyImage(&aux);
    }
    
    // histogram matching
    if (!skip_hist_matching && ref_mri) {
        aux = out;
        out = iftMatchHistogram(out, mask, ref_mri, ref_mask);
        iftDestroyImage(&aux);
    }
    
    return out;
}

