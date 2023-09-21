#include "ift/medical/brain/MSP.h"

#include "ift/imgproc/basic/Histogram.h"
#include "iftCommon.h"
#include "iftGeometric.h"
#include "iftPlane.h"


/********************** PRIVATE FUNCTIONS *************************/
/**
 * @brief Params for the MSP searching algorithm.
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
typedef struct ift_msp_search_params {
    /** Factor that shrinks the original brain image for speeding the search up */
    float downsampling_factor;
    /** Histogram percentage threshold used to binarize the gradient image in the algorithm. */
    float hist_thres;
    /** Step considering during the msp search (main for-loop). */
    int step;
    /** (Closed) Interval used to vary the first point of the MSP during the search (main for-loop). */
    iftInterval xslice_range_A;
    /** (Closed) Interval used to vary the second point of the MSP during the search (main for-loop). */
    iftInterval xslice_range_B;
    /** (Closed) Interval used to vary the second third of the MSP during the search (main for-loop). */
    iftInterval xslice_range_C;
} _iftMSPSearchParams;


/**
 * @brief Gets the parameters for MSP Searching for stage/scale #1.
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
_iftMSPSearchParams _iftGetMSPSearchParamsStage1(const iftImage *img) {
    _iftMSPSearchParams msp_params;
    msp_params.downsampling_factor = 0.25;
    msp_params.hist_thres = 0.92;
    msp_params.step = 2;
    msp_params.xslice_range_A.begin = 0;
    
    int xsize_interp = fabs(msp_params.downsampling_factor * img->xsize);
    
    msp_params.xslice_range_A.end = iftRound(xsize_interp - 1);
    msp_params.xslice_range_B.begin = 0;
    msp_params.xslice_range_B.end = iftRound(xsize_interp - 1);
    msp_params.xslice_range_C.begin = 0;
    msp_params.xslice_range_C.end = iftRound(xsize_interp - 1);
    
    return msp_params;
}


/**
 * @brief Gets the parameters for MSP Searching for stage/scale #2.
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
_iftMSPSearchParams _iftGetMSPSearchParamsStage2(const iftImage *img, iftPlane *msp_new_scale,
                                                 float downsampling_factor, float scaling_factor2,
                                                 _iftMSPSearchParams msp_params1) {
    _iftMSPSearchParams msp_params2;
    msp_params2.downsampling_factor = downsampling_factor;
    msp_params2.hist_thres = 0.95;
    msp_params2.step = 1;
    
    // Transforms the previous MSP to the new scale
    int ysize_interp = fabs(msp_params2.downsampling_factor * img->ysize);
    int zsize_interp = fabs(msp_params2.downsampling_factor * img->zsize);
    
    // extreme points from the image axes/edges on the plane
    iftPoint A = iftGetPointFromYZPlaneProjection(msp_new_scale, 0, 0);
    iftPoint B = iftGetPointFromYZPlaneProjection(msp_new_scale, ysize_interp - 1, 0);
    iftPoint C = iftGetPointFromYZPlaneProjection(msp_new_scale, 0, zsize_interp - 1);
    int range = msp_params1.step * scaling_factor2;
    
    msp_params2.xslice_range_A.begin = A.x - range;
    msp_params2.xslice_range_A.end = A.x + range;
    msp_params2.xslice_range_B.begin = B.x - range;
    msp_params2.xslice_range_B.end = B.x + range;
    msp_params2.xslice_range_C.begin = C.x - range;
    msp_params2.xslice_range_C.end = C.x + range;
    
    return msp_params2;
}


/**
 * @brief Gets the parameters for MSP Searching for stage/scale #3.
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
_iftMSPSearchParams _iftGetMSPSearchParamsStage3(const iftImage *img, iftPlane *msp_new_scale,
                                                 float downsampling_factor, float scaling_factor3,
                                                 _iftMSPSearchParams msp_params2) {
    _iftMSPSearchParams msp_params3;
    msp_params3.downsampling_factor = downsampling_factor;
    msp_params3.hist_thres = 0.95;
    msp_params3.step = 1;
    
    // Transforms the previous MSP to the new scale
    int ysize_interp = fabs(msp_params3.downsampling_factor * img->ysize);
    int zsize_interp = fabs(msp_params3.downsampling_factor * img->zsize);
    
    // extreme points from the image axes/edges on the plane
    iftPoint A = iftGetPointFromYZPlaneProjection(msp_new_scale, 0, 0);
    iftPoint B = iftGetPointFromYZPlaneProjection(msp_new_scale, ysize_interp - 1, 0);
    iftPoint C = iftGetPointFromYZPlaneProjection(msp_new_scale, 0, zsize_interp - 1);
    int range = msp_params2.step * scaling_factor3;
    
    msp_params3.xslice_range_A.begin = A.x - range;
    msp_params3.xslice_range_A.end = A.x + range;
    msp_params3.xslice_range_B.begin = B.x - range;
    msp_params3.xslice_range_B.end = B.x + range;
    msp_params3.xslice_range_C.begin = C.x - range;
    msp_params3.xslice_range_C.end = C.x + range;
    
    return msp_params3;
}


// WE CONSIDER THANT THE NORM OF THE PLANE IS A UNIT VECTOR

/**
 * @brief Computes a score of symmetry based on a given MSP for a set of object voxels <obj_voxels> from a binary image.
 *
 * @attention The function considers that the normal vector of the MSP a unit vector.
 *
 * For each object voxel from <obj_voxels>, it find its correspoding symmetric voxel on the side of the MSP, and check
 * if it also an object voxel.
 *
 * @param bin_mask Binary image.
 * @param obj_voxels Array of selected object voxels used during score computation.
 * @param msp Mid-Sagittal Plane used to find the corresponding symmetric voxel for each object voxel from the set.
 * @return Resulting computed score.
 *
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
float _iftSymmetricScoreBin(const iftImage *bin_mask, const iftIntArray *obj_voxels, const iftPlane *msp) {
    iftVector unit_norm = msp->normal; // it considers that the normal vector is a unit vector
    
    // general plane equation
    // a*x + b*y + c*z + d = 0, where norm = (a, b, c), and d is the translation of the plane with the origin (bias)
    // If bias is not computed from a unit norm it does make sense to use it for the distance computations
    float d = - (unit_norm.x * msp->pos.x + unit_norm.y * msp->pos.y + unit_norm.z * msp->pos.z);
    
    float n_symmetric_obj_voxels = 0;

    #pragma omp parallel for reduction(+:n_symmetric_obj_voxels)
    for (int i = 0; i < obj_voxels->n; i++) {
        iftVoxel u = iftGetVoxelCoord(bin_mask, obj_voxels->val[i]);
        
        // signed distance from u to the plane = dot product between u and the unit norm of the plane
        //                                       + the bias d
        float signed_dist = iftVectorInnerProduct(u, unit_norm) + d;
        float dist = fabsf(signed_dist);
        int sign = (signed_dist > 0) ? -1 : 1;
        
        // if the point is not on the plane
        if (dist >= 1) {
            // Find the symmetric point v to the point u
            iftVoxel v;
            v.x = iftRound(u.x + (sign * 2 * dist * unit_norm.x));
            v.y = iftRound(u.y + (sign * 2 * dist * unit_norm.y));
            v.z = iftRound(u.z + (sign * 2 * dist * unit_norm.z));
            
            n_symmetric_obj_voxels += (iftValidVoxel(bin_mask, v) && iftImgVoxelVal(bin_mask, v));
        }
    }
    
    return (n_symmetric_obj_voxels / (float) obj_voxels->n);
}


/**
 * @brief Finds the MSP for the input brain Image by considering a set of parameters for the Search.
 *
 * @param img Input Brain images.
 * @param msp_params Parameters considered for MSP searching.
 * @return Resulting MSP.
 *
 * @author Samuel Martins
 * @date Nov 22, 2018
 */
iftPlane *_iftFindMSP(const iftImage *img, _iftMSPSearchParams msp_params) {
    //    iftImage *grad = iftTextGradient(img);  // gradient used in the original paper, but it results in more voxels to compute the matching score
    iftImage *grad = iftImageBasins(img, NULL);
    
    int max_val = iftMaximumValue(grad);
    int nbins = max_val + 1;

    iftHist *hist_norm = iftCalcGrayImageHist(grad, NULL, nbins, max_val, true);
    iftHist *acc_hist_norm = iftCalcAccHist(hist_norm);
    int thres = iftHistArgGreaterThan(acc_hist_norm, msp_params.hist_thres, false);
    iftDestroyHist(&hist_norm);
    iftDestroyHist(&acc_hist_norm);
    
    iftImage *bin_mask = iftThreshold(grad, thres, IFT_INFINITY_INT, 1);
    iftIntArray *obj_voxels = iftGetObjectVoxels(bin_mask);
    iftDestroyImage(&grad);
    
    iftVoxel img_center = iftImageCenter(img);
    float min_dist_between_centers = img->xsize / 6.0;
    iftPlane *best_pl = iftCreatePlane();
    float best_score = 0;
    
    
    for (int xa = msp_params.xslice_range_A.begin; xa <= msp_params.xslice_range_A.end; xa += msp_params.step) {
        for (int xb = msp_params.xslice_range_B.begin; xb <= msp_params.xslice_range_B.end; xb += msp_params.step) {
            for (int xc = msp_params.xslice_range_C.begin; xc <= msp_params.xslice_range_C.end; xc += msp_params.step) {
                iftPoint A = {xa, 0, 0};
                iftPoint B = {xb, img->ysize - 1, 0};
                iftPoint C = {xc, 0, img->zsize - 1};
                
                iftPlane *pl = iftGetPlaneFrom3Points(A, B, C, true);
                
                // consider only planes close to the image center
                iftPoint plane_center = iftGetPointFromYZPlaneProjection(pl, img_center.y, img_center.z);
                float dist_centers = fabs(img_center.x - plane_center.x);
                if (dist_centers > min_dist_between_centers) {
                    iftDestroyPlane(&pl);
                    continue;
                }
                
                float score = _iftSymmetricScoreBin(bin_mask, obj_voxels, pl);
                
                if (best_score < score) {
                    best_score = score;
                    best_pl->normal = pl->normal;
                    best_pl->pos = pl->pos;
                }
                
                iftDestroyPlane(&pl);
            }
        }
    }
    iftDestroyImage(&bin_mask);
    iftDestroyIntArray(&obj_voxels);
    
    return best_pl;
}






/********************** PUBLIC FUNCTIONS *************************/
iftImage *iftAlignBrainByMSP(const iftImage *img, iftPlane **msp_out) {
    puts("\n\n### MSP Extraction ###");
    puts("# FIRST STAGE");
    _iftMSPSearchParams msp_params1 = _iftGetMSPSearchParamsStage1(img);
    iftImage *interp_img = iftInterp(img, msp_params1.downsampling_factor, msp_params1.downsampling_factor,
                                     msp_params1.downsampling_factor);
    iftPlane *msp1 = _iftFindMSP(interp_img, msp_params1);
    iftDestroyImage(&interp_img);
    
    
    puts("# SECOND STAGE");
    float downsampling_factor2 = 0.5;
    float scaling_factor2 = downsampling_factor2 / msp_params1.downsampling_factor;
    iftPlane *msp1_scaled = iftScalePlane(msp1, scaling_factor2);
    
    _iftMSPSearchParams msp_params2 = _iftGetMSPSearchParamsStage2(img, msp1_scaled, downsampling_factor2,
                                                                   scaling_factor2, msp_params1);
    interp_img = iftInterp(img, msp_params2.downsampling_factor, msp_params2.downsampling_factor,
                           msp_params2.downsampling_factor);
    iftPlane *msp2 = _iftFindMSP(interp_img, msp_params2);
    iftDestroyImage(&interp_img);
    
    
    puts("# THIRD STAGE");
    float downsampling_factor3 = 1.0;
    float scaling_factor3 = downsampling_factor3 / msp_params2.downsampling_factor;
    iftPlane *msp2_scaled = iftScalePlane(msp2, scaling_factor3);
    
    _iftMSPSearchParams msp_params3 = _iftGetMSPSearchParamsStage3(img, msp2_scaled, downsampling_factor3,
                                                                   scaling_factor3, msp_params2);
    iftPlane *msp = _iftFindMSP(img, msp_params3);
    //    iftWriteImageByExt(iftMSPAsImage(msp, img), "tmp/msp.nii.gz");
    
    puts("# ALIGNING IMAGE BY MSP");
    iftImage *rotate_img = iftRotateImageToMSP(img, msp, IFT_LINEAR_INTERPOLATION);
    
    iftDestroyPlane(&msp1);
    iftDestroyPlane(&msp1_scaled);
    iftDestroyPlane(&msp2);
    iftDestroyPlane(&msp2_scaled);
    
    if (msp_out)
        *msp_out = msp;
    else iftDestroyPlane(&msp);
    puts("###################\n");
    
    return rotate_img;
}




iftImage *iftRotateImageToMSP(const iftImage *img, const iftPlane *msp, iftInterpolationType interp_type) {
    iftInterpPointFunc imageValueAtPoint = iftImageValueAtPoint;
    if (interp_type == IFT_NEAREST_NEIGHBOR_INTERPOLATION)
        imageValueAtPoint = iftImageValueAtPointNearestNeighbor;
    else if (interp_type == IFT_LINEAR_INTERPOLATION)
        imageValueAtPoint = iftImageValueAtPoint;
    
    
    // general equantion of the plane: ax + by + cz + d = 0, where norm vector of the plane = (a, b, c)
    float a = msp->normal.x;
    float b = msp->normal.y;
    float c = msp->normal.z;
    
    // angles for rotation
    float theta_y = atan(c / a) * (180.0 / IFT_PI);
    float theta_z = - atan(b / a) * (180.0 / IFT_PI);
    
    iftPoint msp_center = iftGetPointFromYZPlaneProjection(msp, img->ysize / 2.0, img->zsize / 2.0);
    iftVector trans_center = {.x = - msp_center.x, .y = - msp_center.y, .z = - msp_center.z};
    float disp_x = (img->xsize / 2.0) - msp_center.x; // puts the MSP on the xsize / 2 slice
    iftVector trans_plane = {.x = msp_center.x + disp_x, .y = msp_center.y, .z = msp_center.z};
    
    iftMatrix *Tc = iftTranslationMatrix(trans_center);
    iftMatrix *Ry = iftRotationMatrix(IFT_AXIS_Y, theta_y);
    iftMatrix *Rz = iftRotationMatrix(IFT_AXIS_Z, theta_z);
    iftMatrix *Tp = iftTranslationMatrix(trans_plane);
    
    iftMatrix *A = iftMultMatricesChain(4, Tp, Rz, Ry, Tc);
    iftMatrix *Ainv = iftInvertMatrix(A);
    
    iftImage *rot_img = iftCreateImageFromImage(img);

    #pragma omp parallel for
    for (int q = 0; q < rot_img->n; q++) {
        iftVoxel v = iftGetVoxelCoord(rot_img, q);
        iftPoint Pv = iftVoxelToPoint(v);
        iftPoint Pu = iftTransformPoint(Ainv, Pv);
        
        iftImgVoxelVal(rot_img, v) = imageValueAtPoint(img, Pu);
    }
    
    iftDestroyMatrix(&Tc);
    iftDestroyMatrix(&Tp);
    iftDestroyMatrix(&Ry);
    iftDestroyMatrix(&Rz);
    iftDestroyMatrix(&A);
    iftDestroyMatrix(&Ainv);
    
    return rot_img;
}


iftImage *iftMSPAsImage(const iftPlane *msp, const iftImage *ref_img) {
    iftImage *msp_img = iftCreateImageFromImage(ref_img);

    #pragma omp parallel for
    for (int z = 0; z < msp_img->zsize; z++)
        for (int y = 0; y < msp_img->ysize; y++) {
            iftVoxel v = iftPointToVoxel(iftGetPointFromYZPlaneProjection(msp, y, z));
            if (iftValidVoxel(msp_img, v))
                iftImgVoxelVal(msp_img, v) = 1;
        }
    
    return msp_img;
}

