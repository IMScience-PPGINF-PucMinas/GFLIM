#include "iftImage.h"

#include "iftPng.h"
#include "jpeglib.h"
#include "nifti1.h" // official definitions of NIfTI header

//#include  "tiffio.h"
#include "iftCompression.h"
#include "iftInterpolation.h"

#include "ift/core/dtypes/IntMatrix.h"
#include "ift/core/io/Dir.h"
#include "ift/core/io/Json.h"
#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/Numerical.h"
#include "ift/core/tools/Regex.h"
#include "ift/core/tools/String.h"


/**
@file
@brief Image manipulation methods using the iftImage structure.
*/

/********************** PRIVATE FUNCTIONS *************************/
float IntSwap (int f)                       //Convert int data to big-endian, the default
{                                           //binary format required by ParaView
    union
    {
        int f;
        unsigned char b[4];
    } dat1, dat2;

    dat1.f = f;
    dat2.b[0] = dat1.b[3];
    dat2.b[1] = dat1.b[2];
    dat2.b[2] = dat1.b[1];
    dat2.b[3] = dat1.b[0];
    return dat2.f;
}


// Official definition of the nifti1 header.  Written by Bob Cox, SSCC, NIMH.

/**
 * @brief Alias for the original nifti1 header struct (file nifti1.h),
 *        written by Bob Cox, SSCC, NIMH.
 * 
 * @author Samuka Martins
 * @date Apr 20, 2018
 */
typedef struct nifti_1_header NIfTI1Header;



NIfTI1Header *_iftReadNIfTI1HeaderAndData(const char *path, char **data_out, iftImageDomain *dom_out, iftVoxelSize *voxel_sizes_out) {
    NIfTI1Header *header = iftAlloc(1, sizeof(NIfTI1Header));

    FILE *fp = fopen(path, "rb");
    if (fread(header, sizeof(char), 348, fp) != 348)
        iftError("Invalid NifTI1 (or ANALYZE 7.5 hdr) file. It can't read the 348 bytes from header", "_iftReadNIfTI1HeaderAndData");

    if ((header->datatype == NIFTI_TYPE_COMPLEX64) || (header->datatype == NIFTI_TYPE_FLOAT64) ||
        (header->datatype == NIFTI_TYPE_RGB24)  ||
        (header->dim[0] <= 0) || (header->dim[0] > 4)) {
        char msg[512] = "Error: Data format not supported or header is corrupted.\n" \
                        "Data that is NOT supported: Temporal series, Statistical images, or multiple volumes in a single file.";
        iftError(msg, "_iftReadNIfTI1HeaderAndData");
    }


    if ((header->dim[0] == 4) && (header->dim[4] != 1))
        iftError("Unsuported Image as Time Series", "_iftReadNIfTI1HeaderAndData");

    iftImageDomain dom = {1, 1, 1};
    iftVoxelSize voxel_sizes = {1.0, 1.0, 1.0};

    if (header->dim[0] >= 1) {
        dom.xsize = header->dim[1];
        voxel_sizes.dx = header->pixdim[1];
    }
    if (header->dim[0] >= 2) {
        dom.ysize = header->dim[2];
        voxel_sizes.dy = header->pixdim[2];
    }
    if (header->dim[0] >= 3) {
        dom.zsize = header->dim[3];
        voxel_sizes.dz = header->pixdim[3];
    }
    int nvoxels = dom.xsize * dom.ysize * dom.zsize;


    if (iftEndsWith(path, ".nii")) {
        // reading the raw image data - from byte #352
        if (header->vox_offset != 352)
            iftError("Invalid NifTI1 file. Voxel offset %d should be 352", "_iftReadNIfTI1HeaderAndData", header->vox_offset);
        fseek(fp, 4, SEEK_CUR); // we previously read 348 bytes, then move the file pointer in 4 bytes ahead
    }
    else if (iftEndsWith(path, ".hdr")) {
        char *hdr_base = iftBasename(path);
        char *img_path = iftConcatStrings(2, hdr_base, ".img");

        if (!iftFileExists(img_path))
            iftError("Image (data) file %s does not exist", "_iftReadNIfTI1HeaderAndData", img_path);
        fclose(fp);

        // opening the image (data) file
        fp = fopen(img_path, "rb");

        iftFree(hdr_base);
        iftFree(img_path);
    }
    else iftError("Invalid file format: %s. Try .nii or .hdr", "_iftReadNIfTI1HeaderAndData", iftFileExt(path));


    int n_bytes_per_voxel = header->bitpix / 8;
    int total_bytes = nvoxels * n_bytes_per_voxel;

    char *data = iftAlloc(total_bytes, sizeof(char)); // void elements has 1 byte
    if (fread(data, sizeof(char), total_bytes, fp) != total_bytes)
        iftError("Invalid total number of image voxel data in NifTI1 file", "_iftReadNIfTI1HeaderAndData");

    fclose(fp);

    *data_out = data;
    *dom_out = dom;
    *voxel_sizes_out = voxel_sizes;

    return header;
}


NIfTI1Header *_iftReadNIfTI1HeaderAndDataFromGZip(const char *path, char **data_out, iftImageDomain *dom_out, iftVoxelSize *voxel_sizes_out) {
    NIfTI1Header *header = iftAlloc(1, sizeof(NIfTI1Header));

    iftGZipFile gzip = iftGZipOpen(path, "rb", true);

    if (iftGZipRead(header, sizeof(char), 348, gzip) != 348)
        iftError("Invalid NifTI1 (or ANALYZE 7.5 hdr) file. It can't read the 348 bytes from header", "_iftReadNIfTI1HeaderAndDataFromGZip");

    if ((header->datatype == NIFTI_TYPE_COMPLEX64) || (header->datatype == NIFTI_TYPE_RGB24) ||
        (header->datatype == NIFTI_TYPE_RGB24) ||
        (header->dim[0] <= 0) || (header->dim[0] > 4)) {
        char msg[512] = "Error: Data format not supported, or header is corrupted.\n" \
                        "Data that is NOT supported: complex, RGB, temporal series and statistical images.";
        iftError(msg, "_iftReadNIfTI1HeaderAndDataFromGZip");
    }

    if ((header->dim[0] == 4) && (header->dim[4] != 1))
        iftError("Unsuported Image as Time Series", "_iftReadNIfTI1HeaderAndData");

    iftImageDomain dom = {1, 1, 1};
    iftVoxelSize voxel_sizes = {1.0, 1.0, 1.0};

    if (header->dim[0] >= 1) {
        dom.xsize = header->dim[1];
        voxel_sizes.dx = header->pixdim[1];
    }
    if (header->dim[0] >= 2) {
        dom.ysize = header->dim[2];
        voxel_sizes.dy = header->pixdim[2];
    }
    if (header->dim[0] >= 3) {
        dom.zsize = header->dim[3];
        voxel_sizes.dz = header->pixdim[3];
    }
    int nvoxels = dom.xsize * dom.ysize * dom.zsize;

    // reading the raw image data - from byte #352
    if (header->vox_offset != 352)
        iftError("Invalid NifTI1 file. Voxel offset %d should be 352", "_iftReadNIfTI1HeaderAndDataFromGZip", header->vox_offset);
    iftGZipSeek(gzip, 4, SEEK_CUR); // we previously read 348 bytes, then move the file pointer in 4 bytes ahead

    int n_bytes_per_voxel = header->bitpix / 8;
    int total_bytes = nvoxels * n_bytes_per_voxel;

    char *data = iftAlloc(nvoxels, n_bytes_per_voxel); // void elements has 1 byte
    if (iftGZipRead(data, sizeof(char), total_bytes, gzip) != total_bytes)
        iftError("Invalid total number of image voxel data in NifTI1 file", "_iftReadNIfTI1HeaderAndDataFromGZip");

    iftGZipClose(&gzip);

    *data_out = data;
    *dom_out = dom;
    *voxel_sizes_out = voxel_sizes;

    return header;
}


NIfTI1Header *_iftCreateNIfTI1HeaderFromImage(const iftImage *img, bool is_analyze) {
    NIfTI1Header *header = iftAlloc(1, sizeof(NIfTI1Header)); // put zeros in all fields
    header->sizeof_hdr = 348;
    header->scl_slope = 1.0;
    header->scl_inter = 0.0;

    // for 2D images
    header->dim[0] = 2;
    header->dim[1] = img->xsize;
    header->dim[2] = img->ysize;
    header->pixdim[1] = img->dx;
    header->pixdim[2] = img->dy;

    if (iftIs3DImage(img)) {
        header->dim[0] = 3;
        header->dim[3] = img->zsize;
        header->pixdim[3] = img->dz;
    }

    int depth   = iftImageDepth(img);
    int img_min = iftMinimumValue(img);

    if (depth <= 8) {
        header->bitpix = 8;
        header->datatype = (img_min < 0) ? NIFTI_TYPE_INT8 : NIFTI_TYPE_UINT8;
    }
    else if (depth <= 16) {
        header->bitpix = 16;
        header->datatype = (img_min < 0) ? NIFTI_TYPE_INT16 : NIFTI_TYPE_UINT16;
    }
    else if (depth <= 32) {
        // libIFT does not support unsigned integer NIFTI_TYPE_UINT32
        header->bitpix = 32;
        header->datatype = NIFTI_TYPE_INT32;
    }

    // in order to use the method 3 (full affine matrix for image orientation) as shown on https://brainder.org/2012/09/23/the-nifti-file-format
    header->qform_code = 0;
    header->sform_code = 1;
    header->scl_slope = 1.0;
    header->scl_inter = 0.0;

    // Rotation matrix to transform the Voxel Coordinates VC (i,j,k) System to the World Coordinate WC (x,y,z) system of NIfTI format
    // The orientation on WC is always RAS+, ie, +x Left to Right, +y Posterior to Anterior, +z Inferior to Superior
    // The default orientation on VC of libIFT is LPS+, then to transform LPS+ to RAS+ we need the
    // rotation matrix:
    //   i   j  k
    // [-1,  0, 0, 0] 
    // [ 0, -1, 0, 0] 
    // [ 0,  0, 1, 0]
    header->srow_x[0] = -1;
    header->srow_x[1] = 0;
    header->srow_x[2] = 0;
    header->srow_x[3] = 0;

    header->srow_y[0] = 0;
    header->srow_y[1] = -1;
    header->srow_y[2] = 0;
    header->srow_y[3] = 0;

    header->srow_z[0] = 0;
    header->srow_z[1] = 0;
    header->srow_z[2] = 1;
    header->srow_z[3] = 0;

    if (is_analyze) {
        header->vox_offset = 0; // analyze offset
        sprintf(header->magic, "ni1");    
    }
    else {
        header->vox_offset = 352; // nifti offset
        sprintf(header->magic, "n+1");
    }

    return header;
}


/**
 * @brief Computes the determinant of a 3X3 float matrix.
 * @author Samuel Martins
 * @date Apr 22, 2018
 */
float _iftMat33Determinant(float R[3][3]) {
    double r11, r12, r13, r21, r22, r23, r31, r32, r33 ;
                                                   /*  INPUT MATRIX:  */
    r11 = R[0][0]; r12 = R[0][1]; r13 = R[0][2];  /* [ r11 r12 r13 ] */
    r21 = R[1][0]; r22 = R[1][1]; r23 = R[1][2];  /* [ r21 r22 r23 ] */
    r31 = R[2][0]; r32 = R[2][1]; r33 = R[2][2];  /* [ r31 r32 r33 ] */

    return (r11*r22*r33) - (r11*r32*r23) - (r21*r12*r33) +
           (r21*r32*r13) + (r31*r12*r23) - (r31*r22*r13);
}

/**
 * @brief Multiply (C = AB) the 3X3 float matrix A with B, storing the results in the 3x3 float matrix C.
 * @author Samuel Martins
 * @date Apr 22, 2018
 */
void _iftMat33Multiply(float A[3][3], float B[3][3], float C[3][3]) {
    for (int i = 0; i < 3; i++ )
        for(int j = 0; j < 3; j++)
            C[i][j] =  A[i][0] * B[0][j] +
                       A[i][1] * B[1][j] +
                       A[i][2] * B[2][j];
}



/**
 * @brief Builds a NIfTI Rotation Matrix (its basis vectors) from the quatern coefficients in the nifti header.
 * 
 * @note This code is a slight adaption from the function <nifti_quatern_to_mat44> from
 * the official nifti reader code nifti1_io.c.
 * @note It is the method 2 shown on https://brainder.org/2012/09/23/the-nifti-file-format
 * 
 * @param header NIfTI header.
 * @param vec_i Reference to store the basis vector i from the rotation matrix.
 * @param vec_j Reference to store the basis vector j from the rotation matrix.
 * @param vec_k Reference to store the basis vector k from the rotation matrix.
 * 
 * @author Samuka Martins
 * @date Sep 13, 2018
 */
void _iftQuaternToNIfTIRotationMatrix(const NIfTI1Header *header, iftVector *vec_i, iftVector *vec_j, iftVector *vec_k) {
    double b = iftFixedFloat(header->quatern_b);
    double c = iftFixedFloat(header->quatern_c);
    double d = iftFixedFloat(header->quatern_d);

    double qfac = (header->pixdim[0] < 0.0) ? -1.0 : 1.0 ;  /* left-handedness? */
    double dx = header->pixdim[1];
    double dy = header->pixdim[2];
    double dz = header->pixdim[3];

    // compute 'a' parameter from b, c, d
    double a = 1.0l - (b*b + c*c + d*d) ;
    // special case
    if (a < 1.e-7l) {
        a = 1.0l / sqrt(b*b + c*c + d*d) ;

        // normalize (b,c,d) vector
        b *= a;
        c *= a;
        d *= a;

        a = 0.0l;  // a = 0 ==> 180 degree rotation
    } else a = sqrt(a);  // angle = 2*arccos(a)

    // load rotation matrix, including scaling factors for voxel sizes
    // make sure are positive
    double xd = (dx > 0.0) ? dx : 1.0l;
    double yd = (dy > 0.0) ? dy : 1.0l;
    double zd = (dz > 0.0) ? dz : 1.0l;

    if (qfac < 0.0)
        zd = -zd; // left handedness?

    vec_i->x = (a*a + b*b - c*c - d*d) * xd;
    vec_j->x = 2.0l * (b*c - a*d) * yd;
    vec_k->x = 2.0l * (b*d + a*c) * zd;

    vec_i->y = 2.0l * (b*c + a*d) * xd;
    vec_j->y = (a*a + c*c - b*b - d*d) * yd;
    vec_k->y = 2.0l * (c*d - a*b) * zd;

    vec_i->z = 2.0l * (b*d - a*c) * xd;
    vec_j->z = 2.0l * (c*d + a*b) * yd;
    vec_k->z = (a*a + d*d - c*c - b*b) * zd;

    // we ignore the offset parameters from the homogeneous coordinates, since they are not necessary
    // to find the orientation. They are only used to transform the voxel coordinates into world coordinates.
}


/**
 * @brief Find the orientation in which the image was stored based on Rotation matrix of the nifti header.
 * 
 * NifTI format works with two coordinate spaces: Voxel Coordinates VC (i,j,k) and World Coordinates WC (x,y,z) 
 * By default, WC always has the orientation RAS+, ie.
 * axis x: Left to Right 
 * axis y: Posterior to Anterior 
 * axis z: Inferior to Superior.
 * 
 * However, the image may have been stored on disk with another orientation, so that its orientation on
 * VC is different from WC.
 * NIfTI file stores an affine matrix to map the voxels from VC to WC. This matrix is 4X4 and consists of:
 * [srow_x[0], srow_x[1], srow_x[2], srow_x[3]] 
 * [srow_y[0], srow_y[1], srow_y[2], srow_y[3]] 
 * [srow_z[0], srow_z[1], srow_z[2], srow_z[3]]
 * [    0    ,     0    ,     0    ,     1    ]
 * 
 * If an exception or invalid rotation matrix is found, the default orientation of libIFT (LPS+) is returned.
 *
 * 
 * @param header NIfTI header.
 * @param x_axis_orient [description]
 * @param y_axis_orient [description]
 * @param z_axis_orient [description]
 */
void _iftFindImageOrientationOnNIfTI(const NIfTI1Header *header, ift3DAxisOrientation *x_axis_orient, ift3DAxisOrientation *y_axis_orient, ift3DAxisOrientation *z_axis_orient) {
    // Let VC be the Voxel Coordinate system (i,j,k) and WC be the World Coordinate system (x,y,z)
    // Let A_vc be a vector/matrix on VC system
    // Let A_wc be a vector/matrix on WC system

    // default orientation used by libIFT
    *x_axis_orient = IFT_R2L;
    *y_axis_orient = IFT_A2P;
    *z_axis_orient = IFT_I2S;

    double epsilon = 1e-4;

    // Rotation matrix (with scaling factors) on WC system extracted from the affine transformation matrix.
    // It is the rotation matrix that transforms a point on VC (i,j,k) to WC (x,y,k).
    // It transforms the image orientation of VC to RAS+ on WC.
    // It is equivalent to the basis from VC described on WC
    iftVector vec_i;
    iftVector vec_j;
    iftVector vec_k;

    // considering the full affine matrix from the srow parameters
    if (header->sform_code > 0) {
        vec_i.x = header->srow_x[0];
        vec_i.y = header->srow_y[0];
        vec_i.z = header->srow_z[0];

        vec_j.x = header->srow_x[1];
        vec_j.y = header->srow_y[1];
        vec_j.z = header->srow_z[1];
        
        vec_k.x = header->srow_x[2];
        vec_k.y = header->srow_y[2];
        vec_k.z = header->srow_z[2];
    }
    // considering the rotation matrix from the quatern parameters
    else if (header->qform_code > 0)
        _iftQuaternToNIfTIRotationMatrix(header, &vec_i, &vec_j, &vec_k);
    else {
        iftWarning("There is no Orientation Information (qform_code <= 0 and sform_code <= 0)\n" \
                    "The standard libIFT's orientation (LPS+) will be adopted", "iftReadImageNIfTI");
            return;
    }

    // if all basis vectors are zero-vectors. Nothing to do here
    if (iftIsZeroVector(vec_i) && iftIsZeroVector(vec_j) && iftIsZeroVector(vec_k))
        return;

    // since the basis (rotation matrix) can also have scaling factors incorporated, and we
    // just want to find out its orientation, we will work with the unit vectors of each axis from the basis
    iftVector vec_i_norm = iftNormalizeVector(vec_i);
    iftVector vec_j_norm = iftNormalizeVector(vec_j);

    // It is expected that all basis vectors are orthogonal each other, but we check if it is true and
    // if not, we orthogonalize vec_j_norm on to vec_i_norm and normalize it.
    float inner_product = iftVectorInnerProduct(vec_i_norm, vec_j_norm);    
    if (fabs(inner_product) > epsilon)
        vec_j_norm = iftNormalizeVector(iftOrthogonalizeVector(vec_j_norm, vec_i_norm));
    
    // Normalize the basis vector vec_k. If it is a zero-vector, it will be the cross product vec_i x vec_j
    iftVector vec_k_norm;
    if (iftIsZeroVector(vec_k))
        vec_k_norm = (iftVector) iftVectorCrossProd(vec_i_norm, vec_j_norm); // the cross product between two unit vectors, it is a unit vector orthogonal to i and j
    else {
        vec_k_norm = iftNormalizeVector(vec_k);

        // check if it is necessary to orthogonalize the vector k_norm with vec_i_norm
        inner_product = iftVectorInnerProduct(vec_k_norm, vec_i_norm);    
        if (fabs(inner_product) > epsilon)
            vec_k_norm = iftNormalizeVector(iftOrthogonalizeVector(vec_k_norm, vec_i_norm));

        // check if it is necessary to orthogonalize the vector vec_k_norm with vec_i_norm
        inner_product = iftVectorInnerProduct(vec_k_norm, vec_k_norm);    
        if (fabs(inner_product) > epsilon)
            vec_k_norm = iftNormalizeVector(iftOrthogonalizeVector(vec_k_norm, vec_j_norm));
    }

    // At this point, Bvc_wc is the basis of VC denoted on WC, only using unit vectors, which consists of the rotation matrix
    // that transforms a point (i,j,k) from VC to (x,y,z) on WC system.
    float Bvc_wc[3][3];
    Bvc_wc[0][0] = vec_i_norm.x; Bvc_wc[0][1] = vec_j_norm.x; Bvc_wc[0][2] = vec_k_norm.x;
    Bvc_wc[1][0] = vec_i_norm.y; Bvc_wc[1][1] = vec_j_norm.y; Bvc_wc[1][2] = vec_k_norm.y;
    Bvc_wc[2][0] = vec_i_norm.z; Bvc_wc[2][1] = vec_j_norm.z; Bvc_wc[2][2] = vec_k_norm.z;

    /*
    The orientation of the patient on WC is always RAS+: +x left to Right, +y posterior to Anterior, +z inferior to Superior
    The Identity Matrix on WC corresponds to RAS+ orientation, ie. x-axis is (1,0,0), y-axis (0,1,0), and z-axis (0,0,1): these values on WC system.
    
    To figure out the image orientation on VC, we need to find a Rotation Matrix P_vc on VC that transformed
    to WC system using the basis Bvc_wc (rotation matrix) outputs the Identity Matrix I_wc on WC in the sense of Bvc_wc.
    That is, M_wc = P_vc * Bvc_wc, such as M_wc ~= Iwc_wc (identity matrix on WC system).
    In other words, P_vc corresponds to the inverse matrix of Bvc_wc, since: inv(Bvc_wc) * Bvc_wc = I_wc
    P_wc is the rotation matrix (basis) to transform a point from WC to VC.
    
    We are only interest in rotations of 0, 90, 270, and 360 degrees in each axis, ie, only values -1, 1 in the rotation matrix P_vc.
    We then check which of the possible rotation matrices permuting -1, 1 along the axis is nearest to I_wc.
    Since we are working only with matrix formed by unity vectors, the permutation matrix nearest
    to the I_wc is the one with largest trace, ie whose sum of the main diagonal is maximum, since
    the I_wc has the largest trace.
    */
    
    float max_trace = IFT_INFINITY_FLT_NEG;
    int i_best = 1, j_best = 2, k_best = 3;
    int i_best_val = 1, j_best_val = 1, k_best_val = 1;
    float P_vc[3][3];

    float detR = _iftMat33Determinant(Bvc_wc);

    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            if (i == j) continue;

            for (int k = 1; k <= 3; k++) {
                if ((i == k) || (j == k)) continue;

                // zero the permutation matrix
                P_vc[0][0] = P_vc[0][1] = P_vc[0][2] = P_vc[1][0] = P_vc[1][1] = P_vc[1][2] = P_vc[2][0] = P_vc[2][1] = P_vc[2][2] = 0.0;

                // i_val, j_val, k_val are -1 or 1 and go into rows #1, #2, #3
                for (int i_val = -1; i_val <= 1; i_val += 2) {
                    for (int j_val = -1; j_val <= 1; j_val += 2) {
                        for (int k_val = -1; k_val <= 1; k_val += 2) {
                            // i permutation for the rotation matrix
                            P_vc[0][i-1] = i_val;
                            P_vc[1][j-1] = j_val;
                            P_vc[2][k-1] = k_val;

                            float detP = _iftMat33Determinant(P_vc); // sign of permutation
                            if (detP * detR <= 0.0) continue;  // it doesn't match the sign of Bvc_wc

                            float M_wc[3][3];
                            _iftMat33Multiply(P_vc, Bvc_wc, M_wc);

                            float trace = M_wc[0][0] + M_wc[1][1] + M_wc[2][2]; // trace
                            if (trace > max_trace) {
                                max_trace = trace;
                                i_best = i;
                                j_best = j;
                                k_best = k;
                                i_best_val = i_val;
                                j_best_val = j_val;
                                k_best_val = k_val;
                            }
                        }
                    }
                }
            }
        }
    }

    /*
    At this point, i_best is 0 or 1 or 3; i_best_val is -1 or +1; etc.

    The matrix P corresponds to the best permutation approximation whose multiplication
    with the basis Bvc_wc outputs the nearest matrix to the Identity matrix from WC.
    In other words, P stores the orientation (only using -1 and +1) of the image along the
    axis (i, j, k) on VC that applied the basis/rotation matrix Bvc_wc output the RAS+ orientation on WC,
    which corresponds to the Identity Matrix.

    For example, the column #1 determines the way as the i-axis are on WC axes.
    If i_best is 1, then the value (-1 or 1) from i-axis (column #0) is set in row #2, which 
    indicates that i-axis is along the y-axis of WC.
    Since the orientation of y-axis on WC is +y Posterior to Inferior, the direction of the i-axis on VC
    is P2A (if i_best_val = +1) or A2P (if i_best_val = -1).

    Examples:
    RAS+ Orientation on WC
     R  A  S
     x  y  z
    [1, 0, 0] row #1
    [0, 1, 0] row #2
    [0, 0, 1] row #3


         i  j  k
    P = [1, 0, 0]
        [0, 1, 0]
        [0, 0, 1]
    the image has the same orientation on VC and WC (RAS+). No rotation/transformation is done.

          i   j  k
    P = [-1,  0, 0]
        [ 0, -1, 0]
        [ 0,  0, 1]
    The i-axis is opposite to x-axis, j-axis is opposite to y-axis, k-axis and z-axis are equal. 
    Then, the image has the orientation LPS+ on VC.

          i  j  k
    P = [ 0, 0, 1]
        [-1, 0, 0]
        [ 0, 1, 0]
    i-axis is along the y-axis in the opposite direction, then i-axis is Anterior to Posterior P+.
    j-axis is along the z-axis in the same direction, then j-axis is Inferior to Superior S+.
    k-axis is along x-axis in the same direction, then k-axis is Left to Right R+.
    Therefore, the image orientation on VC is PSR+.
     */
    
    switch (i_best * i_best_val) {
        case  1: *x_axis_orient = IFT_L2R ; break ;
        case -1: *x_axis_orient = IFT_R2L ; break ;
        case  2: *x_axis_orient = IFT_P2A ; break ;
        case -2: *x_axis_orient = IFT_A2P ; break ;
        case  3: *x_axis_orient = IFT_I2S ; break ;
        case -3: *x_axis_orient = IFT_S2I ; break ;
   }

    switch (j_best * j_best_val) {
        case  1: *y_axis_orient = IFT_L2R ; break ;
        case -1: *y_axis_orient = IFT_R2L ; break ;
        case  2: *y_axis_orient = IFT_P2A ; break ;
        case -2: *y_axis_orient = IFT_A2P ; break ;
        case  3: *y_axis_orient = IFT_I2S ; break ;
        case -3: *y_axis_orient = IFT_S2I ; break ;
   }

    switch (k_best * k_best_val) {
        case  1: *z_axis_orient = IFT_L2R ; break ;
        case -1: *z_axis_orient = IFT_R2L ; break ;
        case  2: *z_axis_orient = IFT_P2A ; break ;
        case -2: *z_axis_orient = IFT_A2P ; break ;
        case  3: *z_axis_orient = IFT_I2S ; break ;
        case -3: *z_axis_orient = IFT_S2I ; break ;
   }
}


/**
 * @brief Computes the 4X4 Transformation/Rotation matrix to reslice an image to libIFT default orientation: LPS+
 * @param  x_axis_orient Orientation on x-axis.
 * @param  y_axis_orient Orientation on y-axis.
 * @param  z_axis_orient Orientation on z-axis.
 * @return               4X4 Transformation matrix.
 */
iftMatrix *_iftRotationMatrixToIFTOrientation(ift3DAxisOrientation x_axis_orient,
                                              ift3DAxisOrientation y_axis_orient,
                                              ift3DAxisOrientation z_axis_orient) {
    iftMatrix *R = iftCreateMatrix(4, 4);

    // unit vector (direction) to LPS+ orientation
    // +x Right to Left
    // +y Anterior to Posterior
    // +z Inferior to Superior
    float R2L_trans[3] = { 1,  0,  0};
    float L2R_trans[3] = {-1,  0,  0};
    float A2P_trans[3] = { 0,  1,  0};
    float P2A_trans[3] = { 0, -1,  0};
    float I2S_trans[3] = { 0,  0,  1};
    float S2I_trans[3] = { 0,  0, -1};

    float *x_trans = NULL;
    switch (x_axis_orient) {
        case IFT_R2L: x_trans = R2L_trans; break;
        case IFT_L2R: x_trans = L2R_trans; break;
        case IFT_A2P: x_trans = A2P_trans; break;
        case IFT_P2A: x_trans = P2A_trans; break;
        case IFT_I2S: x_trans = I2S_trans; break;
        case IFT_S2I: x_trans = S2I_trans; break;
    }

    float *y_trans = NULL;
    switch (y_axis_orient) {
        case IFT_R2L: y_trans = R2L_trans; break;
        case IFT_L2R: y_trans = L2R_trans; break;
        case IFT_A2P: y_trans = A2P_trans; break;
        case IFT_P2A: y_trans = P2A_trans; break;
        case IFT_I2S: y_trans = I2S_trans; break;
        case IFT_S2I: y_trans = S2I_trans; break;
    }

    float *z_trans = NULL;
    switch (z_axis_orient) {
        case IFT_R2L: z_trans = R2L_trans; break;
        case IFT_L2R: z_trans = L2R_trans; break;
        case IFT_A2P: z_trans = A2P_trans; break;
        case IFT_P2A: z_trans = P2A_trans; break;
        case IFT_I2S: z_trans = I2S_trans; break;
        case IFT_S2I: z_trans = S2I_trans; break;
    }

    // R(col, row)
    iftMatrixElem(R, 0, 0) = x_trans[0];
    iftMatrixElem(R, 0, 1) = x_trans[1];
    iftMatrixElem(R, 0, 2) = x_trans[2];
    iftMatrixElem(R, 1, 0) = y_trans[0];
    iftMatrixElem(R, 1, 1) = y_trans[1];
    iftMatrixElem(R, 1, 2) = y_trans[2];
    iftMatrixElem(R, 2, 0) = z_trans[0];
    iftMatrixElem(R, 2, 1) = z_trans[1];
    iftMatrixElem(R, 2, 2) = z_trans[2];

    // no translation
    iftMatrixElem(R, 3, 0) = 0;
    iftMatrixElem(R, 3, 1) = 0;
    iftMatrixElem(R, 3, 2) = 0;

    // homogeneous coordinate values
    iftMatrixElem(R, 0, 3) = 0;
    iftMatrixElem(R, 1, 3) = 0;
    iftMatrixElem(R, 2, 3) = 0;
    iftMatrixElem(R, 3, 3) = 1;

    return R;
}


/********************** PUBLIC FUNCTIONS *************************/
void iftCopyGrayImageVoxel(const iftImage *src, iftVoxel u, iftImage *dst, iftVoxel v) {
    iftImgVoxelVal(dst, v) = iftImgVoxelVal(src, u);
}

void iftCopyColorImageVoxel(const iftImage *src, iftVoxel u, iftImage *dst, iftVoxel v) {
    iftImgVoxelVal(dst, v) = iftImgVoxelVal(src, u);
    iftImgVoxelCb(dst, v) = iftImgVoxelCb(src, u);
    iftImgVoxelCr(dst, v) = iftImgVoxelCr(src, u);
}



void iftVerifyImageDomains(const iftImage *img1, const iftImage *img2, const char *function) {
    if ((img1==NULL)||(img2==NULL)||(img1->xsize!=img2->xsize) || (img1->ysize!=img2->ysize) || (img1->zsize!=img2->zsize)) {
        iftError("Images with different domains:\n" \
                "img1 (xsize, ysize, zsize): (%d, %d, %d)\n" \
                "img2 (xsize, ysize, zsize): (%d, %d, %d)\n",
                 function,
                 img1->xsize, img1->ysize, img1->zsize, img2->xsize, img2->ysize, img2->zsize);
    }
}


iftImageDomain iftGetImageDomain(const iftImage *img) {
    iftImageDomain dom;
    dom.xsize = img->xsize;
    dom.ysize = img->ysize;
    dom.zsize = img->zsize;
    return dom;
}


void iftVerifyImages(const iftImage *img1, const iftImage *img2, const char *external_function) {
    iftVerifyImageDomains(img1, img2, external_function);
    if (!iftIsVoxelSizeEqual(img1, img2))
        iftError("Input Images have different Voxel Sizes\n" \
                 "Image 1 (%lf, %lf, %lf)\n" \
                 "Image 2 (%lf, %lf, %lf)", external_function,
                 img1->dx, img1->dy, img1->dz, img2->dx, img2->dy, img2->dz);
}


iftImage  *iftCreateImage(int xsize,int ysize,int zsize) {
    int *val = iftAllocIntArray(xsize*ysize*zsize);

    return iftCreateImageFromBuffer(xsize, ysize, zsize, val);
}


iftImage *iftCreateImageFromImage(const iftImage *src) {
    iftImage *out = NULL;

    if (src != NULL) {
        if (iftIsColorImage(src)) {
            out = iftCreateColorImage(src->xsize, src->ysize, src->zsize, iftImageDepth(src));
        }
        else {
            out = iftCreateImage(src->xsize, src->ysize, src->zsize);
        }
        iftCopyVoxelSize(src, out);
    }

    return out;
}



// Creates an image from a pre-allocated buffer whose size must be xsize*ysize*zsize
iftImage *iftCreateImageFromBuffer(int xsize, int ysize, int zsize, int *val) {
    iftImage *img = NULL;
    int      y, z, xysize;

    img = (iftImage *) iftAlloc(1, sizeof(iftImage));
    if (img == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateImage");
    }

    img->val   = val;
    img->Cb    = img->Cr = NULL;
    img->alpha = NULL;
    img->xsize = xsize;
    img->ysize = ysize;
    img->zsize = zsize;
    img->dx    = 1.0;
    img->dy    = 1.0;
    img->dz    = 1.0;
    img->tby   = iftAllocIntArray(ysize);
    img->tbz   = iftAllocIntArray(zsize);
    img->n     = xsize * ysize * zsize;

    if (img->val == NULL || img->tbz == NULL || img->tby == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateImage");
    }

    img->tby[0] = 0;
    for (y = 1; y < ysize; y++)
        img->tby[y] = img->tby[y - 1] + xsize;

    img->tbz[0] = 0;
    xysize = xsize * ysize;
    for (z = 1; z < zsize; z++)
        img->tbz[z] = img->tbz[z - 1] + xysize;

    return (img);
}

iftImage  *iftCreateColorImage(int xsize,int ysize,int zsize, int depth)
{
    iftImage *img=NULL;
    img = iftCreateImage(xsize, ysize, zsize);

    iftSetCbCr(img, (iftMaxImageRange(depth)+1)/2);

    return(img);
}


iftImageArray *iftCreateImageArray(size_t n) {
    iftImageArray *img_arr = (iftImageArray*) iftAlloc(1, sizeof(iftImageArray));
    img_arr->n             = n;
    img_arr->val           = (iftImage**) iftAlloc(n, sizeof(iftImage*));

    return img_arr;
}


void  iftDestroyPyImage(iftImage *img)
{
    iftDestroyImage(&img);
}

void iftDestroyImage(iftImage **img) {
    if(img != NULL) {
        iftImage *aux = *img;

        if (aux != NULL) {
            if (aux->val != NULL) iftFree(aux->val);
            if (aux->Cb != NULL) iftFree(aux->Cb);
            if (aux->Cr != NULL) iftFree(aux->Cr);
            if (aux->alpha != NULL) iftFree(aux->alpha);
            if (aux->tby != NULL) iftFree(aux->tby);
            if (aux->tbz != NULL) iftFree(aux->tbz);
            iftFree(aux);
            *img = NULL;
        }
    }
}


void iftDestroyImageArray(iftImageArray **img_arr) {

    if (img_arr != NULL && *img_arr != NULL) {
        iftImageArray *aux = *img_arr;

        for (size_t i = 0; i < aux->n; i++)
            iftDestroyImage(&aux->val[i]);
        iftFree(aux);
        *img_arr = NULL;
    }
}



iftImageArray *iftReadImageSet(const iftFileSet *img_set) {
    iftImageArray *img_arr = iftCreateImageArray(img_set->n);

    #pragma omp parallel for
    for (long i = 0; i < img_set->n; i++)
        img_arr->val[i] = iftReadImageByExt(img_set->files[i]->path);

    return img_arr;
}

iftIntMatrix *iftReadGrayImageSetAsIntMatrix(const iftFileSet *img_set) {
    iftImage *img0 = iftReadImageByExt(img_set->files[0]->path);
    iftIntMatrix *mat = iftCreateIntMatrix(img0->n, img_set->n);
    iftDestroyImage(&img0);

    #pragma omp parallel for
    for (int f = 0; f < img_set->n; f++) {
        iftImage *img = iftReadImageByExt(img_set->files[f]->path);
        
        for (int p = 0; p < img->n; p++) {
            iftMatrixElem(mat, p, f) = img->val[p];
        }

        iftDestroyImage(&img);
    }

    return mat;
}



void iftCopyCbCr(const iftImage *src, iftImage *dst) {
    if (iftIsColorImage(src)) {
        if (!iftIsDomainEqual(src, dst))
            iftError("Images must have the same domain", "iftCopyCbCr");

        if (dst->Cb == NULL){
            dst->Cb = iftAllocUShortArray(dst->n);
            dst->Cr = iftAllocUShortArray(dst->n);
        }

        #pragma omp parallel for schedule(auto)
        for (int p = 0; p < dst->n; p++) {
            dst->Cb[p] = src->Cb[p];
            dst->Cr[p] = src->Cr[p];
        }
    }
}


void    iftSetCbCr(iftImage *img, ushort value)
{
    int p;

    if (!iftIsColorImage(img)){
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
    }
    for (p=0; p < img->n; p++) {
        img->Cb[p] = value;
        img->Cr[p] = value;
    }
}

void  iftSetAlpha(iftImage *img, ushort value)
{
    int p;
    if(img->alpha == NULL){
        img->alpha = iftAllocUShortArray(img->n);
    }

    for (p=0; p < img->n; p++) {
        img->alpha[p] = value;
    }
}

void  iftSetCbCrAlpha(iftImage *img, ushort value)
{
    int p;
    if(img->alpha == NULL){
        img->alpha = iftAllocUShortArray(img->n);
    }
    if(img->Cb == NULL){
        img->Cb = iftAllocUShortArray(img->n);
    }
    if(img->Cr == NULL){
        img->Cr = iftAllocUShortArray(img->n);
    }

    for (p=0; p < img->n; p++) {
        img->alpha[p] = value;
        img->Cb[p] = value;
        img->Cr[p] = value;
    }
}



void    iftSetCb(iftImage *img, ushort value)
{
    int p;

    if (img->Cb == NULL){
        img->Cb = iftAllocUShortArray(img->n);
        for (p=0; p < img->n; p++) {
            img->Cb[p] = value;
        }
    }
}

void    iftSetCr(iftImage *img, ushort value)
{
    int p;

    if (img->Cr == NULL){
        img->Cr = iftAllocUShortArray(img->n);
        for (p=0; p < img->n; p++) {
            img->Cr[p] = value;
        }
    }
}


int iftMaximumValueInRegion(const iftImage *img, iftBoundingBox bb) {
    // checkers
    if (img == NULL)
        iftError("Image is NULL", "iftMaximumValueInRegion");
    if (!iftValidVoxel(img, bb.begin))
        iftError("Beginning voxel (%d, %d, %d) from Region (Bound. Box) is not in the Image Domain\n" \
                 "Img (xsize, ysize, zsize): (%d, %d, %d)", "iftMaximumValueInRegion",
                 bb.begin.x, bb.begin.y, bb.begin.z, img->xsize, img->ysize, img->zsize);
    if (!iftValidVoxel(img, bb.end))
        iftError("Ending voxel (%d, %d, %d) from Region (Bound. Box) is not in the Image Domain\n" \
                 "Img (xsize, ysize, zsize): (%d, %d, %d)", "iftMaximumValueInRegion",
                 bb.end.x, bb.end.y, bb.end.z, img->xsize, img->ysize, img->zsize);

    int img_max_val = IFT_INFINITY_INT_NEG;

    iftVoxel v;
    for (v.z = bb.begin.z; v.z <= bb.end.z; v.z++) {
        for (v.y = bb.begin.y; v.y <= bb.end.y; v.y++) {
            for (v.x = bb.begin.x; v.x <= bb.end.x; v.x++) {
                int p = iftGetVoxelIndex(img, v);
                if (img_max_val < img->val[p]) {
                    img_max_val = img->val[p];
                }
            }
        }
    }

    return img_max_val;
}


int iftMaximumValue(const iftImage *img) {
    if (img == NULL)
        iftError("Image is NULL", "iftMaximumValue");

    iftBoundingBox bb;
    bb.begin.x = bb.begin.y = bb.begin.z = 0;
    bb.end.x   = img->xsize-1;
    bb.end.y   = img->ysize-1;
    bb.end.z   = img->zsize-1;

    return iftMaximumValueInRegion(img, bb);
}



int iftMinimumValue(const iftImage *img) {
    int img_min_val = IFT_INFINITY_INT;

    for (int p = 0; p < img->n; p++)
        if (img_min_val > img->val[p])
            img_min_val = img->val[p];

    return img_min_val;
}


void iftMinMaxValueInRegion(const iftImage *img, const iftImage *region, int *min, int *max) {
    *min = IFT_INFINITY_INT;
    *max = IFT_INFINITY_INT_NEG;

    for (int p = 0; p < img->n; p++) {
        if (region->val[p]) {
            if (img->val[p] < *min)
                *min = img->val[p];
            else if (img->val[p] > *max)
                *max = img->val[p];
        }
    }
}

void iftMinMaxValues(const iftImage *img, int *min, int *max) {
    *min = *max = img->val[0];

    for (int p = 1; p < img->n; p++) {
        if (img->val[p] < *min)
            *min = img->val[p];
        else if (img->val[p] > *max)
            *max = img->val[p];
    }
}


int iftMaximumCb(const iftImage *img)
{
    int p, max;

    if (img->Cb == NULL)
        iftError("Image is grayscale", "iftMaximumCb");

    max = IFT_INFINITY_INT_NEG;
    for (p=0; p < img->n; p++)
        if (max < img->Cb[p])
            max = img->Cb[p];

    return(max);
}

int iftMaximumCr(const iftImage *img)
{
    int p, max;

    if (img->Cr == NULL)
        iftError("Image is grayscale", "iftMaximumCr");

    max = IFT_INFINITY_INT_NEG;
    for (p=0; p < img->n; p++)
        if (max < img->Cr[p])
            max = img->Cr[p];

    return(max);
}


iftIntArray *iftMaximumObjectValues(const iftImage *img, const iftImage *label_img, iftIntArray **voxels_idxs_out) {
    iftVerifyImages(img, label_img, "iftMaximumObjectValues");

    int n_objects = iftMaximumValue(label_img);

    iftIntArray *max_obj_vals = iftIntRepeat(IFT_INFINITY_INT_NEG, n_objects + 1);
    iftIntArray *voxels_idxs = iftIntRepeat(-1, n_objects + 1);

    for (int p = 0; p < img->n; p++) {
        int label = label_img->val[p];

        if (img->val[p] > max_obj_vals->val[label]) {
            max_obj_vals->val[label] = img->val[p];
            voxels_idxs->val[label] = p;
        }
    }

    if (voxels_idxs_out == NULL)
        iftDestroyIntArray(&voxels_idxs);
    else *voxels_idxs_out = voxels_idxs;

    return max_obj_vals;
}



int iftMinimumCb(const iftImage *img)
{
    int p, min;

    if (img->Cb == NULL)
        iftError("Image is grayscale", "iftMinimumCb");

    min = IFT_INFINITY_INT;
    for (p=0; p < img->n; p++)
        if (min > img->Cb[p])
            min = img->Cb[p];

    return(min);
}

int iftMinimumCr(const iftImage *img)
{
    int p, min;

    if (img->Cr == NULL)
        iftError("Image is grayscale", "iftMinimumCr");

    min = IFT_INFINITY_INT;
    for (p=0; p < img->n; p++)
        if (min > img->Cr[p])
            min = img->Cr[p];

    return(min);
}


int iftMaximumValueInAdjacency(const iftImage *img, int p, iftAdjRel *A) {
    int i, max_val = IFT_INFINITY_INT_NEG;
    int q;
    iftVoxel u, v;

    u = iftGetVoxelCoord(img, p);

    for(i = 0; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A, u, i);

        if(iftValidVoxel(img, v)) {
            q = iftGetVoxelIndex(img, v);
            max_val = iftMax(max_val, img->val[q]);
        }
    }

    return max_val;
}

int iftMinimumValueInAdjacency(const iftImage *img, int p, iftAdjRel *A) {
    int i, min_val = IFT_INFINITY_INT;
    int q;
    iftVoxel u, v;

    u = iftGetVoxelCoord(img, p);

    for(i = 0; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A, u, i);

        if(iftValidVoxel(img, v)) {
            q = iftGetVoxelIndex(img, v);
            min_val = iftMin(min_val, img->val[q]);
        }
    }

    return min_val;
}

int iftMedianValueInAdjacency(iftImage *img, int p, iftAdjRel *A) {
    int i, i0, median = IFT_INFINITY_INT;
    int q, *values = NULL, *index = NULL;
    iftVoxel u, v;

    u = iftGetVoxelCoord(img, p);
    values = iftAllocIntArray(A->n);
    index = iftAllocIntArray(A->n);

    i0 = 0;
    for(i = 0; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A, u, i);

        if(iftValidVoxel(img, v)) {
            q = iftGetVoxelIndex(img, v);
            values[i0] = img->val[q];
            index[i0]  = i0;
            i0++;
        }
    }

    // Sorting image values
    iftQuickSort(values, index, 0, i0-1, IFT_INCREASING);

    median = values[i0/2];

    iftFree(values);
    iftFree(index);

    return median;
}


char *iftImageExt(const iftImage *img) {
    char *ext = NULL;

    if (iftIs3DImage(img)) {
        ext = iftCopyString(".scn");
    }
    else {
        if (iftIsColorImage(img))
            ext = iftCopyString(".ppm");
        else ext = iftCopyString(".pgm");
    }

    return ext;
}


iftImage *iftReadImage(const char *format, ...) {
    iftImage *img    = NULL;
    FILE     *fp     = NULL;
    uchar    *data8  = NULL;
    ushort   *data16 = NULL;
    int      *data32 = NULL;
    char     type[10];
    int      p, v, xsize, ysize, zsize;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "rb");

    if (fp == NULL) {
        iftError("Cannot open file: \"%s\"", "iftReadImage", filename);
    }
    if (fscanf(fp, "%s\n", type) != 1)
        iftError("Reading error: Image type", "iftReadImage");

    if (iftCompareStrings(type, "SCN")) {

        //iftSkipComments(fp);

        if (fscanf(fp, "%d %d %d\n", &xsize, &ysize, &zsize) != 3)
            iftError("Reading error: Image resolution/size", "iftReadImage");
        img = iftCreateImage(xsize, ysize, zsize);
        if (fscanf(fp, "%f %f %f\n", &img->dx, &img->dy, &img->dz) != 3) {
            iftError("Reading error: Pixel/Voxel size", "iftReadImage");
        }
        if (fscanf(fp, "%d", &v) != 1)
            iftError("Reading error", "iftReadImage");

        while (fgetc(fp) != '\n');

        if (v == 8) {
            data8 = iftAllocUCharArray(img->n);
            if (fread(data8, sizeof(uchar), img->n, fp) != img->n)
                iftError("Reading error", "iftReadImage");
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data8[p];
            iftFree(data8);
        } else if (v == 16) {
            data16 = iftAllocUShortArray(img->n);
            if (fread(data16, sizeof(ushort), img->n, fp) != img->n)
                iftError("Reading error 16 bits", "iftReadImage");
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data16[p];
            iftFree(data16);
        } else if (v == 32) {
            data32 = iftAllocIntArray(img->n);
            if (fread(data32, sizeof(int), img->n, fp) != img->n)
                iftError("Reading error", "iftReadImage");
            for (p = 0; p < img->n; p++)
                img->val[p] = data32[p];
            iftFree(data32);
        } else {
            iftError("Input scene must be 8, 16, or 32 bit", "iftReadImage");
        }
    } else {
        iftError("Invalid file type", "iftReadImage");
    }

    fclose(fp);
    return (img);
}


iftImage *iftReadImageNIfTI(const char *format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    if (!iftRegexMatch(filename, "^.+\\.(nii|nii.gz|hdr)$"))
        iftError("Invalid file format: %s. Try .nii, .nii.gz, or .hdr", "iftReadImageNIfTI", iftFileExt(filename));

    NIfTI1Header *header = NULL;
    iftImageDomain dom;
    iftVoxelSize voxel_sizes;
    char *data = NULL;

    if ((iftEndsWith(filename, ".nii")) || (iftEndsWith(filename, ".hdr")))
        header = _iftReadNIfTI1HeaderAndData(filename, &data, &dom, &voxel_sizes);
    else if (iftEndsWith(filename, ".nii.gz"))
        header = _iftReadNIfTI1HeaderAndDataFromGZip(filename, &data, &dom, &voxel_sizes);
    else iftError("Invalid file format: %s. Try .nii or .hdr", "iftReadImageNIfTI", iftFileExt(filename));

    int nvoxels = dom.xsize * dom.ysize * dom.zsize;

    float *val = iftAllocFloatArray(nvoxels);

    if (header->datatype == NIFTI_TYPE_UINT8) {
        #pragma omp parallel for
        for (int p = 0; p < nvoxels; p++)
            val[p] = (float) ((uchar *) data)[p];
    }
    else if (header->datatype == NIFTI_TYPE_INT8) {
        #pragma omp parallel for
        for (int p = 0; p < nvoxels; p++)
            val[p] = (float) ((char *) data)[p];
    }
    else if (header->datatype == NIFTI_TYPE_UINT16) {
        #pragma omp parallel for
        for (int p = 0; p < nvoxels; p++)
            val[p] = (float) ((ushort *) data)[p];
    }
    else if (header->datatype == NIFTI_TYPE_INT16) {
        #pragma omp parallel for
        for (int p = 0; p < nvoxels; p++)
            val[p] = (float) ((short *) data)[p];
    }
    else if (header->datatype == NIFTI_TYPE_INT32) {
        #pragma omp parallel for
        for (int p = 0; p < nvoxels; p++)
            val[p] = (float) ((int *) data)[p];
    }
    else if (header->datatype == NIFTI_TYPE_UINT32) {
        #pragma omp parallel for
        for (int p = 0; p < nvoxels; p++)
            val[p] = (float) ((uint *) data)[p];
    }
    else if (header->datatype == NIFTI_TYPE_FLOAT32) {
        #pragma omp parallel for
        for (int p = 0; p < nvoxels; p++)
            val[p] = ((float *) data)[p];
    } else if (header->datatype == NIFTI_TYPE_FLOAT64) {
        #pragma omp parallel for
        for (int p = 0; p < nvoxels; p++) {
            double d = ((double *) data)[p];
            val[p] = (float) d;
        }
    }

    // apply the "descompression"
    if (((header->scl_slope != 1.0) && (header->scl_slope != 0.0)) || (header->scl_inter != 0.0)) {
        #pragma omp parallel for
        for (int p = 0; p < nvoxels; p++) {
            val[p] = (val[p] * header->scl_slope) + header->scl_inter;
        }
    }

    float min = IFT_INFINITY_FLT;
    float max = IFT_INFINITY_FLT_NEG;
    for (int p = 0; p < nvoxels; p++) {
        if (val[p] < min)
            min = val[p];
        if (val[p] > max)
            max = val[p];
    }

    // shifting image value range
    if (min < 0.0)
        iftWarning("- Image with negative values within [%f, %f]\n", "iftReadImageNIfTI", min, max);

    // normalization: if the original image were stored as float numbers, we normalize them to 16 bits
    // otherwise (integer data), we only cast the values, assuming that this won't loose much information
    float factor = 1.0;
    float bias = 0.0;
    if ((header->datatype == NIFTI_TYPE_FLOAT32) || (header->datatype == NIFTI_TYPE_FLOAT64)) {
        iftWarning("- The original image data is float or double. Normalizing it to integer data of 12 bits\n", "iftReadImageNIfTI");
        // scaling to fit in range [-2048, 2047]
        if (min < 0.0) {
            if (fabs(min) <= fabs(max))
                factor = 2047 / iftMax(1, max);
            else factor = 2048 / iftMax(1, fabs(min));
        }
        // scaling to fit in range [0, 4095]
        else factor = 4095.0 / iftMax(1, max);
    }

    iftImage *img = iftCreateImage(dom.xsize, dom.ysize, dom.zsize);
    img->dx = voxel_sizes.dx;
    img->dy = voxel_sizes.dy;
    img->dz = voxel_sizes.dz;

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++)
        img->val[p] = iftRound((val[p] * factor) - bias);

    // reorienting/reslacing the image to default libIFT orientation: LPS+
    ift3DAxisOrientation x_axis_orient, y_axis_orient, z_axis_orient;
    _iftFindImageOrientationOnNIfTI(header, &x_axis_orient, &y_axis_orient, &z_axis_orient);

    if ((x_axis_orient != IFT_R2L) || (y_axis_orient != IFT_A2P) || (z_axis_orient != IFT_I2S)) {
        iftMatrix *R = _iftRotationMatrixToIFTOrientation(x_axis_orient, y_axis_orient, z_axis_orient);
        
        iftImage *reslice_img = iftResliceImageByTransMatrix(img, R);
        iftDestroyImage(&img);
        iftDestroyMatrix(&R);

        img = reslice_img;
    }

    iftFree(data);
    iftFree(val);
    iftFree(header);

    return img;
}


iftImage *iftReadImageAnalyze(const char *format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    return iftReadImageNIfTI(filename);
}


iftImage *iftReadImageGZip(const char *format, ...) {
    iftImage *img    = NULL;
    iftGZipFile fp   = NULL;
    uchar    *data8  = NULL;
    ushort   *data16 = NULL;
    int      *data32 = NULL;
    char     type[10], header[IFT_STR_DEFAULT_SIZE];
    int      p, v, xsize, ysize, zsize;
    float    dx, dy, dz;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = iftGZipOpen(filename, "rb", 1);

    if (fp == NULL) {
        iftError("Cannot open file: \"%s\"", "iftReadImageGZip", filename);
    }

    // Reading the header (everything until a \n is found)
    iftGZipGets(header, IFT_STR_DEFAULT_SIZE, fp);

    // Reading all info from the header
    if (sscanf(header, "%s %d %d %d %f %f %f %d", type, &xsize, &ysize, &zsize, &dx, &dy, &dz, &v) != 8)
        iftError("Reading error! Image header does not match what is expected: %s", "iftReadImageGZip", header);

    if (iftCompareStrings(type, "SCN")) {
        size_t nread;

        img = iftCreateImage(xsize, ysize, zsize);
        img->dx = dx;
        img->dy = dy;
        img->dz = dz;

        if (v == 8) {
            data8 = iftAllocUCharArray(img->n);
            if ((nread = iftGZipRead(data8, sizeof(uchar), img->n, fp)) != img->n)
                iftError("Reading error 8 bits. %lu voxels read when %d were expected", "iftReadImageGZip", nread,
                         img->n);
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data8[p];
            iftFree(data8);
        } else if (v == 16) {
            data16 = iftAllocUShortArray(img->n);
            if ((nread = iftGZipRead(data16, sizeof(ushort), img->n, fp)) != img->n)
                iftError("Reading error 16 bits. %lu voxels read when %d were expected", "iftReadImageGZip", nread,
                         img->n);
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data16[p];
            iftFree(data16);
        } else if (v == 32) {
            data32 = iftAllocIntArray(img->n);
            if ((nread = iftGZipRead(data32, sizeof(int), img->n, fp)) != img->n)
                iftError("Reading error 32 bits. %lu voxels read when %d were expected", "iftReadImageGZip", nread,
                         img->n);
            for (p = 0; p < img->n; p++)
                img->val[p] = data32[p];
            iftFree(data32);
        } else {
            iftError("Input scene must be 8, 16, or 32 bit", "iftReadImageGZip");
        }
    } else {
        iftError("Invalid file type", "iftReadImageGZip");
    }

    iftGZipClose(&fp);

    return (img);
}





void iftWritePngImageAux(const char *file_name, png_bytep *row_pointers, int width, int height, int bit_depth,
                         int color_type) {

    /* create file */
    FILE *fp = fopen(file_name, "wb");
    if (!fp)
        iftError("Internal Error: File %s could not be opened for writing", "iftWritePngImageAux", file_name);


    /* initialize stuff */
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
        iftError("Internal Error: png_create_write_struct failed", "iftWriteImagePNG");

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
        iftError("Internal Error: png_create_info_struct failed", "iftWriteImagePNG");

    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during init_io", "iftWriteImagePNG");

    png_init_io(png_ptr, fp);


    /* write header */
    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during writing header", "iftWriteImagePNG");

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);


    /* write bytes */
    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during writing bytes", "iftWriteImagePNG");

    png_write_image(png_ptr, row_pointers);


    /* end write */
    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during end of write", "iftWriteImagePNG");

    png_write_end(png_ptr, NULL);

    png_destroy_write_struct(&png_ptr, &info_ptr);

    /* cleanup heap allocation */
    for (int y=0; y<height; y++)
        iftFree(row_pointers[y]);
    iftFree(row_pointers);

    fclose(fp);
}


uchar iftImageDepth(const iftImage *img) {
    int img_min, img_max;
    iftMinMaxValues(img, &img_min, &img_max);
    
    long max_range;

    if (img_min >= 0)
        max_range = iftNormalizationValue(img_max) + 1;
    else
        max_range = iftNormalizationValue(img_max - img_min) + 1;
    
    return (uchar) iftLog(max_range, 2);
}

bool iftIsValidFormat(const char *filename){
    const char *ext = iftFileExt(filename);
    if (iftCompareStrings(ext, ".zscn") || iftCompareStrings(ext, ".scn.gz") ||
        iftCompareStrings(ext, ".scn") || iftCompareStrings(ext, ".nii") ||
        iftCompareStrings(ext, ".nii.gz") || iftCompareStrings(ext, ".hdr") ||
        iftCompareStrings(ext, ".png") || iftCompareStrings(ext, ".pgm") ||
        iftCompareStrings(ext, ".ppm") || iftCompareStrings(ext, ".mimg"))
        return true;
    return false;
}

bool iftIsValid3DFormat(const char *filename){
    const char *ext = iftFileExt(filename);
    if (iftCompareStrings(ext, ".zscn") || iftCompareStrings(ext, ".scn.gz") ||
        iftCompareStrings(ext, ".scn") || iftCompareStrings(ext, ".nii") ||
        iftCompareStrings(ext, ".nii.gz") || iftCompareStrings(ext, ".hdr"))
        return true;
    return false;
}



void iftWriteImagePNG(const iftImage* img, const char* format, ...) {

    png_bytep *row_pointers;
    int width, height, depth, byteshift;

    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    width = img->xsize;
    height = img->ysize;
    png_byte color_type;
    depth = iftImageDepth(img);

    if(depth<=8) {
        depth = 8;
    } else {
        depth = 16;
    }

    byteshift = depth/8;
    //int offset = depth==16?1:0;//to read the second byte first in cases of 16bit images

    size_t numberOfChannels=1;
    if(iftIsColorImage(img)){
        if(img->alpha == NULL){
            numberOfChannels = 3;//RGB
            color_type = PNG_COLOR_TYPE_RGB;
        }else{
            numberOfChannels = 4;//RGB_ALPHA
            color_type = PNG_COLOR_TYPE_RGB_ALPHA;
        }
    }else{
        if(img->alpha == NULL){
            numberOfChannels = 1;//GRAY
            color_type = PNG_COLOR_TYPE_GRAY;
        }else{
            numberOfChannels = 2;//GRAY_ALPHA
            color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
        }
    }

    //size_t pixel_size = (iftIsColorImage(img)?3:1 ) * byteshift;
    row_pointers = (png_bytep*) iftAlloc(height, sizeof(png_bytep));
    for (int y=0; y<height; y++)
        row_pointers[y] = (png_byte*) iftAlloc(width, numberOfChannels*byteshift);

    if(color_type == PNG_COLOR_TYPE_GRAY){
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                ptr[0] = img->val[p] & 0xFF;//get first byte

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[1] = ptr[0];
                    ptr[0] = (img->val[p]>>8) & 0xFF;//get second byte
                }

                p++;
            }
        }
    }else if(color_type == PNG_COLOR_TYPE_GRAY_ALPHA){
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                if(depth==8){
                    ptr[0] = img->val[p] & 0xFF;//get first byte
                    ptr[1] = img->alpha[p] & 0xFF;//get second byte
                }


                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[0] = img->val[p]>>8;//get first byte
                    ptr[1] = img->val[p] & 0xFF;//get second byte


                    ptr[2] = img->alpha[p]>>8;//get first byte;
                    ptr[3] = img->alpha[p] & 0xFF;;//get second byte
                }
                p++;
            }
        }
    }else if(color_type == PNG_COLOR_TYPE_RGB){
        iftColor rgb, ycbcr;
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                ycbcr.val[0] = img->val[p];
                ycbcr.val[1] = img->Cb[p];
                ycbcr.val[2] = img->Cr[p];

                rgb = iftYCbCrBT2020toRGB(ycbcr, depth, depth);

                ptr[0*byteshift] = rgb.val[0] & 0xFF;//get first byte
                ptr[1*byteshift] = rgb.val[1] & 0xFF;
                ptr[2*byteshift] = rgb.val[2] & 0xFF;

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[(0*byteshift)+1] = ptr[0*byteshift];
                    ptr[(1*byteshift)+1] = ptr[1*byteshift];
                    ptr[(2*byteshift)+1] = ptr[2*byteshift];

                    ptr[0*byteshift] = ((rgb.val[0]>>8) & 0xFF);//get second byte
                    ptr[1*byteshift] = ((rgb.val[1]>>8) & 0xFF);
                    ptr[2*byteshift] = ((rgb.val[2]>>8) & 0xFF);
                }

                p++;
            }
        }

    }else if(color_type == PNG_COLOR_TYPE_RGB_ALPHA){
        iftColor rgb, ycbcr;
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                ycbcr.val[0] = img->val[p];
                ycbcr.val[1] = img->Cb[p];
                ycbcr.val[2] = img->Cr[p];
                ushort alpha = img->alpha[p];

                rgb = iftYCbCrBT2020toRGB(ycbcr, depth, depth);

                ptr[0*byteshift] = rgb.val[0] & 0xFF;//get first byte
                ptr[1*byteshift] = rgb.val[1] & 0xFF;
                ptr[2*byteshift] = rgb.val[2] & 0xFF;
                ptr[3*byteshift] = alpha & 0xFF;

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[(0*byteshift)+1] = ptr[0*byteshift];
                    ptr[(1*byteshift)+1] = ptr[1*byteshift];
                    ptr[(2*byteshift)+1] = ptr[2*byteshift];
                    ptr[(3*byteshift)+1] = ptr[(3*byteshift)];

                    ptr[0*byteshift] = ((rgb.val[0]>>8) & 0xFF);//get second byte
                    ptr[1*byteshift] = ((rgb.val[1]>>8) & 0xFF);
                    ptr[2*byteshift] = ((rgb.val[2]>>8) & 0xFF);
                    ptr[(3*byteshift)] = ((alpha>>8) & 0xFF);
                }
                p++;
            }
        }

    }else{
        iftError("Unknwon color scape", "iftWriteImagePNG");
    };


    iftWritePngImageAux(filename, row_pointers, width, height, depth, color_type);
}

//TODO: suppport to 16-bits channel image
//void iftWriteImageTIFF(const iftImage* image, const char* format, ...){
//    va_list args;
//    char filename[IFT_STR_DEFAULT_SIZE];
//
//    va_start(args, format);
//    vsprintf(filename, format, args);
//    va_end(args);
//
//    TIFF *out= TIFFOpen(filename, "w");
//    int numberOfChannels; // or 3 if there is no alpha channel, you should get a understanding of alpha in class soon.
//
//
//
//    if(iftIsColorImage(image)){
//        if(image->alpha == NULL){
//            numberOfChannels = 3;//RGB
//        }else{
//            numberOfChannels = 4;//RGB_ALPHA
//        }
//    }else{
//        if(image->alpha == NULL){
//            numberOfChannels = 1;//GRAY
//        }else{
//            numberOfChannels = 2;//GRAY_ALPHA
//        }
//    }
//
//    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image->xsize);  // set the width of the image
//    TIFFSetField(out, TIFFTAG_IMAGELENGTH, image->ysize);    // set the height of the image
//    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, numberOfChannels);   // set number of channels per pixel
//    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);    // set the size of the channels
//    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
////   Some other essential fields to set that you do not have to understand for now.
//    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
//    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
//
//    tsize_t linebytes = numberOfChannels * image->xsize;     // length in memory of one row of pixel in the image.
//    unsigned char *buf = NULL;        // buffer used to store the row of pixel information for writing to file
////    Allocating memory to store the pixels of current row
//    buf =(unsigned char *)_TIFFmalloc(linebytes);
//    // We set the strip size of the file to be size of one row of pixels
//    TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, image->xsize*numberOfChannels));
//
//    int k=0;
//    if(numberOfChannels == 1){
//        for (uint32 row = 0; row < (unsigned long)image->ysize; row++)
//        {
//            k=0;
//            for (uint32 col = 0; col < (unsigned long)image->xsize; ++col) {
//                buf[k] = (unsigned char)iftImgVal(image,col,row,0);
//                k++;
//            }
//        }
//    }else if(numberOfChannels == 3){
//        for (uint32 row = 0; row < (unsigned long)image->ysize; row++)
//        {
//            k=0;
//            for (uint32 col = 0; col < (unsigned long)image->xsize; ++col) {
//                buf[k] = (unsigned char)iftImgVal(image,col,row,0);
//                k++;
//                buf[k] = (unsigned char)iftImgCb(image,col,row,0);
//                k++;
//                buf[k] = (unsigned char)iftImgCr(image,col,row,0);
//                k++;
//            }
//        }
//    }else if(numberOfChannels == 4){
//        for (uint32 row = 0; row < (unsigned long)image->ysize; row++)
//        {
//            k=0;
//            for (uint32 col = 0; col < (unsigned long)image->xsize; ++col) {
//                buf[k] = (unsigned char)iftImgVal(image,col,row,0);
//                k++;
//                buf[k] = (unsigned char)iftImgCb(image,col,row,0);
//                k++;
//                buf[k] = (unsigned char)iftImgCr(image,col,row,0);
//                k++;
//                buf[k] = (unsigned char)iftImgAlpha(image,col,row,0);
//            }
//        }
//    }else {
//        iftError("Unknown Color space","iftWriteImageTIFF");
//    }
//
//    (void) TIFFClose(out);
//    if (buf){
//        _TIFFfree(buf);
//    }
//}

void iftWriteImageJPEG(const iftImage* img, const char* format, ...){
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    //code based on external/libjpeg/source/example.c
    /* This struct contains the JPEG compression parameters and pointers to
 * working space (which is allocated as needed by the JPEG library).
 * It is possible to have several such structures, representing multiple
 * compression/decompression processes, in existence at once.  We refer
 * to any one struct (and its associated working data) as a "JPEG object".
 */
    struct jpeg_compress_struct cinfo;

    /* This struct represents a JPEG error handler.  It is declared separately
 * because applications often want to supply a specialized error handler
 * (see the second half of this file for an example).  But here we just
 * take the easy way out and use the standard error handler, which will
 * print a message on stderr and call exit() if compression fails.
 * Note that this struct must live as long as the main JPEG parameter
 * struct, to avoid dangling-pointer problems.
 */
    struct jpeg_error_mgr jerr;

    /* More stuff */
    FILE * outfile;		/* target file */
    JSAMPARRAY buffer;
    /* Step 1: allocate and initialize JPEG compression object */

    /* We have to set up the error handler first, in case the initialization
     * step fails.  (Unlikely, but it could happen if you are out of memory.)
     * This routine fills in the contents of struct jerr, and returns jerr's
     * address which we place into the link field in cinfo.
     */
    cinfo.err = jpeg_std_error(&jerr);

    /* Now we can initialize the JPEG compression object. */
    jpeg_create_compress(&cinfo);
    /* Step 2: specify data destination (eg, a file) */
    /* Note: steps 2 and 3 can be done in either order. */

    /* Here we use the library-supplied code to send compressed data to a
     * stdio stream.  You can also write your own code to do something else.
     * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
     * requires it in order to write binary files.
     */
    if ((outfile = fopen(filename, "wb")) == NULL) {
        fprintf(stderr, "can't open %s\n", filename);
        exit(1);
    }
    jpeg_stdio_dest(&cinfo, outfile);

    /* First we supply a description of the input image.
* Four fields of the cinfo struct must be filled in:
*/


    cinfo.image_width = img->xsize; 	/* image width and height, in pixels */
    cinfo.image_height = img->ysize;
    cinfo.data_precision = iftImageDepth(img);

    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_YCbCr;
    cinfo.jpeg_color_space = JCS_YCbCr;



    /* Now use the library's routine to set default compression parameters.
* (You must set at least cinfo.in_color_space before calling this,
* since the defaults depend on the source color space.)
*/
    jpeg_set_defaults(&cinfo);

    /* Now you can set any non-default parameters you wish to.
* Here we just illustrate the use of quality (quantization table) scaling:
*/
    int quality = 100;
    jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

    /* Step 4: Start compressor */
    /* TRUE ensures that we will write a complete interchange-JPEG file.
     * Pass TRUE unless you are very sure of what you're doing.
     */
    jpeg_start_compress(&cinfo, TRUE);

    /* Step 5: while (scan lines remain to be written) */
    /*           jpeg_write_scanlines(...); */

    /* Here we use the library's state variable cinfo.next_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     * To keep things simple, we pass one scanline per call; you can pass
     * more if you wish, though.
     */
    int row_stride = cinfo.image_width * cinfo.num_components;
    buffer = (*cinfo.mem->alloc_sarray)
            ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    unsigned int imageRow = 0;
    while (cinfo.next_scanline < cinfo.image_height) {
        /* jpeg_write_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could pass
         * more than one scanline at a time if that's more convenient.
         */
        unsigned int imageCol = 0;
        for (unsigned int i = 0; i < (unsigned int)cinfo.image_width; i++) {
            int shift = i*cinfo.num_components;
            buffer[0][(shift+0)] = (unsigned char)iftImgVal(img,imageCol,imageRow,0);
            buffer[0][(shift+1)] = (unsigned char)iftImgCb(img,imageCol,imageRow,0);
            buffer[0][(shift+2)] = (unsigned char)iftImgCr(img,imageCol,imageRow,0);

            imageCol++;
        }
        imageRow++;
        (void) jpeg_write_scanlines(&cinfo, buffer, 1);
    }

    /* Step 6: Finish compression */

    jpeg_finish_compress(&cinfo);
    /* After finish_compress, we can close the output file. */
    fclose(outfile);

    /* Step 7: release JPEG compression object */

    /* This is an important step since it will release a good deal of memory. */
    jpeg_destroy_compress(&cinfo);



}

void iftWriteImagePNGWithAlpha(const iftImage* img, const iftImage *alpha, const char* format, ...) {

    png_bytep *row_pointers;
    int width, height, depth, byteshift, maxVal;

    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    width = img->xsize;
    height = img->ysize;

    depth = iftImageDepth(img);
    byteshift = depth/8;
    // int offset = depth==16?1:0;//to read the second byte first in cases of 16bit images

    maxVal = depth==8?255:65535;

    size_t pixel_size = (iftIsColorImage(img)?4:2 ) * byteshift;

    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (int y=0; y<height; y++)
        row_pointers[y] = (png_byte*) malloc(pixel_size*width);

    if(iftIsColorImage(img)) {
        iftColor rgb, ycbcr;

        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                int a;
                png_byte* ptr = &(row[x*3*byteshift]);

                ycbcr.val[0] = img->val[p];
                ycbcr.val[1] = img->Cb[p];
                ycbcr.val[2] = img->Cr[p];

                rgb = iftYCbCrtoRGB(ycbcr, maxVal);
                a = alpha->val[p];

                ptr[0*byteshift] = rgb.val[0] & 0xFF;//get first byte
                ptr[1*byteshift] = rgb.val[1] & 0xFF;
                ptr[2*byteshift] = rgb.val[2] & 0xFF;
                ptr[3*byteshift] = a & 0xFF;

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[0*byteshift+1] = ptr[0*byteshift];
                    ptr[1*byteshift+1] = ptr[1*byteshift];
                    ptr[2*byteshift+1] = ptr[2*byteshift];
                    ptr[3*byteshift+1] = ptr[3*byteshift];

                    ptr[0*byteshift] = ((rgb.val[0]>>8) & 0xFF);//get second byte
                    ptr[1*byteshift] = ((rgb.val[1]>>8) & 0xFF);
                    ptr[2*byteshift] = ((rgb.val[2]>>8) & 0xFF);
                    ptr[3*byteshift] = ((a>>8) & 0xFF);
                }

                p++;
            }
        }

    }
    else {

        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*2*byteshift]);

                ptr[0] = img->val[p] & 0xFF;//get first byte
                ptr[1] = alpha->val[p] & 0xFF;//get first byte

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[1] = ptr[0];
                    ptr[0] = (img->val[p]>>8) & 0xFF;//get second byte

                    ptr[3] = alpha->val[p] & 0xFF;
                    ptr[2] = (alpha->val[p]>>8) & 0xFF;//get second byte
                }

                p++;
            }
        }

    }

    iftWritePngImageAux(filename, row_pointers, width, height, depth,
                        iftIsColorImage(img) ? PNG_COLOR_TYPE_RGB_ALPHA : PNG_COLOR_TYPE_GRAY_ALPHA);
}

png_bytep* iftReadPngImageAux(const char *file_name, png_structp *png_ptr, png_infop *info_ptr)
{
    png_byte header[8];    // 8 is the maximum size that can be checked

    /* open file and test for it being a png */
    FILE *fp = fopen(file_name, "rb");
    if (!fp)
        iftError("File %s could not be opened for reading", "iftReadPngImageAux", file_name);
    if (fread(header, 1, 8, fp)!=8) iftError("Reading error", "iftReadPngImageAux");
    if (png_sig_cmp(header, 0, 8))
        iftError("File %s is not recognized as a PNG file", "iftReadPngImageAux", file_name);

    int height;
    png_bytep * row_pointers;

    /* initialize stuff */
    png_structp ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!ptr)
        iftError("Internal error: png_create_read_struct failed", "iftReadImagePNG");

    *png_ptr = ptr;
    *info_ptr = png_create_info_struct(*png_ptr);
    int depth = png_get_bit_depth((*png_ptr), (*info_ptr));
    if(depth < 8){
        png_set_expand_gray_1_2_4_to_8(ptr);
    }


    if (!(*info_ptr))
        iftError("Internal error: png_create_info_struct failed", "iftReadImagePNG");

    if (setjmp(png_jmpbuf(*png_ptr)))
        iftError("Internal error: Error during init_io", "iftReadImagePNG");


    png_init_io(*png_ptr, fp);
    png_set_sig_bytes(*png_ptr, 8);

    png_read_info(*png_ptr, *info_ptr);
	// reduces the pixels back down to the original bit depth
	//png_color_8p sig_bit = NULL;
	//if (png_get_sBIT(*png_ptr, *info_ptr, &sig_bit)) {
	//	png_set_shift(*png_ptr, sig_bit);
	//}
    height = png_get_image_height(*png_ptr, *info_ptr);
    png_read_update_info(*png_ptr, *info_ptr);


    /* read file */
    if (setjmp(png_jmpbuf(*png_ptr)))
        iftError("Internal error: Error during read_image", "iftReadImagePNG");

    row_pointers = (png_bytep*) iftAlloc(height, sizeof(png_bytep));
    for (int y=0; y<height; y++)
        row_pointers[y] = (png_byte*) iftAlloc(png_get_rowbytes(*png_ptr, *info_ptr), 1);


    png_read_image(*png_ptr, row_pointers);

    fclose(fp);

    return row_pointers;
}

iftImage* iftReadImagePNG(const char* format, ...) {

    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    png_infop info_ptr;
    png_structp png_ptr;
    png_bytep *row_pointers;

    row_pointers = iftReadPngImageAux(filename, &png_ptr, &info_ptr);

    int width, height, color_type, depth;

    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    color_type = png_get_color_type(png_ptr, info_ptr);
    depth = png_get_bit_depth(png_ptr, info_ptr);
    iftImage* img = iftCreateImage(width, height, 1);
    unsigned int numberChannels = png_get_channels(png_ptr, info_ptr);

    int byteshift = depth/8;

    int x, y;

    int p = 0;

    if(color_type==PNG_COLOR_TYPE_GRAY)//gray image
    {
        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                img->val[p] = ptr[0];
                if(depth==16) {
                    img->val[p] = (img->val[p]<<8)+ptr[1];
                }
                p++;
            }
        }
    }else if(color_type==PNG_COLOR_TYPE_GRAY_ALPHA ){
        if(img->alpha == NULL){
            iftSetAlpha(img,0);
        }
        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                if(depth == 8){
                    img->val[p] = ptr[0];
                    img->alpha[p] = ptr[1];
                }
                else if(depth==16) {
                    img->val[p] = ptr[0];
                    img->val[p] = (img->val[p]<<8)+ptr[1];
                    img->alpha[p] = ptr[2];
                    img->alpha[p] = (img->alpha[p]<<8)+ptr[3];
                }
                p++;
            }
        }
    }
    else if(color_type == PNG_COLOR_TYPE_RGB){//color image

        iftSetCbCr(img, 128);
        iftColor rgb, ycbcr;

        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];

            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                rgb.val[0] = ptr[0*byteshift];
                rgb.val[1] = ptr[1*byteshift];
                rgb.val[2] = ptr[2*byteshift];

                if(depth==16) { //read second byte in case of 16bit images
                    rgb.val[0] = (rgb.val[0]<<8) + ptr[1];
                    rgb.val[1] = (rgb.val[1]<<8) + ptr[3];
                    rgb.val[2] = (rgb.val[2]<<8) + ptr[5];
                }

                ycbcr = iftRGBtoYCbCrBT2020(rgb, depth, depth);

                img->val[p] = ycbcr.val[0];
                img->Cb[p]  = ycbcr.val[1];
                img->Cr[p]  = ycbcr.val[2];

                p++;
            }
        }
    }else if(color_type == PNG_COLOR_TYPE_RGB_ALPHA){
        iftSetCbCr(img, 128);
        iftColor rgb, ycbcr;
        if(img->alpha == NULL){
            iftSetAlpha(img,0);
        }

        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                rgb.val[0] = ptr[0*byteshift];
                rgb.val[1] = ptr[1*byteshift];
                rgb.val[2] = ptr[2*byteshift];
                ushort alpha = ptr[3*byteshift];

                if(depth==16) { //read second byte in case of 16bit images
                    rgb.val[0] = (rgb.val[0]<<8) + ptr[1];
                    rgb.val[1] = (rgb.val[1]<<8) + ptr[3];
                    rgb.val[2] = (rgb.val[2]<<8) + ptr[5];
                    alpha = (alpha<<8) +  ptr[7];
                }

                ycbcr = iftRGBtoYCbCr(rgb, depth==8?255:65535);

                img->val[p] = ycbcr.val[0];
                img->Cb[p] = ycbcr.val[1];
                img->Cr[p] = ycbcr.val[2];
                img->alpha[p] = alpha;

                p++;
            }
        }
    }

    for (y = 0; y < height; ++y) {
        iftFree(row_pointers[y]);
    }

    iftFree(row_pointers);

    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

    img->dz = 0.0;

    return img;
}

void iftConvertRGBImagetoYCbCrImage(iftImage* image_rgb, int normalization_value){
    iftColor yCbCr;
    iftColor rgb;
    int k = 0;
    for (int z = 0; z < image_rgb->zsize; ++z) {
        for (int y = 0; y < image_rgb->ysize; ++y) {
            for (int x = 0; x < image_rgb->xsize; ++x) {
                rgb.val[0] = image_rgb->val[k];
                rgb.val[1] = image_rgb->Cb[k];
                rgb.val[2] = image_rgb->Cr[k];
                yCbCr = iftRGBtoYCbCr(rgb, normalization_value);
                image_rgb->val[k] = yCbCr.val[0];
                image_rgb->Cb[k] = yCbCr.val[1];
                image_rgb->Cr[k] = yCbCr.val[2];
                k++;
            }
        }
    }
}

//TODO: suppport to 16-bits channel image
//iftImage* iftReadImageTIFF(const char* format, ...){
//    va_list args;
//    char filename[IFT_STR_DEFAULT_SIZE];
//
//    va_start(args, format);
//    vsprintf(filename, format, args);
//    va_end(args);
//
//    TIFF *tif=TIFFOpen(filename, "r");
//    uint32 width;
//    uint32 height;
//    uint32 npixels;
//    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
//    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
//
//    npixels = width*height;
//    uint32 * raster=(uint32 *) _TIFFmalloc(npixels *sizeof(uint32));
//
//    int status = TIFFReadRGBAImage(tif, width, height, raster, 0);//read the entire image
//    if(status==0){
//        iftError("An error occured when reading image","iftWriteImageTIFF");
//    }
//
//    iftImage* image = iftCreateImage(width,height,1);
//    iftSetCbCrAlpha(image,0);
//    int k = 0;
//    iftColor rgb, ycbcr;
//    for (int y = height-1; y >=0 ; --y) {
//        for (int x = 0; x < ((int)width); ++x) {
//            //ABGR - 8bits per channel
//            int A = raster[k] >> 24;
//            rgb.val[2] = (raster[k] >> 16) & 0x0000FF;//B
//            rgb.val[1] = (raster[k] >> 8) & 0x0000FF;//G
//            rgb.val[0] =  (raster[k]) & 0x0000FF;//R
//            ycbcr = iftRGBtoYCbCr(rgb,255);
//
//            iftImgVal(image,x,y,0) = ycbcr.val[0];
//            iftImgCb(image,x,y,0) = ycbcr.val[1];
//            iftImgCr(image,x,y,0) = ycbcr.val[2];
//            iftImgAlpha(image,x,y,0) = A;
//            k++;
//        }
//    }
//
//    _TIFFfree(raster);
//    TIFFClose(tif);
//    return image;
//}

iftImage* iftReadImageJPEG(const char* format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    iftImage* image = NULL;
    //code based on externals/libjpeg/source/example.c
    /* This struct contains the JPEG decompression parameters and pointers to
* working space (which is allocated as needed by the JPEG library).
*/
    struct jpeg_decompress_struct cinfo;
    /* We use our private extension JPEG error handler.
* Note that this struct must live as long as the main JPEG parameter
* struct, to avoid dangling-pointer problems.
*/
    struct jpeg_error_mgr jerr;

    /* More stuff */
    FILE * infile;		/* source file */
    JSAMPARRAY buffer;		/* Output row buffer */
    int row_stride;		/* physical row width in output buffer */
    /* In this example we want to open the input file before doing anything else,
 * so that the setjmp() error recovery below can assume the file is open.
 * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
 * requires it in order to read binary files.
 */

    if ((infile = fopen(filename, "rb")) == NULL) {
        printf("[readImageJPEG] can't open %s\n",filename);
        return NULL;
    }

    /* Step 1: allocate and initialize JPEG decompression object */

    /* We set up the normal JPEG error routines, then override error_exit. */
    cinfo.err = jpeg_std_error(&jerr);
    //jerr.pub.error_exit = my_error_exit;

    /* Establish the setjmp return context for my_error_exit to use. */
    jmp_buf setjmp_buffer;
    if (setjmp(setjmp_buffer)) {
        /* If we get here, the JPEG code has signaled an error.
         * We need to clean up the JPEG object, close the input file, and return.
         */
        jpeg_destroy_decompress(&cinfo);
        printf("[readImageJPEG] code has signaled an error\n");
        fclose(infile);
        return NULL;
    }

    /* Now we can initialize the JPEG decompression object. */
    jpeg_create_decompress(&cinfo);
    /* Step 2: specify data source (eg, a file) */
    jpeg_stdio_src(&cinfo, infile);
    /* Step 3: read file parameters with jpeg_read_header() */
    (void) jpeg_read_header(&cinfo, TRUE);
    /* We can ignore the return value from jpeg_read_header since
     *   (a) suspension is not possible with the stdio data source, and
     *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
     * See libjpeg.txt for more info.
     */

    /* Step 4: set parameters for decompression */

    /* In this example, we don't need to change any of the defaults set by
     * jpeg_read_header(), so we do nothing here.
     */

    /* Step 5: Start decompressor */
    (void) jpeg_start_decompress(&cinfo);
    /* We can ignore the return value since suspension is not possible
 * with the stdio data source.
 */

    /* We may need to do some setup of our own at this point before reading
 * the data.  After jpeg_start_decompress() we have the correct scaled
 * output image dimensions available, as well as the output colormap
 * if we asked for color quantization.
 * In this example, we need to make an output work buffer of the right size.
 */
    /* JSAMPLEs per row in output buffer */
    row_stride = cinfo.output_width * cinfo.output_components;
    /* Make a one-row-high sample array that will go away when done with image */
    buffer = (*cinfo.mem->alloc_sarray)
            ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


    /* Step 6: while (scan lines remain to be read) */
    /*           jpeg_read_scanlines(...); */

    /* Here we use the library's state variable cinfo.output_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     */
    image = iftCreateImage(cinfo.output_width,cinfo.output_height,1);

    //0 - JCS_GRAYSCALE
    //1 - JCS_RGB
    //2 - JCS_YCbCr
    //3 - JCS_CMYK
    //4 - JCS_YCCK
    //5 - JCS_BG_RGB
    //6 - JCS_BG_YCC
    unsigned  int imageRow = 0;
    unsigned int imageCol = 0;
    iftColor rgb;
    iftColor YCbCr;
    float scalingFactor = pow(2,cinfo.data_precision)-1;
    switch (cinfo.out_color_space){
        case JCS_GRAYSCALE:

            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    iftImgVal(image,imageCol,imageRow,0) = buffer[0][shift];
                    imageCol++;
                }
                imageRow++;
            }
            break;

        case JCS_RGB:

            iftSetCbCr(image,128);
            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    rgb.val[0] = buffer[0][shift+0];
                    rgb.val[1] = buffer[0][shift+1];
                    rgb.val[2] = buffer[0][shift+2];
                    YCbCr = iftRGBtoYCbCr(rgb,scalingFactor);
                    iftImgVal(image,imageCol,imageRow,0) = YCbCr.val[0];
                    iftImgCb(image,imageCol,imageRow,0) = YCbCr.val[1];
                    iftImgCr(image,imageCol,imageRow,0) = YCbCr.val[2];
                    imageCol++;
                }
                imageRow++;
            }

            break;

        case JCS_YCbCr:

            iftSetCbCr(image,128);
            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    iftImgVal(image,imageCol,imageRow,0) = buffer[0][shift+0];
                    iftImgCb(image,imageCol,imageRow,0) = buffer[0][shift+1];
                    iftImgCr(image,imageCol,imageRow,0) = buffer[0][shift+2];
                    imageCol++;
                }
                imageRow++;
            }

            break;
        case JCS_CMYK:

            iftSetCbCr(image,128);
            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    //convert CMYK to RGB (reference: http://www.rapidtables.com/convert/color/cmyk-to-rgb.htm)
                    rgb.val[0] = 255*(100-buffer[0][shift+0])*(100-buffer[0][shift+3]);
                    rgb.val[1] = 255*(100-buffer[0][shift+1])*(100-buffer[0][shift+3]);;
                    rgb.val[2] = 255*(100-buffer[0][shift+2])*(100-buffer[0][shift+3]);;
                    YCbCr = iftRGBtoYCbCr(rgb,scalingFactor);
                    iftImgVal(image,imageCol,imageRow,0) = YCbCr.val[0];
                    iftImgCb(image,imageCol,imageRow,0) = YCbCr.val[1];
                    iftImgCr(image,imageCol,imageRow,0) = YCbCr.val[2];
                    imageCol++;
                }
                imageRow++;
            }

            break;
        case JCS_YCCK:

            iftSetCbCr(image,128);
            imageRow = 0;
            iftWarning("Image is Y/Cb/Cr/K color space. The channel K is ignored", "iftReadImageJPEG");
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    iftImgVal(image,imageCol,imageRow,0) = buffer[0][shift+0];
                    iftImgCb(image,imageCol,imageRow,0) = buffer[0][shift+1];
                    iftImgCr(image,imageCol,imageRow,0) = buffer[0][shift+2];
                    imageCol++;
                }
                imageRow++;
            }

            break;
        case JCS_BG_RGB:
    
            iftError("Big gamut red/green/blue color space not supported", "iftReadImageJPEG");

            break;
        case JCS_BG_YCC:
    
            iftError("Big gamut red/green/blue color space not supported", "iftReadImageJPEG");

            break;
        default:
    
            iftError("Unkwon color space", "iftReadImageJPEG");

            break;
    }

    /* Step 7: Finish decompression */
    (void) jpeg_finish_decompress(&cinfo);

    /* This is an important step since it will release a good deal of memory. */
    jpeg_destroy_decompress(&cinfo);

    /* After finish_decompress, we can close the input file.
     * Here we postpone it until after no more JPEG errors are possible,
     * so as to simplify the setjmp error logic above.  (Actually, I don't
     * think that jpeg_destroy can do an error exit, but why assume anything...)
     */

    fclose(infile);
    /* At this point you may want to check to see whether any corrupt-data
     * warnings occurred (test whether jerr.pub.num_warnings is nonzero).
     */

    //jerr.num_warnings; //useful to know about corrupted data
    //printf("%ld\n",jerr.num_warnings);

    return image;
}

void iftWriteImage(const iftImage *img, const char *format, ...) {
    FILE   *fp     = NULL;
    int    p;
    uchar  *data8  = NULL;
    ushort *data16 = NULL;
    int    *data32 = NULL;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    int img_min_val = iftMinimumValue(img);
    int img_max_val = iftMaximumValue(img);

    if (img_min_val < 0) {
        char msg[200];
        sprintf(msg, "Shifting image values from [%d,%d] to [%d,%d] on the original image\n",
                img_min_val, img_max_val, 0, img_max_val - img_min_val);
        iftWarning(msg, "iftWriteImage");
        for (p = 0; p < img->n; p++)
            img->val[p] = img->val[p] - img_min_val;
        img_max_val = img_max_val - img_min_val;
    }

    fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteImage", filename);

    fprintf(fp, "SCN\n");
    fprintf(fp, "%d %d %d\n", img->xsize, img->ysize, img->zsize);
    fprintf(fp, "%f %f %f\n", img->dx, img->dy, img->dz);


    if (img_max_val < 256) {
        fprintf(fp, "%d\n", 8);
        data8 = iftAllocUCharArray(img->n);
        for (p = 0; p < img->n; p++)
            data8[p] = (uchar) img->val[p];
        fwrite(data8, sizeof(uchar), img->n, fp);
        iftFree(data8);
    } else if (img_max_val < 65536) {
        fprintf(fp, "%d\n", 16);
        data16 = iftAllocUShortArray(img->n);
        for (p = 0; p < img->n; p++)
            data16[p] = (ushort) img->val[p];
        fwrite(data16, sizeof(ushort), img->n, fp);

        iftFree(data16);
    } else if (img_max_val < IFT_INFINITY_INT) {
        fprintf(fp, "%d\n", 32);
        data32 = iftAllocIntArray(img->n);
        for (p = 0; p < img->n; p++)
            data32[p] = img->val[p];
        fwrite(data32, sizeof(int), img->n, fp);
        iftFree(data32);
    }

    fclose(fp);
}


void iftWriteImageNIfTI(const iftImage *img, const char *format, ...) {
    va_list args;
    char path[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(path, format, args);
    va_end(args);

    if (!iftRegexMatch(path, "^.+\\.(nii|hdr)$"))
        iftError("Invalid file format: %s. Try .nii, or .hdr", "iftWriteImageNIfTI", iftFileExt(path));

    bool is_analyze = iftEndsWith(path, ".hdr");

    NIfTI1Header *header = _iftCreateNIfTI1HeaderFromImage(img, is_analyze);

    FILE *fp = fopen(path, "wb");
    fwrite(header, sizeof(*header), 1, fp); // storing the nifti/analyze header with 348 bytes

    // ANALYZE 7.5 file
    if (is_analyze) {
        char *hdr_base = iftBasename(path);
        char *img_path = iftConcatStrings(2, hdr_base, ".img");
        fclose(fp);

        // opening the image (data) file
        fp = fopen(img_path, "wb");

        iftFree(hdr_base);
        iftFree(img_path);
    }
    else { // nifti
        char nifti_ext[4] = {0, 0, 0, 0};
        fwrite(nifti_ext, 4, 1, fp); // storing the extra 4 bytes of nifti files
    }


    int nbytes = header->bitpix / 8;

    if (header->datatype == NIFTI_TYPE_UINT8) {
        for (int p = 0; p < img->n; p++) {
            uchar val = (uchar) img->val[p];
            fwrite(&val, nbytes, 1, fp); // storing the image data
        }
    }
    else if (header->datatype == NIFTI_TYPE_INT8) {
        for (int p = 0; p < img->n; p++) {
            char val = (char) img->val[p];
            fwrite(&val, nbytes, 1, fp); // storing the image data
        }
    }
    else if (header->datatype == NIFTI_TYPE_UINT16) {
        for (int p = 0; p < img->n; p++) {
            ushort val = (ushort) img->val[p];
            fwrite(&val, nbytes, 1, fp); // storing the image data
        }
    }
    else if (header->datatype == NIFTI_TYPE_INT16) {
        for (int p = 0; p < img->n; p++) {
            short val = (short) img->val[p];
            fwrite(&val, nbytes, 1, fp); // storing the image data
        }
    }
    else if (header->datatype == NIFTI_TYPE_INT32) {
        fwrite(img->val, nbytes, img->n, fp); // storing the image data
    }
    else iftError("Unsupported NIfTI datatype", "iftWriteImageNIfTI");

    fclose(fp);

    iftFree(header);
}


void iftWriteImageNIfTIGZip(const iftImage *img, const char *format, ...) {
    va_list args;
    char path[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(path, format, args);
    va_end(args);

    if (!iftEndsWith(path, ".nii.gz"))
        iftError("Invalid file format: %s. Try .nii.gz", "iftWriteImageNIfTIGZip", iftFileExt(path));

    NIfTI1Header *header = _iftCreateNIfTI1HeaderFromImage(img, false);

    iftGZipFile fp = iftGZipOpen(path, "wb", true);
    iftGZipWrite(header, sizeof(*header), 1, fp); // storing the nifti/analyze header with 348 bytes
    char nifti_ext[4] = {0, 0, 0, 0};
    iftGZipWrite(nifti_ext, 4, 1, fp); // storing the extra 4 bytes of nifti files

    int nbytes = header->bitpix / 8;
    
    if (header->datatype == NIFTI_TYPE_UINT8) {
        for (int p = 0; p < img->n; p++) {
            uchar val = (uchar) img->val[p];
            iftGZipWrite(&val, nbytes, 1, fp); // storing the image data
        }
    }
    else if (header->datatype == NIFTI_TYPE_INT8) {
        for (int p = 0; p < img->n; p++) {
            char val = (char) img->val[p];
            iftGZipWrite(&val, nbytes, 1, fp); // storing the image data
        }
    }
    else if (header->datatype == NIFTI_TYPE_UINT16) {
        for (int p = 0; p < img->n; p++) {
            ushort val = (ushort) img->val[p];
            iftGZipWrite(&val, nbytes, 1, fp); // storing the image data
        }
    }
    else if (header->datatype == NIFTI_TYPE_INT16) {
        for (int p = 0; p < img->n; p++) {
            short val = (short) img->val[p];
            iftGZipWrite(&val, nbytes, 1, fp); // storing the image data
        }
    }
    else if (header->datatype == NIFTI_TYPE_INT32) {
        iftGZipWrite(img->val, nbytes, img->n, fp); // storing the image data
    }
    else iftError("Unsupported NIfTI datatype", "iftWriteImageNIfTIGZip");
    
    iftGZipClose(&fp);

    iftFree(header);
}





void iftWriteImageAnalyze(const iftImage *img, const char *format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    iftWriteImageNIfTI(img, filename);
}

void iftWriteImageGZip(const iftImage *img, const char *format, ...) {
    iftGZipFile fp = NULL;
    int    p;
    uchar  *data8  = NULL;
    ushort *data16 = NULL;
    int    *data32 = NULL;
    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE], header[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    int img_min_val = iftMinimumValue(img);
    int img_max_val = iftMaximumValue(img);

    if (img_min_val < 0) {
        char msg[200];
        sprintf(msg, "Shifting image values from [%d,%d] to [%d,%d] on the original image\n",
                img_min_val, img_max_val, 0, img_max_val - img_min_val);
        iftWarning(msg, "iftWriteImage");
        for (p = 0; p < img->n; p++)
            img->val[p] = img->val[p] - img_min_val;
        img_max_val = img_max_val - img_min_val;
    }

    fp = iftGZipOpen(filename, "wb", 1);
    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteImageGZip", filename);


    if (img_max_val < 256) {
        // The header is saved with a single line because we use iftGZipGets to read it in iftReadImageGZip, and that
        // function reads everything until a '\n' is found. If the gzip file is uncompressed with an external tool,
        // iftReadImage may still be used to read it
        sprintf(header, "SCN %d %d %d %f %f %f %d\n", img->xsize, img->ysize, img->zsize,
                img->dx, img->dy, img->dz, 8);

        // Writes everything including the '\n', which is important for reading the header with iftGZipGets
        iftGZipPuts(header, fp);

        data8 = iftAllocUCharArray(img->n);
        for (p = 0; p < img->n; p++)
            data8[p] = (uchar) img->val[p];
        iftGZipWrite(data8, sizeof(uchar), img->n, fp);
        iftFree(data8);
    } else if (img_max_val < 65536) {
        // The header is saved with a single line because we use iftGZipGets to read it in iftReadImageGZip, and that
        // function reads everything until a '\n' is found. If the gzip file is uncompressed with an external tool,
        // iftReadImage may still be used to read it
        sprintf(header, "SCN %d %d %d %f %f %f %d\n", img->xsize, img->ysize, img->zsize,
                img->dx, img->dy, img->dz, 16);

        // Writes everything including the '\n', which is important for reading the header with iftGZipGets
        iftGZipPuts(header, fp);

        data16 = iftAllocUShortArray(img->n);
        for (p = 0; p < img->n; p++)
            data16[p] = (ushort) img->val[p];
        iftGZipWrite(data16, sizeof(ushort), img->n, fp);

        iftFree(data16);
    } else if (img_max_val < IFT_INFINITY_INT) {
        // The header is saved with a single line because we use iftGZipGets to read it in iftReadImageGZip, and that
        // function reads everything until a '\n' is found. If the gzip file is uncompressed with an external tool,
        // iftReadImage may still be used to read it
        sprintf(header, "SCN %d %d %d %f %f %f %d\n", img->xsize, img->ysize, img->zsize,
                img->dx, img->dy, img->dz, 32);
        // Writes everything including the '\n', which is important for reading the header with iftGZipGets
        iftGZipPuts(header, fp);

        data32 = iftAllocIntArray(img->n);
        for (p = 0; p < img->n; p++)
            data32[p] = img->val[p];
        iftGZipWrite(data32, sizeof(int), img->n, fp);
        iftFree(data32);
    }

    iftGZipClose(&fp);
}

iftImage *iftReadImageP5(const char *format, ...) {
    iftImage *img    = NULL;
    FILE     *fp     = NULL;
    uchar    *data8  = NULL;
    ushort   *data16 = NULL;
    char     type[10];
    int      p, v, xsize, ysize, zsize, hi, lo;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageP5", filename);
    }

    if (fscanf(fp, "%s\n", type) != 1) {
        iftError("Reading error", "iftReadImageP5");
    }

    if (iftCompareStrings(type, "P5")) {

        iftSkipComments(fp);

        if (fscanf(fp, "%d %d\n", &xsize, &ysize) != 2)
            iftError("Reading error", "iftReadImageP5");
        zsize = 1;

        img = iftCreateImage(xsize, ysize, zsize);
        img->dz = 0.0;

        if (fscanf(fp, "%d", &v) != 1)
            iftError("Reading error", "iftReadImageP5");

        while (fgetc(fp) != '\n');

        if ((v <= 255) && (v > 0)) {
            data8 = iftAllocUCharArray(img->n);

            if (fread(data8, sizeof(uchar), img->n, fp) != img->n)
                iftError("Reading error", "iftReadImageP5");

            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data8[p];

            iftFree(data8);

        } else if ((v <= 65535) && (v > 255)) {
            data16 = iftAllocUShortArray(img->n);

            for (p = 0; p < img->n; p++) {
                if ((hi = fgetc(fp)) == EOF)
                    iftError("Reading error", "iftReadImageP5");
                if ((lo = fgetc(fp)) == EOF)
                    iftError("Reading error", "iftReadImageP5");

                data16[p] = (hi << 8) + lo;
            }

            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data16[p];

            iftFree(data16);

        } else {
            iftError("Invalid maximum value", "iftReadImageP5");
        }
    } else {
        iftError("Invalid image type", "iftReadImageP5");
    }

    fclose(fp);
    return (img);
}

void iftWriteImageP5(const iftImage *img, const char *format, ...) {
    FILE   *fp     = NULL;
    int    p, hi, lo;
    uchar  *data8  = NULL;
    ushort *data16 = NULL;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageP5", filename);

    fprintf(fp, "P5\n");
    fprintf(fp, "%d %d\n", img->xsize, img->ysize);

    int img_max_val = iftMaximumValue(img);
    int img_min_val = iftMinimumValue(img);

    if ((img_max_val < 256) && (img_min_val >= 0)) {
        fprintf(fp, "%d\n", 255);
        data8 = iftAllocUCharArray(img->n);
        for (p = 0; p < img->n; p++)
            data8[p] = (uchar) img->val[p];
        fwrite(data8, sizeof(uchar), img->n, fp);
        iftFree(data8);
    } else if (img_max_val < 65536) {
        fprintf(fp, "%d\n", 65535);
        data16 = iftAllocUShortArray(img->n);
        for (p = 0; p < img->n; p++)
            data16[p] = (ushort) img->val[p];

        {
#define HI(num) (((num) & 0x0000FF00) >> 8)
#define LO(num) ((num) & 0x000000FF)
            for (p = 0; p < img->n; p++) {
                hi = HI(data16[p]);
                lo = LO(data16[p]);
                fputc(hi, fp);
                fputc(lo, fp);
            }
        }

        iftFree(data16);
    } else {
        char msg[200];
        sprintf(msg, "Cannot write image as P5 (%d/%d)", img_max_val, img_min_val);
        iftError(msg, "iftWriteImageP5");
    }
    fclose(fp);
}

iftImage *iftReadImageAsP5(const char *format, ...)
{
    char command[400];
    iftImage *img;

    va_list args;
    char filename[300];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    sprintf(command,"convert %s temp.pgm",filename);
    if (system(command)==-1) iftError("Command error", "iftReadImageAsP5");
    img = iftReadImageP5("temp.pgm");
    if (system("rm -f temp.pgm")==-1) iftError("Command error", "iftReadImageAsP5");
    return(img);
}

void iftWriteImageExtFormat(iftImage *img, const char *format, ...)
{
    char command[400];

    va_list args;
    char filename[300];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    if (iftIs3DImage(img))
        iftError("It is not handling 3D images yet", "iftWriteImageExtFormat");

    if (iftIsColorImage(img)){
        iftWriteImageP6(img,"temp.ppm");
        sprintf(command,"convert temp.ppm %s",filename);
        if (system(command)==-1) iftError("Command error", "iftWriteImageExtFormat");
        if (system("rm -f temp.ppm")==-1) iftError("Command error", "iftWriteImageExtFormat");
    }else{
        iftWriteImageP5(img,"temp.pgm");
        sprintf(command,"convert temp.pgm %s",filename);
        if (system(command)==-1) iftError("Command error", "iftWriteImageExtFormat");
        if (system("rm -f temp.pgm")==-1) iftError("Command error", "iftWriteImageExtFormat");
    }
}

iftImage *iftReadImageAsP6(const char *format, ...)
{
    char command[400];
    iftImage *img;

    va_list args;
    char filename[300];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    sprintf(command,"convert %s temp.ppm",filename);
    if(system(command)==-1) iftError("Command error", "iftReadImageAsP6");
    img = iftReadImageP6("temp.ppm");
    if(system("rm -f temp.ppm")==-1) iftError("Command error", "iftReadImageAsP6");;
    return(img);
}

iftImage *iftReadImageP6(const char *format, ...)
{
    iftImage  *img=NULL;
    FILE    *fp=NULL;
    char    type[10];
    int     p,v,xsize,ysize,zsize;
    ushort rgb16[3];
    iftColor RGB,YCbCr;

    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename,"r");
    if (fp == NULL){
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageP6", filename);
    }

    if(fscanf(fp,"%s\n",type)!=1)
        iftError("Reading error", "iftReadImageP6");
    if(iftCompareStrings(type,"P6")){

        iftSkipComments(fp);

        if(fscanf(fp,"%d %d\n",&xsize,&ysize)!=2)
            iftError("Reading error", "iftReadImageP6");

        zsize = 1;
        img = iftCreateImage(xsize,ysize,zsize);
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
        img->dz = 0.0;
        if (fscanf(fp,"%d",&v)!=1)
            iftError("Reading error", "iftReadImageP6");

        while(fgetc(fp) != '\n');

        if (v >= 0 && v < 256) {
            for (p=0; p < img->n; p++) {
                RGB.val[0] = fgetc(fp);
                RGB.val[1] = fgetc(fp);
                RGB.val[2] = fgetc(fp);
                YCbCr      = iftRGBtoYCbCr(RGB,255);
                img->val[p]=YCbCr.val[0];
                img->Cb[p] =(ushort)YCbCr.val[1];
                img->Cr[p] =(ushort)YCbCr.val[2];
            }
        } else if (v >= 256 && v <= 65536) {

            int rgbBitDepth = ceil(iftLog(v, 2));
            int ycbcrBitDepth = rgbBitDepth;

            if(ycbcrBitDepth<10)
                ycbcrBitDepth = 10;
            else if(ycbcrBitDepth < 12)
                ycbcrBitDepth = 12;
            else if(ycbcrBitDepth < 16)
                ycbcrBitDepth = 16;

            for (p=0; p < img->n; p++) {
                // read 6 bytes for each image pixel
                if (fread(rgb16, 2, 3, fp) == 3) {
                    // the PPM format specifies 2-byte integers as big endian,
                    // so we need to swap the bytes if the architecture is little endian
                    RGB.val[0]  = ((rgb16[0] & 0xff) << 8) | ((ushort) rgb16[0] >> 8);
                    RGB.val[1]  = ((rgb16[1] & 0xff) << 8) | ((ushort) rgb16[1] >> 8);
                    RGB.val[2]  = ((rgb16[2] & 0xff) << 8) | ((ushort) rgb16[2] >> 8);
                    YCbCr       = iftRGBtoYCbCrBT2020(RGB, rgbBitDepth, ycbcrBitDepth);
//                    YCbCr = iftRGBtoYCbCr(RGB, v);
                    img->val[p] = YCbCr.val[0];
                    img->Cb[p]  = (ushort)YCbCr.val[1];
                    img->Cr[p]  = (ushort)YCbCr.val[2];
                }
            }
        } else {
            iftError("Invalid maximum value", "iftReadImageP6");
        }
    }else{
        iftError("Invalid image type", "iftReadImageP6");
    }

    fclose(fp);
    return(img);
}


void iftWriteImageP6(const iftImage *img, const char *format, ...) {
    FILE     *fp = NULL;
    int      p;
    ushort   rgb16[3];
    iftColor YCbCr, RGB;

    if (!iftIsColorImage(img))
        iftError("Image is not colored", "iftWriteImageP6");

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "w");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageP6", filename);

    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", img->xsize, img->ysize);

    int img_max_val = iftMaximumValue(img);
    int img_min_val = iftMinimumValue(img);

    if (img_min_val < 0) {
        iftError("Cannot write image as P6", "iftWriteImageP6");
    }
    if (img_max_val < 256) {
        fprintf(fp, "%d\n", 255);
        for (p = 0; p < img->n; p++) {
            YCbCr.val[0] = img->val[p];
            YCbCr.val[1] = img->Cb[p];
            YCbCr.val[2] = img->Cr[p];

            RGB = iftYCbCrtoRGB(YCbCr, 255);

            fputc(((uchar) RGB.val[0]), fp);
            fputc(((uchar) RGB.val[1]), fp);
            fputc(((uchar) RGB.val[2]), fp);
        }
    } else if (img_max_val < 65536) {
//        int rgbBitDepth = 9;
//        // find the bit depth for the maximum value img_max_val
//        while ((1 << rgbBitDepth) <= img_max_val) {
//            rgbBitDepth++;
//        }

        int rgbBitDepth = ceil(iftLog(img_max_val, 2));

        fprintf(fp, "%d\n", (1 << rgbBitDepth) - 1);
        for (p = 0; p < img->n; p++) {
            YCbCr.val[0] = img->val[p];
            YCbCr.val[1] = img->Cb[p];
            YCbCr.val[2] = img->Cr[p];
            RGB = iftYCbCrBT2020toRGB(YCbCr, rgbBitDepth, rgbBitDepth);
//            RGB = iftYCbCrtoRGB(YCbCr, img_max_val);
            // the PPM format specifies 2-byte integers as big endian,
            // so we need to swap the bytes if the architecture is little endian
            rgb16[0] = ((RGB.val[0] & 0xff) << 8) | ((ushort) RGB.val[0] >> 8);
            rgb16[1] = ((RGB.val[1] & 0xff) << 8) | ((ushort) RGB.val[1] >> 8);
            rgb16[2] = ((RGB.val[2] & 0xff) << 8) | ((ushort) RGB.val[2] >> 8);
            // write 6 bytes for each image pixel
            if (fwrite(rgb16, 2, 3, fp) != 3) {
                iftError("Cannot write 16-bit image as P6", "iftWriteImageP6");
            }
        }
    } else {
        iftError("Cannot write image as P6", "iftWriteImageP6");
    }
    fclose(fp);
}


iftImage *iftReadImageP2(const char *format, ...) {
    iftImage *img = NULL;
    FILE     *fp  = NULL;
    char     type[10];
    int      p, v, xsize, ysize, zsize;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "r");

    if (fp == NULL) {
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageP2", filename);
    }

    if (fscanf(fp, "%s\n", type) != 1)
        iftError("Reading error", "iftReadImageP2");

    if (iftCompareStrings(type, "P2")) {

        iftSkipComments(fp);

        if (fscanf(fp, "%d %d\n", &xsize, &ysize) != 2)
            iftError("Reading error", "iftReadImageP2");
        zsize = 1;
        img   = iftCreateImage(xsize, ysize, zsize);
        img->dz = 0.0;
        if (fscanf(fp, "%d", &v) != 1)
            iftError("Reading error", "iftReadImageP2");

        while (fgetc(fp) != '\n');

        for (p = 0; p < img->n; p++)
            if (fscanf(fp, "%d", &img->val[p]) != 1)
                iftError("Reading error", "iftReadImageP2");

    } else {
        iftError("Invalid image type", "iftReadImageP2");
    }

    fclose(fp);
    return (img);
}

void iftWriteImageP2(const iftImage *img, const char *format, ...) {
    FILE *fp = NULL;
    int  p;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];
    int     depth = iftImageDepth(img);
      
    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "w");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageP2", filename);

    fprintf(fp, "P2\n");
    fprintf(fp, "%d %d\n", img->xsize, img->ysize);

    int img_max_val = iftMaxImageRange(depth);
    fprintf(fp, "%d\n", img_max_val);
    for (p = 0; p < img->n; p++) {
        fprintf(fp, "%d ", img->val[p]);
        if (iftGetXCoord(img, p) == (img->xsize - 1)) fprintf(fp, "\n");
    }

    fclose(fp);
}


void iftInsertObject(iftImage *bin, iftImage *label, int obj_code, iftVoxel  pos)
{
    int p,q;
    iftVoxel u,v;

    if (!iftValidVoxel(label,pos))
        iftError("Invalid position for insertion", "iftInsertObject");

    if ((bin->xsize > label->xsize) ||
        (bin->ysize > label->ysize) ||
        (bin->zsize > label->zsize))
        iftError("Object cannot be inserted", "iftInsertObject");

    for (u.z=0; u.z < bin->zsize; u.z++)
        for (u.y=0; u.y < bin->ysize; u.y++)
            for (u.x=0; u.x < bin->xsize; u.x++){
                v.x = u.x + pos.x;
                v.y = u.y + pos.y;
                v.z = u.z + pos.z;
                if (iftValidVoxel(label,v)){
                    p = iftGetVoxelIndex(bin,u);
                    q = iftGetVoxelIndex(label,v);
                    if (bin->val[p])
                        label->val[q] = obj_code;
                }
            }
}


iftImage *iftRemoveLabels(const iftImage *label_img, const iftIntArray *labels) {
    iftImage *out_label_img = iftCreateImageFromImage(label_img);

    #pragma omp parallel for
    for (int p = 0; p < label_img->n; p++) {
        bool is_any_target_label = false;

        for (int i = 0; i < labels->n; i++)
            if (label_img->val[p] == labels->val[i]) {
                is_any_target_label = true;
                break;
            }

        if (!is_any_target_label)
            out_label_img->val[p] = label_img->val[p];
    }

    return out_label_img;
}


iftImage *iftRelabelImage(const iftImage *label_img) {
    iftIntArray *labels = iftGetObjectLabels(label_img);
    if (labels->n == 0)
        iftError("Label Image is Empty", "iftRelabelImage");

    // values of labels in the dict keys and their indices in the dict values
    iftDict *new_labels_dict = iftIntArrayToDict(labels->val, labels->n);
    iftImage *relabel_img = iftCreateImageFromImage(label_img);

    #pragma omp parallel for
    for (int p = 0; p < label_img->n; p++) {
        if (label_img->val[p]) {
            // the first object label of the labels array starts at 0, then we must increment the mapping 
            relabel_img->val[p] = iftGetLongValFromDict(label_img->val[p], new_labels_dict) + 1;
        }
    }

    iftDestroyDict(&new_labels_dict);
    iftDestroyIntArray(&labels);

    return relabel_img;
}


iftImage *iftCopyImage(const iftImage *img) {
    if (img == NULL)
        return NULL;

    iftImage *imgc=iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftCopyImageInplace(img, imgc);

    return(imgc);
}


void iftCopyImageInplace(const iftImage *src, iftImage *dest) {
    int p;

    iftVerifyImageDomains(src, dest, "iftCopyImageInplace");

    iftCopyVoxelSize(src, dest);

    for (p=0; p < src->n; p++)
        dest->val[p]= src->val[p];

    if (src->Cb != NULL) {
        if(dest->Cb == NULL)
            dest->Cb = iftAllocUShortArray(src->n);
        if(dest->Cr == NULL)
            dest->Cr = iftAllocUShortArray(src->n);
        for (p=0; p < src->n; p++) {
            dest->Cb[p]= src->Cb[p];
            dest->Cr[p]= src->Cr[p];
        }
    }
}


void iftCopyImageInsideRegion(const iftImage *src, const iftImage *mask, iftImage *dst) {
    iftVerifyImages(src, mask, "iftCopyImageInsideRegion");
    iftVerifyImages(src, dst, "iftCopyImageInsideRegion");
    
    if (iftIsColorImage(src) && !iftIsColorImage(dst))
        iftError("Both images (src and dst) must be colored or gray", "iftCopyImageInsideRegion");
    
    #pragma omp parallel for
    for (int p = 0; p < src->n; p++) {
        if (mask->val[p])
            dst->val[p] = src->val[p];
    }
}

bool iftIsImageEqual(const iftImage *img1, const iftImage *img2) {
    if (!iftIsVoxelSizeEqual(img1, img2)) {
        printf("Different voxel size (dx, dy, dz)\nimg1: %.2f, %.2f, %.2f\nimg2: %.2f, %.2f, %.2f\n",
               img1->dx, img1->dy, img1->dz, img2->dx, img2->dy, img2->dz);
        return false;
    }
    if (!iftIsDomainEqual(img1, img2)) {
        printf("Different domains (xsize, ysize, zsize)\nimg1: %d, %d, %d\nimg2: %d, %d, %d\n",
               img1->xsize, img1->ysize, img1->zsize, img2->xsize, img2->ysize, img2->zsize);
        return false;
    }

    for (int p = 0; p < img1->n; p++) {
        if (img1->val[p] != img2->val[p]) {
            printf("img1->val[%d] != img2->val[%d] (%d != %d)\n",
                   p, p, img1->val[p], img2->val[p]);
            // return false;
        }
    }

    return true;
}



iftImage *iftCreateCuboid(int xsize, int ysize, int zsize, float perc, int val) {
    if ((perc < 0.0) || (perc > 1.0))
        iftError("Invalid percentage %f... Try [0,1]", "iftCreateCuboid", perc);

    iftImage *img = iftCreateImage(xsize, ysize, zsize);

    iftVoxel uo;
    uo.x = (int) ((1-perc) * xsize) / 2;
    uo.y = (int) ((1-perc) * ysize) / 2;
    uo.z = (int) ((1-perc) * zsize) / 2;

    iftVoxel uf;
    uf.x = uo.x + (ceil((perc * xsize))) - 1;
    uf.y = uo.y + (ceil((perc * ysize))) - 1;
    uf.z = uo.z + (ceil((perc * zsize))) - 1;

    iftVoxel u;

    for (u.z = uo.z; u.z <= uf.z; u.z++) {
        for (u.y = uo.y; u.y <= uf.y; u.y++) {
            for (u.x = uo.x; u.x <= uf.x; u.x++){
                int p       = iftGetVoxelIndex(img, u);
                img->val[p] = val;
            }
        }
    }

    return img;
}


iftImage *iftExtractSlice(const iftImage *vol_img, iftImagePlaneOrientation plane_orientation, long slice) {
    // CHECKERS
    if (vol_img == NULL)
        iftError("Volumetric Image is NULL", "iftExtractSlice");
    if (!iftIs3DImage(vol_img))
        iftError("Image is not Volumetric (3D)", "iftExtractSlice");
    if (slice < 0)
        iftError("Invalid Slice %ld (< 0)", "iftExtractSlice", slice);

    // Gets the Slice
    iftImage *out_img = NULL;
    int q = 0;
    iftVoxel u;
    switch(plane_orientation) {
        case IFT_AXIAL:
            if (slice >= vol_img->zsize)
                iftError("Invalid Slice %ld (> the last Axial Slice)", "iftExtractSlice", slice);

            out_img = iftCreateImage(vol_img->xsize, vol_img->ysize, 1);

            u.z = slice;
            if (iftIsColorImage(vol_img)) {
                iftSetCbCr(out_img, (iftMaxImageRange(iftImageDepth(vol_img))+1)/2);
                for (u.y = 0; u.y < vol_img->ysize; u.y++)
                    for (u.x = 0; u.x < vol_img->xsize; u.x++) {
                        int p = iftGetVoxelIndex(vol_img, u);
                        out_img->val[q] = vol_img->val[p];
                        out_img->Cb[q]  = vol_img->Cb[p];
                        out_img->Cr[q]  = vol_img->Cr[p];
                        q++;
                    }
            }
            else {
                for (u.y = 0; u.y < vol_img->ysize; u.y++)
                    for (u.x = 0; u.x < vol_img->xsize; u.x++) {
                        int p = iftGetVoxelIndex(vol_img, u);
                        out_img->val[q++] = vol_img->val[p];
                    }
            }
            break;
        case IFT_CORONAL:
            if (slice >= vol_img->ysize)
                iftError("Invalid Slice %ld (> the last Coronal Slice)", "iftExtractSlice", slice);

            out_img = iftCreateImage(vol_img->xsize, vol_img->zsize, 1);

            u.y = slice;
            if (iftIsColorImage(vol_img)) {
                iftSetCbCr(out_img, (iftMaxImageRange(iftImageDepth(vol_img))+1)/2);
                for (u.z = 0; u.z < vol_img->zsize; u.z++)
                    for (u.x = 0; u.x < vol_img->xsize; u.x++) {
                        int p = iftGetVoxelIndex(vol_img, u);
                        out_img->val[q] = vol_img->val[p];
                        out_img->Cb[q]  = vol_img->Cb[p];
                        out_img->Cr[q]  = vol_img->Cr[p];
                        q++;
                    }
            }
            else {
                for (u.z = 0; u.z < vol_img->zsize; u.z++)
                    for (u.x = 0; u.x < vol_img->xsize; u.x++) {
                        int p = iftGetVoxelIndex(vol_img, u);
                        out_img->val[q++] = vol_img->val[p];
                    }
            }
            break;
        case IFT_SAGITTAL:
            if (slice >= vol_img->xsize)
                iftError("Invalid Slice %ld (> the last Sagittal Slice)", "iftExtractSlice", slice);

            out_img = iftCreateImage(vol_img->ysize, vol_img->zsize, 1);

            u.x = slice;
            if (iftIsColorImage(vol_img)) {
                iftSetCbCr(out_img, (iftMaxImageRange(iftImageDepth(vol_img))+1)/2);
                for (u.z = 0; u.z < vol_img->zsize; u.z++)
                    for (u.y = 0; u.y < vol_img->ysize; u.y++) {
                        int p = iftGetVoxelIndex(vol_img, u);
                        out_img->val[q] = vol_img->val[p];
                        out_img->Cb[q]  = vol_img->Cb[p];
                        out_img->Cr[q]  = vol_img->Cr[p];
                        q++;
                    }
            }
            else {
                for (u.z = 0; u.z < vol_img->zsize; u.z++)
                    for (u.y = 0; u.y < vol_img->ysize; u.y++) {
                        int p = iftGetVoxelIndex(vol_img, u);
                        out_img->val[q++] = vol_img->val[p];
                    }
            }
            break;
        default:
            iftError("Invalid Image Plane", "iftExtractSlice");
    }

    iftCopyVoxelSize(vol_img, out_img);


    return out_img;
}


char iftAdjacentVoxels(iftImage *img, iftAdjRel *A, iftVoxel u, iftVoxel v)
{
    int i;

    for (i=0; i < A->n; i++) {
        if ((A->dx[i]==(v.x-u.x))&&
            (A->dy[i]==(v.y-u.y))&&
            ((A->dz[i]==(v.z-u.z))))
            return 1;
    }

    return 0;
}

iftImage *iftCSVtoImage(const char *format,...)
{
    char      basename[100],*ext = NULL,newfilename[150];
    iftImage  *img=NULL;
    int       len,p;
    iftVoxel  um,u;
    FILE     *fp=NULL;

    va_list args;
    char filename[150];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    ext = iftLowerString(iftFileExt(filename));

    if (iftCompareStrings(ext,"csv")){
        len = strlen(filename);
        strncpy(basename,filename,len-4);
        basename[len-4]='\0';
        sprintf(newfilename,"%s.scn",basename);
        fp = fopen(filename,"r");
        um.x=um.y=um.z= IFT_INFINITY_INT_NEG;
        while (!feof(fp)){
            if (fscanf(fp,"%d,%d,%d",&u.x,&u.y,&u.z)!=3) iftError("Reading error", "iftCSVtoImage");
            if (u.x > um.x) um.x = u.x;
            if (u.y > um.y) um.y = u.y;
            if (u.z > um.z) um.z = u.z;
        }
        fclose(fp);
        img = iftCreateImage(um.x+10,um.y+10,um.z+10);
        fp = fopen(filename,"r");
        while (!feof(fp)){
            if(fscanf(fp,"%d,%d,%d",&u.x,&u.y,&u.z)!=3) iftError("Reading error", "iftCSVtoImage");;
            u.x=u.x+5; u.y=u.y+5; u.z=u.z+5;
            p = iftGetVoxelIndex(img,u);
            img->val[p]=255;
        }
        fclose(fp);
        iftWriteImage(img,newfilename);
    }else
        iftError(MSG_FILE_OPEN_ERROR, "iftCSVtoImage", filename);

    iftFree(ext);

    return(img);
}


void iftSetFrame(iftImage *img, int sz, int val) {
    iftSetRectangularBoxFrame(img, sz, sz, sz, val);
}


void iftSetRectangularBoxFrame(iftImage *img, int sx, int sy, int sz, int val) {
    iftBoundingBox bb = {.begin = {sx, sy, sz},
                         .end = {img->xsize - sx - 1, img->ysize - sy - 1, img->zsize - sz - 1}};
    
    iftImage *roi = iftExtractROI(img, bb);
    iftSetImage(img, val);
    iftInsertROI(roi, img, bb.begin);
    iftDestroyImage(&roi);
}


iftImage *iftAddFrame(const iftImage *img, int sz, int value)
{
    return iftAddRectangularBoxFrame(img, sz, sz, sz, value);
}

iftImage *iftAddRectangularBoxFrame(const iftImage *img, int dx, int dy, int dz, int value)
{
    iftImage *fimg;
    int p, q, xsize, ysize, zsize;
    iftVoxel u;
    int depth = iftImageDepth(img);

    if (!iftIs3DImage(img))
      dz = 0;
    
    xsize     = img->xsize+2*dx;
    ysize     = img->ysize+2*dy;
    zsize     = img->zsize+2*dz;
      
    if (iftIsColorImage(img)){
      fimg = iftCreateColorImage(xsize,ysize,zsize,depth); 
      iftSetImage(fimg,value);
      p = 0;
      for (u.z=dz; u.z < fimg->zsize-dz; u.z++)
	for (u.y=dy; u.y < fimg->ysize-dy; u.y++)
	  for (u.x=dx; u.x < fimg->xsize-dx; u.x++){
	    q = iftGetVoxelIndex(fimg,u);
	    fimg->val[q] = img->val[p];
	    fimg->Cb[q]  = img->Cb[p];
	    fimg->Cr[q]  = img->Cr[p];
	    p++;
	  }
    }else{
      fimg = iftCreateImage(xsize,ysize,zsize); 
      iftSetImage(fimg,value);
      p = 0;
      for (u.z=dz; u.z < fimg->zsize-dz; u.z++)
	for (u.y=dy; u.y < fimg->ysize-dy; u.y++)
	  for (u.x=dx; u.x < fimg->xsize-dx; u.x++){
	    q = iftGetVoxelIndex(fimg,u);
	    fimg->val[q] = img->val[p];
	    p++;
	  }
    }
    
    iftCopyVoxelSize(img,fimg);
    
    return(fimg);
}

iftImage *iftRemFrame(const iftImage *fimg, int sz) {
    return iftRemRectangularBoxFrame(fimg, sz, sz, sz);
}

iftImage *iftRemRectangularBoxFrame(const iftImage *fimg, int dx, int dy, int dz) {

  iftImage *img;
  int p, q, xsize, ysize, zsize;
  iftVoxel u;
  int depth = iftImageDepth(fimg);

  if (!iftIs3DImage(fimg))
    dz = 0;

  xsize     = fimg->xsize-2*dx;
  ysize     = fimg->ysize-2*dy;
  zsize     = fimg->zsize-2*dz;
      
  if (iftIsColorImage(fimg)){
    img = iftCreateColorImage(xsize,ysize,zsize,depth); 

    q = 0;
    for (u.z=dz; u.z < fimg->zsize-dz; u.z++)
      for (u.y=dy; u.y < fimg->ysize-dy; u.y++)
	for (u.x=dx; u.x < fimg->xsize-dx; u.x++){
	  p = iftGetVoxelIndex(fimg,u);
	  img->val[q] = fimg->val[p];
	  img->Cb[q]  = fimg->Cb[p];
	  img->Cr[q]  = fimg->Cr[p];
	  q++;
	}
    }else{
      img = iftCreateImage(xsize,ysize,zsize); 

      q = 0;
      for (u.z=dz; u.z < fimg->zsize-dz; u.z++)
	for (u.y=dy; u.y < fimg->ysize-dy; u.y++)
	  for (u.x=dx; u.x < fimg->xsize-dx; u.x++){
	    p = iftGetVoxelIndex(fimg,u);
	    img->val[q] = fimg->val[p];
	    q++;
	  }
    }
    
    iftCopyVoxelSize(fimg,img);
    
    return(img);
}

void iftSetImage(iftImage *img, int value) {
    for (int p = 0; p < img->n; p++)
        img->val[p] = value;
}

/**************************************************/
iftImage *iftGetXYSlice(const iftImage *img, int zcoord)
{
    iftImage *slice;
    iftVoxel  u;
    int       p,q;

    if ( (zcoord < 0) || (zcoord >= img->zsize))
        iftError("Invalid z coordinate", "iftGetXYSlice");

    if(iftIsColorImage(img))
        slice = iftCreateColorImage(img->xsize,img->ysize,1, iftImageDepth(img));
    else
        slice = iftCreateImage(img->xsize,img->ysize,1);

    u.z   = zcoord;
    q     = 0;
    for (u.y = 0; u.y < img->ysize; u.y++)
        for (u.x = 0; u.x < img->xsize; u.x++)
        {
            p = iftGetVoxelIndex(img,u);
            slice->val[q] = img->val[p];
            if(iftIsColorImage(img))
            {
                slice->Cb[q] = img->Cb[p];
                slice->Cr[q] = img->Cr[p];
            }
            q++;
        }
    iftCopyVoxelSize(img,slice);

    return(slice);
}

void iftPutXYSlice(iftImage *img, const iftImage *slice, int zcoord)
{
    iftVoxel  u;
    int       p,q;

    if ( (zcoord < 0) || (zcoord >= img->zsize))
        iftError("Invalid z coordinate", "iftPutXYSlice");

    if ( (img->ysize!=slice->ysize)||(img->xsize!=slice->xsize) )
        iftError("Image and slice are incompatibles", "iftPutXYSlice");

    u.z   = zcoord;
    p     = 0;
    for (u.y = 0; u.y < img->ysize; u.y++)
        for (u.x = 0; u.x < img->xsize; u.x++)
        {
            q = iftGetVoxelIndex(img,u);
            img->val[q] = slice->val[p];
            if(iftIsColorImage(img))
            {
                img->Cb[q] = slice->Cb[p];
                img->Cr[q] = slice->Cr[p];
            }
            p++;
        }
}

iftImage *iftGetZXSlice(const iftImage *img, int ycoord)
{
    iftImage *slice;
    iftVoxel  u;
    int       p,q;

    if ( (ycoord < 0) || (ycoord >= img->ysize))
        iftError("Invalid y coordinate", "iftGetZXSlice");

    if(iftIsColorImage(img))
        slice = iftCreateColorImage(img->zsize,img->xsize,1, iftImageDepth(img));
    else
        slice = iftCreateImage(img->zsize,img->xsize,1);

    u.y   = ycoord;
    q = 0;
    for (u.x = 0; u.x < img->xsize; u.x++)
        for (u.z = 0; u.z < img->zsize; u.z++)
        {
            p = iftGetVoxelIndex(img,u);
            slice->val[q] = img->val[p];
            if(iftIsColorImage(img))
            {
                slice->Cb[q] = img->Cb[p];
                slice->Cr[q] = img->Cr[p];
            }
            q++;
        }
    iftCopyVoxelSize(img,slice);

    return(slice);
}

void iftPutZXSlice(const iftImage *img, iftImage *slice, int ycoord)
{
    iftVoxel  u;
    int       p,q;

    if ( (ycoord < 0) || (ycoord >= img->ysize))
        iftError("Invalid y coordinate", "iftPutZXSlice");

    if ( (img->xsize!=slice->ysize)||(img->zsize!=slice->xsize) )
        iftError("Image and slice are incompatibles", "iftPutZXSlice");

    u.y   = ycoord;
    p     = 0;
    for (u.x = 0; u.x < img->xsize; u.x++)
        for (u.z = 0; u.z < img->zsize; u.z++)
        {
            q = iftGetVoxelIndex(img,u);
            img->val[q] = slice->val[p];
            if(iftIsColorImage(img))
            {
                img->Cb[q] = slice->Cb[p];
                img->Cr[q] = slice->Cr[p];
            }
            p++;
        }
}

iftImage *iftGetYZSlice(const iftImage *img, int xcoord)
{
    iftImage *slice;
    iftVoxel  u;
    int       p,q;

    if ( (xcoord < 0) || (xcoord >= img->xsize))
        iftError("Invalid x coordinate", "iftGetYZSlice");

    if(iftIsColorImage(img))
        slice = iftCreateColorImage(img->ysize,img->zsize,1, iftImageDepth(img));
    else
        slice = iftCreateImage(img->ysize,img->zsize,1);

    u.x   = xcoord;
    q     = 0;
    for (u.z = 0; u.z < img->zsize; u.z++)
        for (u.y = 0; u.y < img->ysize; u.y++)
        {
            p = iftGetVoxelIndex(img,u);
            slice->val[q] = img->val[p];
            if(iftIsColorImage(img))
            {
                slice->Cb[q] = img->Cb[p];
                slice->Cr[q] = img->Cr[p];
            }
            q++;
        }
    iftCopyVoxelSize(img,slice);

    return(slice);
}

void iftPutYZSlice(iftImage *img, iftImage *slice, int xcoord)
{
    iftVoxel  u;
    int       p,q;

    if ( (xcoord < 0) || (xcoord >= img->xsize))
        iftError("Invalid x coordinate", "iftPutYZSlice");

    if ( (img->zsize!=slice->ysize)||(img->ysize!=slice->xsize) )
        iftError("Image and slice are incompatibles", "iftPutYZSlice");

    u.x   = xcoord;
    p     = 0;
    for (u.z = 0; u.z < img->zsize; u.z++)
        for (u.y = 0; u.y < img->ysize; u.y++)
        {
            q = iftGetVoxelIndex(img,u);
            img->val[q] = slice->val[p];
            if(iftIsColorImage(img))
            {
                img->Cb[q] = slice->Cb[p];
                img->Cr[q] = slice->Cr[p];
            }
            p++;
        }
}

iftVoxel iftGeometricCenterVoxel(const iftImage *bin_mask) {
    iftPoint c = {0, 0, 0};
    int n = 0;

    for (int p = 0; p < bin_mask->n; p++) {
        if (bin_mask->val[p] != 0) {
            iftVoxel u = iftGetVoxelCoord(bin_mask, p);
            c.x += u.x;
            c.y += u.y;
            c.z += u.z;
            n++;
        }
    }

    if (n == 0)
        iftError("No object inside the image", "iftGeometricCenterVoxel");

    iftVoxel gc;
    gc.x = iftRound(c.x / n);
    gc.y = iftRound(c.y / n);
    gc.z = iftRound(c.z / n);

    // geometric center is out of the object
    // gets the closest object's voxel to it
    if (iftImgVoxelVal(bin_mask, gc) == 0) {
        long min_dist = IFT_INFINITY_LONG;
        iftVoxel closest_voxel = {-1, -1, -1};

        for (int p = 0; p < bin_mask->n; p++) {
            if (bin_mask->val[p] > 0) {
                iftVoxel u = iftGetVoxelCoord(bin_mask, p);
                long dist = iftSquaredVoxelDistance(u, gc);

                if (dist < min_dist) {
                    closest_voxel = u;
                    min_dist = dist;
                }
            }
        }
        gc = closest_voxel;
    }

    return gc;
}


iftVoxelArray *iftGeometricCenterVoxelsLabelImage(const iftImage *label_img) {
    int n_labels = iftMaximumValue(label_img);

    iftVoxelArray *gcs = iftCreateVoxelArray(n_labels + 1);
    iftIntArray *n_voxels = iftCreateIntArray(n_labels + 1);
    iftImage *gc_img = iftCreateImageFromImage(label_img);

    for (int p = 0; p < label_img->n; p++) {
        iftVoxel u = iftGetVoxelCoord(label_img, p);
        int label = label_img->val[p];

        gcs->val[label].x += u.x;
        gcs->val[label].y += u.y;
        gcs->val[label].z += u.z;
        n_voxels->val[label]++;
    }

    for (int label = 1; label <= n_labels; label++) {
        if (n_voxels->val[label] > 0) {
            gcs->val[label].x = iftRound(gcs->val[label].x / (float) n_voxels->val[label]);
            gcs->val[label].y = iftRound(gcs->val[label].y / (float) n_voxels->val[label]);
            gcs->val[label].z = iftRound(gcs->val[label].z / (float) n_voxels->val[label]);

            // geometric center is out of the object, then it gets the closest voxel to it from the object
            // ps: this code is not the fastest one. An IFT could be used instead.
            if (iftImgVoxelVal(label_img, gcs->val[label]) != label) {
                int min_dist = IFT_INFINITY_INT;
                iftVoxel closest_voxel;
            
                for (int p = 0; p < label_img->n; p++) {
                    if (label_img->val[p] == label) {
                        iftVoxel u = iftGetVoxelCoord(label_img, p);
                        int dist   = iftSquaredVoxelDistance(u, gcs->val[label]);

                        if (dist < min_dist) {
                            closest_voxel = iftGetVoxelCoord(label_img, p);
                            min_dist = dist;
                        }
                    }
                }
                gcs->val[label] = closest_voxel;
            }
            gc_img->val[iftGetVoxelIndex(label_img, gcs->val[label])] = label;
        }
    }
    iftDestroyIntArray(&n_voxels);

    return gcs;
}


int iftObjectDiagonal(iftImage *obj)
{
    int p,xsize,ysize,zsize,diag;
    iftVoxel u,min,max;

    min.x = min.y = min.z = IFT_INFINITY_INT;
    max.x = max.y = max.z = IFT_INFINITY_INT_NEG;
    for (u.z=0; u.z < obj->zsize; u.z++)
        for (u.y=0; u.y < obj->ysize; u.y++)
            for (u.x=0; u.x < obj->xsize; u.x++){
                p = iftGetVoxelIndex(obj,u);
                if (obj->val[p]!=0){
                    if (u.x < min.x) min.x = u.x;
                    if (u.y < min.y) min.y = u.y;
                    if (u.z < min.z) min.z = u.z;
                    if (u.x > max.x) max.x = u.x;
                    if (u.y > max.y) max.y = u.y;
                    if (u.z > max.z) max.z = u.z;
                }
            }
    xsize = max.x - min.x + 1;
    ysize = max.y - min.y + 1;
    zsize = max.z - min.z + 1;
    diag  = (int)(sqrtf(xsize*xsize + ysize*ysize + zsize*zsize)+0.5);

    return(diag);
}


iftPoint iftGeometricCenter(iftImage *obj)
{
    iftVoxel u;
    iftPoint c;
    int p;
    unsigned long n=0;

    c.x = c.y = c.z = 0.0;
    for (u.z=0; u.z < obj->zsize; u.z++)
        for (u.y=0; u.y < obj->ysize; u.y++)
            for (u.x=0; u.x < obj->xsize; u.x++) {
                p = iftGetVoxelIndex(obj,u);
                if (obj->val[p]!=0){
                    c.x += u.x;
                    c.y += u.y;
                    c.z += u.z;
                    n++;
                }
            }
    if (n==0)
        iftError("Empty image", "iftGeometricCenter");

    c.x /= n;
    c.y /= n;
    c.z /= n;

    return(c);
}


iftImage *iftImageGradientMagnitude(const iftImage *img, iftAdjRel *Ain)
{
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(img))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }
    else A = Ain;

    float   dist,gx,gy,gz, g, gmax;
    float   gxCb , gyCb , gzCb, gxCr , gyCr , gzCr;
    iftVoxel   u,v;
    int i;
    float     *mag =iftAllocFloatArray(A->n), *weight = iftAllocFloatArray(A->n);
    float     _2sigma2;
    int        dx, dy, dz;
    iftImage  *grad=iftCreateImage(img->xsize,img->ysize,img->zsize);

    iftCopyVoxelSize(img,grad);

    iftMaxAdjShifts(A, &dx, &dy, &dz);
    _2sigma2 = 2.0*(dx*dx+dy*dy+dz*dz)/9.0;
    for (i=0; i < A->n; i++){
        mag[i]=sqrtf(A->dx[i]*A->dx[i]+A->dy[i]*A->dy[i]+A->dz[i]*A->dz[i]);
        weight[i]=exp(-mag[i]/_2sigma2)/mag[i];
    }

    if ( (img->Cb == NULL) && (img->Cr == NULL) ){
        for (int z=0; z < img->zsize; z++)
            for (int y=0; y < img->ysize; y++)
#pragma omp parallel for shared(z,y,A,img,weight,grad) private(u,v,gx,gy,gz,i,dist)
                    for (int x=0; x < img->xsize; x++) {
                        u.x=x; u.y=y; u.z=z;
                        int p = iftGetVoxelIndex(img,u);
                        gx = gy = gz = 0.0;
                        for (i=1; i < A->n; i++) {
                            v.x = u.x + A->dx[i];
                            v.y = u.y + A->dy[i];
                            v.z = u.z + A->dz[i];
                            if (iftValidVoxel(img,v)){
                                int q = iftGetVoxelIndex(img,v);
                                dist = img->val[q]-img->val[p];
                                gx  += dist*A->dx[i]*weight[i];
                                gy  += dist*A->dy[i]*weight[i];
                                gz  += dist*A->dz[i]*weight[i];
                            }
                        }
                        grad->val[p]=(int)sqrtf(gx*gx + gy*gy + gz*gz);
                    }
    }else{ // colored image
        for (int z=0; z < img->zsize; z++)
            for (int y=0; y < img->ysize; y++)
#pragma omp parallel for shared(z,y,A,img,weight,grad) private(u,v,gx,gy,gz,i,dist,gxCb,gyCb,gzCb,gxCr,gyCr,gzCr,gmax,g)
                    for (int x=0; x < img->xsize; x++) {
                        u.x=x; u.y=y; u.z=z;
                        int p = iftGetVoxelIndex(img,u);
                        gx = gy = gz = 0.0;
                        gxCb = gyCb = gzCb = 0.0;
                        gxCr = gyCr = gzCr = 0.0;
                        for (i=1; i < A->n; i++) {
                            v.x = u.x + A->dx[i];
                            v.y = u.y + A->dy[i];
                            v.z = u.z + A->dz[i];
                            if (iftValidVoxel(img,v)){
                                int q = iftGetVoxelIndex(img,v);
                                dist = img->val[q]-img->val[p];
                                gx  += dist*A->dx[i]*weight[i];
                                gy  += dist*A->dy[i]*weight[i];
                                gz  += dist*A->dz[i]*weight[i];
                                dist = img->Cb[q]-img->Cb[p];
                                gxCb  += dist*A->dx[i]*weight[i];
                                gyCb  += dist*A->dy[i]*weight[i];
                                gzCb  += dist*A->dz[i]*weight[i];
                                dist = img->Cr[q]-img->Cr[p];
                                gxCr  += dist*A->dx[i]*weight[i];
                                gyCr  += dist*A->dy[i]*weight[i];
                                gzCr  += dist*A->dz[i]*weight[i];
                            }
                        }
                        gmax = sqrtf(gx*gx + gy*gy + gz*gz);
                        g    = sqrtf(gxCb*gxCb + gyCb*gyCb + gzCb*gzCb);
                        if (g > gmax)
                            gmax = g;
                        g    = sqrtf(gxCr*gxCr + gyCr*gyCr + gzCr*gzCr);
                        if (g > gmax)
                            gmax = g;
                        grad->val[p] = (int)gmax;
                    }
    }

    iftFree(mag);
    iftFree(weight);

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }

    return(grad);
}


void iftGetDisplayRange(iftImage *img, int *lower, int *higher)
{
    int p, i, n;
    double *hist;

    int img_min_val = iftMinimumValue(img);
    int img_max_val = iftMaximumValue(img);
    n=img_max_val-img_min_val+1;

    hist = iftAllocDoubleArray(n);

    /* Compute histogram */
    for (p=0; p < img->n; p++)
        hist[img->val[p]-img_min_val]++;

    /* Compute normalized and accumulated histogram */

    hist[0] /= img->n;
    for (i=1; i < n; i++) {
        hist[i] = hist[i]/img->n + hist[i-1];
    }

    /* Compute lower value */
    for (i=0; i < n; i++)
        if (hist[i] > 0.020){
            *lower = (i+img_min_val);
            break;
        }

    for (i=n-1; i >= 0; i--)
        if (hist[i] < 0.998){
            *higher = (i+img_min_val);
            break;
        }

    iftFree(hist);
}


inline iftVoxel iftGetVoxelCoord(const iftImage *img, int p)
{
    /* old
     * u.x = (((p) % (((img)->xsize)*((img)->ysize))) % (img)->xsize)
     * u.y = (((p) % (((img)->xsize)*((img)->ysize))) / (img)->xsize)
     * u.z = ((p) / (((img)->xsize)*((img)->ysize)))
     */
    iftVoxel u;
    div_t res1 = div(p, img->xsize * img->ysize);
    div_t res2 = div(res1.rem, img->xsize);

    u.x = res2.rem;
    u.y = res2.quot;
    u.z = res1.quot;
    u.t = 0;
    
    return u;
}

iftImage *iftLuminance(const iftImage *img)
{
    iftImage *Y=iftCreateImage(img->xsize,img->ysize,img->zsize);
    int p;

    for (p=0; p < img->n; p++)
        Y->val[p] = img->val[p];


    iftCopyVoxelSize(img,Y);

    return(Y);
}

iftImage *iftImageCb(const iftImage *img)
{
    iftImage *Cb=iftCreateImage(img->xsize,img->ysize,img->zsize);
    int p;

    if (img->Cb == NULL)
        iftError("There is no color component", "iftImageCb");

    for (p=0; p < img->n; p++)
        Cb->val[p] = img->Cb[p];


    iftCopyVoxelSize(img,Cb);

    return(Cb);
}

iftImage *iftImageCr(const iftImage *img)
{
    iftImage *Cr=iftCreateImage(img->xsize,img->ysize,img->zsize);
    int p;

    if (img->Cr == NULL)
        iftError("There is no color component", "iftImageCb");

    for (p=0; p < img->n; p++)
        Cr->val[p] = img->Cr[p];

    iftCopyVoxelSize(img,Cr);

    return(Cr);
}

iftImage *iftImageRed(iftImage *img)
{
    iftImage *red=iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftColor RGB,YCbCr;
    int p;
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));

    if ((img->Cr == NULL)||(img->Cb == NULL))
        iftError("There are no color components", "iftImageRed");

    for (p=0; p < img->n; p++) {
        YCbCr.val[0] = img->val[p];
        YCbCr.val[1] = img->Cb[p];
        YCbCr.val[2] = img->Cr[p];
        RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
        red->val[p] = RGB.val[0];
    }

    iftCopyVoxelSize(img,red);

    return(red);
}

iftImage *iftImageGreen(iftImage *img)
{
    iftImage *green=iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftColor RGB,YCbCr;
    int p;
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));

    if ((img->Cr == NULL)||(img->Cb == NULL))
        iftError("There are no color components", "iftImageGreen");

    for (p=0; p < img->n; p++) {
        YCbCr.val[0] = img->val[p];
        YCbCr.val[1] = img->Cb[p];
        YCbCr.val[2] = img->Cr[p];
        RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
        green->val[p] = RGB.val[1];
    }

    iftCopyVoxelSize(img,green);

    return(green);
}

iftImage *iftImageGray(iftImage *img)
{
    iftImage *gray=iftCreateImage(img->xsize,img->ysize,img->zsize);
    int p;


    for (p=0; p < img->n; p++) {
        gray->val[p] = img->val[p];
    }

    iftCopyVoxelSize(img,gray);

    return(gray);
}

iftImage *iftImageBlue(iftImage *img)
{
    iftImage *blue=iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftColor RGB,YCbCr;
    int p;
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));

    if ((img->Cr == NULL)||(img->Cb == NULL))
        iftError("There are no color components", "iftImageBlue");

    for (p=0; p < img->n; p++) {
        YCbCr.val[0] = img->val[p];
        YCbCr.val[1] = img->Cb[p];
        YCbCr.val[2] = img->Cr[p];
        RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
        blue->val[p] = RGB.val[2];
    }

    iftCopyVoxelSize(img,blue);

    return(blue);
}

iftImage *iftImageHue(iftImage *img)
{
    iftImage *hue = iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftColor RGB,YCbCr,HSV;
    int p;
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));

    if ((img->Cr == NULL)||(img->Cb == NULL))
        iftError("There are no color components", "iftImageHue");

    for (p=0; p < img->n; p++) {
        YCbCr.val[0] = img->val[p];
        YCbCr.val[1] = img->Cb[p];
        YCbCr.val[2] = img->Cr[p];
        RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
        HSV = iftRGBtoHSV(RGB,normalization_value);
        hue->val[p] = HSV.val[0];
    }

    iftCopyVoxelSize(img,hue);

    return(hue);
}

iftImage *iftImageSaturation(iftImage *img)
{
    iftImage *saturation=iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftColor RGB,YCbCr,HSV;
    int p;
    int normalization_value = iftNormalizationValue(iftMaximumValue(img));

    if ((img->Cr == NULL)||(img->Cb == NULL))
        iftError("There are no color components", "iftImageSaturation");

    for (p=0; p < img->n; p++) {
        YCbCr.val[0] = img->val[p];
        YCbCr.val[1] = img->Cb[p];
        YCbCr.val[2] = img->Cr[p];
        RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
        HSV = iftRGBtoHSV(RGB,normalization_value);
        saturation->val[p] = HSV.val[1];
    }

    iftCopyVoxelSize(img,saturation);

    return(saturation);
}

iftImage *iftCreateGaussian(int xsize, int ysize, int zsize, iftVoxel mean, float stdev, int maxval)
{
    iftImage *img   = iftCreateImage(xsize,ysize,zsize);
    iftVoxel u;
    int p;
    float K, variance=stdev*stdev;

    if (variance <= 0.0)
        iftError("Invalid variance", "iftCreateGaussian");

    K=2.0*variance;

    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++){
                p = iftGetVoxelIndex(img,u);
                img->val[p]=(int)(maxval*expf(-(float)iftSquaredVoxelDistance(u,mean)/K));
            }

    return(img);
}

iftImage *iftRegionBorders(const iftImage *label_img) {
    iftImage *borders = iftCreateImageFromImage(label_img);

    iftAdjRel *A = (iftIs3DImage(label_img)) ? iftSpheric(1.0) : iftCircular(1.0);

    iftVoxel u;
    for (u.z = 0; u.z < label_img->zsize; u.z++)
        for (u.y = 0; u.y < label_img->ysize; u.y++)
            for (u.x = 0; u.x < label_img->xsize; u.x++){
                int p = iftGetVoxelIndex(label_img, u);

                // if p is an object voxel
                if (label_img->val[p]) {
                    int label = label_img->val[p];

                    for (int i = 1; i < A->n; i++) {
                        iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                        if (iftValidVoxel(label_img, v)) {
                            int q = iftGetVoxelIndex(label_img,v);

                            if (label_img->val[q] != label_img->val[p]) {
                                borders->val[p] = label;
                                break;
                            }
                        }
                    }
                }
            }
    iftDestroyAdjRel(&A);

    return borders;
}

iftImage *iftCreateImageWithGaussianHistogram(int xsize, int ysize, int zsize, float mean, float stdev, int maxval)
{
    iftImage *img=iftCreateImage(xsize,ysize,zsize);
    int       p, i, j;
    float    *gauss= iftAllocFloatArray(maxval+1), area, variance=stdev*stdev;


    if ((mean > maxval)||(mean < 0))
        iftError("The mean value must be positive and lower than the maximum value",
                 "iftCreateImageWithGaussianHistogram");

    /* Create a Gaussian histogram with area = number of voxels */

    for (i=0, area=0.0; i <= maxval; i++){
        gauss[i] = exp(-(i-mean)*(i-mean)/(2*variance));
        area    += gauss[i];
    }
    for (i=0; i <= maxval; i++){
        gauss[i] = img->n*gauss[i]/area;
    }

    /* Create image with Gaussian histogram */

    p = 0;
    for (i=0; (i <= maxval)&&(p < img->n); i++){
        for (j=0; (j < gauss[i])&&(p < img->n); j++){
            img->val[p]=i; p++;
        }
    }


    iftFree(gauss);
    return(img);
}


iftImage *iftCreateImageWithTwoGaussiansHistogram(int xsize, int ysize, int zsize, float mean1, float stdev1, float mean2, float stdev2, int maxval)
{
    iftImage *img=iftCreateImage(xsize,ysize,zsize);
    int       p, i, j;
    float    *gauss=iftAllocFloatArray(maxval+1), area;
    float     variance1=stdev1*stdev1,variance2=stdev2*stdev2;


    if ((mean1 > maxval)||(mean1 < 0)||(mean2 > maxval)||(mean2 < 0))
        iftError("The mean values must be positive and lower than the maximum value",
                 "iftCreateImageWithTwoGaussiansHistogram");

    /* Create a two-Gaussians histogram with area = number of voxels */

    for (i=0, area=0.0; i <= maxval; i++){
        gauss[i] = exp(-(i-mean1)*(i-mean1)/(2*variance1)) + exp(-(i-mean2)*(i-mean2)/(2*variance2));
        area    += gauss[i];
    }
    for (i=0; i <= maxval; i++){
        gauss[i] = img->n*gauss[i]/area;
    }

    /* Create image with Gaussian histogram */

    p = 0;
    for (i=0; (i <= maxval)&&(p < img->n); i++){
        for (j=0; (j < gauss[i])&&(p < img->n); j++){
            img->val[p]=i; p++;
        }
    }


    iftFree(gauss);
    return(img);
}



iftImage *iftReadRawSlices(char *basename, int first, int last, int xsize, int ysize, int bits_per_voxel)
{
    FILE     *fp;
    int       p, offset, n = xsize*ysize, i, zsize = last - first + 1;
    iftImage *img;
    char      filename[200];
    uchar    *data8;
    ushort   *data16;
    int    *data32;

    img = iftCreateImage(xsize,ysize,zsize);

    switch(bits_per_voxel) {

        case 8:
            data8 = iftAllocUCharArray(n);
            for (i=first; i <= last; i++) {
                offset = n*(i-first);
                sprintf(filename,"%s%03d.raw",basename,i);
                fp  = fopen(filename,"rb");
                if(fread(data8,sizeof(uchar),n,fp)!=n)
                    iftError("Reading error", "iftReadRawSlices");
                for (p=0; p < n; p++) {
                    img->val[p+offset] = data8[p];
                }
                fclose(fp);
            }
            iftFree(data8);
            break;

        case 16:
            data16 = iftAllocUShortArray(n);
            for (i=first; i <= last; i++) {
                offset = n*(i-first);
                sprintf(filename,"%s%03d.raw",basename,i);
                fp  = fopen(filename,"rb");
                if(fread(data16,sizeof(ushort),n,fp)!=n)
                    iftError("Reading error", "iftReadRawSlices");

                for (p=0; p < n; p++) {
                    img->val[p+offset] = data16[p];

                }

                fclose(fp);
            }
            iftFree(data16);
            break;

        case 32:
            data32 = iftAllocIntArray(n);
            for (i=first; i <= last; i++) {
                offset = n*(i-first);
                sprintf(filename,"%s%03d.raw",basename,i);
                fp  = fopen(filename,"rb");
                if(fread(data32,sizeof(int),n,fp)!=n)
                    iftError("Reading error", "iftReadRawSlices");
                for (p=0; p < n; p++) {
                    img->val[p+offset] = data32[p];
                }
                fclose(fp);
            }
            iftFree(data32);
            break;

        default:
            iftError("Invalid number of bits per voxel", "iftReadRawSlices");
    }


    return(img);
}

void iftWriteRawSlices(iftImage *img, char *basename) {
    FILE *fp;
    int  p, offset, n = img->xsize * img->ysize, i;
    char filename[200];

    int img_max_val = iftMaximumValue(img);

    if (img_max_val <= 255) {
        uchar *data8 = iftAllocUCharArray(n);
        for (i = 0; i < img->zsize; i++) {
            offset = n * i;
            for (p = 0; p < n; p++)
                data8[p] = (uchar) img->val[p + offset];
            sprintf(filename, "%s%03d.raw", basename, i);
            fp = fopen(filename, "wb");
            fwrite(data8, n, sizeof(uchar), fp);
            fclose(fp);
        }
        iftFree(data8);
    } else {
        if (img_max_val <= 4095) {
            ushort *data16 = iftAllocUShortArray(n);
            for (i = 0; i < img->zsize; i++) {
                offset = n * i;
                for (p = 0; p < n; p++)
                    data16[p] = (ushort) img->val[p + offset];
                sprintf(filename, "%s%03d.raw", basename, i);
                fp = fopen(filename, "wb");
                fwrite(data16, n, sizeof(ushort), fp);
                fclose(fp);
            }
            iftFree(data16);
        } else {
            int *data32 = iftAllocIntArray(n);
            for (i = 0; i < img->zsize; i++) {
                offset = n * i;
                for (p = 0; p < n; p++)
                    data32[p] = (int) img->val[p + offset];
                sprintf(filename, "%s%03d.raw", basename, i);
                fp = fopen(filename, "wb");
                fwrite(data32, n, sizeof(int), fp);
                fclose(fp);
            }
            iftFree(data32);
        }
    }
}

void iftWriteRawScene(iftImage *img, const char *format, ...) {
    FILE *fp;
    int  p, offset, n = img->xsize * img->ysize, i;

    int img_max_val = iftMaximumValue(img);
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    if (img_max_val <= 255) {
        uchar *data8 = iftAllocUCharArray(n);
        fp = fopen(filename, "wb");

        for (i = 0; i < img->zsize; i++) {
            offset = n * i;
            for (p = 0; p < n; p++)
                data8[p] = (uchar) img->val[p + offset];

            fwrite(data8, n, sizeof(uchar), fp);
        }
        free(data8);
        fclose(fp);
    } else {
        if (img_max_val <= 65535) {
            ushort *data16 = iftAllocUShortArray(n);
            fp = fopen(filename, "wb");

            for (i = 0; i < img->zsize; i++) {
                offset = n * i;
                for (p = 0; p < n; p++)
                    data16[p] = (ushort) img->val[p + offset];
                fwrite(data16, n, sizeof(ushort), fp);
            }
            free(data16);
            fclose(fp);

        } else {
            int *data32 = iftAllocIntArray(n);
            fp = fopen(filename, "wb");

            for (i = 0; i < img->zsize; i++) {
                offset = n * i;
                for (p = 0; p < n; p++)
                    data32[p] = (int) img->val[p + offset];
                fwrite(data32, n, sizeof(int), fp);
            }
            free(data32);
            fclose(fp);
        }
    }
}

iftImage *iftReadRawScene(int xsize, int ysize, int zsize, int bits_per_voxel, const char *format, ...)
{
    FILE     *fp;
    int       p, n;
    iftImage *img;
    uchar    *data8;
    ushort   *data16;
    int      *data32;
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    img = iftCreateImage(xsize,ysize,zsize);
    n = img->n;

    switch(bits_per_voxel) {
        case 8:
            data8 = iftAllocUCharArray(n);
            fp  = fopen(filename,"rb");
            if(fread(data8,sizeof(uchar),n,fp)!=n)
                iftError("Reading error", "iftReadRawScene");
            for (p=0; p < n; p++) {
                img->val[p] = data8[p];
            }
            fclose(fp);
            iftFree(data8);
            break;

        case 16:
            data16 = iftAllocUShortArray(n);
            fp  = fopen(filename,"rb");
            if(fread(data16,sizeof(ushort),n,fp)!=n)
                iftError("Reading error", "iftReadRawScene");
            for (p=0; p < n; p++) {
                img->val[p] = data16[p];
            }
            fclose(fp);
            iftFree(data16);
            break;

        case 32:
            data32 = iftAllocIntArray(n);
            fp  = fopen(filename,"rb");
            if(fread(data32,sizeof(int),n,fp)!=n)
                iftError("Reading error", "iftReadRawScene");
            for (p=0; p < n; p++) {
                img->val[p] = data32[p];
            }
            fclose(fp);

            iftFree(data32);
            break;

        default:
            iftError("Invalid number of bits per voxel", "iftReadRawScene");
    }


    return(img);
}

void iftWriteRawSceneWithInfoOnFilename(iftImage *img, const char *format, ...) {
    int nbits;
    char filename[IFT_STR_DEFAULT_SIZE];
    va_list args;

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    char *bname = iftBasename(filename);
    const char *ext = iftFileExt(filename);
    char suffix[IFT_STR_DEFAULT_SIZE];
    int img_max_val = iftMaximumValue(img);
    int img_min_val = iftMinimumValue(img);

    if (img_min_val >=0 && img_max_val <= 255) {
        nbits = 8;
    } else if (img_min_val >=0 && img_max_val <= 65535) {
        nbits = 16;
    } else {
        nbits = 32;
    }

    sprintf(suffix, "sfx_%d_%d_%d_%d_%f_%f_%f", img->xsize, img->ysize, img->zsize, nbits,
            img->dx, img->dy, img->dz);

    iftWriteRawScene(img, "%s_%s%s", bname, suffix, ext);

    iftFree(bname);
}

iftImage *iftReadRawSceneWithInfoOnFilename(const char *format, ...) {
    char filename[IFT_STR_DEFAULT_SIZE];
    va_list args;

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    char *suffix = iftSplitStringAt(filename, "sfx_", -1);
    int xsize, ysize, zsize, nbits;
    float dx, dy, dz;

    if(sscanf(suffix, "%d_%d_%d_%d_%f_%f_%f", &xsize, &ysize, &zsize, &nbits,
              &dx, &dy, &dz) != 7)
        iftError("Expected raw image suffix of \"sfx_<xsize>_<ysize>_<zsize>_<nbits>_<dx>_<dy>_<dz>\" not found!",
                 "iftReadRawSceneWithInfoOnFilename");

    iftImage *img = iftReadRawScene(xsize, ysize, zsize, nbits, filename);

    img->dx = dx;
    img->dy = dy;
    img->dz = dz;

    iftFree(suffix);

    return img;
}


iftImage *iftReadSlices(const iftFileSet *files)
{
    iftImage *img, *slice = NULL;

    if(!iftIsImageFile(files->files[0]->path))
        iftError("File %s is not an image!", "iftReadSlices", files->files[0]->path);

    slice = iftReadImageByExt(files->files[0]->path);

    if(iftIs3DImage(slice))
        iftError("Cannot use 3D image %s as a slice for computing another 3D volume!", "iftReadSlices",
                 files->files[0]->path);

    img = iftCreateImage(slice->xsize,slice->ysize,files->n);

    iftPutXYSlice(img, slice, 0);

    iftDestroyImage(&slice);

    for(size_t i = 1; i < files->n; i++) {
        if(!iftIsImageFile(files->files[i]->path))
            iftError("File %s is not an image!", "iftReadSlices", files->files[i]->path);

        slice = iftReadImageByExt(files->files[i]->path);

        if(iftIs3DImage(slice))
            iftError("Cannot use 3D image %s as a slice for computing another 3D volume!", "iftReadSlices",
                     files->files[i]->path);

        iftPutXYSlice(img, slice, i);
        iftDestroyImage(&slice);
    }

    return(img);
}

iftImage *iftExtractGreyObject(iftImage *image)
{
    int p,q;
    iftVoxel uo,uf,u;
    iftImage *result=NULL;

    uo.x = uo.y = uo.z = IFT_INFINITY_INT;
    uf.x = uf.y = uf.z = IFT_INFINITY_INT_NEG;
    for (u.z=0; u.z < image->zsize; u.z++)
        for (u.y=0; u.y < image->ysize; u.y++)
            for (u.x=0; u.x < image->xsize; u.x++){
                p = iftGetVoxelIndex(image,u);
                if (image->val[p] > 0){
                    if (u.x < uo.x) uo.x = u.x;
                    if (u.y < uo.y) uo.y = u.y;
                    if (u.z < uo.z) uo.z = u.z;
                    if (u.x > uf.x) uf.x = u.x;
                    if (u.y > uf.y) uf.y = u.y;
                    if (u.z > uf.z) uf.z = u.z;
                }
            }

    result = iftCreateImage(uf.x-uo.x+1,uf.y-uo.y+1,uf.z-uo.z+1);

    q = 0;
    for (u.z=uo.z; u.z <= uf.z; u.z++)
        for (u.y=uo.y; u.y <= uf.y; u.y++)
            for (u.x=uo.x; u.x <= uf.x; u.x++){
                p = iftGetVoxelIndex(image,u);
                result->val[q]=image->val[p];
                q++;
            }

    iftCopyVoxelSize(image,result);

    return(result);
}

iftImage *iftExtractGreyObjectPos(iftImage *image, iftVoxel *pos)
{
    int p,q;
    iftVoxel uo,uf,u;
    iftImage *result=NULL;

    uo.x = uo.y = uo.z = IFT_INFINITY_INT;
    uf.x = uf.y = uf.z = IFT_INFINITY_INT_NEG;
    for (u.z=0; u.z < image->zsize; u.z++)
        for (u.y=0; u.y < image->ysize; u.y++)
            for (u.x=0; u.x < image->xsize; u.x++){
                p = iftGetVoxelIndex(image,u);
                if (image->val[p] > 0){
                    if (u.x < uo.x) uo.x = u.x;
                    if (u.y < uo.y) uo.y = u.y;
                    if (u.z < uo.z) uo.z = u.z;
                    if (u.x > uf.x) uf.x = u.x;
                    if (u.y > uf.y) uf.y = u.y;
                    if (u.z > uf.z) uf.z = u.z;
                }
            }

    result = iftCreateImage(uf.x-uo.x+1,uf.y-uo.y+1,uf.z-uo.z+1);

    q = 0;
    for (u.z=uo.z; u.z <= uf.z; u.z++)
        for (u.y=uo.y; u.y <= uf.y; u.y++)
            for (u.x=uo.x; u.x <= uf.x; u.x++){
                p = iftGetVoxelIndex(image,u);
                result->val[q]=image->val[p];
                q++;
            }

    iftCopyVoxelSize(image,result);
    *pos = uo;

    return(result);
}



iftImage  *iftCrop2DImageByBorder(iftImage *img, int border)
{
    iftImage *imgc=iftCreateImage(img->xsize - 2*border,img->ysize - 2*border,img->zsize);
    int pcrop;
    int x1,y1,z1,x2,y2,z2;

    x1 = border;
    y1 = border;
    z1 = 0;
    x2 = img->xsize - border;
    y2 = img->ysize - border;
    z2 = 1;

    iftCopyVoxelSize(img,imgc);
    pcrop = 0;
    for (int y = y1; y < y2; y++)
        for (int x = x1; x < x2; x++)
            for (int z = z1; z < z2 ; z++){
                int p = x + img->tby[y] + img->tbz[z];
                imgc->val[pcrop]=img->val[p];
                pcrop++;
            }

    if (img->Cb != NULL) {
        imgc->Cb = iftAllocUShortArray(imgc->n);
        imgc->Cr = iftAllocUShortArray(imgc->n);
        pcrop = 0;
        for (int y = y1; y < y2; y++)
            for (int x = x1; x < x2; x++)
                for (int z = z1; z < z2 ; z++){
                    int p = x + img->tby[y] + img->tbz[z];
                    imgc->Cb[pcrop]=img->Cb[p];
                    imgc->Cr[pcrop]=img->Cr[p];
                    pcrop++;
                }
    }
    return(imgc);
}

iftImage  *iftCropImage(iftImage *img, int dx, int dy, int dz)
{
    iftImage *imgc;
    int p,q;
    int xsize,ysize,zsize;
    iftVoxel uo, uf, u, v;

    xsize = img->xsize - 2*dx;
    ysize = img->ysize - 2*dy;
    zsize = img->zsize - 2*dz;

    uo.x  = dx;
    uo.y  = dy;
    uo.z  = dz;
    uf.x  = img->xsize - dx - 1;
    uf.y  = img->ysize - dy - 1;
    uf.z  = img->zsize - dz - 1;

    if (!iftIs3DImage(img)) { /* 2D image */
        if (iftIsColorImage(img)) {
            imgc = iftCreateColorImage(xsize,ysize,0, iftImageDepth(img));
            q = 0;
            u.z=v.z=0;
            for (u.y = uo.y, v.y=0; u.y <= uf.y; u.y++,v.y++)
                for (u.x = uo.x,v.x=0; u.x <= uf.x; u.x++,v.x++){
                    p = iftGetVoxelIndex(img,u);
                    q = iftGetVoxelIndex(imgc,v);
                    imgc->val[q] = img->val[p];
                    imgc->Cb[q]  = img->Cb[p];
                    imgc->Cr[q]  = img->Cr[p];
                }
        } else {
            imgc = iftCreateImage(xsize,ysize,0);
            q = 0;
            u.z=v.z=0;
            for (u.y = uo.y, v.y=0; u.y <= uf.y; u.y++,v.y++)
                for (u.x = uo.x,v.x=0; u.x <= uf.x; u.x++,v.x++){
                    p = iftGetVoxelIndex(img,u);
                    q = iftGetVoxelIndex(imgc,v);
                    imgc->val[q] = img->val[p];
                }
        }
    } else { /* 3D image */
        imgc = iftCreateImage(xsize,ysize,zsize);
        for (u.z = uo.z, v.z=0; u.z <= uf.z; u.z++, v.z++)
            for (u.y = uo.y, v.y=0; u.y <= uf.y; u.y++,v.y++)
                for (u.x = uo.x,v.x=0; u.x <= uf.x; u.x++,v.x++){
                    p = iftGetVoxelIndex(img,u);
                    q = iftGetVoxelIndex(imgc,v);
                    imgc->val[q] = img->val[p];
                }
    }

    iftCopyVoxelSize(img,imgc);
    return(imgc);
}

iftImage *iftConvertColorSpace(iftImage* img, char colSpaceIn, char colSpaceOut) {
    int normValue = iftNormalizationValue(iftMaximumValue(img));

    /* validate the color spaces */
    bool isColorImgIn = iftIsColorImage(img);

    if((isColorImgIn && (colSpaceIn == GRAY_CSPACE || colSpaceIn == GRAYNorm_CSPACE)) ||
        (!isColorImgIn && (colSpaceIn != GRAY_CSPACE && colSpaceIn != GRAYNorm_CSPACE))) {
        iftError("The number of image channels (%d) and the specified input color space does not match",
            "iftConvertColorSpace", isColorImgIn ? 3 : 1);
    }

    if(!isColorImgIn && colSpaceOut != GRAY_CSPACE && colSpaceOut != GRAYNorm_CSPACE) {
        iftError("Color image representation can only be used with color images. The given image is grayscale.",
            "iftConvertColorSpace");
    }

    if(colSpaceIn == YCbCrNorm_CSPACE || colSpaceIn == RGBNorm_CSPACE || colSpaceIn == GRAYNorm_CSPACE ||
        colSpaceIn == LAB_CSPACE || colSpaceIn == LABNorm_CSPACE || colSpaceIn == LABNorm2_CSPACE) {
        iftError("The input color space is float and this method receives an integer image (iftImage). Use the method iftMConvertColorSpace() instead.",
            "iftConvertColorSpace");
    }

    if(colSpaceOut == YCbCrNorm_CSPACE || colSpaceOut == RGBNorm_CSPACE || colSpaceOut == GRAYNorm_CSPACE ||
        colSpaceOut == LAB_CSPACE || colSpaceOut == LABNorm_CSPACE || colSpaceOut == LABNorm2_CSPACE) {
        iftError("The chosen output color space is float and this method returns an integer image (iftImage). Use the method iftMConvertColorSpace() instead.",
            "iftConvertColorSpace");
    }

    if(colSpaceIn == colSpaceOut)
        return iftCopyImage(img);

    /* create a new image whose number of bands depends on the chosen color space */
    iftImage *imgOut = NULL;
    if(colSpaceOut == GRAY_CSPACE || colSpaceOut == GRAYNorm_CSPACE)
        imgOut = iftCreateImage(img->xsize, img->ysize, img->zsize);
    else
        imgOut = iftCreateColorImage(img->xsize, img->ysize, img->zsize, iftImageDepth(img));
    
    bool isColorImgOut = iftIsColorImage(imgOut);

    /* convert the color space of each pixel */
    for(int p = 0; p < img->n; p++) {
        iftFColor colorIn, colorOut;
        if(isColorImgIn) {
            colorIn.val[0] = img->val[p];
            colorIn.val[1] = img->Cb[p];
            colorIn.val[2] = img->Cr[p];
        }
        else {
            colorIn.val[0] = img->val[p];
        }

        colorOut = iftConvertPixelColorSpace(colorIn, colSpaceIn, colSpaceOut, normValue);

        if(isColorImgOut) {
            imgOut->val[p] = colorOut.val[0];
            imgOut->Cb[p] = colorOut.val[1];
            imgOut->Cr[p] = colorOut.val[2];
        }
        else {
            imgOut->val[p] = colorOut.val[0];
        }
    }

    return imgOut;
}

iftImage *iftLooseImage(iftImage *image, int xsize, int ysize, int zsize)
{
    int p, q;
    iftVoxel u, v, c, cResult;
    iftImage *result;

    c = iftGeometricCenterVoxel(image);
    result = iftCreateImage(xsize,ysize,zsize);
    cResult.x = iftRound((result->xsize - 1) / 2);
    cResult.y = iftRound((result->ysize - 1) / 2);
    cResult.z = iftRound((result->zsize - 1) / 2);

    for (u.z=0, v.z= cResult.z - c.z; u.z < image->zsize; u.z++,v.z++)
        for (u.y=0, v.y= cResult.y - c.y; u.y < image->ysize; u.y++,v.y++)
            for (u.x=0, v.x= cResult.x - c.x; u.x < image->xsize; u.x++,v.x++)
                if (iftValidVoxel(result,v))
                {
                    p = iftGetVoxelIndex(image,u);
                    q = iftGetVoxelIndex(result,v);
                    result->val[q] = image->val[p];
                }

    iftCopyVoxelSize(image, result);

    return(result);
}

void iftCenterImages(iftImage *image1, iftImage *image2, iftImage **centeredImage1, iftImage **centeredImage2)
{
    iftImage *image1Fit, *image2Fit;
    int maxX, maxY, maxZ, tolerance;

    tolerance = 5;
    image1Fit = iftExtractGreyObject(image1);
    image2Fit = iftExtractGreyObject(image2);
    maxX = iftMax(image1Fit->xsize, image2Fit->xsize) + tolerance;
    maxY = iftMax(image1Fit->ysize, image2Fit->ysize) + tolerance;
    maxZ = iftMax(image1Fit->zsize, image2Fit->zsize) + tolerance;

    *centeredImage1 = iftLooseImage(image1Fit, maxX, maxY, maxZ);
    *centeredImage2 = iftLooseImage(image2Fit, maxX, maxY, maxZ);

    iftDestroyImage(&image1Fit);
    iftDestroyImage(&image2Fit);
}



iftImage *iftColorTableToImage(iftColorTable *ct, int xsize, int ysize)
{
    int delta_y = (int)((float)ysize/ct->ncolors);
    iftImage *img = iftCreateImage(xsize,ysize,1);
    iftVoxel u;
    iftColor RGB, YCbCr;

    img->Cb = iftAllocUShortArray(img->n);
    img->Cr = iftAllocUShortArray(img->n);
    u.z     = 0;

    if (delta_y==0){ /* there are more colors than the image height */
        int deltaC = (int)((float)ct->ncolors/ysize);
        int i;
        for (u.y=0, i=0; u.y < img->ysize; u.y++, i += deltaC){
            RGB.val[0]  = ct->color[i].val[0];
            RGB.val[1]  = ct->color[i].val[1];
            RGB.val[2]  = ct->color[i].val[2];
            YCbCr       = iftRGBtoYCbCr(RGB,255);
            for (u.x=0; u.x < img->xsize; u.x++){
                int p           = iftGetVoxelIndex(img,u);
                img->val[p] = YCbCr.val[0];
                img->Cb[p]  = YCbCr.val[1];
                img->Cr[p]  = YCbCr.val[2];
            }
        }
    }else{ /* each color must be repeated by delta_y times along the vertical
        of the image */

        u.y = 0;
        for (int i=0; i < ct->ncolors; i++){

            RGB.val[0]  = ct->color[i].val[0];
            RGB.val[1]  = ct->color[i].val[1];
            RGB.val[2]  = ct->color[i].val[2];
            YCbCr       = iftRGBtoYCbCr(RGB,255);

            for (int y=0; (y <= delta_y)&&(u.y < img->ysize); y++, u.y++){
                for (u.x=0; u.x < img->xsize; u.x++){
                    int p           = iftGetVoxelIndex(img,u);
                    img->val[p] = YCbCr.val[0];
                    img->Cb[p]  = YCbCr.val[1];
                    img->Cr[p]  = YCbCr.val[2];
                }
            }
        }
    }

    return(img);
}



void iftTickColorTableImage(iftImage *img, float minval, float maxval, int nticks, const char *filename)
{
    FILE *fp=fopen(filename,"w");
    float value = minval, delta = (maxval-minval)/nticks, delta_y = (float)(img->ysize-1)/(float)nticks;
    int   xmax = (int)(0.1*img->xsize);
    iftVoxel u;
    int   i;

    u.z = 0;
    for (i=0, u.y = 0; (i <= nticks)&&(u.y < img->ysize); i++, value += delta, u.y = (int)(u.y + delta_y)) {
        fprintf(fp,"%f\n",value);
        for (u.x=0; u.x <= xmax; u.x++) {
            int p       = iftGetVoxelIndex(img,u);
            img->val[p] = 16;
            img->Cb[p]  = 128;
            img->Cr[p]  = 128;
        }
    }

    fclose(fp);

}


iftVoxelArray *iftGeometricCentersFromLabelImage(const iftImage *label_img) {
    int max_label = iftMaximumValue(label_img);

    iftVoxelArray *gcs = iftCreateVoxelArray(max_label + 1);
    iftIntArray *volumes = iftCreateIntArray(max_label + 1);

    for (int p = 0; p < label_img->n; p++) {
        int label = label_img->val[p];

        iftVoxel u = iftGetVoxelCoord(label_img, p);
        gcs->val[label].x += u.x;
        gcs->val[label].y += u.y;
        gcs->val[label].z += u.z;

        volumes->val[label]++;
    }

    #pragma omp parallel for
    for (int label = 1; label <= max_label; label++) {
        gcs->val[label].x /= volumes->val[label];
        gcs->val[label].y /= volumes->val[label];
        gcs->val[label].z /= volumes->val[label];

        if (iftImgVoxelVal(label_img, gcs->val[label]) != label) {
            int min_dist = IFT_INFINITY_INT;
            iftVoxel closest_voxel;
                
            for (int p = 0; p < label_img->n; p++) {
                if (label_img->val[p] == label) {
                    iftVoxel u = iftGetVoxelCoord(label_img, p);
                    int dist   = iftSquaredVoxelDistance(u, gcs->val[label]);

                    if (dist < min_dist) {
                        min_dist = dist;
                        closest_voxel = u;
                    }
                }
            }
            gcs->val[label] = closest_voxel;
        }
    }

    return gcs;
}



iftImage *iftReadImageByExt(const char *format, ...) {
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    if (!iftFileExists(filename))
        iftError("Image %s does not exist", "iftReadImageByExt", filename);


    iftImage *img = NULL;
    char *ext     = iftLowerString(iftFileExt(filename));

    if(iftCompareStrings(ext, ".png")) {
        img = iftReadImagePNG(filename);
    }
    else if (iftCompareStrings(ext, ".pgm")){
        FILE *fp = fopen(filename,"r");
        char type[10];
        if(fscanf(fp,"%s",type)!=1) iftError("Reading Error", "iftReadImageByExt");
        if (iftCompareStrings(type,"P5")){
            fclose(fp);
            img   = iftReadImageP5(filename);
        } else {
            fclose(fp);
            img   = iftReadImageP2(filename);
        }
    } else if (iftCompareStrings(ext, ".ppm")){
        img   = iftReadImageP6(filename);
    } else if (iftCompareStrings(ext, ".scn")){
        img   = iftReadImage(filename);
    } else if (iftCompareStrings(ext, ".zscn") || iftCompareStrings(ext, ".scn.gz")) {
        img   = iftReadImageGZip(filename);
    } else if (iftCompareStrings(ext, ".jpg") || iftCompareStrings(ext, ".jpeg")){
        img = iftReadImageJPEG(filename);
    } else if (iftCompareStrings(ext, ".tif") || iftCompareStrings(ext, ".tiff") || iftCompareStrings(ext, ".TIFF") ){
        iftError("ift do not support (yet) tiff format", "iftReadImageByExt");
    } else if (iftCompareStrings(ext, ".hdr")) {
        img = iftReadImageAnalyze(filename);
    } else if (iftCompareStrings(ext, ".nii") || iftCompareStrings(ext, ".nii.gz")) {
        img = iftReadImageNIfTI(filename);
//    } else if (iftCompareStrings(ext, ".npy")) {
//        img = iftReadNumPyAsImage(filename);
    }
    else {
        iftError("Invalid image format: \"%s\" - Try .scn, .zscn, .scn.gz, .pgm, .pgm, .png, .nii, .nii.gz, .hdr",
                 "iftReadImageByExt", ext);
    }

    iftFree(ext);
    return(img);
}

void iftWriteImageByExt(const iftImage *img, const char *format, ...) {
    if (img == NULL)
        iftWarning("Image is NULL... Nothing to write", "iftWriteImageByExt");
    else {
        char command[400];

        va_list args;
        char filename[300];

        va_start(args, format);
        vsprintf(filename, format, args);
        va_end(args);

        char *parent_dir = iftParentDir(filename);
        if (!iftDirExists(parent_dir))
            iftMakeDir(parent_dir);
        iftFree(parent_dir);

        char *ext = iftLowerString(iftFileExt(filename));

        if(iftCompareStrings(ext, ".png")) {
            iftWriteImagePNG(img,filename);
        } else if (iftCompareStrings(ext, ".scn")) {
            iftWriteImage(img, filename);
        } else if (iftCompareStrings(ext, ".scn.gz") || iftCompareStrings(ext, ".zscn")) {
            iftWriteImageGZip(img, filename);
        }else if (iftCompareStrings(ext, ".pgm")) {
            if (iftMaximumValue(img)>255)
                iftWriteImageP2(img,filename);
            else
                iftWriteImageP5(img,filename);
        } else if (iftCompareStrings(ext, ".ppm")){
            iftWriteImageP6(img,filename);
        } else if (iftIsColorImage(img)){
            iftWriteImageP6(img,"temp.ppm");
            sprintf(command,"convert temp.ppm %s",filename);
            if (system(command)==-1)
                iftError("Program convert failed or is not installed", "iftWriteImageByExt");
            if (system("rm -f temp.ppm")==-1)
                iftError("Cannot remore temp.ppm", "iftWriteImageByExt");
        } else if(iftCompareStrings(ext, ".jpg") || iftCompareStrings(ext, ".jpeg")) {
            iftWriteImageJPEG(img,filename);
        } else if(iftCompareStrings(ext, ".hdr")) {
            iftWriteImageAnalyze(img, filename);
        } else if(iftCompareStrings(ext, ".tif") || iftCompareStrings(ext, ".tiff") || iftCompareStrings(ext, ".TIFF")){
            iftError("ift do not support (yet) tiff format", "iftReadImageByExt");
        } else if (iftCompareStrings(ext, ".nii")) {
            iftWriteImageNIfTI(img, filename);
        } else if (iftCompareStrings(ext, ".nii.gz")) {
            iftWriteImageNIfTIGZip(img, filename);
//        } else if (iftCompareStrings(ext, ".npy")) {
//            iftWriteImageAsNumPy(img, filename);
        }
        else {
            printf("Invalid image format: %s. Please select among the accepted ones: .scn, .zscn, .scn.gz, .ppm, .pgm, .png, .nii, .hdr, .img\n",ext);
            exit(-1);
        }


        iftFree(ext);
    }
}


int iftNumberOfElements(const iftImage *mask)
{
    int p, nnodes = 0;

    for (p=0; p < mask->n; p++)
        if (mask->val[p]>0){
            nnodes++;
        }
    return(nnodes);
}


iftImage *iftSelectImageDomain(int xsize, int ysize, int zsize)
{
    iftImage *mask=iftCreateImage(xsize,ysize,zsize);

    iftSetImage(mask,1);
    return(mask);
}


iftImage *iftSelectRegionOfInterest(int xsize, int ysize, int zsize, iftVoxel uo, iftVoxel uf)
{
    iftImage *mask = iftCreateImage(xsize,ysize,zsize);
    int  p;
    iftVoxel u;

#pragma omp parallel for private(u,p)
    for (int x=uo.x; x <= uf.x; x++)
        for (int y=uo.y; y <= uf.y; y++)
            for (int z=uo.z; z <= uf.z; z++){
                u.x=x;u.y=y;u.z=z;
                p = iftGetVoxelIndex(mask,u);
                mask->val[p]=1;
            }

    return(mask);
}

iftImage *iftSwitchXByZ(iftImage *img)
{
    iftImage *nimg = iftCreateImage(img->zsize,img->ysize,img->xsize);
    iftVoxel  u, v;
    int       p, q;

    nimg->dx = img->dz;
    nimg->dy = img->dy;
    nimg->dz = img->dx;

    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++){
                v.x = u.z;
                v.y = img->ysize-1-u.y;
                v.z = u.x;
                p = iftGetVoxelIndex(img,u);
                q = iftGetVoxelIndex(nimg,v);
                nimg->val[q] = img->val[p];
            }

    return(nimg);
}

iftImage *iftSwitchYByZ(iftImage *img)
{
    iftImage *nimg = iftCreateImage(img->xsize,img->zsize,img->ysize);
    iftVoxel  u, v;
    int       p, q;

    nimg->dx = img->dx;
    nimg->dy = img->dz;
    nimg->dz = img->dy;

    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++){
                v.x = img->xsize-1-u.x;
                v.y = u.z;
                v.z = u.y;
                p = iftGetVoxelIndex(img,u);
                q = iftGetVoxelIndex(nimg,v);
                nimg->val[q] = img->val[p];
            }

    return(nimg);
}

iftImage *iftSwitchXByY(iftImage *img)
{
    iftImage *nimg = iftCreateImage(img->ysize,img->xsize,img->zsize);
    iftVoxel  u, v;
    int       p, q;

    nimg->dx = img->dy;
    nimg->dy = img->dx;
    nimg->dz = img->dz;

    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++){
                v.x = img->ysize-1-u.y;
                v.y = u.x;
                v.z = u.z;
                p = iftGetVoxelIndex(img,u);
                q = iftGetVoxelIndex(nimg,v);
                nimg->val[q] = img->val[p];
            }

    return(nimg);
}

iftImage *iftInvertX(iftImage *img)
{
    iftImage *nimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftVoxel  u, v;
    int       p, q;

    nimg->dx = img->dx;
    nimg->dy = img->dy;
    nimg->dz = img->dz;

    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++){
                v.x = img->xsize-1-u.x;
                v.y = u.y;
                v.z = u.z;
                p = iftGetVoxelIndex(img,u);
                q = iftGetVoxelIndex(nimg,v);
                nimg->val[q] = img->val[p];
            }

    return(nimg);
}

iftImage *iftInvertY(iftImage *img)
{
    iftImage *nimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftVoxel  u, v;
    int       p, q;

    nimg->dx = img->dx;
    nimg->dy = img->dy;
    nimg->dz = img->dz;

    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++){
                v.x = u.x;
                v.y = img->ysize-1-u.y;
                v.z = u.z;
                p = iftGetVoxelIndex(img,u);
                q = iftGetVoxelIndex(nimg,v);
                nimg->val[q] = img->val[p];
            }

    return(nimg);
}

iftImage *iftInvertZ(iftImage *img)
{
    iftImage *nimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftVoxel  u, v;
    int       p, q;

    nimg->dx = img->dx;
    nimg->dy = img->dy;
    nimg->dz = img->dz;

    for (u.z=0; u.z < img->zsize; u.z++)
        for (u.y=0; u.y < img->ysize; u.y++)
            for (u.x=0; u.x < img->xsize; u.x++){
                v.x = u.x;
                v.y = u.y;
                v.z = img->zsize-1-u.z;
                p = iftGetVoxelIndex(img,u);
                q = iftGetVoxelIndex(nimg,v);
                nimg->val[q] = img->val[p];
            }

    return(nimg);
}

iftVector iftObjectAxesVariance(iftImage *img) {
    iftVector result;

    unsigned long npoints = 0;
    iftPoint center = iftGeometricCenter(img);

    result.x = result.y = result.z = result.t = 0;
    for (int p = 0; p < img->n; ++p) {

        if(img->val[p]<=0)
            continue;

        iftVoxel v = iftGetVoxelCoord(img, p);
        result.x += (v.x - center.x)*(v.x - center.x);
        result.y += (v.y - center.y)*(v.y - center.y);
        result.z += (v.z - center.z)*(v.z - center.z);

        npoints++;
    }

    result = (iftVector)iftVectorScalarProd(result, 1.0/npoints);

    return result;
}


iftVoxelArray *iftFindClosestObjectVoxels(const iftImage *label_img, const iftVoxelArray *voxel_arr) {
    iftImage *root_img = NULL;
    iftImage *edt = iftEuclDistTrans(label_img, NULL, IFT_EXTERIOR, &root_img, NULL, NULL);
    
    iftVoxelArray *closest_obj_voxels = iftCreateVoxelArray(voxel_arr->n);

    #pragma omp parallel for
    for (int i = 0; i < voxel_arr->n; i++) {
        int root_idx = iftImgVoxelVal(root_img, voxel_arr->val[i]);
        closest_obj_voxels->val[i] = iftGetVoxelCoord(label_img, root_idx);
    }

    iftDestroyImage(&root_img);
    iftDestroyImage(&edt);

    return closest_obj_voxels;
}


void iftFillBoundingBoxInImage(iftImage *img, iftBoundingBox bb, int value) {
    if (bb.begin.x < 0)
        bb.begin.x = 0;
    if (bb.begin.y < 0)
        bb.begin.y = 0;
    if (bb.begin.z < 0)
        bb.begin.z = 0;
    if (bb.end.x >= img->xsize) {
        iftWarning("Bounding Box: end.x = %d out of Image Domain. Changing to %d\n", "iftFillBoundingBoxInImage",
                   bb.end.x, img->xsize-1);
        bb.end.x = img->xsize-1;
    }
    if (bb.end.y >= img->ysize) {
        iftWarning("Bounding Box: end.y = %d out of Image Domain. Changing to %d\n", "iftFillBoundingBoxInImage",
                   bb.end.y, img->ysize-1);
        bb.end.y = img->ysize-1;
    }
    if (bb.end.z >= img->zsize) {
        iftWarning("Bounding Box: end.z = %d out of Image Domain. Changing to %d\n", "iftFillBoundingBoxInImage",
                   bb.end.z, img->zsize-1);
        bb.end.z = img->zsize-1;
    }

    #pragma omp parallel for
    for (int z = bb.begin.z; z <= bb.end.z; z++)
        for (int y = bb.begin.y; y <= bb.end.y; y++)
            for (int x = bb.begin.x; x <= bb.end.x; x++)
                iftImgVal(img, x, y, z) = value;
}

void iftFillBoundingBoxInColorImage(iftImage *img, iftBoundingBox bb, iftColor YCbCr) {
    if (!iftIsColorImage(img))
        iftError("Input image is not a Color Image", "iftFillBoundingBoxInColorImage");

    if (bb.begin.x < 0)
        bb.begin.x = 0;
    if (bb.begin.y < 0)
        bb.begin.y = 0;
    if (bb.begin.z < 0)
        bb.begin.z = 0;
    if (bb.end.x >= img->xsize) {
        iftWarning("Bounding Box: end.x = %d out of Image Domain. Changing to %d\n", "iftFillBoundingBoxInColorImage",
                   bb.end.x, img->xsize-1);
        bb.end.x = img->xsize-1;
    }
    if (bb.end.y >= img->ysize) {
        iftWarning("Bounding Box: end.y = %d out of Image Domain. Changing to %d\n", "iftFillBoundingBoxInColorImage",
                   bb.end.y, img->ysize-1);
        bb.end.y = img->ysize-1;
    }
    if (bb.end.z >= img->zsize) {
        iftWarning("Bounding Box: end.z = %d out of Image Domain. Changing to %d\n", "iftFillBoundingBoxInColorImage",
                   bb.end.z, img->zsize-1);
        bb.end.z = img->zsize-1;
    }

    #pragma omp parallel for
    for (int z = bb.begin.z; z <= bb.end.z; z++)
        for (int y = bb.begin.y; y <= bb.end.y; y++)
            for (int x = bb.begin.x; x <= bb.end.x; x++) {
                iftImgVal(img, x, y, z) = YCbCr.val[0];
                iftImgCb(img, x, y, z) = YCbCr.val[1];
                iftImgCr(img, x, y, z) = YCbCr.val[2];
            }
}


iftBoundingBox iftMinBoundingBox(const iftImage *img, iftVoxel *gc_out) {
    if (img == NULL)
        iftError("Image is NULL", "iftMinBoundingBox");

    long n = 0; // number of spels non-background (non-zero)
    iftVoxel gc = {0.0, 0.0, 0.0};
    iftBoundingBox mbb;
    mbb.begin.x = mbb.begin.y = mbb.begin.z = IFT_INFINITY_INT;
    mbb.end.x = mbb.end.y = mbb.end.z = IFT_INFINITY_INT_NEG;

    for (long p = 0; p < img->n; p++) {
        if (img->val[p] != 0) {
            iftVoxel v = iftGetVoxelCoord(img, p);

            mbb.begin.x = iftMin(mbb.begin.x, v.x);
            mbb.begin.y = iftMin(mbb.begin.y, v.y);
            mbb.begin.z = iftMin(mbb.begin.z, v.z);

            mbb.end.x = iftMax(mbb.end.x, v.x);
            mbb.end.y = iftMax(mbb.end.y, v.y);
            mbb.end.z = iftMax(mbb.end.z, v.z);

            gc.x += v.x;
            gc.y += v.y;
            gc.z += v.z;
            n++;
        }
    }

    if (mbb.begin.x == IFT_INFINITY_INT) {
        mbb.begin.x = mbb.begin.y = mbb.begin.z = -1;
        mbb.end.x   = mbb.end.y   = mbb.end.z   = -1;
        gc.x        = gc.y        = gc.z        = -1.0;
    } else {
        gc.x /= n;
        gc.y /= n;
        gc.z /= n;
    }

    if (gc_out != NULL)
        *gc_out = gc;

    return mbb;
}


iftBoundingBox iftMinObjectBoundingBox(const iftImage *img, int obj_label, iftVoxel *gc_out) {
    if (img == NULL)
        iftError("Image is NULL", "iftMinObjectBoundingBox");

    long n = 0; // number of spels from the object <obj_label>
    iftVoxel gc = {0.0, 0.0, 0.0};
    iftBoundingBox mbb;
    mbb.begin.x = mbb.begin.y = mbb.begin.z = IFT_INFINITY_INT;
    mbb.end.x = mbb.end.y = mbb.end.z = IFT_INFINITY_INT_NEG;


    for (long p = 0; p < img->n; p++) {
        if (img->val[p] == obj_label) {
            iftVoxel v = iftGetVoxelCoord(img, p);

            mbb.begin.x = iftMin(mbb.begin.x, v.x);
            mbb.begin.y = iftMin(mbb.begin.y, v.y);
            mbb.begin.z = iftMin(mbb.begin.z, v.z);

            mbb.end.x = iftMax(mbb.end.x, v.x);
            mbb.end.y = iftMax(mbb.end.y, v.y);
            mbb.end.z = iftMax(mbb.end.z, v.z);

            gc.x += v.x;
            gc.y += v.y;
            gc.z += v.z;
            n++;
        }
    }
    if (mbb.begin.x == IFT_INFINITY_INT) {
        mbb.begin.x = mbb.begin.y = mbb.begin.z = -1;
        mbb.end.x   = mbb.end.y   = mbb.end.z   = -1;
        gc.x        = gc.y        = gc.z        = -1.0;
    } else {
        gc.x /= n;
        gc.y /= n;
        gc.z /= n;
    }

    if (gc_out != NULL)
        *gc_out = gc;

    return mbb;
}


iftBoundingBox *iftMinLabelsBoundingBox(const iftImage *img, const iftIntArray *labels, iftVoxelArray **gcs_out) {
    int n_objs = labels->n;

    int *n_voxels        = iftAllocIntArray(n_objs); // used to compute the geo center of each bounding box
    iftVoxelArray *gcs = iftCreateVoxelArray(n_objs);
    iftBoundingBox *mbbs = (iftBoundingBox*) iftAlloc(n_objs+1, sizeof(iftBoundingBox));

    iftDict *labels_dict = iftCreateDict(); // indexes the labels into a dict

    for (int o = 0; o < n_objs; o++) {
        iftInsertIntoDict(labels->val[o], o, labels_dict);

        gcs->val[o].x   = gcs->val[o].y   = gcs->val[o].z   = 0.0;
        mbbs[o].begin.x = mbbs[o].begin.y = mbbs[o].begin.z = IFT_INFINITY_INT;
        mbbs[o].end.x   = mbbs[o].end.y   = mbbs[o].end.z   = IFT_INFINITY_INT_NEG;
        n_voxels[o]     = 0;
    }

    // finds the min bounding boxes
    for (long p = 0; p < img->n; p++) {
        int label = img->val[p];

        // comparisons extremely faster than iftDictContainKey
        if ((label != 0) && (iftIntArrayContainValue(labels->val, labels->n, label))) {
            int o      = iftGetLongValFromDict(label, labels_dict); // object index in the label array
            iftVoxel v = iftGetVoxelCoord(img, p);

            mbbs[o].begin.x = iftMin(mbbs[o].begin.x, v.x);
            mbbs[o].begin.y = iftMin(mbbs[o].begin.y, v.y);
            mbbs[o].begin.z = iftMin(mbbs[o].begin.z, v.z);

            mbbs[o].end.x = iftMax(mbbs[o].end.x, v.x);
            mbbs[o].end.y = iftMax(mbbs[o].end.y, v.y);
            mbbs[o].end.z = iftMax(mbbs[o].end.z, v.z);

            gcs->val[o].x += v.x;
            gcs->val[o].y += v.y;
            gcs->val[o].z += v.z;
            n_voxels[o]++;
        }
    }
    
    for (int o = 0; o < n_objs; o++) {
        // if there is not the object o
        if (mbbs[o].begin.x == IFT_INFINITY_INT) {
            mbbs[o].begin.x = mbbs[o].begin.y = mbbs[o].begin.z = -1;
            mbbs[o].end.x   = mbbs[o].end.y   = mbbs[o].end.z   = -1;
            gcs->val[o].x   = gcs->val[o].y   = gcs->val[o].z   = -1;
        }
        else {
            gcs->val[o].x = iftRound(gcs->val[o].x / n_voxels[o]);
            gcs->val[o].y = iftRound(gcs->val[o].y / n_voxels[o]);
            gcs->val[o].z = iftRound(gcs->val[o].z / n_voxels[o]);
        }
    }

    // computes the geometric center
    if (gcs_out == NULL) { iftFree(gcs); }
    else {
        *gcs_out = iftFindClosestObjectVoxels(img, gcs);
        iftDestroyVoxelArray(&gcs);
    }

    iftDestroyDict(&labels_dict);
    iftFree(n_voxels);

    return mbbs;
}


iftBoundingBox *iftAllMinBoundingBox(const iftFileSet *img_paths, iftVoxel **gcs_out) {
    iftBoundingBox *mbbs = (iftBoundingBox*) iftAlloc(img_paths->n, sizeof(iftBoundingBox));
    iftVoxel *gcs         = (iftVoxel*) iftAlloc(img_paths->n, sizeof(iftVoxel));

    #pragma omp parallel for
    for (size_t i = 0; i < img_paths->n; i++) {
        iftImage *img = iftReadImageByExt(img_paths->files[i]->path);
        mbbs[i]       = iftMinBoundingBox(img, &gcs[i]);
        iftDestroyImage(&img);
    }

    if (gcs_out == NULL)
        iftFree(gcs);
    else
        *gcs_out = gcs;

    return mbbs;
}


iftBoundingBox *iftAllMinObjectBoundingBox(const iftFileSet *img_paths, int obj_label, iftVoxel **gcs_out) {
    iftBoundingBox *mbbs  = (iftBoundingBox*) iftAlloc(img_paths->n, sizeof(iftBoundingBox));
    iftVoxel *gcs         = (iftVoxel*) iftAlloc(img_paths->n, sizeof(iftVoxel));

    for (size_t i = 0; i < img_paths->n; i++) {
        iftImage *img = iftReadImageByExt(img_paths->files[i]->path);
        mbbs[i]       = iftMinObjectBoundingBox(img, obj_label, &gcs[i]);
        iftDestroyImage(&img);
    }

    if (gcs_out == NULL)
        iftFree(gcs);
    else
        *gcs_out = gcs;

    return mbbs;
}


iftBoundingBox **iftAllMinLabelsBoundingBox(const iftFileSet *label_paths, const iftIntArray *labels,
                                            iftVoxel ***gcs_out) {
    iftBoundingBox **mbbs = (iftBoundingBox**) iftAlloc(labels->n, sizeof(iftBoundingBox*));
    iftVoxel **gcs        = (iftVoxel**) iftAlloc(labels->n, sizeof(iftVoxel*));
    for (int o = 0; o < labels->n; o++) {
        mbbs[o] = (iftBoundingBox*) iftAlloc(label_paths->n, sizeof(iftBoundingBox));
        gcs[o]  = (iftVoxel*) iftAlloc(label_paths->n, sizeof(iftVoxel));
    }

    for (size_t i = 0; i < label_paths->n; i++) {
        iftImage *label_img = iftReadImageByExt(label_paths->files[i]->path);

        iftVoxelArray *gcs_arr_in_img = NULL;
        iftBoundingBox *mbbs_img = iftMinLabelsBoundingBox(label_img, labels, &gcs_arr_in_img);

        for (int o = 0; o < labels->n; o++) {
            mbbs[o][i] = mbbs_img[o];
            gcs[o][i] = gcs_arr_in_img->val[o];
        }

        iftDestroyImage(&label_img);
        iftDestroyVoxelArray(&gcs_arr_in_img);
    }

    if (gcs_out == NULL) {
        for (int o = 0; o < labels->n; o++)
            iftFree(gcs[o]);
        iftFree(gcs);
    }
    else *gcs_out = gcs;

    return mbbs;
}


iftImage *iftExtractROI(const iftImage *img, iftBoundingBox bb) {
    if (img == NULL)
        iftError("Image is NULL", "iftExtractROI");

    iftVoxel uo = bb.begin;
    iftVoxel uf = bb.end;

    if (!iftValidVoxel(img, uo))
        iftError(
                "Initial Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for image with size (%d, %d, %d)",
                "iftExtractROI", uo.x, uo.y, uo.z, img->xsize, img->ysize, img->zsize);
    if (!iftValidVoxel(img, uf)) {
        iftWarning("The ROI Image DOES NOT fit entirely inside the Image.\n" \
                   "It will copied/inserted what it is possible", "iftExtractROI");
        // gets the valid ending voxel
        uf.x = iftMin(uf.x, img->xsize - 1);
        uf.y = iftMin(uf.y, img->ysize - 1);
        uf.z = iftMin(uf.z, img->zsize - 1);
    }


    iftImage *roi = iftCreateImage(uf.x - uo.x + 1, uf.y - uo.y + 1, uf.z - uo.z + 1);
    iftCopyVoxelSize(img, roi);

    if (iftIsColorImage(img)) {
        iftSetCbCr(roi, (iftMaxImageRange(iftImageDepth(img))+1)/2);
        int q = 0;
        iftVoxel u;

        for (u.z = uo.z; u.z <= uf.z; u.z++)
            for (u.y = uo.y; u.y <= uf.y; u.y++)
                for (u.x = uo.x; u.x <= uf.x; u.x++) {
                    int p = iftGetVoxelIndex(img, u);
                    roi->val[q] = img->val[p];
                    roi->Cb[q]  = img->Cb[p];
                    roi->Cr[q]  = img->Cr[p];
                    q++;
                }
    }
    else {
        int q = 0;
        iftVoxel u;

        for (u.z = uo.z; u.z <= uf.z; u.z++)
            for (u.y = uo.y; u.y <= uf.y; u.y++)
                for (u.x = uo.x; u.x <= uf.x; u.x++) {
                    int p = iftGetVoxelIndex(img, u);
                    roi->val[q] = img->val[p];
                    q++;
                }
    }


    return roi;
}


iftImage *iftExtractObjectInsideROI(const iftImage *src, iftBoundingBox bb, int obj_label) {
    if (src == NULL)
        iftError("Source Image is NULL", "iftExtractObjectInsideROI");

    iftImage *roi = iftExtractROI(src, bb);

    #pragma omp parallel for
    for (int p = 0; p < roi->n; p++) {
        if (roi->val[p] != obj_label)
            roi->val[p] = 0;
    }

    return roi;
}


iftImage *iftExtractLabelsInsideROI(const iftImage *src_img, iftBoundingBox bb, const iftIntArray *labels) {
    iftImage *roi = iftExtractROI(src_img, bb);

    for (int p = 0; p < roi->n; p++) {
        int label = 0;

        for (int o = 0; o < labels->n; o++) {
            // if the voxel belongs to a required object, it label is kept
            if (roi->val[p] == labels->val[o]) {
                label = roi->val[p];
                break;
            }
        }
        roi->val[p] = label;

    }

    return roi;
}


void iftInsertROI(const iftImage *roi, iftImage *target, iftVoxel begin) {
    if (roi == NULL)
        iftError("ROI Image is NULL", "iftInsertROI");
    if (target == NULL)
        iftError("Source Image is NULL", "iftInsertROI");

    iftVoxel uo = begin;
    iftVoxel uf;
    uf.x = uo.x + roi->xsize - 1;
    uf.y = uo.y + roi->ysize - 1;
    uf.z = uo.z + roi->zsize - 1;

    if (!iftValidVoxel(target, uo))
        iftError(
                "Initial Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for image with size (%d, %d, %d)",
                "iftInsertROI", uo.x, uo.y, uo.z, target->xsize, target->ysize, target->zsize);
    if (!iftValidVoxel(target, uf)) {
        iftWarning("The ROI Image DOES NOT fit entirely inside the Target Image.\n" \
                   "It will copied/inserted what it is possible", "iftInsertROI");
        // gets the valid ending voxel
        uf.x = iftMin(uf.x, target->xsize - 1);
        uf.y = iftMin(uf.y, target->ysize - 1);
        uf.z = iftMin(uf.z, target->zsize - 1);
    }


    if (iftIsColorImage(target)) {
        iftVoxel u;

        for (u.z = uo.z; u.z <= uf.z; u.z++)
            for (u.y = uo.y; u.y <= uf.y; u.y++)
                for (u.x = uo.x; u.x <= uf.x; u.x++) {
                    iftVoxel v;
                    v.x = u.x - uo.x;
                    v.y = u.y - uo.y;
                    v.z = u.z - uo.z;

                    int p = iftGetVoxelIndex(target, u);
                    int q = iftGetVoxelIndex(roi, v);

                    target->val[p] = roi->val[q];
                    target->Cb[p]  = roi->Cb[q];
                    target->Cr[p]  = roi->Cr[q];
                }
    } else {
        iftVoxel u;

        for (u.z = uo.z; u.z <= uf.z; u.z++)
            for (u.y = uo.y; u.y <= uf.y; u.y++)
                for (u.x = uo.x; u.x <= uf.x; u.x++) {
                    iftVoxel v;
                    v.x = u.x - uo.x;
                    v.y = u.y - uo.y;
                    v.z = u.z - uo.z;

                    int p       = iftGetVoxelIndex(target, u);
                    int q       = iftGetVoxelIndex(roi, v);

                    target->val[p] = roi->val[q];
                }
    }
}


void iftInsertROIByCenter(const iftImage *roi, iftImage *target) {
    if (roi == NULL)
        iftError("ROI Image is NULL", "iftInsertROIByCenter");
    if (target == NULL)
        iftError("Source Image is NULL", "iftInsertROIByCenter");

    iftVoxel roi_center;
    roi_center.x = roi->xsize / 2;
    roi_center.y = roi->ysize / 2;
    roi_center.z = roi->zsize / 2;

    iftVoxel target_center;
    target_center.x = target->xsize / 2;
    target_center.y = target->ysize / 2;
    target_center.z = target->zsize / 2;

    iftInsertROIByPosition(roi, roi_center, target, target_center);
/*
//    iftVoxel uo; // initial point = displacement between the centers
//    uo.x = target_center.x - roi_center.x;
//    uo.y = target_center.y - roi_center.y;
//    uo.z = target_center.z - roi_center.z;
//
//    if (!iftValidVoxel(target, uo)) {
//        char msg[512];
//        sprintf(msg, "Initial Position to Insert Roi doesn't is valid: (%d, %d, %d)\n", uo.x, uo.y, uo.z);
//
//        // since the image centers are aligned, the initial point of Roi in the Target Image will never be
//        // larger than the center point in Target Image... Then, if it is invalid is because
//        // some coordinate is less than 0
//        uo.x = iftMax(uo.x, 0);
//        uo.y = iftMax(uo.y, 0);
//        uo.z = iftMax(uo.z, 0);
//
//        sprintf(msg, "New Initial Point: (%d, %d, %d)", uo.x, uo.y, uo.z);
//        iftWarning(msg, "iftInsertROIByCenter");
//    }
//
//    iftVoxel uf;
//    uf.x = uo.x + roi->xsize - 1;
//    uf.y = uo.y + roi->ysize - 1;
//    uf.z = uo.z + roi->zsize - 1;
//
//
//    if (!iftValidVoxel(target, uf)) {
//        iftWarning("The ROI Image DOES NOT fit entirely inside the Target Image.\n" \
//                   "It will copied/inserted what it is possible", "iftInsertROIByCenter");
//        // gets the valid ending voxel
//        uf.x = iftMin(uf.x, target->xsize - 1);
//        uf.y = iftMin(uf.y, target->ysize - 1);
//        uf.z = iftMin(uf.z, target->zsize - 1);
//    }
//
//
//    if (iftIsColorImage(target)) {
//        iftVoxel u;
//
//        for (u.z = uo.z; u.z <= uf.z; u.z++)
//            for (u.y = uo.y; u.y <= uf.y; u.y++)
//                for (u.x = uo.x; u.x <= uf.x; u.x++) {
//                    iftVoxel v;
//                    v.x = u.x - uo.x;
//                    v.y = u.y - uo.y;
//                    v.z = u.z - uo.z;
//
//                    int p = iftGetVoxelIndex(target, u);
//                    int q = iftGetVoxelIndex(roi, v);
//
//                    target->val[p] = roi->val[q];
//                    target->Cb[p]  = roi->Cb[q];
//                    target->Cr[p]  = roi->Cr[q];
//                }
//    } else {
//        iftVoxel u;
//
//        for (u.z = uo.z; u.z <= uf.z; u.z++)
//            for (u.y = uo.y; u.y <= uf.y; u.y++)
//                for (u.x = uo.x; u.x <= uf.x; u.x++) {
//                    iftVoxel v;
//                    v.x = u.x - uo.x;
//                    v.y = u.y - uo.y;
//                    v.z = u.z - uo.z;
//
//                    int p       = iftGetVoxelIndex(target, u);
//                    int q       = iftGetVoxelIndex(roi, v);
//
//                    target->val[p] = roi->val[q];
//                }
//    }
*/
}


void iftInsertROIByPosition(const iftImage *roi, const iftVoxel roi_pos, iftImage *target, const iftVoxel target_pos) {
    iftVoxel uo, uf;

    if (roi == NULL)
        iftError("ROI Image is NULL", "iftInsertROIByPosition");
    if (target == NULL)
        iftError("Source Image is NULL", "iftInsertROIByPosition");

    // Computing the first valid position in the target image that intersects the ROI
    uo.x = iftMax(0, target_pos.x - roi_pos.x);
    uo.y = iftMax(0, target_pos.y - roi_pos.y);
    uo.z = iftMax(0, target_pos.z - roi_pos.z);

    // Computing the last valid position in the target image that intersects the ROI
    uf.x = iftMin(target->xsize - 1, uo.x + roi->xsize - 1);
    uf.y = iftMin(target->ysize - 1, uo.y + roi->ysize - 1);
    uf.z = iftMin(target->zsize - 1, uo.z + roi->zsize - 1);

    if (iftIsColorImage(target)) {
        iftVoxel u;

        // Iterating over the target image and copying the values from the ROI
        for (u.z = uo.z; u.z <= uf.z; u.z++)
            for (u.y = uo.y; u.y <= uf.y; u.y++)
                for (u.x = uo.x; u.x <= uf.x; u.x++) {
                    iftVoxel v;

                    // NOTE: this is almost the same as subtracting from uo
                    v.x = (u.x - target_pos.x) + roi_pos.x;
                    v.y = (u.y - target_pos.y) + roi_pos.y;
                    v.z = (u.z - target_pos.z) + roi_pos.z;

                    int p = iftGetVoxelIndex(target, u);
                    int q = iftGetVoxelIndex(roi, v);

                    target->val[p] = roi->val[q];
                    target->Cb[p]  = roi->Cb[q];
                    target->Cr[p]  = roi->Cr[q];
                }
    } else {
        iftVoxel u;

        // Iterating over the target image and copying the values from the ROI
        for (u.z = uo.z; u.z <= uf.z; u.z++)
            for (u.y = uo.y; u.y <= uf.y; u.y++)
                for (u.x = uo.x; u.x <= uf.x; u.x++) {
                    iftVoxel v;

                    v.x = (u.x - target_pos.x) + roi_pos.x;
                    v.y = (u.y - target_pos.y) + roi_pos.y;
                    v.z = (u.z - target_pos.z) + roi_pos.z;

                    int p       = iftGetVoxelIndex(target, u);
                    int q       = iftGetVoxelIndex(roi, v);

                    target->val[p] = roi->val[q];
                }
    }
}


iftImage *iftCopyImageROI(const iftImage *src, iftBoundingBox bb) {
    if (src == NULL)
        iftError("Source Image is NULL", "iftCopyImageROI");

    iftImage *copy = iftCreateImage(src->xsize, src->ysize, src->zsize);
    iftCopyVoxelSize(src, copy);

    iftImage *roi  = iftExtractROI(src, bb);
    iftInsertROI(roi, copy, bb.begin);

    iftDestroyImage(&roi);

    return copy;
}


iftImage *iftExtractObject(const iftImage *src, int obj_label) {
    if (src == NULL)
        iftError("Source Image is NULL", "iftExtractObject");

    iftBoundingBox bb;
    bb.begin = iftGetVoxelCoord(src, 0);
    bb.end   = iftGetVoxelCoord(src, src->n - 1);

    return iftExtractObjectInsideROI(src, bb, obj_label);
}


iftImage *iftExtractLabels(const iftImage *src_img, const iftIntArray *labels) {
    iftBoundingBox bb;
    bb.begin = iftGetVoxelCoord(src_img, 0);
    bb.end   = iftGetVoxelCoord(src_img, src_img->n - 1);

    return iftExtractLabelsInsideROI(src_img, bb, labels);
}


int iftCountObjectSpels(const iftImage *label, int obj_label) {
    if (label == NULL)
        iftError("Label Image is NULL", "iftCountObjectSpels");

    int n_obj_spels = 0;

    #pragma omp parallel for reduction(+:n_obj_spels)
    for (int p = 0; p < label->n; p++) {
        if (label->val[p] == obj_label)
            n_obj_spels++;
    }

    return n_obj_spels;
}

int iftCountBackgroundSpels(const iftImage *label_img) {
    return iftCountObjectSpels(label_img, 0);
}


int iftCountTotalObjectSpels(const iftImage *label_img) {
    return (label_img->n - iftCountBackgroundSpels(label_img));
}


int iftCountObjectSpelsFromBoundingBox(const iftImage *label_img, int obj_label, iftBoundingBox bb) {
    if (label_img == NULL)
        iftError("Label Image is NULL", "iftCountObjectVoxels");
    if (!iftValidVoxel(label_img, bb.begin))
        iftError("Begin Voxel (%d, %d, %d) from Bounding Box does not belong to Label Image Domain",
                 "iftCountObjectVoxelsFromBoundingBox", bb.begin.x, bb.begin.y, bb.begin.z);
    if (!iftValidVoxel(label_img, bb.end))
        iftError("Ending Voxel (%d, %d, %d) from Bounding Box does not belong to Label Image Domain",
                 "iftCountObjectVoxelsFromBoundingBox", bb.end.x, bb.end.y, bb.end.z);

    int n_obj_spels = 0;

    iftVoxel u;
    for (u.z = bb.begin.z; u.z <= bb.end.z; u.z++)
        for (u.y = bb.begin.y; u.y <= bb.end.y; u.y++)
            for (u.x = bb.begin.x; u.x <= bb.end.x; u.x++) {
                int p = iftGetVoxelIndex(label_img, u);
                n_obj_spels += (label_img->val[p] == obj_label);
            }

    return n_obj_spels;
}


iftIntArray *iftCountLabelSpels(const iftImage *label) {
    if (label == NULL)
        iftError("Label Image is NULL", "iftCountObjectVoxels");

    int n_labels               = iftMaximumValue(label);
    iftIntArray *n_label_spels = iftCreateIntArray(n_labels + 1); // indices in [0, n_labels]

    for (int p = 0; p < label->n; p++) {
        n_label_spels->val[label->val[p]]++;
    }

    return n_label_spels;
}


float iftAreaVolumeOfObject(const iftImage *label, int obj_label) {
    if (label == NULL)
        iftError("Label Image is NULL", "iftAreaVolumeOfObject");

    int n_labels = iftMaximumValue(label);
    if ((obj_label < 0) || (obj_label > n_labels))
        iftError("Invalid Object Label: %d\nTry [0..%d]", "iftCountObjectSpels", obj_label, n_labels);

    float dv;
    if (iftIs3DImage(label))
        dv = label->dx * label->dy * label->dz; // volume (in mm³) from a single spel (pixel/voxel)
    else
        dv = label->dx * label->dy; // area (in mm²) from a single spel (pixel/voxel)

    int n_spels  = iftCountObjectSpels(label, obj_label);
    float volume = n_spels * dv;

    return volume;
}


double iftAreaVolumeOfObjectFromBoundingBox(const iftImage *label_img, int obj_label, iftBoundingBox bb) {
    if (label_img == NULL)
        iftError("Label Image is NULL", "iftAreaVolumeOfObject");

    double dv;
    if (iftIs3DImage(label_img))
        dv = label_img->dx * label_img->dy * label_img->dz; // volume (in mm³) from a single spel (pixel/voxel)
    else
        dv = label_img->dx * label_img->dy; // area (in mm²) from a single spel (pixel/voxel)

    int n_spels   = iftCountObjectSpelsFromBoundingBox(label_img, obj_label, bb);
    double volume = n_spels * dv;

    return volume;
}


iftFloatArray *iftAreaVolumeOfLabels(const iftImage *label) {
    if (label == NULL)
        iftError("Label Image is NULL", "iftAreaVolumeOfLabels");

    float dv;
    if (iftIs3DImage(label))
        dv = label->dx * label->dy * label->dz; // volume (in mm³) from a single spel (pixel/voxel)
    else
        dv = label->dx * label->dy; // area (in mm²) from a single spel (pixel/voxel)

    iftIntArray *n_label_spels = iftCountLabelSpels(label); // [0] to [n_labels]
    iftFloatArray *vols        = iftCreateFloatArray(n_label_spels->n); // [0] to [n_labels]
    int n_labels               = n_label_spels->n-1; // exclude the position 0 (background)

    for (int o = 1; o <= n_labels; o++)
        vols->val[o] = n_label_spels->val[o] * dv;

    iftDestroyIntArray(&n_label_spels);

    return vols;
}


iftImageTiles *iftCreateImageTiles(iftBoundingBox *tile_coords, int ntiles_x, int ntiles_y, int ntiles_z,
                                   iftBoundingBox bb) {
    iftImageTiles *tiles = NULL;

    tiles = (iftImageTiles*)iftAlloc(1, sizeof(iftImageTiles));

    tiles->ntiles_x    = ntiles_x;
    tiles->ntiles_y    = ntiles_y;
    tiles->ntiles_z    = ntiles_z;
    tiles->ntiles      = ntiles_x*ntiles_y*ntiles_z;
    tiles->bb          = bb;
    tiles->tile_coords = tile_coords;

    return tiles;
}

void iftDestroyImageTiles(iftImageTiles **tiles)
{
    if(tiles != NULL && *tiles != NULL)
    {
        if((*tiles)->tile_coords != NULL) iftFree((*tiles)->tile_coords);
        iftFree(*tiles);

        *tiles = NULL;
    }
}


iftImageTiles* iftReadImageTiles(const char *filename)
{
    int ntiles_x, ntiles_y, ntiles_z, n, ndim = 6, i;
    iftBoundingBox bb;
    iftBoundingBox *tile_coords = NULL;
    iftDict *json = NULL;
    iftIntArray *read_tile_coords = NULL;

    if(!iftFileExists(filename)) {
        return NULL;
    }

    json = iftReadJson(filename);

    ntiles_x = iftGetLongValFromDict("ntiles_x", json);
    ntiles_y = iftGetLongValFromDict("ntiles_y", json);
    ntiles_z = iftGetLongValFromDict("ntiles_z", json);

    bb.begin.x = iftGetLongValFromDict("bb_begin_x", json);
    bb.begin.y = iftGetLongValFromDict("bb_begin_y", json);
    bb.begin.z = iftGetLongValFromDict("bb_begin_z", json);

    bb.end.x = iftGetLongValFromDict("bb_end_x", json);
    bb.end.y = iftGetLongValFromDict("bb_end_y", json);
    bb.end.z = iftGetLongValFromDict("bb_end_z", json);

    read_tile_coords = iftGetIntArrayFromDict("tile_coords", json);

    // Loading coordinates
    n = ntiles_x*ntiles_y*ntiles_z;

    tile_coords = (iftBoundingBox*)iftAlloc(n, sizeof(iftBoundingBox));

    for(i = 0; i < n; i++) {
        tile_coords[i].begin.x = read_tile_coords->val[ndim * i];
        tile_coords[i].begin.y = read_tile_coords->val[ndim * i + 1];
        tile_coords[i].begin.z = read_tile_coords->val[ndim * i + 2];
        tile_coords[i].end.x = read_tile_coords->val[ndim * i + 3];
        tile_coords[i].end.y = read_tile_coords->val[ndim * i + 4];
        tile_coords[i].end.z = read_tile_coords->val[ndim * i + 5];
    }

    iftDestroyDict(&json);
    iftDestroyIntArray(&read_tile_coords);

    return iftCreateImageTiles(tile_coords, ntiles_x, ntiles_y, ntiles_z, bb);
}

void iftWriteImageTiles(iftImageTiles *tiles, const char *filename)
{
    int ndim = 6, i;
    iftDict *json = NULL;
    iftIntArray *tile_coords = NULL;

    json = iftCreateDict();

    iftInsertIntoDict("ntiles_x", tiles->ntiles_x, json);
    iftInsertIntoDict("ntiles_y", tiles->ntiles_y, json);
    iftInsertIntoDict("ntiles_z", tiles->ntiles_z, json);

    iftInsertIntoDict("bb_begin_x", tiles->bb.begin.x, json);
    iftInsertIntoDict("bb_begin_y", tiles->bb.begin.y, json);
    iftInsertIntoDict("bb_begin_z", tiles->bb.begin.z, json);

    iftInsertIntoDict("bb_end_x", tiles->bb.end.x, json);
    iftInsertIntoDict("bb_end_y", tiles->bb.end.y, json);
    iftInsertIntoDict("bb_end_z", tiles->bb.end.z, json);

    tile_coords = iftCreateIntArray(tiles->ntiles*ndim);
    for(i = 0; i < tiles->ntiles; i++)
    {
        tile_coords->val[ndim*i]     = tiles->tile_coords[i].begin.x;
        tile_coords->val[ndim*i + 1] = tiles->tile_coords[i].begin.y;
        tile_coords->val[ndim*i + 2] = tiles->tile_coords[i].begin.z;
        tile_coords->val[ndim*i + 3] = tiles->tile_coords[i].end.x;
        tile_coords->val[ndim*i + 4] = tiles->tile_coords[i].end.y;
        tile_coords->val[ndim*i + 5] = tiles->tile_coords[i].end.z;
    }

    iftInsertIntoDict("tile_coords", tile_coords, json);
    iftWriteJson(json, filename);

    iftDestroyDict(&json);
    iftDestroyIntArray(&tile_coords);
}

iftImageTiles *iftCopyImageTiles(const iftImageTiles *tiles) {
    if (tiles == NULL)
        iftError("Image Tiles is NULL", "iftCopyImageTiles");

    iftBoundingBox *tile_coords = (iftBoundingBox*) iftAlloc(tiles->ntiles, sizeof(iftBoundingBox));
    for (int t = 0; t < tiles->ntiles; t++) {
        tile_coords[t].begin = tiles->tile_coords[t].begin;
        tile_coords[t].end   = tiles->tile_coords[t].end;
    }

    return iftCreateImageTiles(tile_coords, tiles->ntiles_x, tiles->ntiles_y, tiles->ntiles_z, tiles->bb);
}


void iftNumberOfEquallyDimensionedTilesWithGivenMaximumSize(int xsize, int ysize, int zsize, unsigned long tile_number_of_voxels,
                                                            int *ntiles_x, int *ntiles_y, int *ntiles_z) {
    unsigned long n;
    double f;
    double xsize_tile, ysize_tile, zsize_tile;

    n = xsize*ysize*zsize;

    if(zsize > 1) {
        // computing the downscaling factor as considering the cubic root of the
        // fraction between tile_number_of_voxels and n (i.e., the volume diagonal). Then
        // the scaling factor is used to reduce the sides of the volume
        f = pow(tile_number_of_voxels / (double)n, 1.0 / 3.0);
        xsize_tile = xsize*f;
        ysize_tile = ysize*f;
        zsize_tile = zsize*f;
    } else {
        // Same reasoning as for the 3D case, this time considering the diagonal of
        // the rectangle
        f = sqrt(tile_number_of_voxels / (double)n);
        xsize_tile = xsize*f;
        ysize_tile = ysize*f;
        zsize_tile = 1;
    }

    // Determining the number of vertical and horizontal tiles
    *ntiles_x = iftMax(iftRound(ceil(xsize / xsize_tile)), 1);
    *ntiles_y = iftMax(iftRound(ceil(ysize / ysize_tile)), 1);
    *ntiles_z = iftMax(iftRound(ceil(zsize / zsize_tile)), 1);
}

iftImageTiles *iftImageTilesByEquallyDividingAxes(const iftImage *img, int ntiles_x, int ntiles_y, int ntiles_z) {
    iftBoundingBox bb;

    bb.begin.x = bb.begin.y = bb.begin.z = 0;
    bb.end.x = img->xsize - 1;
    bb.end.y = img->ysize - 1;
    bb.end.z = img->zsize - 1;

    return iftBoundingBoxImageTilesByEquallyDividingAxes(bb, ntiles_x, ntiles_y, ntiles_z);
}

iftImageTiles *iftBoundingBoxImageTilesByEquallyDividingAxes(iftBoundingBox bb, int ntiles_x, int ntiles_y, int
ntiles_z)
{
    int z, y, x, z_tile, y_tile, x_tile, i, n;
    iftBoundingBox *tile_coords = NULL;
    iftImageTiles *tiles = NULL;
    int xsize, ysize, zsize;
    iftVoxel start, end;

    if(ntiles_x <= 0 || ntiles_y <= 0 || ntiles_z <= 0)
        iftError("The number of tiles must be positive and not (%d, %d, %d) !",
                 "iftBoundingBoxImageTilesByEquallyDividingAxes",
                 ntiles_x, ntiles_y, ntiles_z);

    xsize = bb.end.x - bb.begin.x + 1;
    ysize = bb.end.y - bb.begin.y + 1;
    zsize = bb.end.z - bb.begin.z + 1;

    (ntiles_x) = ((ntiles_x) >= xsize) ? xsize : (ntiles_x);
    (ntiles_y) = ((ntiles_y) >= ysize) ? ysize : (ntiles_y);
    (ntiles_z) = ((ntiles_z) >= zsize) ? zsize : (ntiles_z);

    x_tile = xsize/ntiles_x;
    y_tile = ysize/ntiles_y;
    z_tile = zsize/ntiles_z;

    x_tile += ceill((xsize-x_tile*ntiles_x) / (double) ntiles_x);
    y_tile += ceill((ysize-y_tile*ntiles_y) / (double) ntiles_y);
    z_tile += ceill((zsize-z_tile*ntiles_z) / (double) ntiles_z);

    // Since the tile size is rounded up, the number of tiles may actually decrease. Hence, we update the number of tiles
    // for each axis
    ntiles_x = ceill(xsize / (double) x_tile);
    ntiles_y = ceill(ysize / (double) y_tile);
    ntiles_z = ceill(zsize / (double) z_tile);

    n = (ntiles_x)*(ntiles_y)*(ntiles_z);

    tile_coords = (iftBoundingBox*)iftAlloc(n, sizeof(iftBoundingBox));

    i = 0;
    start.z = 0;
    for(z = 0; z < ntiles_z; z++) {
        start.y = 0;
        end.z = iftMin(start.z + z_tile - 1, zsize - 1);

        for (y = 0; y < ntiles_y; y++) {
            start.x = 0;
            end.y = iftMin(start.y + y_tile - 1, ysize - 1);

            for (x = 0; x < ntiles_x; x++) {
                end.x = iftMin(start.x + x_tile - 1, xsize - 1);

                tile_coords[i].begin = start;
                tile_coords[i].end = end;

                start.x = end.x + 1;
                i++;
            }

            start.y = end.y + 1;
        }

        start.z = end.z + 1;
    }

    tiles = iftCreateImageTiles(tile_coords, ntiles_x, ntiles_y, ntiles_z, bb);

    return tiles;
}


iftImageTiles *iftImageTilesByStriding(iftImage *img, int tile_xsize, int tile_ysize, int tile_zsize,
                                       int xstride, int ystride, int zstride) {
    iftBoundingBox bb;

    bb.begin.x = bb.begin.y = bb.begin.z = 0;
    bb.end.x = img->xsize - 1;
    bb.end.y = img->ysize - 1;
    bb.end.z = img->zsize - 1;

    // Forcing the zstride and tile_zsize to be 1 if the image is 2D
    if(!iftIs3DImage(img)) {
        zstride = 1;
        tile_zsize = 1;
    }

    return iftBoundingBoxImageTilesByStriding(bb, tile_xsize, tile_ysize, tile_zsize, xstride, ystride, zstride);


}

iftImageTiles *iftBoundingBoxImageTilesByStriding(iftBoundingBox bb, int tile_xsize, int tile_ysize,
                                                  int tile_zsize, int xstride, int ystride, int zstride) {
    iftImageTiles *tiles = NULL;
    int i = 0, ntiles_x, ntiles, ntiles_y, ntiles_z, xsize, ysize, zsize;
    iftVoxel start, end;

    iftBoundingBox *tile_coords = NULL;

    // Forcing the zstride and tile_zsize to be 1 if the image/bounding box is 2D
    if(bb.end.z-bb.begin.z+1 == 1) {
        zstride = 1;
        tile_zsize = 1;
    }

    if(xstride <= 0 || ystride <= 0 || zstride <= 0) {
        iftError("The stride used for computing image tiles must be positive, and not (%d, %d, %d) as provided!",
                 "iftBoundingBoxImageTilesByStriding", xstride, ystride, zstride);
    }

    xsize = bb.end.x - bb.begin.x + 1;
    ysize = bb.end.y - bb.begin.y + 1;
    zsize = bb.end.z - bb.begin.z + 1;

    ntiles_x = ceill(xsize / (double) xstride);
    ntiles_y = ceill(ysize / (double) ystride);
    ntiles_z = ceill(zsize / (double) zstride);

    ntiles = ntiles_x * ntiles_y * ntiles_z;

    tile_coords = (iftBoundingBox*)iftAlloc(ntiles, sizeof(iftBoundingBox));

    for (start.z = 0; start.z < zsize; start.z += zstride) {
        for (start.y = 0; start.y < ysize; start.y += ystride) {
            for (start.x = 0; start.x < xsize; start.x += xstride) {
                end.x = iftMin(start.x + tile_xsize - 1, xsize - 1);
                end.y = iftMin(start.y + tile_ysize - 1, ysize - 1);
                end.z = iftMin(start.z + tile_zsize - 1, zsize - 1);

                tile_coords[i].begin = start;
                tile_coords[i].end = end;
                i++;
            }
        }
    }

    tiles = iftCreateImageTiles(tile_coords, ntiles_x, ntiles_y, ntiles_z, bb);

    return tiles;
}

iftSet *iftGetIndicesFromTilesIntersectingVoxel(iftImageTiles *tiles, iftVoxel v, bool get_only_first_tile) {
    int cur_tile;
    bool stop = false;
    iftVoxel uo, uf, u_tile;
    iftSet *tile_indices = NULL;

    if(v.x >= tiles->bb.begin.x && v.x <= tiles->bb.end.x && v.y >= tiles->bb.begin.y && v.y <= tiles->bb.end.y
       && v.z >= tiles->bb.begin.z && v.z <= tiles->bb.end.z) {
        // Computing the tile index for the current voxel if it is within the dimensions
        // of the bounding box that encloses the tiles
        for (u_tile.z = 0; u_tile.z < tiles->ntiles_z && !stop; u_tile.z++) {
            for (u_tile.y = 0; u_tile.y < tiles->ntiles_y && !stop; u_tile.y++) {
                for (u_tile.x = 0; u_tile.x < tiles->ntiles_x && !stop; u_tile.x++) {
                    // Getting the tiles->bb.begin voxel from the current tile
                    cur_tile = iftGetTileIndex(tiles, u_tile.x, u_tile.y, u_tile.z);
                    uo = tiles->tile_coords[cur_tile].begin;
                    uf = tiles->tile_coords[cur_tile].end;

                    // First voxel in the coordinates of the original image
                    uo.x += tiles->bb.begin.x;
                    uo.y += tiles->bb.begin.y;
                    uo.z += tiles->bb.begin.z;

                    // Last voxel in the coordinates of the original image
                    uf.x += tiles->bb.begin.x;
                    uf.y += tiles->bb.begin.y;
                    uf.z += tiles->bb.begin.z;

                    // If the voxel is within the current tile then we have found the proper index
                    if (v.x >= uo.x && v.x <= uf.x && v.y >= uo.y && v.y <= uf.y && v.z >= uo.z && v.z <= uf.z) {
                        iftInsertSet(&tile_indices, cur_tile);
                        stop = (get_only_first_tile) ? true : false;
                    }
                }
            }
        }
    }

    return tile_indices;
}


int iftGetIndexFromFirstTileIntersectingVoxel(iftImageTiles *tiles, iftVoxel v)
{
    int tile = IFT_NIL;
    iftSet *intersecting_tiles = NULL;

    intersecting_tiles = iftGetIndicesFromTilesIntersectingVoxel(tiles, v, true);

    if(intersecting_tiles != NULL)
        tile = intersecting_tiles->elem;

    iftDestroySet(&intersecting_tiles);

    return tile;
}

iftSet* iftGetIndicesFromAllTilesIntersectingVoxel(iftImageTiles *tiles, iftVoxel v)
{
    return iftGetIndicesFromTilesIntersectingVoxel(tiles, v, false);
}

iftImage *iftExtractTile(iftImage *img, iftImageTiles *tiles, int tile)
{
    iftImage *tile_img = NULL;
    iftBoundingBox bb;

    bb = iftGetTileInOriginalImageCoordinates(tiles, tile);

    if(!iftValidVoxel(img, bb.begin) || !iftValidVoxel(img, bb.end))
        iftError("The image tile structure does not refer to the image from which tile %d is being extracted!" \
                 "The sizes differ: image (%d, %d, %d), tile structure bounding box (%d, %d, %d)-(%d, %d, %d), "\
                 "tile bounding box (%d, %d, %d)-(%d, %d, %d).", "iftExtractTile", tile, img->xsize, img->ysize,
                 img->zsize,
                 tiles->bb.begin.x, tiles->bb.begin.y, tiles->bb.begin.z, tiles->bb.end.x, tiles->bb.end.y,
                 tiles->bb.end.z,
                 bb.begin.x, bb.begin.y, bb.begin.z, bb.end.x, bb.end.y, bb.end.z);

    tile_img = iftExtractROI(img, bb);

    return tile_img;
}

iftBoundingBox iftGetTileInOriginalImageCoordinates(iftImageTiles *tiles, int tile)
{
    iftBoundingBox bb;

    bb.begin.x = tiles->tile_coords[tile].begin.x + tiles->bb.begin.x;
    bb.begin.y = tiles->tile_coords[tile].begin.y + tiles->bb.begin.y;
    bb.begin.z = tiles->tile_coords[tile].begin.z + tiles->bb.begin.z;

    bb.end.x = tiles->tile_coords[tile].end.x + tiles->bb.begin.x;
    bb.end.y = tiles->tile_coords[tile].end.y + tiles->bb.begin.y;
    bb.end.z = tiles->tile_coords[tile].end.z + tiles->bb.begin.z;

    return bb;
}

iftImage **iftSplitImageIntoTilesByEquallyDividingAxes(iftImage *img, int ntiles_x, int ntiles_y, int ntiles_z,
                                                       iftImageTiles **tiles) {
    iftImage **img_tiles = NULL;

    iftDestroyImageTiles(tiles);

    if(tiles == NULL)
        iftError("Pointer <tiles> may not be NULL since it returns the data structure\n"
                 "with the actual number of computed tiles!", "iftSplitImageIntoTilesByEquallyDividingAxes");

    *tiles = iftImageTilesByEquallyDividingAxes(img, ntiles_x, ntiles_y, ntiles_z);

    img_tiles = (iftImage**)iftAlloc((*tiles)->ntiles, sizeof(iftImage*));

    for(int i = 0; i < (*tiles)->ntiles; i++){
        img_tiles[i] = iftExtractTile(img, (*tiles), i);
    }

    return img_tiles;
}

iftImage **iftSplitImageIntoTilesByStriding(iftImage *img, int tile_xsize, int tile_ysize, int tile_zsize, int xstride,
                                            int ystride, int zstride, iftImageTiles **tiles) {
    iftImage **img_tiles = NULL;

    iftDestroyImageTiles(tiles);

    if(tiles == NULL)
        iftError("Pointer <tiles> may not be NULL since it returns the data structure\n"
                 "with the actual number of computed tiles!", "iftSplitImageIntoTilesByStriding");

    *tiles = iftImageTilesByStriding(img, tile_xsize, tile_ysize, tile_zsize, xstride, ystride, zstride);

    img_tiles = (iftImage**)iftAlloc((*tiles)->ntiles, sizeof(iftImage*));

    for(int i = 0; i < (*tiles)->ntiles; i++) {
        img_tiles[i] = iftExtractTile(img, (*tiles), i);
    }

    return img_tiles;
}
iftImage *iftRecomposeImageFromTiles(iftImage **img_tiles, iftImageTiles *tiles)
{
    int i, n;
    int xsize, ysize, zsize;
    iftImage *img = NULL;

    n = tiles->ntiles;
    // Determining the recomposed image's size from the last tile's coordinates
    // and size
    xsize = tiles->tile_coords[n - 1].begin.x + img_tiles[n - 1]->xsize;
    ysize = tiles->tile_coords[n - 1].begin.y + img_tiles[n - 1]->ysize;
    zsize = tiles->tile_coords[n - 1].begin.z + img_tiles[n - 1]->zsize;

    img = iftCreateImage(xsize, ysize, zsize);
    // We assume that all tiles have the same pixel size
    iftCopyVoxelSize(img_tiles[0], img);

    for(i = 0; i < n; i++) {
        iftInsertROI(img_tiles[i], img, tiles->tile_coords[i].begin);
    }

    return img;
}


// Similarly to iftAddFrame, iftAddPadding pads the image with the given value
iftImage *iftAddPadding(iftImage *img, int dx, int dy, int dz, int value) {
  return iftAddRectangularBoxFrame(img, dx, dy, dz, value);
}

// Similarly to iftRemFrame, iftRemPadding removes image padding
iftImage *iftRemPadding(iftImage *fimg, int dx, int dy, int dz) {
    return iftRemRectangularBoxFrame(fimg, dx, dy, dz);
}

bool iftIsLabelImage(const iftImage *label_img, int n_objects) {
    if (label_img == NULL)
        iftError("Label Image is NULL", "iftIsLabelImage");
    if (n_objects <= 0)
        iftError("Number of Objects %d is <= 0", "iftIsLabelImage", n_objects);

    iftBMap *contain_label = iftCreateBMap(n_objects+1); // [0] is the background

    for (int p = 0; p < label_img->n; p++)
        if ((label_img->val[p] < 0) || (label_img->val[p] > n_objects)) {
            puts("\n*****");
            iftWarning("Label %d is out of the range [0..%d]", "iftIsLabelImage",
                       label_img->val[p], n_objects);
            puts("*****");

            iftDestroyBMap(&contain_label);
            return false;
        }
        else {
            iftBMapSet1(contain_label, label_img->val[p]);
        }


    for (int o = 1; o <= n_objects; o++)
        if (!iftBMapValue(contain_label, o)) {
            puts("\n*****");
            iftWarning("Label %d is missing!", "iftIsLabelImage", o);
            puts("*****");
            return false;
        }

    return true;
}


int iftGetNumberOfObjectsFromImage(const iftImage *label_img) {
    if (label_img == NULL)
        iftError("Label Image is NULL", "iftGetNumberOfObjectsFromImage");

    return iftMaximumValue(label_img);
}


int iftGetNumberOfObjectsFromImagePathname(const char *label_img_pathname) {
    iftImage *label_img = iftReadImageByExt(label_img_pathname);
    int n_objects       = iftGetNumberOfObjectsFromImage(label_img);
    iftDestroyImage(&label_img);

    return n_objects;
}


iftImage *iftBMapToBinImage(const iftBMap *bmap, int xsize, int ysize, int zsize) {
    if (bmap == NULL)
        iftError("Bin Map is NULL", "iftBMapToBinImage");

    int n = xsize * ysize * zsize;

    if (bmap->n != n)
        iftError("Image Domain is != of the Bin Map Size\n" \
                 "Img Domain: (%d, %d, %d) = %d spels\nBin Map: %d elems", "iftBMapToBinImage",
                 xsize, ysize, zsize, n, bmap->n);

    iftImage *bin = iftCreateImage(xsize, ysize, zsize);

    for (int p = 0; p < bmap->n; p++)
        bin->val[p] = iftBMapValue(bmap, p);

    return bin;
}


iftBMap *iftBinImageToBMap(const iftImage *bin_img) {
    if (bin_img == NULL)
        iftError("Bin Image is NULL", "iftBinImageToBMap");

    iftBMap *bmap = iftCreateBMap(bin_img->n);

    for (int i = 0; i < bin_img->n; i++)
        if (bin_img->val[i])
            iftBMapSet1(bmap, i);

    return bmap;
}


bool iftIsBinaryImage(const iftImage *img) {  
  iftIntArray *labels = iftGetObjectLabels(img);
  bool is_binary = (labels->n == 1);
  iftDestroyIntArray(&labels);
  return is_binary;
}


iftIntArray *iftGetObjectLabels(const iftImage *label_img) {
    if (label_img == NULL)
        iftError("Label Image is NULL", "iftGetObjectLabels");

    int max_label       = iftMaximumValue(label_img);
    iftBMap *label_bmap = iftCreateBMap(max_label+1);

    for (int p = 0; p < label_img->n; p++)
        iftBMapSet1(label_bmap, label_img->val[p]);

    int n_labels = 0;
    for (int i = 1; i <= max_label; i++)
        if (iftBMapValue(label_bmap, i))
            n_labels++;

    iftIntArray *labels = iftCreateIntArray(n_labels);
    int idx = 0;
    for (int i = 1; i <= max_label; i++)
        if (iftBMapValue(label_bmap, i)) {
            labels->val[idx] = i;
            idx++;
        }


    iftDestroyBMap(&label_bmap);

    return labels;
}

iftIntArray *iftFindObjectLabels(const iftImage *label_img) {
    int max_label = iftMaximumValue(label_img);
    iftIntArray *map = iftCreateIntArray(max_label + 1);

    #pragma omp parallel for
    for (int p = 0; p < label_img->n; p++) {
        int label = label_img->val[p];
        
        if (label_img->val[p]) {
            map->val[label] = true;
        }
    }

    return map;
}


iftIntArray *iftGetObjectVoxels(const iftImage *label_img) {
    iftSet *obj_set = NULL;
    
    long n = 0;
    
    for (int p = 0; p < label_img->n; p++)
        if (label_img->val[p]) {
            iftInsertSet(&obj_set, p);
            n++;
        }
        
    iftIntArray *obj_arr = iftCreateIntArray(n);
    
    long i = 0;
    while (obj_set != NULL)
        obj_arr->val[i++] = iftRemoveSet(&obj_set);
    
    return obj_arr;
}



int iftGetLabelByMajorityVoting(const iftImage *label_img, bool ignore_zero) {
    if (label_img == NULL)
        iftError("Label Image is NULL", "iftGetLabelByMajorityVoting");

    int max_label = iftMaximumValue(label_img);

    if ((max_label == 0) && (ignore_zero)) {
        iftWarning("Maximum value from the Image is 0 and ignore_zero is True. The zero will be returned",
                   "iftGetLabelByMajorityVoting");
        return 0;
    }


    long *occurs = iftAllocLongIntArray(max_label+1);

    // counts the label occurrences
    for (int p = 0; p < label_img->n; p++)
        occurs[label_img->val[p]]++;

    // sets the initial label for the majority voting
    int start_label = 0;
    if (ignore_zero)
        start_label = 1;


    int most_frequent_label = 0;
    long max_occur          = IFT_INFINITY_LONG_NEG;

    // majority voting
    for (int label = start_label; label <= max_label; label++) {
        if (max_occur < occurs[label]) {
            most_frequent_label = label;
            max_occur           = occurs[label];
        }
    }

    iftFree(occurs);

    return most_frequent_label;
}


iftIntArray *iftGetSupervoxelTrueLabelByMajorityVoting(const iftImage *label_img, const iftImage *super_img) {
    if (label_img == NULL)
        iftError("Label Image is NULL", "iftGetSupervoxelTrueLabelByMajorityVoting");
    if (super_img == NULL)
        iftError("Superpixel Image is NULL", "iftGetSupervoxelTrueLabelByMajorityVoting");
    iftVerifyImageDomains(label_img, super_img, "iftGetSupervoxelTrueLabelByMajorityVoting");

    int max_true_label = iftMaximumValue(label_img);
    int n_cells        = iftMaximumValue(super_img);

    if (n_cells <= 0) {
        iftError("Invalid Number of Superpixels: %d\n", "iftGetSupervoxelTrueLabelByMajorityVoting",
                 n_cells);
    }


    // matrix (n_supervoxel, n_true_labels+1) of the true label occurrence for each superpixel
    iftIntMatrix *occurs = iftCreateIntMatrix(max_true_label+1, n_cells);

    // counts the true label occurrences for each superpixel
    for (int p = 0; p < super_img->n; p++) {
        const int cell       = super_img->val[p] - 1; // ex: the supervoxel with label 2 has cell (idx) 1
        const int true_label = label_img->val[p];
        if (true_label >= 0) {
            int idx = iftGetMatrixIndex(occurs, true_label, cell);
            occurs->val[idx]++;
        }
    }


    iftIntArray *most_frequent_labels = iftCreateIntArray(n_cells);

    // majority voting
    for (int cell = 0; cell < n_cells; cell++)
    {
        long max_occur = 0;
        most_frequent_labels->val[cell] = -1;
        for (int true_label = 0; true_label <= max_true_label; true_label++) {
            int idx = iftGetMatrixIndex(occurs, true_label, cell);

            if (max_occur < occurs->val[idx]) {
                most_frequent_labels->val[cell] = true_label;
                max_occur                       = occurs->val[idx];
            }
        }
    }

    iftDestroyIntMatrix(&occurs);

    return most_frequent_labels;
}


iftImage *iftAdjRelToImage(const iftAdjRel *adj, int xsize, int ysize, int zsize, iftVoxel ref_voxel) {
    // checkers
    if (adj == NULL)
        iftError("Adjacency is NULL", "iftAdjRelToImage");
    if (xsize <= 0)
        iftError("X Size %d <= 0", "iftAdjRelToImage", xsize);
    if (ysize <= 0)
        iftError("Y Size %d <= 0", "iftAdjRelToImage", ysize);
    if (zsize <= 0)
        iftError("Z Size %d <= 0", "iftAdjRelToImage", zsize);

    iftImage *bin_img = iftCreateImage(xsize, ysize, zsize);

    if (!iftValidVoxel(bin_img, ref_voxel))
        iftError("Reference Voxel (%d, %d, %d) does not belong to the Resulting Bin Image Domain:\n" \
                 "xsize = %d\nysize = %d\nzsize = %d\n",
                 "iftValidVoxel", ref_voxel.x, ref_voxel.y, ref_voxel.z,
                 bin_img->xsize, bin_img->ysize, bin_img->zsize);

#pragma omp parallel for
    for (int i = 0; i < adj->n; i++) {
        iftVoxel v;
        v.x = ref_voxel.x + adj->dx[i];
        v.y = ref_voxel.y + adj->dy[i];
        v.z = ref_voxel.z + adj->dz[i];

        if (iftValidVoxel(bin_img, v)) {
            int p = iftGetVoxelIndex(bin_img, v);
            bin_img->val[p] = 1;
        }
    }

    return bin_img;
}

/* Assigns colors to the luminance values of an input image according
   to the blue to red (rainbow) color map */

iftImage *iftColorCoding(iftImage* img)
{
    float   maxval = iftMaximumValue(img);
    float   minval = iftMinimumValue(img);
    float   v;
    iftImage *cimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
    int     normalization_value;
    if (iftMaximumValue(img)<=255){
        normalization_value=255;
        iftSetCbCr(cimg,128);
    } else {
        normalization_value=65535;
        iftSetCbCr(cimg,(normalization_value+1)/2);
    }

    for (int p=0; p < img->n; p++) {
        v = img->val[p];
        v = (v-minval)/(maxval-minval);
        float r,g,b;
        iftHeatColorMapping(v, &r, &g, &b);
        iftColor RGB,YCbCr;
        RGB.val[0]   = r*normalization_value;
        RGB.val[1]   = g*normalization_value;
        RGB.val[2]   = b*normalization_value;
        YCbCr        = iftRGBtoYCbCr(RGB,normalization_value);
        cimg->val[p] = (int)YCbCr.val[0];
        cimg->Cb[p]  = (unsigned short)YCbCr.val[1];
        cimg->Cr[p]  = (unsigned short)YCbCr.val[2];
    }

    return(cimg);
}



iftImagePlaneOrientation iftGetImagePlaneOrientation(const char *orientation) {
    char *upper_orientation = iftUpperString(orientation);

    if (iftCompareStrings(upper_orientation, "AXIAL")) {
        return IFT_AXIAL;
    }
    else if (iftCompareStrings(upper_orientation, "CORONAL")) {
        return IFT_CORONAL;
    }
    else if (iftCompareStrings(upper_orientation, "SAGITTAL")) {
        return IFT_SAGITTAL;
    }
    else {
        iftError("Invalid string for orientation: \"%s\"\n" \
                 "Try: AXIAL, CORONAL, or SAGITTAL", "iftGetImagePlaneOrientation", orientation);
        return false; // just to avoid warning
    }

    iftFree(upper_orientation);
}

void iftWriteVTKImage(iftImage *img, char *filename)        //Write volume as VTK binary integer volume (for visualization)
{                                                           //Data is saved as cell-data, for ease of use in ParaView
    FILE *fp=NULL;
    int   p;

    int minval = iftMinimumValue(img);
    int maxval = iftMaximumValue(img);

    if (minval < 0)
    {
        char msg[200];
        sprintf(msg,"Shifting image values from [%d,%d] to [%d,%d] on the original image\n",minval,maxval,0,maxval-minval);
        iftWarning(msg,"iftWriteImage");
        for (p=0; p < img->n; p++)
            img->val[p] = img->val[p] - minval;
        maxval = maxval - minval;
    }

    fp = fopen(filename,"wb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, filename, "iftWriteVTKImage");

    fprintf(fp,"# vtk DataFile Version 3.0\n");
    fprintf(fp,"#generated by iftSkel\n");
    fprintf(fp,"BINARY\n");
    fprintf(fp,"DATASET STRUCTURED_POINTS\n");
    fprintf(fp,"DIMENSIONS %d %d %d\n",img->xsize+1,img->ysize+1,img->zsize+1);
    fprintf(fp,"ORIGIN 0 0 0\n");                           //ALEX: add +1 to point-dimensions since we're saving data as cells
    fprintf(fp,"SPACING 1 1 1\n");                          //ALEX: May need to be img->dx,img->dy,img->dz
    fprintf(fp,"CELL_DATA %d\n",img->xsize*img->ysize*img->zsize);
    fprintf(fp,"SCALARS voxel_data int\n");
    fprintf(fp,"LOOKUP_TABLE default\n");

    for(int i=0;i<img->n;++i)
    {
        int v = img->val[i];
        v = IntSwap(v);
        fwrite(&v,1,sizeof(v),fp);
    }

    fclose(fp);
}

iftImage *iftReadVTKImage(const char *filename)         //Load binary unsigned-char VTK volumes. These follow the
{                                                       //very particular (and apparently formally wrong) format
    iftImage  *img=NULL;                                //produced by binvox.
    FILE    *fp=NULL;
    uchar   *data8=NULL;
    int     p,xsize,ysize,zsize,ox,oy,oz,dx,dy,dz,dummy;
    char    buffer[128];

    fp = fopen(filename,"rb");

    if (fp == NULL){
        iftError(MSG_FILE_OPEN_ERROR, filename, "iftReadVTKImage");
    }

    while(!feof(fp))                                     //Skip comments
    {
        if (fgets(buffer, 128, fp)==NULL)
            iftError("Reading error", "iftReadVTKImage");
        if(buffer[0]=='#') continue;
        break;
    }

    if(strcmp(buffer, "BINARY\n"))
        iftError("cannot load ASCII files", "iftReadVTKImage");

    if (fgets(buffer, 128, fp)==NULL)
        iftError("Reading error", "iftReadVTKImage");
    if(strcmp(buffer, "DATASET STRUCTURED_POINTS\n"))
        iftError("can load only VTK structured points", "iftReadVTKImage");

    if (fgets(buffer, 128, fp)==NULL)
        iftError("Reading error", "iftReadVTKImage");
    sscanf(buffer, "DIMENSIONS %d %d %d\n", &zsize, &xsize, &ysize);
    if (fgets(buffer, 128, fp)==NULL)
        iftError("Reading error", "iftReadVTKImage");
    //ALEX: 1st deviation: dimensions come as z,x,y
    sscanf(buffer, "ORIGIN %d %d %d\n", &ox, &oy, &oz);

    if (fgets(buffer, 128, fp)==NULL)
        iftError("Reading error", "iftReadVTKImage");
    sscanf(buffer, "SPACING %d %d %d\n", &dx, &dy, &dz);

    if (fgets(buffer, 128, fp)==NULL)
        iftError("Reading error", "iftReadVTKImage");
    sscanf(buffer, "POINT_DATA %d\n", &dummy);           //ALEX: 2nd deviation: data is saved as points, not cells

    if (fgets(buffer, 128, fp)!=NULL) {};                              //Skip SCALARS info
    if (fgets(buffer, 128, fp)!=NULL) {};                              //Skip LOOKUP_TABLE info

    img = iftCreateImage(xsize,ysize,zsize);
    img->dx = dx;
    img->dy = dy;
    img->dz = dz;

    data8  = iftAllocUCharArray(img->n);
    if (fread(data8,sizeof(uchar),img->n,fp)!=img->n)   //ALEX: 3rd deviation: data is in native endian format,
        iftError("Reading error", "iftReadVTKImage");    //      and not the VTK mandatory big endian
    for (p=0; p < img->n; p++)
        img->val[p] = (int)(data8[p]!=0);               //Threshold the input to binary (0,1) image
    iftFree(data8);

    int minval = iftMinimumValue(img);
    int maxval = iftMaximumValue(img);

    printf("Input: size %d %d %d; data range %f %f\n",img->xsize,img->ysize,img->zsize,(float)minval,(float)maxval);

    fclose(fp);
    return(img);
}

int iftMaximumValueInMask(const iftImage *img, const iftImage *mask)
{
    if (mask == NULL)
        return iftMaximumValue(img);
    else {
        iftVerifyImageDomains(img, mask, "iftMaximumValueInMask");

        int maxval = IFT_INFINITY_INT_NEG;
    
        for (int p = 0; p < img->n; p++) {
            if (mask->val[p]!=0){
                if (img->val[p] > maxval)
                    maxval = img->val[p];
            }
        }
    
        return maxval;
    }
}

int iftMinimumValueInMask(const iftImage *img, iftImage *mask)
{
    char   no_mask=0;
    int    minval=IFT_INFINITY_INT;

    if (mask == NULL) {
        mask = iftCreateImage(img->xsize,img->ysize,img->zsize);
        iftSetImage(mask,1);
        no_mask = 1;
    } else {
        iftVerifyImageDomains(img,mask,"iftMinimumValueInMask");
    }

    for (int p = 0; p < img->n; p++) {
        if (mask->val[p]!=0){
            if (img->val[p] < minval)
                minval = img->val[p];
        }
    }

    if (no_mask)
        iftDestroyImage(&mask);

    return(minval);
}

void iftMinimumValuesForLabels(const iftImage *img, const iftImage *label_img,
                               iftIntArray **min_vals_out, iftIntArray **max_vals_out) {
    int n_objs = iftMaximumValue(label_img);
    iftIntArray *min_vals = iftCreateIntArray(n_objs + 1);
    iftIntArray *max_vals = iftCreateIntArray(n_objs + 1);

    for (int label = 0; label <= n_objs; label++) {
        min_vals->val[label] = IFT_INFINITY_INT;
        max_vals->val[label] = IFT_INFINITY_INT_NEG;
    }


    for (int p = 0; p < img->n; p++) {
        int label = label_img->val[p];

        if (img->val[p] < min_vals->val[label]) {
            min_vals->val[label] = img->val[p];
        }
        if (img->val[p] > max_vals->val[label]) {
            max_vals->val[label] = img->val[p];
        }
    }


    if (min_vals_out)
        *min_vals_out = min_vals;
    else iftDestroyIntArray(&min_vals);

    if (max_vals_out)
        *max_vals_out = max_vals;
    else
        iftDestroyIntArray(&max_vals);
}


iftImage *iftTranslateImageContent(const iftImage *img, iftVector disp_vec) {
    if (img == NULL)
        iftError("Image is NULL", "iftShiftImage");
    
    iftCopyImageVoxelFunc copyImageVoxelFunc;
    if (iftIsColorImage(img))
        copyImageVoxelFunc = iftCopyColorImageVoxel;
    else copyImageVoxelFunc = iftCopyGrayImageVoxel;
    
    iftImage *trans_img = iftCreateImageFromImage(img);

    #pragma omp parallel for
    for (int p = 0; p < img->n; p++) {
        iftVoxel u = iftGetVoxelCoord(img, p);
        iftVoxel v = {.x = u.x + disp_vec.x, .y = u.y + disp_vec.y, .z = u.z + disp_vec.z};
        
        if (iftValidVoxel(img, v))
            copyImageVoxelFunc(img, u, trans_img, v);
    }
    
    return trans_img;
}


void iftSetImageBoundingBoxValue(iftImage *img, iftBoundingBox bb, int value )
{
    if (img == NULL)
        iftError("Image is NULL", "iftSetImageBoundingBoxValue");

    iftVoxel upper_left_voxel = bb.begin;
    iftVoxel bottom_right_voxel = bb.end;

    if (!iftValidVoxel(img, upper_left_voxel))
        iftError(
                "Initial Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for image with size (%d, %d, %d)",
                "iftSetImageBoundingBoxValue", upper_left_voxel.x, upper_left_voxel.y, upper_left_voxel.z, img->xsize,
                img->ysize, img->zsize);

    if (!iftValidVoxel(img, bottom_right_voxel)) {
        iftWarning("Final Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for image with size (%d, %d, %d).\n It will copied/inserted what it is possible", "iftSetImageBoundingBoxValue", bottom_right_voxel.x, bottom_right_voxel.y, bottom_right_voxel.z, img->xsize, img->ysize, img->zsize);
        // gets the valid ending voxel
        bottom_right_voxel.x = iftMin(bottom_right_voxel.x, img->xsize - 1);
        bottom_right_voxel.y = iftMin(bottom_right_voxel.y, img->ysize - 1);
        bottom_right_voxel.z = iftMin(bottom_right_voxel.z, img->zsize - 1);
    }

    iftVoxel u;
    for (u.y = upper_left_voxel.y; u.y <= bottom_right_voxel.y; u.y++)
        for (u.z = upper_left_voxel.z; u.z <= bottom_right_voxel.z; u.z++)
            for (u.x = upper_left_voxel.x; u.x <= bottom_right_voxel.x; u.x++) {
                int p = iftGetVoxelIndex(img, u);
                img->val[p]=value;
            }
}

iftImage *iftRelabelGrayScaleImage(const iftImage *img,bool increment_values){

    if (img == NULL || iftIsColorImage(img))
        iftError("Img must be a grayscale image", "iftRelabelGrayScaleImage");

    iftImage *out_img=iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftDict *dict=iftCreateDict();

    int current_val=0+increment_values;

    for(int p = 0; p < img->n; p++) {
        if (!iftDictContainKey(img->val[p],dict,NULL)){
            iftInsertIntoDict(img->val[p],current_val,dict);
            current_val++;
        }
    }

#pragma omp parallel for
    for(int p = 0; p < img->n; p++) {
        out_img->val[p] = (int)iftGetLongValFromDict(img->val[p],dict);
    }

    iftDestroyDict(&dict);

    return out_img;
}


iftImage *iftExtractLargestObjectInLabelImage(iftImage *label_img){

    iftImage *rlabel,*nlabel;
    iftAdjRel *A;
    if (iftIs3DImage(label_img))
        A = iftSpheric(sqrtf(1.0));
    else
        A = iftCircular(sqrtf(1.0));

    rlabel = iftRelabelRegions(label_img,A);

    int max_label=iftMaximumValue(rlabel);
    int *obj_count=iftAllocIntArray(max_label+1);

    for (int p=0;p<label_img->n;p++){
        if (label_img->val[p] > 0)
            obj_count[rlabel->val[p]]++;
    }

    /*get the biggest object component in the image*/
    int bigger_size=0;
    int bigger_object=-1;
    for (int p=0;p<=max_label;p++)
        if (obj_count[p] > bigger_size){
            bigger_size=obj_count[p];
            bigger_object=p;
        }

    nlabel=iftCreateImage(label_img->xsize,label_img->ysize,label_img->zsize);
    for (int p=0;p<label_img->n;p++)
        if (rlabel->val[p]==bigger_object)
            nlabel->val[p]=1;
        else
            nlabel->val[p]=0;

    iftDestroyAdjRel(&A);
    iftDestroyImage(&rlabel);
    iftFree(obj_count);

    return nlabel;

}

iftImage *iftTranslateObject(iftImage *label, iftVoxel dst, iftVoxel *objCenter)
{
    iftImage *centered = NULL;
    iftBoundingBox bb;
    iftVoxel begin,end,u,v,T;
    int p,q,x,y,z;

    begin.x = label->xsize;
    begin.y = label->ysize;
    begin.z = label->zsize;

    end.x = 0;
    end.y = 0;
    end.z = 0;


    for (int p = 0; p < label->n; p++){
        if (label->val[p] != 0) {
            u = iftGetVoxelCoord(label, p);
            if (u.x < begin.x)
                begin.x = u.x;
            if (u.y < begin.y)
                begin.y = u.y;
            if (u.z < begin.z)
                begin.z = u.z;
            if (u.x > end.x)
                end.x = u.x;
            if (u.y > end.y)
                end.y = u.y;
            if (u.z > end.z)
                end.z = u.z;
        }
    }

    bb.begin = begin;
    bb.end = end;

    objCenter->x = (bb.end.x + bb.begin.x) / 2;
    objCenter->y = (bb.end.y + bb.begin.y) / 2;
    objCenter->z = (bb.end.z + bb.begin.z) / 2;

    T.x = dst.x - objCenter->x;
    T.y = dst.y - objCenter->y;
    T.z = dst.z - objCenter->z;

    centered = iftCreateImageFromImage(label);

    for (x = bb.begin.x; x < bb.end.x; x++)
        for (y = bb.begin.y; y < bb.end.y; y++)
            for (z = bb.begin.z; z < bb.end.z; z++){
                v.x = x;
                v.y = y;
                v.z = z;
                u.x = x + T.x;
                u.y = y + T.y;
                u.z = z + T.z;
                if ((iftValidVoxel(centered,u)) && (iftValidVoxel(label,v))) {
                    p = iftGetVoxelIndex(label, v);
                    q = iftGetVoxelIndex(centered, u);
                    if (label->val[p] != 0) {
                        centered->val[q] = label->val[p];
                    }
                }

            }

    return centered;
}

/*****************************************************************************************************************************/

iftImage  *iftImageDomes(iftImage *img, iftAdjRel *A)
{
  iftImage   *domes, *basins;

  basins = iftImageBasins(img, A);
  domes  = iftComplement(basins);
  iftDestroyImage(&basins);

  return (domes);
}



iftImage *iftImageBasins(const iftImage *img, iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL)
        A = (iftIs3DImage(img)) ? iftSpheric(1.0) : iftCircular(1.5);
    else A = iftCopyAdjacency(Ain);

    int dx, dy, dz;
    iftMaxAdjShifts(A, &dx, &dy, &dz);

    iftImage *basins = iftCreateImage(img->xsize,img->ysize,img->zsize);
    float K = sqrtf(dx*dx + dy*dy + dz*dz);     
    float *w = iftAllocFloatArray(A->n); // weight vectors
    float wt = 0.0; // total weight

    
    for (int i = 1; i < A->n; i++) {
      w[i] = K / sqrtf(iftPowerOfTwo(A->dx[i]) + iftPowerOfTwo(A->dy[i]) + iftPowerOfTwo(A->dz[i]));
      wt  += w[i];
    }

    // normalize the weights (influence) of the adjacent pixels

    for (int i = 1; i < A->n; i++) {
      w[i] = (w[i] / wt);
    }

    if(iftIsColorImage(img)) {
#pragma omp parallel for
        for (int p = 0; p < img->n; p++) {
            iftVoxel u = iftGetVoxelCoord(img, p);
            double acc_dist = 0.0f;

            for (int i = 1; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                if (iftValidVoxel(img, v)) {
                    int q = iftGetVoxelIndex(img, v);
                    acc_dist += w[i] * (abs(img->val[q] - img->val[p])+
					abs(img->Cb[q] - img->Cb[p]) + 
					abs(img->Cr[q] - img->Cr[p]));
                }
            }
            basins->val[p] = iftRound(acc_dist);	    
        }
    } else {
#pragma omp parallel for
        for (int p = 0; p < img->n; p++) {
            iftVoxel u = iftGetVoxelCoord(img, p);
            float acc_dist = 0.0f;

            for (int i = 1; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                if (iftValidVoxel(img, v)) {
                    int q = iftGetVoxelIndex(img, v);
                    acc_dist += w[i] * abs(img->val[q] - img->val[p]);
                }
            }
            basins->val[p] = iftRound(acc_dist);
        }
    }

    iftDestroyAdjRel(&A);
    iftFree(w);

    return basins;
}


int iftGetSuperpixelSize(iftImage *spixLabels, int label)
{
    int s = 0;
    for(int p = 0; p < spixLabels->n; p++)
        if(spixLabels->val[p] == label)
            s++;

    return s;
}


void iftGetLabeledPathSides(const iftImage *label, const iftSet *path, iftSet **obj, iftSet **bkg)
{
    if (!path) return;
    *obj = NULL, *bkg = NULL;

    iftVoxel left, right;
    iftVoxel prev_u = iftGetVoxelCoord(label, path->elem);
    for (const iftSet *s = path->next; s; s = s->next)
    {
        iftVoxel u = iftGetVoxelCoord(label, s->elem);
        iftSidePixelsWithJitter(&prev_u, &u, M_SQRT2, 5, &left, &right);

        if (iftValidVoxel(label, left)) {
            int l = iftGetVoxelIndex(label, left);
            if (label->val[l])
                iftInsertSet(obj, l);
            else
                iftInsertSet(bkg, l);
        }

        if (iftValidVoxel(label, right)) {
            int r = iftGetVoxelIndex(label, right);
            if (label->val[r])
                iftInsertSet(obj, r);
            else
                iftInsertSet(bkg, r);
        }
        prev_u = u;
    }
}

iftImage *iftGrayImageToColorImage(const iftImage *img, iftColorTable *ctb){
  iftImage *colored_image;
  float Imax = iftNormalizationValue(iftMaximumValue(img));

  colored_image = iftCreateColorImage(img->xsize, img->ysize, img->zsize, Imax);
				      
  for(int p = 0; p < img->n; p++){
    
    iftColor ycbcr_color = ctb->color[(int)(((float)img->val[p]/Imax)*(ctb->ncolors-1))];

    colored_image->val[p] = ycbcr_color.val[0];
    colored_image->Cb[p]  = ycbcr_color.val[1];
    colored_image->Cr[p]  = ycbcr_color.val[2];
  }
  
  return(colored_image);
}

iftHessianImage *iftCreateHessianImage()
{
    iftHessianImage *H = (iftHessianImage*)iftAlloc(1,sizeof(iftHessianImage));
    return H;
}

void iftDestroyHessianImage(iftHessianImage **Dimgs)
{
    iftDestroyImage(&(*Dimgs)->Dxx);
    iftDestroyImage(&(*Dimgs)->Dyy);
    iftDestroyImage(&(*Dimgs)->Dzz);
    iftDestroyImage(&(*Dimgs)->Dxy);
    iftDestroyImage(&(*Dimgs)->Dxz);
    iftDestroyImage(&(*Dimgs)->Dyz);
    iftFree(*Dimgs);
    *Dimgs = NULL;
}

/*
 * Computes the Hessian matrix of each spel on the image. The Hessian matrix
 * contains the second order derivatives of a given function f. It is written
 * as:
 *        _                        _
 *       |                          |            _               _
 *       | Df/Dxx   Df/Dxy   Df/Dxz |           |                 |
 *       |                          |           | Df/Dxx   Df/Dxy |
 *  H =  | Df/Dyx   Df/Dyy   Df/Dyz |   or H =  |                 |,
 *       |                          |           | Df/Dyx   Df/Dyy |
 *       | Df/Dzx   Df/Dzy   Df/Dzz |           |_               _|
 *       |_                        _|
 *
 * for 3D and 2D images respectively. The Hessian matrix is interesting to
 * identify tubular-shaped structures.
 */
iftHessianImage *iftHessianImageBySobel(iftImage *img)
{
    iftHessianImage *hessian = iftCreateHessianImage();

    if (iftIs3DImage(img)){

        iftKernel *Xk = iftSobelXKernel();
        iftKernel *Yk = iftSobelYKernel();
        iftKernel *Zk = iftSobelZKernel();

        iftImage *Dx = iftLinearFilter(img,Xk);
        iftImage *Dy = iftLinearFilter(img,Yk);
        iftImage *Dz = iftLinearFilter(img,Zk);

        hessian->Dxx = iftLinearFilter(Dx,Xk); //Dxx
        hessian->Dyy = iftLinearFilter(Dy,Yk); //Dyy
        hessian->Dzz = iftLinearFilter(Dz,Zk); //Dzz
        hessian->Dxy = iftLinearFilter(Dx,Yk); //Dxy
        hessian->Dyz = iftLinearFilter(Dy,Zk); //Dyz
        hessian->Dxz = iftLinearFilter(Dx,Zk); //Dxz

        iftDestroyImage(&Dx);
        iftDestroyImage(&Dy);
        iftDestroyImage(&Dz);

        iftDestroyKernel(&Xk);
        iftDestroyKernel(&Yk);
        iftDestroyKernel(&Zk);

    } else {

        iftKernel *Xk = iftSobelXKernel();
        iftKernel *Yk = iftSobelYKernel();

        iftImage *Dx = iftLinearFilter(img,Xk);
        iftImage *Dy = iftLinearFilter(img,Yk);

        hessian->Dxx = iftLinearFilter(Dx,Xk); //Dxx
        hessian->Dxx = iftLinearFilter(Dy,Yk); //Dyy
        hessian->Dxy = iftLinearFilter(Dx,Yk); //Dxy

        iftDestroyImage(&Dx);
        iftDestroyImage(&Dy);

        iftDestroyKernel(&Xk);
        iftDestroyKernel(&Yk);

    }

    return hessian;
}

iftHessianImage *iftHessianImageByGaussian(iftImage *img, float radius, float stdev)
{
    iftHessianImage *H = iftCreateHessianImage();

    iftHessianKernel *K = iftGaussianHessianKernels(radius,stdev);

    if (iftIs3DImage(img)){

        H->Dxx = iftLinearFilter(img,K->Kxx);
        H->Dxy = iftLinearFilter(img,K->Kxy);
        H->Dxz = iftLinearFilter(img,K->Kxz);
        H->Dyy = iftLinearFilter(img,K->Kyy);
        H->Dyz = iftLinearFilter(img,K->Kyz);
        H->Dzz = iftLinearFilter(img,K->Kzz);

    } else {

        H->Dxx = iftLinearFilter(img,K->Kxx);
        H->Dxy = iftLinearFilter(img,K->Kxy);
        H->Dyy = iftLinearFilter(img,K->Kyy);

    }

    iftDestroyHessianKernel(&K);

    return H;
}
