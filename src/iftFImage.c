#include "iftFImage.h"

#include "ift/core/io/Dir.h"
#include "ift/core/io/File.h"
#include "ift/core/io/NumPy.h"
#include "ift/core/io/Stream.h"
#include "ift/core/tools/String.h"



/* ------------------------ Private functions ------------------------- */

float FloatSwap (float f)                   //Convert float data to big-endian, the default
{                                           //binary format required by ParaView
    union
    {
        float f;
        unsigned char b[4];
    } dat1, dat2;

    dat1.f = f;
    dat2.b[0] = dat1.b[3];
    dat2.b[1] = dat1.b[2];
    dat2.b[2] = dat1.b[1];
    dat2.b[3] = dat1.b[0];
    return dat2.f;
}

/* ----------------- Public functions ----------------------------------- */

void iftVerifyFImageDomains(iftFImage *img1, iftFImage *img2, char *function) {
    if ((img1->xsize != img2->xsize) ||
        (img1->ysize != img2->ysize) ||
        (img1->zsize != img2->zsize)) {
        iftError("FImages must have the same domain", function);
    }
}

/**
@file
@brief A description is missing here
*/
char iftIs3DFImage(const iftFImage *img) {
    if (img->zsize > 1)
        return (1);
    else
        return (0);
}

int iftFXSize(iftFImage *img) {
    return (img->xsize);
}

int iftFYSize(iftFImage *img) {
    return (img->ysize);
}

int iftFZSize(iftFImage *img) {
    return (img->zsize);
}

iftVoxel iftFGetVoxelCoord(const iftFImage *img, int p) {
    iftVoxel u;

    u.x = iftFGetXCoord(img, p);
    u.y = iftFGetYCoord(img, p);
    u.z = iftFGetZCoord(img, p);

    return (u);
}


iftFImage *iftCreateFImage(int xsize, int ysize, int zsize) {
    iftFImage *img = NULL;
    int       y, z, xysize;


    img = (iftFImage *) iftAlloc(1, sizeof(iftFImage));
    if (img == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateFImage");
    }

    img->val   = iftAllocFloatArray(xsize * ysize * zsize);
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
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateFImage");
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

iftFImage *iftCreateFImageFromFImage(const iftFImage *src) {
    iftFImage *dst = iftCreateFImage(src->xsize, src->ysize, src->zsize);
    iftCopyVoxelSize(src, dst);

    return dst;
}


void iftDestroyFImage(iftFImage **img) {
    iftFImage *aux;

    aux = *img;
    if (aux != NULL) {
        if (aux->val != NULL) iftFree(aux->val);
        if (aux->tby != NULL) iftFree(aux->tby);
        if (aux->tbz != NULL) iftFree(aux->tbz);
        iftFree(aux);
        *img = NULL;
    }
}

iftFImage *iftCreateConstantFImage(int xsize, int ysize, int zsize, float val) {
    if (xsize <= 0) {
        iftError("Invalid xsize: %d <= 0", "iftCreateConstantFImage", xsize);
    }
    if (ysize <= 0) {
        iftError("Invalid ysize: %d <= 0", "iftCreateConstantFImage", ysize);
    }
    if (zsize <= 0) {
        iftError("Invalid zsize: %d <= 0", "iftCreateConstantFImage", zsize);
    }

    iftFImage *fimg = iftCreateFImage(xsize, ysize, zsize);

    #pragma omp parallel for
    for (int p = 0; p < fimg->n; p++) {
        fimg->val[p] = val;
    }

    return fimg;
}

char iftFValidVoxel(const iftFImage *img, iftVoxel v) {
    if ((v.x >= 0) && (v.x < img->xsize) &&
        (v.y >= 0) && (v.y < img->ysize) &&
        (v.z >= 0) && (v.z < img->zsize))
        return (1);
    else
        return (0);
}


char iftFValidPoint(iftFImage *img, iftPoint P) {
    if ((P.x >= 0) && (P.x < img->xsize) &&
        (P.y >= 0) && (P.y < img->ysize) &&
        (P.z >= 0) && (P.z < img->zsize))
        return (1);
    else
        return (0);
}


void iftFSetImage(iftFImage *img, float value) {
    int p;

    for (p = 0; p < img->n; p++)
        img->val[p] = value;

}

void iftFCopyVoxelSize(const iftFImage *img1, iftFImage *img2) {
    img2->dx = img1->dx;
    img2->dy = img1->dy;
    img2->dz = img1->dz;
}

void iftFCopyVoxelSizeFromImage(const iftImage *img1, iftFImage *img2) {
    img2->dx = img1->dx;
    img2->dy = img1->dy;
    img2->dz = img1->dz;
}

void iftFCopyVoxelSizeToImage(iftFImage *img1, iftImage *img2) {
    img2->dx = img1->dx;
    img2->dy = img1->dy;
    img2->dz = img1->dz;
}


float iftFMaximumValue(const iftFImage *img) {
    float img_max_val = IFT_INFINITY_FLT_NEG;

    for (int p = 0; p < img->n; p++)
        if (img_max_val < img->val[p])
            img_max_val = img->val[p];

    return (img_max_val);
}


float iftFMinimumValue(const iftFImage *img) {
    float img_min_val = IFT_INFINITY_FLT;

    for (int p = 0; p < img->n; p++)
        if (img_min_val > img->val[p])
            img_min_val = img->val[p];

    return (img_min_val);
}

iftFImage *iftFNormalizeImageLocally(iftFImage *img, iftAdjRel *A) {
    iftFImage *fimg = iftCreateFImage(img->xsize, img->ysize, img->zsize);
    float     sum;
    int       p, q, i;
    iftVoxel  u, v;

    for (p = 0; p < img->n; p++) {
        u   = iftFGetVoxelCoord(img, p);
        sum = 0.0;
        for (i = 0; i < A->n; i++) {
            v = iftGetAdjacentVoxel(A, u, i);
            if (iftFValidVoxel(img, v)) {
                q = iftFGetVoxelIndex(img, v);
                sum += (img->val[q] * img->val[q]);
            }
        }
        if (sum > 0.0) {
            sum = sqrtf(sum);
            fimg->val[p] = (float) img->val[p] / sum;
        }
    }

    return (fimg);
}

iftFImage *iftImageToFImage(const iftImage *img1) {
    iftFImage *img2 = iftCreateFImage(img1->xsize, img1->ysize, img1->zsize);
    int       p;

    for (p = 0; p < img2->n; p++) {
        img2->val[p] = (float) img1->val[p];
    }

    img2->dx = img1->dx;
    img2->dy = img1->dy;
    img2->dz = img1->dz;

    return (img2);
}


iftImage *iftFImageToImage(const iftFImage *img1, int Imax) {
    iftImage *img2 = iftCreateImage(img1->xsize, img1->ysize, img1->zsize);
    int      p;

    float img1_max_val = iftFMaximumValue(img1);
    float img1_min_val = iftFMinimumValue(img1);

    if (img1_max_val > img1_min_val) {
        for (p = 0; p < img2->n; p++) {
            img2->val[p] = (int) (Imax * ((img1->val[p] - img1_min_val) /
                                  (img1_max_val - img1_min_val)));
        }
    } else {
        iftWarning("Image is empty", "iftFImageToImage");
    }

    img2->dx = img1->dx;
    img2->dy = img1->dy;
    img2->dz = img1->dz;

    return (img2);
}

iftFImage *iftImageToFImageMaxVal(const iftImage *img, float fmax) {
    iftFImage *fimg = iftCreateFImage(img->xsize, img->ysize, img->zsize);
    iftCopyVoxelSize(img, fimg);

    int img_min_val = iftMinimumValue(img);
    int img_max_val = iftMaximumValue(img);

    if (img_max_val > img_min_val) {
        for (int p = 0; p < fimg->n; p++) {
            fimg->val[p] = (fmax * (img->val[p] - img_min_val) /
                            (img_max_val - img_min_val) * 1.0);
        }
    } else {
        iftWarning("Image is empty", "iftFImageToImage");
    }

    return (fimg);
}


iftFImage *iftFCopyImage(const iftFImage *fimg) {
    iftFImage *copy = NULL;

    if (fimg != NULL) {
        copy = iftCreateFImage(fimg->xsize, fimg->ysize, fimg->zsize);
        iftFCopyVoxelSize(fimg, copy);

        #pragma omp parallel for
        for (int p = 0; p < fimg->n; p++)
            copy->val[p] = fimg->val[p];
    }

    return copy;
}


iftFImage *iftReadFImage(const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);
    
    iftNumPyHeader *header = NULL;
    void *data = iftReadNumPy(npy_path, &header);
    iftCDataType cdtype = iftNumPyDTypeAsCDataType(header->dtype);
    
    // xsize = ncols, ysize = nrows, zsize = nslices
    int xsize = 1, ysize = 1, zsize = 1;
    if (header->n_dims == 1)
        xsize = (int) header->shape[0];
    else if (header->n_dims == 2) {
        ysize = (int) header->shape[0];
        xsize = (int) header->shape[1];
    }
    else if (header->n_dims == 3) {
        zsize = (int) header->shape[0];
        ysize = (int) header->shape[1];
        xsize = (int) header->shape[2];
    }
    else iftError("Numpy array shape is not (1-D or 2-D or 3-D)... it is actually %d-D", "iftReadFImage", header->n_dims);
    
    
    iftFloatArray *arr = NULL;
    if (cdtype == IFT_FLT_TYPE) {
        arr = iftCreateFloatArray(header->size);
        iftFree(arr->val);
        arr->val = data;
    }
    else {
        arr = iftConvertToFloatArray(data, header->size, cdtype);
        iftFree(data);
    }
    
    iftFImage *fimg = iftCreateFImage(xsize, ysize, zsize);
    
    // fortran order: column-major - the matrix is stored inverted as (ncols, nrows)
    if (header->fortran_order) {
        iftFImage *fortran_fimg = iftCreateFImage(zsize, ysize, xsize);
        iftFree(fortran_fimg->val);
        fortran_fimg->val = arr->val;
        arr->val = NULL;
    
        for (int z = 0; z < zsize; z++)
            for (int y = 0; y < ysize; y++)
                for (int x = 0; x < xsize; x++)
                    iftImgVal(fimg, x, y, z) = iftImgVal(fortran_fimg, z, y, x);
    
        iftDestroyFImage(&fortran_fimg);
    }
    else {
        iftFree(fimg->val);
        fimg->val = arr->val;
        arr->val = NULL;
    }
    
    iftDestroyFloatArray(&arr);
    iftDestroyNumPyHeader(&header);
    
    return fimg;
}


void iftWriteFImage(const iftFImage *fimg, const char *format, ...) {
    va_list args;
    char npy_path[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(npy_path, format, args);
    va_end(args);

    char *parent_dir = iftDirname(npy_path);

    if (!iftDirExists(parent_dir)) {
        iftMakeDir(parent_dir);
    }
    iftFree(parent_dir);
    
    long shape[3];
    size_t n_dims;
    if (iftIs3DFImage(fimg)) {
        n_dims = 3;
        shape[0] = fimg->zsize;
        shape[1] = fimg->ysize;
        shape[2] = fimg->xsize;
    }
    else {
        n_dims = 2;
        shape[0] = fimg->ysize;
        shape[1] = fimg->xsize;
    }
    
    iftNumPyHeader *header = iftCreateNumPyHeader(IFT_FLT_TYPE, shape, n_dims);
    iftWriteNumPy(header, fimg->val, npy_path);
    iftDestroyNumPyHeader(&header);
}


iftFImage *iftFGetXYSlice(iftFImage *img, int zcoord) {
    iftFImage *slice;
    iftVoxel  u;
    int       p, q;

    if ((zcoord < 0) || (zcoord >= img->zsize))
        iftError("Invalid z coordinate", "iftFGetXYSlice");

    slice = iftCreateFImage(img->xsize, img->ysize, 1);
    u.z = zcoord;
    q = 0;
    for (u.y = 0; u.y < img->ysize; u.y++)
        for (u.x = 0; u.x < img->xsize; u.x++) {
            p = iftFGetVoxelIndex(img, u);
            slice->val[q] = img->val[p];
            q++;
        }
    iftFCopyVoxelSize(img, slice);

    return (slice);
}

void iftFPutXYSlice(iftFImage *img, iftFImage *slice, int zcoord) {
    iftVoxel u;
    int      p, q;

    if ((zcoord < 0) || (zcoord >= img->zsize))
        iftError("Invalid z coordinate", "iftFPutXYSlice");

    if ((img->ysize != slice->ysize) || (img->xsize != slice->xsize))
        iftError("Image and slice are incompatibles", "iftFPutXYSlice");

    u.z = zcoord;
    p = 0;
    for (u.y = 0; u.y < img->ysize; u.y++)
        for (u.x = 0; u.x < img->xsize; u.x++) {
            q = iftFGetVoxelIndex(img, u);
            img->val[q] = slice->val[p];
            p++;
        }
}

iftFImage *iftFGetZXSlice(iftFImage *img, int ycoord) {
    iftFImage *slice;
    iftVoxel  u;
    int       p, q;

    if ((ycoord < 0) || (ycoord >= img->ysize))
        iftError("Invalid y coordinate", "iftFGetZXSlice");

    slice = iftCreateFImage(img->zsize, img->xsize, 1);
    u.y = ycoord;
    q = 0;
    for (u.x = 0; u.x < img->xsize; u.x++)
        for (u.z = 0; u.z < img->zsize; u.z++) {
            p = iftFGetVoxelIndex(img, u);
            slice->val[q] = img->val[p];
            q++;
        }
    slice->dx = img->dz;
    slice->dy = img->dx;
    slice->dz = img->dy;

    return (slice);
}

void iftFPutZXSlice(iftFImage *img, iftFImage *slice, int ycoord) {
    iftVoxel u;
    int      p, q;

    if ((ycoord < 0) || (ycoord >= img->ysize))
        iftError("Invalid y coordinate", "iftFPutZXSlice");

    if ((img->xsize != slice->ysize) || (img->zsize != slice->xsize))
        iftError("Image and slice are incompatibles", "iftFPutZXSlice");

    u.y = ycoord;
    p = 0;
    for (u.x = 0; u.x < img->xsize; u.x++)
        for (u.z = 0; u.z < img->zsize; u.z++) {
            q = iftFGetVoxelIndex(img, u);
            img->val[q] = slice->val[p];
            p++;
        }
}

iftFImage *iftFGetYZSlice(iftFImage *img, int xcoord) {
    iftFImage *slice;
    iftVoxel  u;
    int       p, q;

    if ((xcoord < 0) || (xcoord >= img->xsize))
        iftError("Invalid x coordinate", "iftGetYZSlice");

    slice = iftCreateFImage(img->ysize, img->zsize, 1);
    u.x = xcoord;
    q = 0;
    for (u.z = 0; u.z < img->zsize; u.z++)
        for (u.y = 0; u.y < img->ysize; u.y++) {
            p = iftFGetVoxelIndex(img, u);
            slice->val[q] = img->val[p];
            q++;
        }
    slice->dx = img->dy;
    slice->dy = img->dz;
    slice->dz = img->dx;

    return (slice);
}

void iftFPutYZSlice(iftFImage *img, iftFImage *slice, int xcoord) {
    iftVoxel u;
    int      p, q;

    if ((xcoord < 0) || (xcoord >= img->xsize))
        iftError("Invalid x coordinate", "iftFPutYZSlice");

    if ((img->zsize != slice->ysize) || (img->ysize != slice->xsize))
        iftError("Image and slice are incompatibles", "iftFPutYZSlice");

    u.x = xcoord;
    p = 0;
    for (u.z = 0; u.z < img->zsize; u.z++)
        for (u.y = 0; u.y < img->ysize; u.y++) {
            q = iftFGetVoxelIndex(img, u);
            img->val[q] = slice->val[p];
            p++;
        }
}


iftFImage *iftFReadRawSlices(char *basename, int first, int last, int xsize, int ysize) {
    FILE      *fp;
    int       p, offset, n = xsize * ysize, i, zsize = last - first + 1;
    iftFImage *img;
    char      filename[200];
    float     *data        = iftAllocFloatArray(n);

    img = iftCreateFImage(xsize, ysize, zsize);
    for (i = first; i <= last; i++) {
        offset = n * (i - first);
        sprintf(filename, "%s%03d.raw", basename, i);
        fp = fopen(filename, "rb");
        if (fread(data, sizeof(float), n, fp) != n) iftError("Reading error", "iftFReadRawSlices");
        for (p = 0; p < n; p++) {
            img->val[p + offset] = data[p];
        }

        fclose(fp);
    }

    iftFree(data);

    return (img);
}

void iftFWriteRawSlices(iftFImage *img, char *basename) {
    FILE  *fp;
    int   p, offset, n = img->xsize * img->ysize, i;
    char  filename[200];
    float *data        = iftAllocFloatArray(n);

    for (i = 0; i < img->zsize; i++) {
        offset = n * i;
        for (p = 0; p < n; p++)
            data[p] = img->val[p + offset];
        sprintf(filename, "%s%03d.raw", basename, i);
        fp = fopen(filename, "wb");
        fwrite(data, n, sizeof(float), fp);
        fclose(fp);
    }

    iftFree(data);

}


float iftFImageValueAtPoint(const iftFImage *img, iftPoint P)
{
    iftVoxel u[8];
    int p[8], i;
    float dx,dy,dz, val[6], value;

    u[0].x = (int)P.x;      u[0].y = (int)P.y;       u[0].z = (int)P.z;
    u[1].x = u[0].x+1;      u[1].y = u[0].y;         u[1].z = u[0].z;
    u[2].x = u[0].x;        u[2].y = u[0].y + 1;     u[2].z = u[0].z;
    u[3].x = u[0].x+1;      u[3].y = u[0].y + 1;     u[3].z = u[0].z;
    u[4].x = u[0].x;        u[4].y = u[0].y;         u[4].z = u[0].z + 1;
    u[5].x = u[0].x+1;      u[5].y = u[0].y;         u[5].z = u[0].z + 1;
    u[6].x = u[0].x;        u[6].y = u[0].y + 1;     u[6].z = u[0].z + 1;
    u[7].x = u[0].x+1;      u[7].y = u[0].y + 1;     u[7].z = u[0].z + 1;

    for (i=0; i < 8; i++) {
        if (iftValidVoxel(img,u[i])){
            p[i] = iftGetVoxelIndex(img,u[i]);
        }else{
            p[0]   = iftGetVoxelIndex(img,u[0]);
            return(img->val[p[0]]);
        }
    }

    dx = dy = dz = 1.0;

    val[0] =(float)img->val[p[1]]*dx+(float)img->val[p[0]]*(1.0-dx);
    val[1] =(float)img->val[p[3]]*dx+(float)img->val[p[2]]*(1.0-dx);
    val[2] =(float)img->val[p[5]]*dx+(float)img->val[p[4]]*(1.0-dx);
    val[3] =(float)img->val[p[7]]*dx+(float)img->val[p[6]]*(1.0-dx);
    val[4] = val[1]*dy + val[0]*(1.0-dy);
    val[5] = val[3]*dy + val[2]*(1.0-dy);
    value  = (val[5]*dz + val[4]*(1.0-dz));

    return(value);
}


float iftFImageValueAtPoint2D(const iftFImage *img, iftPoint P)
{
    iftVoxel u[4];
    int      p[4], i;
    float    dx, dy, val[2], value;

    u[0].x = (int)P.x;         u[0].y = (int)P.y;          u[0].z = 0;
    u[1].x = u[0].x + 1;       u[1].y = u[0].y;            u[1].z = 0;
    u[2].x = u[0].x;           u[2].y = u[0].y + 1;        u[2].z = 0;
    u[3].x = u[0].x + 1;       u[3].y = u[0].y + 1;        u[3].z = 0;

    for (i=0; i < 4; i++) {
        if (iftValidVoxel(img,u[i])){
            p[i] = iftGetVoxelIndex(img,u[i]);
        }else{
            p[0]   = iftGetVoxelIndex(img,u[0]);
            return(img->val[p[0]]);
        }
    }

    dx    = P.x-u[0].x;
    dy    = P.y-u[0].y;

    val[0] =(float)img->val[p[1]]*dx+(float)img->val[p[0]]*(1.0-dx);
    val[1] =(float)img->val[p[3]]*dx+(float)img->val[p[2]]*(1.0-dx);
    value  = (val[1]*dy + val[0]*(1.0-dy));

    return(value);
}


/* Converts CT images in attenuation coefficient unit to Hounsfield
   Unit (HU) */

iftImage *iftAttCoefToHU(iftFImage *attcoef, double mean_of_water) {
    iftImage *hu = iftCreateImage(attcoef->xsize, attcoef->ysize, attcoef->zsize);

    if (fabs(mean_of_water) < IFT_EPSILON)
        iftError("Invalid mean value of water", "iftAttCoefToHU");

    for (int p = 0; p < hu->n; p++)
        hu->val[p] = (int) (1000.0 * (attcoef->val[p] - mean_of_water) / mean_of_water);

    return (hu);
}


iftFImage *iftFExtractROI(const iftFImage *fimg, iftBoundingBox bb) {
    iftVoxel uo = bb.begin;
    iftVoxel uf = bb.end;

    if (!iftValidVoxel(fimg, uo))
        iftError(
                "Initial Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for fimage with size (%d, %d, %d)",
                "iftFExtractROI", uo.x, uo.y, uo.z, fimg->xsize, fimg->ysize, fimg->zsize);
    if (!iftValidVoxel(fimg, uf)) {
        iftWarning("The ROI Image DOES NOT fit entirely inside the FImage.\n" \
                   "It will copied/inserted what it is possible", "iftFExtractROI");
        // gets the valid ending voxel
        uf.x = iftMin(uf.x, fimg->xsize - 1);
        uf.y = iftMin(uf.y, fimg->ysize - 1);
        uf.z = iftMin(uf.z, fimg->zsize - 1);
    }

    iftFImage *froi = iftCreateFImage(uf.x - uo.x + 1, uf.y - uo.y + 1, uf.z - uo.z + 1);
    iftCopyVoxelSize(fimg, froi);

    int      q = 0;
    iftVoxel u;

    for (u.z = uo.z; u.z <= uf.z; u.z++)
        for (u.y = uo.y; u.y <= uf.y; u.y++)
            for (u.x = uo.x; u.x <= uf.x; u.x++) {
                int p = iftGetVoxelIndex(fimg, u);
                froi->val[q] = fimg->val[p];
                q++;
            }

    return froi;
}


// Inserts the roi image into img. Note that only the portion of the
// roi image that lies within img will be copied.
void iftFInsertROI(const iftFImage *roi, iftFImage *target, iftVoxel begin) {
    iftVoxel uo = begin;
    iftVoxel uf;
    uf.x = uo.x + roi->xsize - 1;
    uf.y = uo.y + roi->ysize - 1;
    uf.z = uo.z + roi->zsize - 1;

    if (!iftValidVoxel(target, uo))
        iftError(
                "Initial Point (%d, %d, %d) of the Bound. Box is not a valid pixel/voxel for image with size (%d, %d, %d)",
                "iftFInsertROI", uo.x, uo.y, uo.z, target->xsize, target->ysize, target->zsize);
    if (!iftValidVoxel(target, uf)) {
        iftWarning("The ROI Image DOES NOT fit entirely inside the Target Image.\n" \
                   "It will copied/inserted what it is possible", "iftFInsertROI");
        // gets the valid ending voxel
        uf.x = iftMin(uf.x, target->xsize - 1);
        uf.y = iftMin(uf.y, target->ysize - 1);
        uf.z = iftMin(uf.z, target->zsize - 1);
    }

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
            }
}


void iftFInsertROIByPosition(const iftFImage *roi, const iftVoxel roi_pos, iftFImage *target, const iftVoxel target_pos) {
    iftVoxel uo, uf;

    if (roi == NULL)
        iftError("ROI Image is NULL", "iftFInsertROIByPosition");
    if (target == NULL)
        iftError("Source Image is NULL", "iftFInsertROIByPosition");

    // Computing the first valid position in the target image that intersects the ROI
    uo.x = iftMax(0, target_pos.x - roi_pos.x);
    uo.y = iftMax(0, target_pos.y - roi_pos.y);
    uo.z = iftMax(0, target_pos.z - roi_pos.z);

    // Computing the last valid position in the target image that intersects the ROI
    uf.x = iftMin(target->xsize - 1, uo.x + roi->xsize - 1);
    uf.y = iftMin(target->ysize - 1, uo.y + roi->ysize - 1);
    uf.z = iftMin(target->zsize - 1, uo.z + roi->zsize - 1);

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


iftFImage *iftFAddFrame(const iftFImage *img, int sz, float value) {
    return iftFAddRectangularBoxFrame(img, sz, sz, sz, value);
}

/*! @brief Generalization of function iftFAddFrame for adding frames
  with different sizes for every axis
 */
iftFImage *iftFAddRectangularBoxFrame(const iftFImage *img, int sx, int sy, int sz, float value) {
    iftFImage *fimg;
    int       p, q;
    iftVoxel  u;

    if (iftIs3DFImage(img)) {
        fimg = iftCreateFImage(img->xsize + (2 * sx), img->ysize + (2 * sy), img->zsize + (2 * sz));
        iftFCopyVoxelSize(img, fimg);
        iftFSetImage(fimg, value);

        p = 0;
        for (u.z = sz; u.z < fimg->zsize - sz; u.z++)
            for (u.y = sy; u.y < fimg->ysize - sy; u.y++)
                for (u.x = sx; u.x < fimg->xsize - sx; u.x++) {
                    q = iftFGetVoxelIndex(fimg, u);
                    fimg->val[q] = img->val[p];
                    p++;
                }
    } else {
        fimg = iftCreateFImage(img->xsize + (2 * sx), img->ysize + (2 * sy), 1);
        iftFCopyVoxelSize(img, fimg);
        iftFSetImage(fimg, value);

        p = 0;
        u.z = 0;
        for (u.y = sy; u.y < fimg->ysize - sy; u.y++)
            for (u.x = sx; u.x < fimg->xsize - sx; u.x++) {
                q = iftFGetVoxelIndex(fimg, u);
                fimg->val[q] = img->val[p];
                p++;
            }
    }


    return (fimg);
}

iftFImage *iftFRemFrame(const iftFImage *fimg, int sz) {
    return iftFRemRectangularBoxFrame(fimg, sz, sz, sz);
}

/*! @brief
  Generalization of function iftFRemFrame for removing frames with
  different sizes for every axis
*/
iftFImage *iftFRemRectangularBoxFrame(const iftFImage *fimg, int sx, int sy, int sz) {
    iftFImage *img;
    int       p, q;
    iftVoxel  u;

    if (iftIs3DFImage(fimg)) {
        img = iftCreateFImage(fimg->xsize - (2 * sx), fimg->ysize - (2 * sy),
                              fimg->zsize - (2 * sz));
        iftFCopyVoxelSize(fimg, img);

        p = 0;
        for (u.z = sz; u.z < fimg->zsize - sz; u.z++)
            for (u.y = sy; u.y < fimg->ysize - sy; u.y++)
                for (u.x = sx; u.x < fimg->xsize - sx; u.x++) {
                    q = iftFGetVoxelIndex(fimg, u);
                    img->val[p] = fimg->val[q];
                    p++;
                }
    } else {
        img = iftCreateFImage(fimg->xsize - (2 * sx), fimg->ysize - (2 * sy), 1);
        iftFCopyVoxelSize(fimg, img);

        p = 0;
        u.z = 0;
        for (u.y = sy; u.y < fimg->ysize - sy; u.y++)
            for (u.x = sx; u.x < fimg->xsize - sx; u.x++) {
                q = iftFGetVoxelIndex(fimg, u);
                img->val[p] = fimg->val[q];
                p++;
            }
    }


    return (img);
}


iftBoundingBox iftFMinBoundingBox(const iftFImage *fimg, iftPoint *gc_out) {
    if (fimg == NULL)
        iftError("fImage is NULL", "iftFMinBoundingBox");

    long           n  = 0; // number of spels non-background (non-zero)
    iftPoint       gc = {0.0, 0.0, 0.0};
    iftBoundingBox mbb;
    mbb.begin.x = mbb.begin.y = mbb.begin.z = IFT_INFINITY_INT;
    mbb.end.x   = mbb.end.y   = mbb.end.z   = IFT_INFINITY_INT_NEG;

    for (long p = 0; p < fimg->n; p++) {
        if (fimg->val[p] > 0.0) if (!iftAlmostZero(fimg->val[p])) {
            iftVoxel v = iftFGetVoxelCoord(fimg, p);

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


iftFImage *iftFShiftImage(const iftFImage *fimg, int dx, int dy, int dz) {
    if (fimg == NULL)
        iftError("fImage is NULL", "iftFShiftImage");

    // no shifting
    if ((dx == 0) && (dy == 0) && (dz == 0)) {
        return iftFCopyImage(fimg);
    }

    iftFImage *shift_img = iftCreateFImage(fimg->xsize, fimg->ysize, fimg->zsize);

    #pragma omp parallel for
    for (int p = 0; p < fimg->n; p++) {
        iftVoxel v = iftFGetVoxelCoord(fimg, p);
        v.x += dx;
        v.y += dy;
        v.z += dz;

        if (iftValidVoxel(fimg, v)) {
            int q             = iftFGetVoxelIndex(shift_img, v);
            shift_img->val[q] = fimg->val[p];
        }
    }

    return shift_img;
}


iftFImage *iftFExtractSlice(const iftFImage *vol_fimg, iftImagePlaneOrientation plane_orientation, long slice) {
    // CHECKERS
    if (vol_fimg == NULL)
        iftError("Volumetric FImage is NULL", "iftFExtractSlice");
    if (!iftIs3DFImage(vol_fimg))
        iftError("FImage is not Volumetric (3D)", "iftFExtractSlice");
    if (slice < 0)
        iftError("Invalid Slice %ld (< 0)", "iftFExtractSlice", slice);

    // Gets the Slice
    iftFImage *out_fimg = NULL;
    int q = 0;
    iftVoxel u;
    switch(plane_orientation) {
        case IFT_AXIAL:
            if (slice >= vol_fimg->zsize)
                iftError("Invalid Slice %ld (> the last Axial Slice)", "iftFExtractSlice", slice);

            out_fimg = iftCreateFImage(vol_fimg->xsize, vol_fimg->ysize, 1);

            u.z = slice;
            for (u.y = 0; u.y < vol_fimg->ysize; u.y++)
                for (u.x = 0; u.x < vol_fimg->xsize; u.x++) {
                    int p = iftGetVoxelIndex(vol_fimg, u);
                    out_fimg->val[q++] = vol_fimg->val[p];
                }
            break;
        case IFT_CORONAL:
            if (slice >= vol_fimg->ysize)
                iftError("Invalid Slice %ld (> the last Coronal Slice)", "iftFExtractSlice", slice);

            out_fimg = iftCreateFImage(vol_fimg->xsize, vol_fimg->zsize, 1);

            u.y = slice;
            for (u.z = 0; u.z < vol_fimg->zsize; u.z++)
                for (u.x = 0; u.x < vol_fimg->xsize; u.x++) {
                    int p = iftGetVoxelIndex(vol_fimg, u);
                    out_fimg->val[q++] = vol_fimg->val[p];
                }
           break;
        case IFT_SAGITTAL:
            if (slice >= vol_fimg->xsize)
                iftError("Invalid Slice %ld (> the last Sagittal Slice)", "iftFExtractSlice", slice);

            out_fimg = iftCreateFImage(vol_fimg->ysize, vol_fimg->zsize, 1);

            u.x = slice;
            for (u.z = 0; u.z < vol_fimg->zsize; u.z++)
                for (u.y = 0; u.y < vol_fimg->ysize; u.y++) {
                    int p = iftGetVoxelIndex(vol_fimg, u);
                    out_fimg->val[q++] = vol_fimg->val[p];
                }
            break;
        default:
            iftError("Invalid Image Plane", "iftFExtractSlice");
    }

    iftCopyVoxelSize(vol_fimg, out_fimg);


    return out_fimg;
}

void iftFWriteVTKImage(iftFImage *img, char *filename)      //Write volume as VTK binary float volume (for visualization)
{                                                           //Data is saved as cell-data, for ease of use in ParaView
    FILE       *fp=NULL;

    fp = fopen(filename,"wb");
    if (fp == NULL){
        iftError(MSG_FILE_OPEN_ERROR, filename, "iftFWriteVTKImage");
    }
    float minval = iftFMinimumValue(img);
    float maxval = iftFMaximumValue(img);

    fprintf(fp,"# vtk DataFile Version 3.0\n");
    fprintf(fp,"#generated by iftSkel\n");
    fprintf(fp,"BINARY\n");
    fprintf(fp,"DATASET STRUCTURED_POINTS\n");
    fprintf(fp,"DIMENSIONS %d %d %d\n",img->xsize+1,img->ysize+1,img->zsize+1);
    fprintf(fp,"ORIGIN 0 0 0\n");                           //ALEX: add +1 to point-dimensions since we're saving data as cells
    fprintf(fp,"SPACING 1 1 1\n");                          //ALEX: May need to be img->dx,img->dy,img->dz
    fprintf(fp,"CELL_DATA %d\n",img->xsize*img->ysize*img->zsize);
    fprintf(fp,"SCALARS voxel_data float\n");
    fprintf(fp,"LOOKUP_TABLE default\n");

    for(int i=0;i<img->n;++i)
    {
        float v = FloatSwap(img->val[i]);
        fwrite(&v,1,sizeof(v),fp);
    }

    fclose(fp);

    printf("Written skeleton: range %f %f\n",(float)minval,(float)maxval);
}

void iftFWriteVTKImageVector(iftFImage **img, char *filename)
{
    FILE       *fp=NULL;

    fp = fopen(filename,"wb");
    if (fp == NULL){
        iftError(MSG_FILE_OPEN_ERROR, filename, "iftFWriteVTKImageVector");
    }

    fprintf(fp,"# vtk DataFile Version 3.0\n");
    fprintf(fp,"#generated by iftSkel\n");
    fprintf(fp,"BINARY\n");
    fprintf(fp,"DATASET STRUCTURED_POINTS\n");
    fprintf(fp,"DIMENSIONS %d %d %d\n",img[0]->xsize+1,img[0]->ysize+1,img[0]->zsize+1);
    fprintf(fp,"ORIGIN 0 0 0\n");                           //ALEX: add +1 to point-dimensions since we're saving data as cells
    fprintf(fp,"SPACING 1 1 1\n");                          //ALEX: May need to be img->dx,img->dy,img->dz
    fprintf(fp,"CELL_DATA %d\n",img[0]->n);
    fprintf(fp,"VECTORS voxel_data float\n");
    fprintf(fp,"LOOKUP_TABLE default\n");

    for(int i=0;i<img[0]->n;++i)
    {
      for (int j=0; j < 3; j++) {
        float v = FloatSwap(img[j]->val[i]);
        fwrite(&v,1,sizeof(v),fp);
      }
    }
    fclose(fp);
}


iftFImage *iftFMask(const iftFImage *fimg, const iftImage *mask) {
    if ((fimg->xsize != mask->xsize) || (fimg->ysize != mask->ysize) ||
        (fimg->zsize != mask->zsize))
        iftError("FImage and Mask have different domains\n" \
                "FImage: (%d, %d, %d)\n" \
                "Mask: (%d, %d, %d)", "iftFMask",
                fimg->xsize, fimg->ysize, fimg->zsize,
                mask->xsize, mask->ysize, mask->zsize);

    iftFImage *out = iftCreateFImage(fimg->xsize, fimg->ysize, fimg->zsize);
    for (int p = 0; p < fimg->n; p++)
        if (mask->val[p])
            out->val[p] = fimg->val[p];

    return out;
}

int iftFImageDepth(iftFImage* fimg) {
    int max = iftNormalizationValue(iftFMaximumValue(fimg)) + 1;
    double depth = iftLog(max, 2);
    return depth;
}



iftFImage *iftFlipFImage(const iftFImage *fimg, char axis) {
    iftFImage *flip_fimg = iftCreateFImageFromFImage(fimg);
    iftVoxel u;

    if (axis == IFT_AXIS_X) {
        for (u.z = 0; u.z < fimg->zsize; u.z++)  
            for (u.y = 0; u.y < fimg->ysize; u.y++)  
                for (u.x = 0; u.x < fimg->xsize; u.x++) {
                    iftVoxel v;
                    v.x = fimg->xsize - 1 - u.x;
                    v.y = u.y; v.z = u.z;
                    iftImgVoxelVal(flip_fimg, v) = iftImgVoxelVal(fimg, u);
                }
    }
    else if (axis == IFT_AXIS_Y) {
        for (u.z = 0; u.z < fimg->zsize; u.z++)  
            for (u.y = 0; u.y < fimg->ysize; u.y++)  
                for (u.x = 0; u.x < fimg->xsize; u.x++) {
                    iftVoxel v;
                    v.y = fimg->ysize - 1 - u.y;
                    v.x = u.x; v.z = u.z;
                    iftImgVoxelVal(flip_fimg, v) = iftImgVoxelVal(fimg, u);
                }
    }
    else if (axis == IFT_AXIS_Z) {
        for (u.z = 0; u.z < fimg->zsize; u.z++)  
            for (u.y = 0; u.y < fimg->ysize; u.y++)  
                for (u.x = 0; u.x < fimg->xsize; u.x++) {
                    iftVoxel v;
                    v.z = fimg->zsize - 1 - u.z;
                    v.x = u.x; v.y = u.y;
                    iftImgVoxelVal(flip_fimg, v) = iftImgVoxelVal(fimg, u);
                }
    }
    else iftError("You must inform IFT_AXIS_X, IFT_AXIS_Y, or IFT_AXIS_Z", "iftFlipFImage");
    
    return flip_fimg;
}


void iftFImageSetValue(iftFImage *fimg, const iftSet *S, float value)
{
    for (const iftSet *s = S; s; s = s->next)
        fimg->val[s->elem] = value;
}
