%include "iftAdjacency.i"
%include "iftCommon.i"
%include "iftSort.i"



typedef struct {
    /** Brightness pixels array. */
    int *val;
    /** Blue component pixels array */
    ushort *Cb;
    /** Red component pixels array */
    ushort *Cr;
    /** alpha component pixels array */
    ushort *alpha;

    /** X axis size. */
    int xsize;
    /** Y axis size. */
    int ysize;
    /** Z axis size. */
    int zsize;

    /** X axis voxel size. */
    float dx;
    /** Y axis voxel size. */
    float dy;
    /** Z axis voxel size. */
    float dz;

    /** speed-up voxel access tables */
    int *tby, *tbz;
    /** Number of pixels. */
    int n; // number of voxels
} iftImage;

%extend iftImage {

	%newobject __getitem__;
	PyObject* __getitem__(int pixel)
	{
	    iftImage *img = ($self);
	
	    if(iftIsColorImage(img))
	        return Py_BuildValue("(iii)", img->val[pixel], img->Cb[pixel], img->Cr[pixel]);
	
	    return Py_BuildValue("i", img->val[pixel]);
	}
	
	void __setitem__(int pixel, PyObject* color) {
	    iftImage *img = ($self);
	
	    img->val[pixel] = (int) PyInt_AsLong(PyTuple_GetItem(color, 0));
	    if(iftIsColorImage(img)){
	        img->Cb[pixel] = (ushort) PyInt_AsLong(PyTuple_GetItem(color, 1));
	        img->Cr[pixel] = (ushort) PyInt_AsLong(PyTuple_GetItem(color, 2));
	    }
	}
	
	int GetVoxelIndex(iftVoxel* voxel) {
	    iftImage *s = ($self);
	    iftVoxel v = *voxel;
	    return iftGetVoxelIndex(s,v);
	}
	
	
	int GetCoordIndex(int x, int y, int z)
	{
	    iftImage *s = ($self);
	    iftVoxel v = {.x = x, .y = y, .z = z};
	    return iftGetVoxelIndex(s, v);
	}
	
	
	iftVoxel GetVoxel(int p) {
	    iftVoxel v = iftGetVoxelCoord(($self), p);
	    return v;
	}
	
	PyObject *GetVoxelSizes(void) {
	    iftImage *img = ($self);
	    npy_intp dims[1];
	
	    PyArrayObject *result = NULL;
	    if (iftIs3DImage(img)) {
	        dims[0] = 3;
	        result = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	        float *ptr = (float*) array_data(result);
	
	        ptr[0] = img->dz;
	        ptr[1] = img->dy;
	        ptr[2] = img->dx;
	    }
	    else {
	        dims[0] = 2;
	        result = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
	        float *ptr = (float*) array_data(result);
	
	        ptr[0] = img->dy;
	        ptr[1] = img->dx;
	    }
	
	    return PyArray_Return(result);
	}
	
	void SetVoxelSizes(PyObject *input) {
	    if (PyArray_TYPE(input) != NPY_FLOAT32) {
	        SWIG_Error(0, "Numpy array of voxel sizes should be float32.");
	        return;
	    }
	
	    PyArrayObject *npy = obj_to_array_no_conversion(input, NPY_FLOAT32);
	
	    iftImage *img = ($self);
	
	    int ndim = PyArray_NDIM(npy);
	    if (ndim != 1) {
	        SWIG_Error(0, "Numpy array is not unidimensional.");
	        return;
	    }
	
	    int *dims = PyArray_SHAPE(npy);
	    float *ptr = (float*) array_data(npy);
	
	    if ((dims[0] != 2) && (dims[0] != 3)) {
	        SWIG_Error(0, "Invalid numpy array of voxel sizes. Try an array with 2 or 3 dimensions.");
	        return;
	    }
	
	    if (iftIs3DImage(img)) {
	        if (dims[0] != 3) {
	            SWIG_Error(0, "Image is 3D but the array of voxel sizes does not have 3 elements");
	            return;
	        }
	        img->dx = ptr[2];
	        img->dy = ptr[1];
	        img->dz = ptr[0];
	    }
	    else {
	        if (dims[0] != 2) {
	            SWIG_Error(0, "Image is 2D but the array of voxel sizes does not have 2 elements");
	            return;
	        }
	        img->dx = ptr[1];
	        img->dy = ptr[0];
	        img->dz = 0.0;
	    }
	}
	
	
	PyObject* AsNumPy(void) {
	    iftImage *img = ($self);
	    npy_intp *dims = NULL;
	    int n_dims;
	    PyArrayObject *result = NULL;
	
	    if(iftIsColorImage(img)){        
	        if (iftIs3DImage(img)) {
	            // (zsize, ysize, xsize, channels)
	            n_dims = 4;
	            dims = calloc(n_dims, sizeof(npy_intp));
	            dims[0] = img->zsize;
	            dims[1] = img->ysize;
	            dims[2] = img->xsize;
	            dims[3] = 3;
	        }
	        else {
	            // (ysize, xsize, channels)
	            n_dims = 3;
	            dims = calloc(n_dims, sizeof(npy_intp));
	            dims[0] = img->ysize;
	            dims[1] = img->xsize;
	            dims[2] = 3;
	        }
	        result = PyArray_SimpleNew(n_dims, dims, NPY_INT32);
	
	        int *ptr = (int*) array_data(result);
	
	        for(int p = 0; p < img->n; p++){
	            ptr[3 * p]       = (int) img->val[p];
	            ptr[(3 * p) + 1] = (int) img->Cb[p];
	            ptr[(3 * p) + 2] = (int) img->Cr[p];
	        }
	    } else {
	        if (iftIs3DImage(img)) {
	            // (zsize, ysize, xsize)
	            n_dims = 3;
	            dims = calloc(n_dims, sizeof(npy_intp));
	            dims[0] = img->zsize;
	            dims[1] = img->ysize;
	            dims[2] = img->xsize;
	        }
	        else {
	            // (ysize, xsize)
	            n_dims = 2;
	            dims = calloc(n_dims, sizeof(npy_intp));
	            dims[0] = img->ysize;
	            dims[1] = img->xsize;
	        }
	
	        result = PyArray_SimpleNew(n_dims, dims, NPY_INT32);
	        int *ptr = (int*) array_data(result);
	
	        for(int i = 0; i < img->n; i++)
	            ptr[i] = (int) img->val[i];
	    }
	
	    free(dims);
	    return PyArray_Return(result);
	}
	
	void FromNumPy(PyObject* input, bool is3D = false) {
	
	    PyArrayObject* ary = obj_to_array_no_conversion(input, NPY_INT32);
	
	    if (ary == NULL) return;
	
	    if (!require_contiguous(ary))
	        SWIG_Error(0, "Input numpy array is not contiguous");
	
	    iftImage* img = ($self);
	    int* ptr = array_data(ary);
	
	    if (is3D) {
	
	        int n_dims = array_numdims(ary); 
	        int ary_xsize = array_size(ary, 2);
	        int ary_ysize = array_size(ary, 1);
	        int ary_zsize = array_size(ary, 0);    
	        int ary_n = ary_xsize * ary_ysize * ary_zsize;
	
	        if (img->n != ary_n) {
	            // if image is of different size, realloc
	            img->val = iftRealloc(img->val, ary_n * sizeof (*img->val));
	            
	            if (!img->Cb) {
	                img->Cb = iftRealloc(img->Cb, ary_n * sizeof (*img->Cb));
	                if (!img->Cb) iftError(MSG_MEMORY_ALLOC_ERROR, "FromNumPy");
	            }
	            if (!img->Cr) {
	                img->Cr = iftRealloc(img->Cr, ary_n * sizeof (*img->Cr));
	                if (!img->Cr) iftError(MSG_MEMORY_ALLOC_ERROR, "FromNumPy");
	            }
	
	            img->alpha = NULL;
	            img->xsize = ary_xsize;
	            img->ysize = ary_ysize;
	            img->zsize = ary_zsize;
	            img->dx    = 1.0;
	            img->dy    = 1.0;
	            img->dz    = 1.0;
	            img->tby   = iftAllocIntArray(ary_ysize);
	            img->tbz   = iftAllocIntArray(ary_zsize);
	            img->n     = ary_n;
	
	            if (img->val == NULL || img->tbz == NULL || img->tby == NULL) {
	                iftError(MSG_MEMORY_ALLOC_ERROR, "FromNumPy");
	            }
	
	            img->tby[0] = 0;
	            for (int y = 1; y < ary_ysize; y++)
	                img->tby[y] = img->tby[y - 1] + ary_xsize;
	
	            img->tbz[0] = 0;
	            int ary_xysize = ary_xsize * ary_ysize;
	            for (int z = 1; z < ary_zsize; z++)
	                img->tbz[z] = img->tbz[z - 1] + ary_xysize;           
	        }
	
	        if (n_dims == 4 && array_size(ary, 3) == 3) {
	            // color 3d image
	            if (!img->Cb) img->Cb = iftAlloc(ary_n, sizeof (*img->Cb));
	            if (!img->Cr) img->Cr = iftAlloc(ary_n, sizeof (*img->Cr));
	
	            int ary_2n = ary_n * 2;
	            for (int i = 0; i < img->n; i++) {
	                img->val[i] = (int) ptr[i];
	                img->Cb[i] = (ushort) ptr[i + ary_n];
	                img->Cr[i] = (ushort) ptr[i + ary_2n];
	            }
	
	        } else if (n_dims == 3) {
	            // grayscale 3d image
	            for (int i = 0; i < img->n; i++) {
	                img->val[i] = (int) ptr[i];
	            }
	
	        } else {
	            SWIG_Error(0, "Input does not respect IFT standard for 3D image");
	        }
	
	    } else {
	        // is 2d image
	
	        int n_dims = array_numdims(ary); 
	        int ary_xsize = array_size(ary, 1);
	        int ary_ysize = array_size(ary, 0);    
	        int ary_n = ary_xsize * ary_ysize;
	
	        if (img->n != ary_n) {
	            // if image is of different size, realloc
	            img->val = iftRealloc(img->val, ary_n * sizeof (*img->val));
	            
	            if (!img->Cb) {
	                img->Cb = iftRealloc(img->Cb, ary_n * sizeof (*img->Cb));
	                if (!img->Cb) iftError(MSG_MEMORY_ALLOC_ERROR, "FromNumPy");
	            }
	            if (!img->Cr) {
	                img->Cr = iftRealloc(img->Cr, ary_n * sizeof (*img->Cr));
	                if (!img->Cr) iftError(MSG_MEMORY_ALLOC_ERROR, "FromNumPy");
	            }
	
	            img->alpha = NULL;
	            img->xsize = ary_xsize;
	            img->ysize = ary_ysize;
	            img->zsize = 1;
	            img->dx    = 1.0;
	            img->dy    = 1.0;
	            img->dz    = 1.0;
	            img->tby   = iftAllocIntArray(ary_ysize);
	            img->tbz   = iftAllocIntArray(1);
	            img->n     = ary_n;
	
	            if (img->val == NULL || img->tbz == NULL || img->tby == NULL) {
	                iftError(MSG_MEMORY_ALLOC_ERROR, "FromNumPy");
	            }
	
	            img->tby[0] = 0;
	            for (int y = 1; y < ary_ysize; y++)
	                img->tby[y] = img->tby[y - 1] + ary_xsize;
	
	            img->tbz[0] = 0;
	        }
	
	        if (n_dims == 3 && array_size(ary, 2) == 3) {
	            // color 2d image
	            if (!img->Cb) img->Cb = iftAlloc(ary_n, sizeof (*img->Cb));
	            if (!img->Cr) img->Cr = iftAlloc(ary_n, sizeof (*img->Cr));
	
	            int ary_2n = ary_n * 2;
	            for (int i = 0; i < img->n; i++) {
	                img->val[i] = (int) ptr[i];
	                img->Cb[i] = (int) ptr[i + ary_n];
	                img->Cr[i] = (int) ptr[i + ary_2n];
	            }
	
	        } else if (n_dims == 2) {
	            // grayscale 2d image
	            for (int i = 0; i < img->n; i++) {
	                img->val[i] = (int) ptr[i];
	            }
	
	        } else {
	            SWIG_Error(0, "Input does not respect IFT standard for 2D image");
	        }
	
	    }
	}
	
	PyObject* ToPlot(void) {
	    iftImage* img = ($self);
	
	    if(iftIsColorImage(img))
	    {
	        int normalization_value = iftNormalizationValue(iftMaximumValue(img));
	        iftColor RGB;
	        npy_intp dims[3]= {img->ysize, img->xsize, 3};
	        PyArrayObject* result = PyArray_SimpleNew(3, dims, NPY_FLOAT32);
	        float* p = (float*) array_data(result);
	
	        for(int i = 0, j = 0; i < img->n; i++, j += 3){
	            RGB = iftGetRGB(img, i, normalization_value);
	            p[j]     = (float) RGB.val[0] / normalization_value;
	            p[j + 1] = (float) RGB.val[1] / normalization_value;
	            p[j + 2] = (float) RGB.val[2] / normalization_value;
	        }
	
	        return PyArray_Return(result);
	
	    } else {
	        int max = iftMaximumValue(img);
	        npy_intp dims[3] = {img->ysize, img->xsize, 3};
	        PyArrayObject* result = PyArray_SimpleNew(3, dims, NPY_FLOAT32);
	        float* p = (float*) array_data(result);
	
	        for(int i = 0, j = 0; i < img->n; i++, j += 3){
	            float v = (float) img->val[i] / max;
	            p[j] = v;
	            p[j + 1] = v;
	            p[j + 2] = v;
	        }
	
	        return PyArray_Return(result);
	    }
	}
	
	void Write(char* filename){
	    iftImage* img = ($self);
	    iftWriteImageByExt(img, filename);
	}
	
	
	int ImgVal(long x, long y, long z) {
	    iftImage *img = ($self);
	    return iftImgVal(img, x, y, z);
	}
	
	
	int ImgVoxelVal(iftVoxel v) {
	    iftImage *img = ($self);
	    return iftImgVoxelVal(img, v);
	}
	
	
	void SetList(PyObject *list, int value)
	{
	    iftImage *img = ($self);
	    for (int i = 0; i < PyList_Size(list); i++) {
	        long idx = PyLong_AsLong(PyList_GetItem(list, i));
	        img->val[idx] = value;
	    }
	}
	
	
	void SetValue(int value)
	{
	    iftImage *img = ($self);
	    for (int i = 0; i < img->n; i++)
	        img->val[i] = value;
	}
	
	
	long Sum()
	{
	    iftImage *img = ($self);
	    long sum = 0;
	    for (int i = 0; i < img->n; i++)
	        sum += img->val[i];
	    return sum;
	}
	

	~iftImage() {
		iftImage* ptr = ($self);
		iftDestroyImage(&ptr);
	}
};

%newobject iftGetImageDomain;
%feature("autodoc", "2");
iftImageDomain iftGetImageDomain(const iftImage *img);

%feature("autodoc", "2");
static inline bool iftIsColorImage(const iftImage *img) ;

%feature("autodoc", "2");
static inline bool iftIs3DImage(const iftImage *img) ;

%newobject iftCreateImage;
%feature("autodoc", "2");
iftImage *iftCreateImage(int xsize, int ysize, int zsize);

%newobject iftCreateImageFromImage;
%feature("autodoc", "2");
iftImage *iftCreateImageFromImage(const iftImage *src);

%newobject iftReadGrayImageSetAsIntMatrix;
%feature("autodoc", "2");
iftIntMatrix *iftReadGrayImageSetAsIntMatrix(const iftFileSet *img_set);

%feature("autodoc", "2");
void iftSetCbCr(iftImage *img, ushort value);

%feature("autodoc", "2");
int iftMaximumValueInRegion(const iftImage *img, iftBoundingBox bb);

%feature("autodoc", "2");
int iftMaximumValueInMask(const iftImage *img, const iftImage *mask);

%feature("autodoc", "2");
int iftMinimumValueInMask(const iftImage *img, iftImage *mask);

%feature("autodoc", "2");
void iftMinimumValuesForLabels(const iftImage *img, const iftImage *label_img,
                               iftIntArray **min_vals_out, iftIntArray **max_vals_out);

%feature("autodoc", "2");
int iftMaximumValue(const iftImage *img);

%newobject iftReadImageByExt;
%feature("autodoc", "2");
iftImage *iftReadImageByExt(const char *filename, ...);

%feature("autodoc", "2");
void iftWriteImageByExt(const iftImage *img, const char *filename, ...);

%newobject iftCopyImage;
%feature("autodoc", "2");
iftImage *iftCopyImage(const iftImage *img);

%feature("autodoc", "2");
void iftCopyImageInsideRegion(const iftImage *src, const iftImage *mask, iftImage *dst);

%newobject iftImageGradientMagnitude;
%feature("autodoc", "2");
iftImage *iftImageGradientMagnitude(const iftImage *img, iftAdjRel *Ain);

%feature("autodoc", "2");
void iftSetImage(iftImage *img, int value);

%newobject iftGetXYSlice;
%feature("autodoc", "2");
iftImage *iftGetXYSlice(const iftImage *img, int zcoord);

%newobject iftGetZXSlice;
%feature("autodoc", "2");
iftImage *iftGetZXSlice(const iftImage *img, int ycoord);

%newobject iftGetYZSlice;
%feature("autodoc", "2");
iftImage *iftGetYZSlice(const iftImage *img, int xcoord);

%newobject iftPutXYSlice;
%feature("autodoc", "2");
void iftPutXYSlice(iftImage *img, const iftImage *slice, int zcoord);

%newobject iftImageCb;
%feature("autodoc", "2");
iftImage *iftImageCb(const iftImage *img);

%newobject iftImageCr;
%feature("autodoc", "2");
iftImage *iftImageCr(const iftImage *img);

%newobject iftImageGray;
%feature("autodoc", "2");
iftImage *iftImageGray(iftImage *img);

%newobject iftRegionBorders;
%feature("autodoc", "2");
iftImage *iftRegionBorders(const iftImage *label_img);

%newobject iftNumberOfElements;
%feature("autodoc", "2");
int iftNumberOfElements(const iftImage *mask);

%newobject iftSelectImageDomain;
%feature("autodoc", "2");
iftImage *iftSelectImageDomain(int xsize, int ysize, int zsize);

%feature("autodoc", "2");
void iftFillBoundingBoxInImage(iftImage *img, iftBoundingBox bb, int value);

%newobject iftMinBoundingBox;
%feature("autodoc", "2");
iftBoundingBox iftMinBoundingBox(const iftImage *img, iftVoxel *gc_out);

%newobject iftMinObjectBoundingBox;
%feature("autodoc", "2");
iftBoundingBox iftMinObjectBoundingBox(const iftImage *img, int obj_label, iftVoxel *gc_out);

%newobject iftExtractROI;
%feature("autodoc", "2");
iftImage *iftExtractROI(const iftImage *img, iftBoundingBox bb);

%newobject iftInsertROI;
%feature("autodoc", "2");
void iftInsertROI(const iftImage *roi, iftImage *target, iftVoxel begin);

%newobject iftExtractObject;
%feature("autodoc", "2");
iftImage *iftExtractObject(const iftImage *src, int obj_label);

%newobject iftExtractLabels;
%feature("autodoc", "2");
iftImage *iftExtractLabels(const iftImage *src_img, const iftIntArray *labels);

%newobject iftRelabelImage;
%feature("autodoc", "2");
iftImage *iftRelabelImage(const iftImage *label_img);

%newobject iftGetObjectLabels;
%feature("autodoc", "2");
iftIntArray *iftGetObjectLabels(const iftImage *label_img);

%newobject iftFindObjectLabels;
%feature("autodoc", "2");
iftIntArray *iftFindObjectLabels(const iftImage *label_img);

%newobject iftGetObjectVoxels;
%feature("autodoc", "2");
iftIntArray *iftGetObjectVoxels(const iftImage *label_img);

%newobject iftImageBasins;
%feature("autodoc", "2");
iftImage *iftImageBasins(const iftImage *img, iftAdjRel *Ain);

%feature("autodoc", "2");
void iftGetLabeledPathSides(const iftImage *label, const iftSet *path, iftSet **obj, iftSet **bkg);

%newobject iftGrayImageToColorImage;
%feature("autodoc", "2");
iftImage *iftGrayImageToColorImage(const iftImage *img, iftColorTable *ctb);

%inline %{

iftImage *CreateImageFromNumPy(PyObject *input, bool is3D) {
    PyArrayObject *ary = obj_to_array_no_conversion(input, NPY_INT32);

    if (ary == NULL)
        return NULL;

    if (!require_contiguous(ary))
        SWIG_Error(0, "Input numpy array is not contiguous");

    // supported shapes
    // 2D gray image: (ysize, xsize)
    // 2D color image: (ysize, xsize, 3)
    // 3D gray image: (zsize, ysize, xsize)
    int n_dims = array_numdims(ary);


    if ((n_dims != 2) && (n_dims != 3) && (n_dims != 4)) {
        char error[512];
        sprintf(error, "Image Numpy Array should have 2 or 3 or 4 dimensions. Found: %d", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    // 3D gray image is actually a 2D gray image
    if ((is3D) && (n_dims == 2))
        is3D = false;
    
    bool is_color_image = false;

    int xsize, ysize, zsize;

    if (is3D) {
        xsize = array_size(ary, 2);
        ysize = array_size(ary, 1);
        zsize = array_size(ary, 0);    

        // color image
        is_color_image = (n_dims == 4);
    }
    else {
        xsize = array_size(ary, 1);
        ysize = array_size(ary, 0);
        zsize = 1;

        // color image
        is_color_image = (n_dims == 3) && (array_size(ary, 2) == 3);
    }


    int *ptr = array_data(ary);
    iftImage *img = NULL;

    if (is_color_image) {
        img = iftCreateColorImage(xsize, ysize, zsize, 8);

        for (int p = 0; p < img->n; p++) {
            img->val[p] = ptr[3 * p];
            img->Cb[p] = ptr[(3 * p) + 1];
            img->Cr[p] = ptr[(3 * p) + 2];
        }
    }
    else {
        img = iftCreateImage(xsize, ysize, zsize);
        for (int p = 0; p < img->n; p++)
            img->val[p] = ptr[p];
    }

    return img;
}

iftImage *NDArrayToLabelImage(PyObject *input)
{
    int is_new = 0;
    PyArrayObject *ary = obj_to_array_allow_conversion(input, NPY_DOUBLE, &is_new);

    if (ary == NULL)
        return NULL;

    if (!require_contiguous(ary))
        SWIG_Error(0, "Input numpy array is not contiguous");

    int n_dims = array_numdims(ary);

    if (n_dims != 3) {
        char error[256];
        sprintf(error, "Image Numpy Array must have 3 dimensions, %d found.", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    int n_labels = array_size(ary, 0);
    int xsize = array_size(ary, 2);
    int ysize = array_size(ary, 1);
    double *ptr = array_data(ary);
    iftImage *img = iftCreateImage(xsize, ysize, 1);

    for (int p = 0; p < img->n; p++) {
        double max_val = -1.0;
        int max_lb = -1;
        for (int b = 0; b < n_labels; b++)
        {
            int idx = img->n * b + p;
            if (ptr[idx] > max_val) {
                max_lb = b;
                max_val = ptr[idx];
            }
        }
        img->val[p] = max_lb;
    }

    return img;
}

iftImage *ImageFromNumPyLabelIndex(PyObject *index, PyObject *label,
                              int xsize, int ysize, int zsize, bool is_whole_array)
{
    int is_new = 0;
    PyArrayObject *idx = obj_to_array_allow_conversion(index, NPY_INT64, &is_new);
    PyArrayObject *lbl = obj_to_array_allow_conversion(label, NPY_INT64, &is_new);

    if (!idx || !lbl)
        return NULL;

    if (!require_contiguous(lbl) || !require_contiguous(idx))
        SWIG_Error(0, "Inputs must be contiguous array");

    if (array_numdims(idx) != 1 || array_numdims(lbl) != 1)
        SWIG_Error(0, "Inputs must be one dimensional arrays");

    int size = xsize * ysize * zsize;
    if (is_whole_array && (size != array_size(idx, 0) || size != array_size(lbl, 0)))
        SWIG_Error(0, "Array size must be equal to image frame");

    iftImage *img = iftCreateImage(xsize, ysize, zsize);

    long *d_idx = array_data(idx);
    long *d_lbl = array_data(lbl);

    for (int i = 0; i < img->n; i++) {
        img->val[ d_idx[i] ] = d_lbl[i];
    }

    return img;
}


iftImage *ImageFromNumPyMapping(const iftImage *components, PyObject *reference, PyObject *label)
{
    int is_new = 0;
    PyArrayObject *ref = obj_to_array_allow_conversion(reference, NPY_INT64, &is_new);
    PyArrayObject *lbl = obj_to_array_allow_conversion(label, NPY_INT64, &is_new);

    if (!lbl || !ref)
        return NULL;

    if (!require_contiguous(lbl) || !require_contiguous(ref))
        SWIG_Error(0, "Inputs must be contiguous array");

    if (array_numdims(ref) != 1 || array_numdims(lbl) != 1)
        SWIG_Error(0, "Inputs must be one dimensional arrays");

    iftImage *img = iftCreateImage(components->xsize, components->ysize, components->zsize);

    long *d_ref = array_data(ref);
    long *d_lbl = array_data(lbl);

    for (int i = 0; i < img->n; i++) {
        int idx = d_ref[ components->val[i] - 1 ];
        img->val[i] = d_lbl[idx];
    }

    return img;
}



PyObject *iftPathValue(PyObject *in_image, PyObject *in_markers, float penalization)
{
    int is_new_image = 0, is_new_markers = 0;
    PyArrayObject *image = obj_to_array_allow_conversion(in_image, NPY_FLOAT, &is_new_image);
    if (!image || !require_contiguous(image))
        Py_RETURN_NONE;

    PyArrayObject *markers = obj_to_array_allow_conversion(in_markers, NPY_INT32, &is_new_markers);
    if (!markers || !require_contiguous(markers))
        Py_RETURN_NONE;

    if (!require_dimensions(image, 3) || !require_dimensions(markers, 2))
        Py_RETURN_NONE;

    npy_intp *image_shape = PyArray_SHAPE(image);
    if (!require_size(markers, image_shape + 1, 2))
        Py_RETURN_NONE;

    PyArrayObject *dist = PyArray_SimpleNew(2, image_shape + 1, NPY_FLOAT);
    if (!dist || !require_contiguous(dist))
        Py_RETURN_NONE;
    int height = image_shape[1], width = image_shape[2], size = width * height;
    float *dist_ptr = array_data(dist);
    iftFHeap *Q = iftCreateFHeap(size, dist_ptr);

    int *markers_ptr = array_data(markers);
    for (int i = 0; i < size; i++) {
        if (markers_ptr[i]) {
            dist_ptr[i] = 0;
            iftInsertFHeap(Q, i);
        } else {
            dist_ptr[i] = IFT_INFINITY_FLT;
        }
    }

    iftAdjRel *A = iftCircular(1.0f);

    float *image_ptr = array_data(image);
    while (!iftEmptyFHeap(Q))
    {
        int p = iftRemoveFHeap(Q);
        div_t res = div(p, width);
        iftVoxel u = {.x = res.rem, .y = res.quot, .z = 0};

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (v.x < 0 || v.x >= width || v.y < 0 || v.y >= height)
                continue;
            int q = v.y * width + v.x;
            if (dist_ptr[q] < dist_ptr[p])
                continue;
            float norm = 0;
            for (int j = 0, row = 0; j < image_shape[0]; j++, row += size) {
                float dif = image_ptr[row + p] - image_ptr[row + q];
                norm += dif * dif;
            }
            norm = sqrtf(norm);
            norm = iftMax(norm, dist_ptr[p] + penalization);
            if (norm < dist_ptr[q]) {
                dist_ptr[q] = norm;
                if (Q->color[q] == IFT_GRAY)
                    iftGoUpFHeap(Q, Q->pos[q]);
                else
                    iftInsertFHeap(Q, q);
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftDestroyFHeap(&Q);
    if (is_new_image)
        Py_DECREF(image);
    if (is_new_markers)
        Py_DECREF(markers);

    return PyArray_Return(dist);
}

static inline float squared(float x) {
    return x * x;
}

PyObject *iftGaussian(PyObject *in_image, PyObject *in_markers, float sigma)
{
    int is_new_image = 0, is_new_markers = 0;
    PyArrayObject *image = obj_to_array_allow_conversion(in_image, NPY_FLOAT, &is_new_image);
    if (!image || !require_contiguous(image))
        Py_RETURN_NONE;

    PyArrayObject *markers = obj_to_array_allow_conversion(in_markers, NPY_INT32, &is_new_markers);
    if (!markers || !require_contiguous(markers))
        Py_RETURN_NONE;

    if (!require_dimensions(image, 3) || !require_dimensions(markers, 2))
        Py_RETURN_NONE;

    npy_intp *image_shape = PyArray_SHAPE(image);
    if (!require_size(markers, image_shape + 1, 2))
        Py_RETURN_NONE;

    PyArrayObject *dist = PyArray_SimpleNew(2, image_shape + 1, NPY_FLOAT);
    if (!dist || !require_contiguous(dist))
        Py_RETURN_NONE;
    int height = image_shape[1], width = image_shape[2], size = width * height;
    float *dist_ptr = array_data(dist);
    iftFHeap *Q = iftCreateFHeap(size, dist_ptr);
    iftSetRemovalPolicyFHeap(Q, IFT_MAXVALUE);

    int *markers_ptr = array_data(markers);
    iftVoxel *root = calloc(size, sizeof (*root));
    for (int i = 0; i < size; i++) {
        if (markers_ptr[i]) {
            dist_ptr[i] = 1.0f;
            div_t res = div(i, width);
            root[i].x = res.rem;
            root[i].y = res.quot;
            iftInsertFHeap(Q, i);
        } else {
            dist_ptr[i] = 0.0f;
        }
    }

    iftAdjRel *A = iftCircular(1.0f);

    float *image_ptr = array_data(image);
    while (!iftEmptyFHeap(Q))
    {
        int p = iftRemoveFHeap(Q);
        div_t res = div(p, width);
        iftVoxel u = {.x = res.rem, .y = res.quot, .z = 0};

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (v.x < 0 || v.x >= width || v.y < 0 || v.y >= height)
                continue;
            int q = v.y * width + v.x;
            if (dist_ptr[q] > dist_ptr[p])
                continue;
            float gauss_dist = expf( - (squared(v.x-root[p].x) + squared(v.y-root[p].y)) / sigma);
            if (gauss_dist > dist_ptr[q]) {
                dist_ptr[q] = gauss_dist;
                root[q].x = root[p].x; root[q].y = root[p].y;
                if (Q->color[q] == IFT_GRAY)
                    iftGoUpFHeap(Q, Q->pos[q]);
                else
                    iftInsertFHeap(Q, q);
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftDestroyFHeap(&Q);
    if (is_new_image)
        Py_DECREF(image);
    if (is_new_markers)
        Py_DECREF(markers);
    free(root);

    return PyArray_Return(dist);
}



%}

%newobject CreateImageFromNumPy;
%feature("autodoc", "2");
iftImage *CreateImageFromNumPy(PyObject *input, bool is3D) ;

%newobject NDArrayToLabelImage;
%feature("autodoc", "2");
iftImage *NDArrayToLabelImage(PyObject *input)
;

%newobject ImageFromNumPyLabelIndex;
%feature("autodoc", "2");
iftImage *ImageFromNumPyLabelIndex(PyObject *index, PyObject *label,
                              int xsize, int ysize, int zsize, bool is_whole_array)
;

%newobject ImageFromNumPyMapping;
%feature("autodoc", "2");
iftImage *ImageFromNumPyMapping(const iftImage *components, PyObject *reference, PyObject *label)
;

