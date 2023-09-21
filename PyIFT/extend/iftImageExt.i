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
