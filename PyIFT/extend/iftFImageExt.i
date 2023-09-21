%newobject __getitem__;
PyObject* __getitem__(int pixel) {
    PyObject *color = PyFloat_FromDouble(($self)->val[pixel]);
    return color;
}

void __setitem__(int pixel, PyObject* color){
    ($self)->val[pixel] = (float) PyFloat_AsDouble(color);
}

void Add(float value) {
    iftFImage *fimg = ($self);
    for (int i = 0; i < fimg->n; i++)
        fimg->val[i] += value;
}

void Sub(float value) {
    iftFImage *fimg = ($self);
    for (int i = 0; i < fimg->n; i++)
        fimg->val[i] -= value;
}

void Mult(float value) {
    iftFImage *fimg = ($self);
    for (int i = 0; i < fimg->n; i++)
        fimg->val[i] *= value;
}

void Div(float value) {
    iftFImage *fimg = ($self);
    for (int i = 0; i < fimg->n; i++)
        fimg->val[i] /= value;
}

int GetVoxelIndex(iftVoxel* voxel) {
    iftImage *s = ($self);
    iftVoxel v = *voxel;
    return iftGetVoxelIndex(s,v);
}

iftVoxel GetVoxel(int p) {
    iftVoxel v = iftGetVoxelCoord(($self),p);
    return v;
}

PyObject* AsNumPy() {
    iftFImage* img = ($self);
    PyArrayObject* npy = NULL, *ref = NULL;

    if (iftIs3DFImage(img)) {
        
        npy_intp dims[3] = {img->zsize, img->ysize, img->xsize};
        ref = PyArray_SimpleNewFromData(3, dims, NPY_FLOAT32, (float*) img->val);
        npy = PyArray_SimpleNew(3, dims, NPY_FLOAT32);
        PyArray_CopyInto(npy, ref);
    } else {

        npy_intp dims[2] = {img->ysize, img->xsize};
        ref = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, (float*) img->val);
        npy = PyArray_SimpleNew(2, dims, NPY_FLOAT32);
        PyArray_CopyInto(npy, ref);
    }

    return PyArray_Return(npy);
}

void FromNumPy(PyObject* input) {

    PyArrayObject* ary = obj_to_array_no_conversion(input, NPY_FLOAT32);

    if (ary == NULL) return;

    if (!require_contiguous(ary))
        SWIG_Error(0, "Input numpy array is not contiguous");

    iftFImage* img = ($self);

    int n_dims = array_numdims(ary); 

    if (n_dims == 3) {
        
        int ary_xsize = array_size(ary, 2);
        int ary_ysize = array_size(ary, 1);
        int ary_zsize = array_size(ary, 0);    
        int ary_n = ary_xsize * ary_ysize * ary_zsize;

        if (img->n != ary_n) {
            // if image is of different size, realloc
            img->val = iftRealloc(img->val, ary_n * sizeof (*img->val));
            
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

    } else if (n_dims == 2) {
        
        int ary_xsize = array_size(ary, 1);
        int ary_ysize = array_size(ary, 0);    
        int ary_n = ary_xsize * ary_ysize;

        if (img->n != ary_n) {

            // if image is of different size, realloc
            img->val = iftRealloc(img->val, ary_n * sizeof (*img->val));

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

    } else {
        img->n = 0;
        SWIG_Error(0, "Input does not respect IFT standard for float image");
    }

    float* ptr = (float*) array_data(ary);
    memcpy(img->val, ptr, sizeof (*img->val) * img->n);
}

PyObject* ToPlot()
{
    iftFImage *img = ($self);
    float normalization_value = iftFMaximumValue(img);

    npy_intp dims[3] = {img->ysize, img->xsize, 3};
    PyArrayObject* result = PyArray_SimpleNew(3, dims, NPY_FLOAT32);
    float* p = (float*) array_data(result);

    for(int i = 0, j = 0; i < img->n; i++, j += 3){
        p[j] = (float) img->val[i] / normalization_value;
        p[j + 1] = (float) img->val[i] /  normalization_value;
        p[j + 2] = (float) img->val[i] / normalization_value;
    }

    return PyArray_Return(result);
}

void Write(char *filename) {
    iftFImage *img = ($self);
    iftWriteImageByExt(iftFImageToImage(img, 255), filename);
}

void SetMask(const iftImage *mask, float value)
{
    iftFImage *img = ($self);
    if (mask->xsize != img->xsize ||
        mask->ysize != img->ysize ||
        mask->zsize != img->zsize)
        SWIG_Error(0, "Float image and mask have different domains");

    for(int i = 0; i < img->n; i++) {
        if (mask->val[i])
            img->val[i] = value;
    }
}

void SetInfinity(const iftImage *mask)
{
    iftFImage *img = ($self);
    if (mask != NULL) {
        if (mask->xsize != img->xsize ||
            mask->ysize != img->ysize ||
            mask->zsize != img->zsize)
            SWIG_Error(0, "Float image and mask have different domains");

        for (int i = 0; i < img->n; i++) {
            if (mask->val[i])
                img->val[i] = IFT_INFINITY_FLT;
        }
    } else {
        for (int i = 0; i < img->n; i++) {
            img->val[i] = IFT_INFINITY_FLT;
        }
    }
}

void SetInfinityNeg(const iftImage *mask)
{
    iftFImage *img = ($self);
    if (mask != NULL) {
        if (mask->xsize != img->xsize ||
            mask->ysize != img->ysize ||
            mask->zsize != img->zsize)
            SWIG_Error(0, "Float image and mask have different domains");

        for (int i = 0; i < img->n; i++) {
            if (mask->val[i])
                img->val[i] = IFT_INFINITY_FLT_NEG;
        }
    } else {
        for (int i = 0; i < img->n; i++) {
            img->val[i] = IFT_INFINITY_FLT_NEG;
        }
    }
}

void SetList(PyObject *list, float value)
{
    iftFImage *img = ($self);
    for (int i = 0; i < PyList_Size(list); i++) {
        long idx = PyLong_AsLong(PyList_GetItem(list, i));
        img->val[idx] = value;
    }
}