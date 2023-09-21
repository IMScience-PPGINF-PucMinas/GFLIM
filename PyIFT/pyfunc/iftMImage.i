//! swig(newobject, stable)
iftMImage *CreateMImageFromNumPy(PyObject *input) {
    int is_new = 0;
    PyArrayObject *ary = obj_to_array_allow_conversion(input, NPY_FLOAT, &is_new);

    if (ary == NULL)
        return NULL;

    if (!require_contiguous(ary)) {
        SWIG_Error(0, "Input numpy array is not contiguous");
        return NULL;
    }

    int n_dims = array_numdims(ary);

    // (ysize, xsyze, nbands) or
    // (zsize, ysize, xsyze, nbands)
    int xsize, ysize, zsize, n_bands; 

    if (n_dims == 3) {
        zsize = 1;
        ysize = array_size(ary, 0);
        xsize = array_size(ary, 1);
        n_bands = array_size(ary, 2);
    }
    else if (n_dims == 4) {
        zsize = array_size(ary, 0);
        ysize = array_size(ary, 1);
        xsize = array_size(ary, 2);
        n_bands = array_size(ary, 3);
    }
    else {
        char error[256];
        sprintf(error, "Image Numpy Array must have 3 or 4 dimensions, %d found.", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    iftMImage *mimg = iftCreateMImage(xsize, ysize, zsize, n_bands);

    float *ptr = array_data(ary);

    for (int p = 0; p < mimg->n; p++) {
        int start_voxel_idx = p * n_bands;

        for (int b = 0; b < n_bands; b++) {
            mimg->val[p][b] = ptr[start_voxel_idx + b];
        }
    }
    
    return mimg;
}
