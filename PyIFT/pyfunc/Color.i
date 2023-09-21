//! swig(newobject, stable)
iftColorTable *CreateColorTableFromNumPy(PyObject *input) {
    PyArrayObject *arr = obj_to_array_no_conversion(input, NPY_INT32);
    if (!require_contiguous(arr))
        SWIG_Error(0, "Input numpy array is not contiguous");

    // (n_colors, 3) ==> 3 channels: [:, 0] = Y, [:, 0] = Cb, [:, 0] = Cr
    int n_dims = array_numdims(arr);
    if (n_dims != 2) {
        char error[512];
        sprintf(error, "Color Numpy Array must have 2 dimensions. Found: %d", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    npy_intp *shape = array_dimensions(arr); // (n_colors, 3)
    if (shape[1] != 3) {
        char error[512];
        sprintf(error, "Number of Channels in the Color Number Array must have 3. Found: %d", shape[1]);
        SWIG_Error(0, error);
        return NULL;
    }

    int *ptr = array_data(arr);

    iftColorTable *cmap = iftCreateColorTable(shape[0]);
    
    for (int c = 0; c < cmap->ncolors; c++) {
        int i = c * 3;
        cmap->color[c].val[0] = ptr[i];   // Y
        cmap->color[c].val[1] = ptr[i+1]; // Cb
        cmap->color[c].val[2] = ptr[i+2]; // Cr
    }

    return cmap;
}


