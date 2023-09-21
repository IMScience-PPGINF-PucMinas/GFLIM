//! swig(newobject, stable)
iftIntArray *CreateIntArrayFromNumPy(PyObject *input) {
    PyArrayObject *arr = obj_to_array_no_conversion(input, NPY_INT32);
    if (!require_contiguous(arr))
        SWIG_Error(0, "Input numpy array is not contiguous");
    
    int n_dims = array_numdims(input);
    if (n_dims != 1) {
        SWIG_Error(0, "Input Integer Array must be 1-dimensional");
        return NULL;
    }
    
    npy_intp *shape = array_dimensions(input);
    
    iftIntArray *iarr = iftCreateIntArray(shape[0]);
    
    int *ptr = array_data(input);
    memcpy(iarr->val, ptr, sizeof (*ptr) * iarr->n);
    
    return iarr;
}