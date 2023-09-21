//! swig(newobject, stable)
iftDblArray *CreateDblArrayFromNumPy(PyObject *input) {
    PyArrayObject *arr = obj_to_array_no_conversion(input, NPY_DOUBLE);
    if (!require_contiguous(arr))
        SWIG_Error(0, "Input numpy array is not contiguous");
    
    int n_dims = array_numdims(input);
    if (n_dims != 1) {
        SWIG_Error(0, "Input Array must be 1-dimensional");
        return NULL;
    }
    
    npy_intp *shape = array_dimensions(input);
    
    iftDblArray *darr = iftCreateDblArray(shape[0]);
    
    double *ptr = array_data(input);
    memcpy(darr->val, ptr, sizeof (*ptr) * darr->n);
    
    return darr;
}