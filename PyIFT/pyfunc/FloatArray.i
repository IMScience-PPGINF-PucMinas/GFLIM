//! swig(newobject, stable)
iftFloatArray *CreateFloatArrayFromNumPy(PyObject *input) {
    PyArrayObject *arr = obj_to_array_no_conversion(input, NPY_FLOAT32);
    if (!require_contiguous(arr))
        SWIG_Error(0, "Input numpy array is not contiguous");
    
    int n_dims = array_numdims(input);
    if (n_dims != 1) {
        SWIG_Error(0, "Input Array must be 1-dimensional");
        return NULL;
    }
    
    npy_intp *shape = array_dimensions(input);
    
    iftFloatArray *farr = iftCreateFloatArray(shape[0]);
    
    float *ptr = array_data(input);
    memcpy(farr->val, ptr, sizeof (*ptr) * farr->n);
    
    return farr;
}