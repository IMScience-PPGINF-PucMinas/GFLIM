PyObject *AsNumPy(void) {
    iftIntMatrix *m = ($self);
    npy_intp dims[2] = {m->nrows, m->ncols};
    PyArrayObject* result = PyArray_SimpleNew(2, dims, NPY_INT32);

    int *data = array_data(result);
    memcpy(data, m->val, m->ncols * m->nrows * sizeof (*data));

    return PyArray_Return(result);
}

void FromNumPy(PyObject* input) {
    PyArrayObject *ary = obj_to_array_no_conversion(input, NPY_INT32);

    if (ary == NULL) return;

    if (!require_contiguous(ary))
        SWIG_Error(0, "Input numpy array is not contiguous");

    iftIntMatrix *m = ($self);
    float *ptr = array_data(ary);

    int n_dims = array_numdims(ary);
    if (n_dims != 2)
        SWIG_Error(0, "Input is not a 2D array");

    int ary_xsize = array_size(ary, 1);
    int ary_ysize = array_size(ary, 0);

    if (ary_ysize != m->nrows || ary_xsize == m->ncols) {
        m->nrows = ary_ysize;
        m->ncols = ary_xsize;
        m->n = ary_xsize * ary_ysize;
        m->val = iftRealloc(m->val, m->n * sizeof *m->val);
        if (m->val == NULL)
            iftError(MSG_MEMORY_ALLOC_ERROR, "FromNumpy");
    }
    memcpy(m->val, ptr, m->n * sizeof *m->val);
}

