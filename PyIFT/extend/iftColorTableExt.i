
PyObject *AsNumPy(void) {
    iftColorTable *ctb = ($self);

    npy_intp dims[2] = {ctb->ncolors, 3};
    PyArrayObject *npy = PyArray_SimpleNew(2, dims, NPY_INT32);

    int *ptr = (int *) array_data(npy);

    for (int c = 0; c < ctb->ncolors; c++) {
        int i = c * 3;
        ptr[i]   = ctb->color[c].val[0]; // R
        ptr[i+1] = ctb->color[c].val[1]; // G
        ptr[i+2] = ctb->color[c].val[2]; // B
    }

    return PyArray_Return(npy);
}

void __setitem__(int idx, const char *hex)
{
    iftColorTable *ctb = ($self);
    
    if (idx >= ctb->ncolors)
        SWIG_Error(12, "Index out of iftColorTable range");
    else
        ctb->color[idx] = iftRGBtoYCbCr(iftHexaColorToRGB(hex), 255);
}
