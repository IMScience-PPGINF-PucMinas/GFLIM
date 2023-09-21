iftImage* __getitem__(int band){
    iftMImage* mimg = ($self);
    iftImage* img;

    img = iftMImageToImage(mimg, 255, band);
    return img;
}

void __setitem__(int band, iftImage* img){
    iftMImage* mimg = ($self);

    if(band >= mimg->m){
        SWIG_Error(12, "Index out of MImage bands");
    }

    for(int i = 0; i < img->n; i++){
        mimg->val[i][band] = (float) img->val[i];
    }
}

/* bugado
void AddImage(iftImage* img){
    iftMImage* mimg = ($self);
    mimg->m++;
    iftRealloc(mimg->band, sizeof(iftBand) * mimg->m);

    for(int i = 0; i < img->n; i++){
        mimg->band[mimg->m - 1].val[i] = (float) img->val[i];
    }
}
*/

PyObject* AsNumPy(void) {

    iftMImage *mimg = ($self);

    int n_dims = 4;
    npy_intp *dims = calloc(n_dims, sizeof (*dims));
    dims[0] = mimg->zsize;
    dims[1] = mimg->ysize;
    dims[2] = mimg->xsize;
    dims[3] = mimg->m;

    int n_bands = mimg->m;

    PyArrayObject *result = PyArray_SimpleNew(n_dims, dims, NPY_FLOAT);
    float *data = array_data(result);

    for (int p = 0; p < mimg->n; p++) {
        int start_voxel_idx = p * n_bands;
        for (int b = 0; b < n_bands; b++) {
            data[start_voxel_idx + b] = mimg->val[p][b];
        }
    }

    free(dims);
    return PyArray_Return(result);
}

void WriteAsNumPy(const char *filename) {
    iftMImage *mimg = ($self);
    iftWriteMImageAsNumPy(mimg, filename);
}
