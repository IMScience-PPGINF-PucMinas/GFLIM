
%newobject InsideMask;
iftSet *InsideMask(iftImage *mask)
{
    iftSet *s = ($self);
    iftSet *out = NULL;
    for (iftSet *cur = s; cur; cur = cur->next)
    {
        int p = cur->elem;
        if (mask->val[p])
            iftInsertSet(&out, p);
    }
    return out;
}


%newobject __add__;
iftSet *__add__(iftSet *other)
{
    iftSet *s = ($self);
    return iftSetConcat(s, other);
}


%newobject AsVoxelArray;
iftVoxelArray *AsVoxelArray(iftImage *image)
{
    iftSet *s = ($self);
    iftVoxelArray *array = iftCreateVoxelArray(iftSetSize(s));
    int i = 0;
    for (iftSet *S = s; S; S = S->next) {
        iftVoxel v = iftGetVoxelCoord(image, S->elem);
        array->val[i].x = v.x; array->val[i].y = v.y; array->val[i].z = v.z;
        i++;
    }
    return array;
}


PyObject *AsNumPy(void)
{
    iftSet *s = ($self);
    npy_intp dims[1] = {iftSetSize(s)};
    PyArrayObject *result = PyArray_SimpleNew(1, dims, NPY_INT32);

    int *data = array_data(result);
    int i = 0;
    for (iftSet *t = s; t; t = t->next) {
        data[i] = t->elem;
        i++;
    }

    return PyArray_Return(result);

}

PyObject *AsNumPy(const iftMImage *mimg)
{
    iftSet *s = ($self);
    npy_intp dims[2] = {iftSetSize(s), mimg->m};
    PyArrayObject *result = PyArray_SimpleNew(2, dims, NPY_FLOAT);

    float *data = array_data(result);

    int i = 0;
    for (iftSet *t = s; t; t = t->next)
    {
        int row = i * mimg->m;
        for (int j = 0; j < mimg->m; j++) {
            data[row + j] = mimg->val[t->elem][j];
        }
        i++;
    }

    return PyArray_Return(result);
}

PyObject *AsNumPy(const iftFImage *img)
{
    iftSet *s = ($self);
    npy_intp dims[1] = {iftSetSize(s)};
    PyArrayObject *result = PyArray_SimpleNew(1, dims, NPY_FLOAT);

    float *data = array_data(result);

    int i = 0;
    for (iftSet *t = s; t; t = t->next) {
        data[i] = img->val[t->elem];
        i++;
    }

    return PyArray_Return(result);
}


PyObject *AsNumPy(const iftImage *img)
{
    iftSet *s = ($self);
    npy_intp dims[1] = {iftSetSize(s)};
    PyArrayObject *result = PyArray_SimpleNew(1, dims, NPY_INT32);

    int *data = array_data(result);

    int i = 0;
    for (iftSet *t = s; t; t = t->next) {
        data[i] = img->val[t->elem];
        i++;
    }

    return PyArray_Return(result);
}


PyObject *CumulativeMean(const iftMImage *mimg)
{
    iftSet *s = ($self);
    npy_intp dims[2] = {iftSetSize(s), mimg->m};
    PyArrayObject *result = PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    double *data = array_data(result);

    int i = 0;
    for (iftSet *t = s; t; t = t->next)
    {
        int row = i * mimg->m;
        for (int j = 0; j < mimg->m; j++) {
            if (i == 0) {
                data[row + j] = (double) mimg->val[t->elem][j];
            } else {
                data[row + j] = data[row - mimg->m + j] + mimg->val[t->elem][j];
            }
        }
        i++;
    }

    for (int i = 0; i < dims[0]; i++) {
        int row = i * mimg->m;
        for (int j = 0; j < dims[1]; j++) {
            data[row + j] /= (i + 1);
        }
    }

    return PyArray_Return(result);
}
