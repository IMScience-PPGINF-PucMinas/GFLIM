PyObject *AsList(void)
{
    iftVoxelArray *v_arry = ($self);

    PyObject *list = PyList_New(v_arry->n);
    for (int i = 0; i < v_arry->n; i++) {
        PyObject *voxel = PyList_New(3);
        PyList_SetItem(voxel, 0, PyLong_FromLong((long) v_arry->val[i].x));
        PyList_SetItem(voxel, 1, PyLong_FromLong((long) v_arry->val[i].y));
        PyList_SetItem(voxel, 2, PyLong_FromLong((long) v_arry->val[i].z));
        PyList_SetItem(list, i, voxel);
    }
    return list;
}


%newobject __add__;
iftVoxelArray *__add__(iftVoxelArray *other)
{
    iftVoxelArray *current = ($self);
    iftVoxelArray *out = iftCreateVoxelArray(other->n + current->n);
    int i = 0;
    for (; i < current->n; i++) {
        iftCopyVoxel(&out->val[i], &current->val[i]);
    }

    for (int j = 0; j < other->n; i++, j++) {
        iftCopyVoxel(&out->val[i], &current->val[j]);
    }
    return out;
}


PyObject *__getitem__(int index)
{
    iftVoxelArray *arr = ($self);
    if (index < 0 || index >= arr->n)
    {
        char error[512];
        sprintf(error, "Index %d, out of range [0, %d)\n", index, arr->n);
        SWIG_Error(0, error);
    }


    PyObject *list = PyList_New(3);
    PyList_SetItem(list, 0, PyLong_FromLong((long) arr->val[index].x));
    PyList_SetItem(list, 1, PyLong_FromLong((long) arr->val[index].y));
    PyList_SetItem(list, 2, PyLong_FromLong((long) arr->val[index].z));

    return list;
}


iftVoxelArray *__setitem__(int index, PyObject *list)
{
    if (!PyList_Check(list) || PyList_Size(list) != 3)
        return PyErr_Format(PyExc_RuntimeError, NULL);

    iftVoxelArray *arr = ($self);
    if (index < 0 || index >= arr->n)
    {
        char error[512];
        sprintf(error, "Index %d, out of range [0, %d)\n", index, arr->n);
        SWIG_Error(0, error);
    }


    arr->val[index].x = PyLong_AsLong(PyList_GetItem(list, 0));
    arr->val[index].y = PyLong_AsLong(PyList_GetItem(list, 1));
    arr->val[index].z = PyLong_AsLong(PyList_GetItem(list, 2));

    return arr;
}


PyObject *__len__(void)
{
    return PyLong_FromLong((long) ($self)->n);
}


iftVoxelArray *Insert(int index, PyObject *list)
{
    if (!PyList_Check(list) || PyList_Size(list) != 3)
        SWIG_Error(0, "Expected PyList");

    iftVoxelArray *arr = ($self);
    if (index < 0 || index > arr->n)
        SWIG_Error(0, "Index out of range");

    arr->n += 1;
    arr->val = iftRealloc(arr->val, arr->n * (sizeof *arr->val));
    if (!arr->val)
        SWIG_Error(0, "Could not allocate more memory");

    for (int i = arr->n - 1; i > index; i--) {
        iftCopyVoxel(&arr->val[i - 1], &arr->val[i]);
    }

    arr->val[index].x = PyLong_AsLong(PyList_GetItem(list, 0));
    arr->val[index].y = PyLong_AsLong(PyList_GetItem(list, 1));
    arr->val[index].z = PyLong_AsLong(PyList_GetItem(list, 2));
    return arr;
}


iftVoxelArray *Remove(int index)
{
    iftVoxelArray *v_arr = ($self);
    if (index > v_arr->n) {
        SWIG_Error(0, "Index out of range");
    }

    v_arr->n--;
    for (int i = index; i < v_arr->n; i++) {
        iftCopyVoxel(&v_arr->val[i+1], &v_arr->val[i]);
    }

    return v_arr;
}


PyObject *AsNumPy(const iftMImage *mimg, const iftAdjRel *A)
{
    iftVoxelArray *v_arr = ($self);

    npy_intp dims[2] = {v_arr->n, A->n * mimg->m};
    PyArrayObject *npy = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *data = (double*) array_data(npy);

    for (int i = 0; i < v_arr->n; i++)
    {
        iftVoxel *u = &v_arr->val[i];
        if (!iftValidVoxel(mimg, (*u))) {
            char buf[128];
            sprintf(buf, "Voxel (%d, %d, %d), on index %d is not valid for this iftMImage",
                    u->x, u->y, u->z, i);
            SWIG_Error(0, buf);
        }
        int p = iftGetVoxelIndex(mimg, (*u));
        int row = i * mimg->m * A->n;
        for (int j = 0, k = 0; j < A->n; j++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, (*u), j);
            // if is not valid keep the previous p
            if (iftValidVoxel(mimg, v))
                p = iftGetVoxelIndex(mimg, v);

            for (int b = 0; b < mimg->m; b++, k++) {
                data[row + k] = (double) mimg->val[p][b];
            }
        }
    }

    return PyArray_Return(npy);
}


%newobject AsSet;
iftSet *AsSet(const iftImage *image)
{
    iftVoxelArray *arr = ($self);
    iftSet *s = NULL;

    for (int i = 0; i < arr->n; i++)
        iftInsertSet(&s, iftGetVoxelIndex(image, arr->val[i]));
    return s;
}
