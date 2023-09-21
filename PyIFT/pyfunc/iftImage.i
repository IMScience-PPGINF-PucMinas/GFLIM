//! swig(newobject, stable)
iftImage *CreateImageFromNumPy(PyObject *input, bool is3D) {
    PyArrayObject *ary = obj_to_array_no_conversion(input, NPY_INT32);

    if (ary == NULL)
        return NULL;

    if (!require_contiguous(ary))
        SWIG_Error(0, "Input numpy array is not contiguous");

    // supported shapes
    // 2D gray image: (ysize, xsize)
    // 2D color image: (ysize, xsize, 3)
    // 3D gray image: (zsize, ysize, xsize)
    int n_dims = array_numdims(ary);


    if ((n_dims != 2) && (n_dims != 3) && (n_dims != 4)) {
        char error[512];
        sprintf(error, "Image Numpy Array should have 2 or 3 or 4 dimensions. Found: %d", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    // 3D gray image is actually a 2D gray image
    if ((is3D) && (n_dims == 2))
        is3D = false;
    
    bool is_color_image = false;

    int xsize, ysize, zsize;

    if (is3D) {
        xsize = array_size(ary, 2);
        ysize = array_size(ary, 1);
        zsize = array_size(ary, 0);    

        // color image
        is_color_image = (n_dims == 4);
    }
    else {
        xsize = array_size(ary, 1);
        ysize = array_size(ary, 0);
        zsize = 1;

        // color image
        is_color_image = (n_dims == 3) && (array_size(ary, 2) == 3);
    }


    int *ptr = array_data(ary);
    iftImage *img = NULL;

    if (is_color_image) {
        img = iftCreateColorImage(xsize, ysize, zsize, 8);

        for (int p = 0; p < img->n; p++) {
            img->val[p] = ptr[3 * p];
            img->Cb[p] = ptr[(3 * p) + 1];
            img->Cr[p] = ptr[(3 * p) + 2];
        }
    }
    else {
        img = iftCreateImage(xsize, ysize, zsize);
        for (int p = 0; p < img->n; p++)
            img->val[p] = ptr[p];
    }

    return img;
}

//! swig(newobject)
iftImage *NDArrayToLabelImage(PyObject *input)
{
    int is_new = 0;
    PyArrayObject *ary = obj_to_array_allow_conversion(input, NPY_DOUBLE, &is_new);

    if (ary == NULL)
        return NULL;

    if (!require_contiguous(ary))
        SWIG_Error(0, "Input numpy array is not contiguous");

    int n_dims = array_numdims(ary);

    if (n_dims != 3) {
        char error[256];
        sprintf(error, "Image Numpy Array must have 3 dimensions, %d found.", n_dims);
        SWIG_Error(0, error);
        return NULL;
    }

    int n_labels = array_size(ary, 0);
    int xsize = array_size(ary, 2);
    int ysize = array_size(ary, 1);
    double *ptr = array_data(ary);
    iftImage *img = iftCreateImage(xsize, ysize, 1);

    for (int p = 0; p < img->n; p++) {
        double max_val = -1.0;
        int max_lb = -1;
        for (int b = 0; b < n_labels; b++)
        {
            int idx = img->n * b + p;
            if (ptr[idx] > max_val) {
                max_lb = b;
                max_val = ptr[idx];
            }
        }
        img->val[p] = max_lb;
    }

    return img;
}

//! swig(newobject)
iftImage *ImageFromNumPyLabelIndex(PyObject *index, PyObject *label,
                              int xsize, int ysize, int zsize, bool is_whole_array)
{
    int is_new = 0;
    PyArrayObject *idx = obj_to_array_allow_conversion(index, NPY_INT64, &is_new);
    PyArrayObject *lbl = obj_to_array_allow_conversion(label, NPY_INT64, &is_new);

    if (!idx || !lbl)
        return NULL;

    if (!require_contiguous(lbl) || !require_contiguous(idx))
        SWIG_Error(0, "Inputs must be contiguous array");

    if (array_numdims(idx) != 1 || array_numdims(lbl) != 1)
        SWIG_Error(0, "Inputs must be one dimensional arrays");

    int size = xsize * ysize * zsize;
    if (is_whole_array && (size != array_size(idx, 0) || size != array_size(lbl, 0)))
        SWIG_Error(0, "Array size must be equal to image frame");

    iftImage *img = iftCreateImage(xsize, ysize, zsize);

    long *d_idx = array_data(idx);
    long *d_lbl = array_data(lbl);

    for (int i = 0; i < img->n; i++) {
        img->val[ d_idx[i] ] = d_lbl[i];
    }

    return img;
}


//! swig(newobject)
iftImage *ImageFromNumPyMapping(const iftImage *components, PyObject *reference, PyObject *label)
{
    int is_new = 0;
    PyArrayObject *ref = obj_to_array_allow_conversion(reference, NPY_INT64, &is_new);
    PyArrayObject *lbl = obj_to_array_allow_conversion(label, NPY_INT64, &is_new);

    if (!lbl || !ref)
        return NULL;

    if (!require_contiguous(lbl) || !require_contiguous(ref))
        SWIG_Error(0, "Inputs must be contiguous array");

    if (array_numdims(ref) != 1 || array_numdims(lbl) != 1)
        SWIG_Error(0, "Inputs must be one dimensional arrays");

    iftImage *img = iftCreateImage(components->xsize, components->ysize, components->zsize);

    long *d_ref = array_data(ref);
    long *d_lbl = array_data(lbl);

    for (int i = 0; i < img->n; i++) {
        int idx = d_ref[ components->val[i] - 1 ];
        img->val[i] = d_lbl[idx];
    }

    return img;
}



PyObject *iftPathValue(PyObject *in_image, PyObject *in_markers, float penalization)
{
    int is_new_image = 0, is_new_markers = 0;
    PyArrayObject *image = obj_to_array_allow_conversion(in_image, NPY_FLOAT, &is_new_image);
    if (!image || !require_contiguous(image))
        Py_RETURN_NONE;

    PyArrayObject *markers = obj_to_array_allow_conversion(in_markers, NPY_INT32, &is_new_markers);
    if (!markers || !require_contiguous(markers))
        Py_RETURN_NONE;

    if (!require_dimensions(image, 3) || !require_dimensions(markers, 2))
        Py_RETURN_NONE;

    npy_intp *image_shape = PyArray_SHAPE(image);
    if (!require_size(markers, image_shape + 1, 2))
        Py_RETURN_NONE;

    PyArrayObject *dist = PyArray_SimpleNew(2, image_shape + 1, NPY_FLOAT);
    if (!dist || !require_contiguous(dist))
        Py_RETURN_NONE;
    int height = image_shape[1], width = image_shape[2], size = width * height;
    float *dist_ptr = array_data(dist);
    iftFHeap *Q = iftCreateFHeap(size, dist_ptr);

    int *markers_ptr = array_data(markers);
    for (int i = 0; i < size; i++) {
        if (markers_ptr[i]) {
            dist_ptr[i] = 0;
            iftInsertFHeap(Q, i);
        } else {
            dist_ptr[i] = IFT_INFINITY_FLT;
        }
    }

    iftAdjRel *A = iftCircular(1.0f);

    float *image_ptr = array_data(image);
    while (!iftEmptyFHeap(Q))
    {
        int p = iftRemoveFHeap(Q);
        div_t res = div(p, width);
        iftVoxel u = {.x = res.rem, .y = res.quot, .z = 0};

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (v.x < 0 || v.x >= width || v.y < 0 || v.y >= height)
                continue;
            int q = v.y * width + v.x;
            if (dist_ptr[q] < dist_ptr[p])
                continue;
            float norm = 0;
            for (int j = 0, row = 0; j < image_shape[0]; j++, row += size) {
                float dif = image_ptr[row + p] - image_ptr[row + q];
                norm += dif * dif;
            }
            norm = sqrtf(norm);
            norm = iftMax(norm, dist_ptr[p] + penalization);
            if (norm < dist_ptr[q]) {
                dist_ptr[q] = norm;
                if (Q->color[q] == IFT_GRAY)
                    iftGoUpFHeap(Q, Q->pos[q]);
                else
                    iftInsertFHeap(Q, q);
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftDestroyFHeap(&Q);
    if (is_new_image)
        Py_DECREF(image);
    if (is_new_markers)
        Py_DECREF(markers);

    return PyArray_Return(dist);
}

static inline float squared(float x) {
    return x * x;
}

PyObject *iftGaussian(PyObject *in_image, PyObject *in_markers, float sigma)
{
    int is_new_image = 0, is_new_markers = 0;
    PyArrayObject *image = obj_to_array_allow_conversion(in_image, NPY_FLOAT, &is_new_image);
    if (!image || !require_contiguous(image))
        Py_RETURN_NONE;

    PyArrayObject *markers = obj_to_array_allow_conversion(in_markers, NPY_INT32, &is_new_markers);
    if (!markers || !require_contiguous(markers))
        Py_RETURN_NONE;

    if (!require_dimensions(image, 3) || !require_dimensions(markers, 2))
        Py_RETURN_NONE;

    npy_intp *image_shape = PyArray_SHAPE(image);
    if (!require_size(markers, image_shape + 1, 2))
        Py_RETURN_NONE;

    PyArrayObject *dist = PyArray_SimpleNew(2, image_shape + 1, NPY_FLOAT);
    if (!dist || !require_contiguous(dist))
        Py_RETURN_NONE;
    int height = image_shape[1], width = image_shape[2], size = width * height;
    float *dist_ptr = array_data(dist);
    iftFHeap *Q = iftCreateFHeap(size, dist_ptr);
    iftSetRemovalPolicyFHeap(Q, IFT_MAXVALUE);

    int *markers_ptr = array_data(markers);
    iftVoxel *root = calloc(size, sizeof (*root));
    for (int i = 0; i < size; i++) {
        if (markers_ptr[i]) {
            dist_ptr[i] = 1.0f;
            div_t res = div(i, width);
            root[i].x = res.rem;
            root[i].y = res.quot;
            iftInsertFHeap(Q, i);
        } else {
            dist_ptr[i] = 0.0f;
        }
    }

    iftAdjRel *A = iftCircular(1.0f);

    float *image_ptr = array_data(image);
    while (!iftEmptyFHeap(Q))
    {
        int p = iftRemoveFHeap(Q);
        div_t res = div(p, width);
        iftVoxel u = {.x = res.rem, .y = res.quot, .z = 0};

        for (int i = 1; i < A->n; i++)
        {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (v.x < 0 || v.x >= width || v.y < 0 || v.y >= height)
                continue;
            int q = v.y * width + v.x;
            if (dist_ptr[q] > dist_ptr[p])
                continue;
            float gauss_dist = expf( - (squared(v.x-root[p].x) + squared(v.y-root[p].y)) / sigma);
            if (gauss_dist > dist_ptr[q]) {
                dist_ptr[q] = gauss_dist;
                root[q].x = root[p].x; root[q].y = root[p].y;
                if (Q->color[q] == IFT_GRAY)
                    iftGoUpFHeap(Q, Q->pos[q]);
                else
                    iftInsertFHeap(Q, q);
            }
        }
    }

    iftDestroyAdjRel(&A);
    iftDestroyFHeap(&Q);
    if (is_new_image)
        Py_DECREF(image);
    if (is_new_markers)
        Py_DECREF(markers);
    free(root);

    return PyArray_Return(dist);
}

