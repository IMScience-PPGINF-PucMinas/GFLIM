PyObject *GetData(void) {
    iftDataSet *Z = ($self);

    if (Z->data == NULL)
        return Py_BuildValue(""); // returning a PyObject as None
    else {
        npy_intp dims[2] = {Z->nsamples, Z->nfeats};
        PyArrayObject *ref = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, (float*) Z->data->val);
        PyArrayObject *npy = PyArray_SimpleNew(2, dims, NPY_FLOAT32);
        PyArray_CopyInto(npy, ref); // (dst, src)

        return PyArray_Return(npy);
    }
}

PyObject *GetProjection(void) {
    iftDataSet *Z = ($self);

    if (Z->projection == NULL)
        return Py_BuildValue(""); // returning a PyObject as None
    else {
        npy_intp dims[2] = {Z->projection->nrows, Z->projection->ncols};
        PyArrayObject *ref = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT64, (double*) Z->projection->val);
        PyArrayObject *npy = PyArray_SimpleNew(2, dims, NPY_FLOAT64);
        PyArray_CopyInto(npy, ref); // (dst, src)

        return PyArray_Return(npy);

    }
}


PyObject *GetTrueLabels(void) {
    iftDataSet *Z = ($self);

    npy_intp dims[1] = {Z->nsamples};
    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_INT32);

    int *ptr = (int *) array_data(npy);
    for (int s = 0; s < Z->nsamples; s++)
        ptr[s] = Z->sample[s].truelabel;

    return PyArray_Return(npy);
}

PyObject *GetLabels(void) {
    iftDataSet *Z = ($self);

    npy_intp dims[1] = {Z->nsamples};
    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_INT32);

    int *ptr = (int *) array_data(npy);
    for (int s = 0; s < Z->nsamples; s++)
        ptr[s] = Z->sample[s].label;

    return PyArray_Return(npy);
}

PyObject *GetIds(void) {
    iftDataSet *Z = ($self);

    npy_intp dims[1] = {Z->nsamples};
    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_INT32);

    int *ptr = (int *) array_data(npy);
    for (int s = 0; s < Z->nsamples; s++)
        ptr[s] = Z->sample[s].id;

    return PyArray_Return(npy);
}

PyObject *GetGroups(void) {
    iftDataSet *Z = ($self);

    npy_intp dims[1] = {Z->nsamples};
    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_INT32);

    int *ptr = (int *) array_data(npy);
    for (int s = 0; s < Z->nsamples; s++)
        ptr[s] = Z->sample[s].group;

    return PyArray_Return(npy);
}

PyObject *GetIsSupervised(void) {
    iftDataSet *Z = ($self);

    npy_intp dims[1] = {Z->nsamples};
    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_BOOL);

    bool *ptr = (bool *) array_data(npy);
    for (int s = 0; s < Z->nsamples; s++)
        ptr[s] = Z->sample[s].isSupervised;

    return PyArray_Return(npy);
}

PyObject *GetIsLabelPropagated(void) {
    iftDataSet *Z = ($self);

    npy_intp dims[1] = {Z->nsamples};
    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_BOOL);

    bool *ptr = (bool *) array_data(npy);
    for (int s = 0; s < Z->nsamples; s++)
        ptr[s] = Z->sample[s].isLabelPropagated;

    return PyArray_Return(npy);
}

PyObject *GetWeights(void) {
    iftDataSet *Z = ($self);

    npy_intp dims[1] = {Z->nsamples};
    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_FLOAT32);

    float *ptr = (float*) array_data(npy);
    for (int s = 0; s < Z->nsamples; s++)
        ptr[s] = Z->sample[s].weight;

    return PyArray_Return(npy);
}

PyObject *GetStatus(void) {
    iftDataSet *Z = ($self);

    npy_intp dims[1] = {Z->nsamples};
    PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_UINT32);

    unsigned int *ptr = (unsigned int *) array_data(npy);
    for (int s = 0; s < Z->nsamples; s++)
        ptr[s] = Z->sample[s].status;

    return PyArray_Return(npy);
}

PyObject *GetAlphas(void) {
    iftDataSet *Z = ($self);

    if (Z->alpha == NULL)
        return Py_BuildValue(""); // returning a PyObject as None
    else {
        npy_intp dims[1] = {Z->nfeats};
        PyArrayObject *ref = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT32, (float*) Z->alpha);
        PyArrayObject *npy = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
        PyArray_CopyInto(npy, ref); // (dst, src)

        return PyArray_Return(npy);
    }
}

PyObject *GetRefData(void) {
    iftDataSet *Z = ($self);

    if (Z->ref_data == NULL) {
        return PyList_New(0); // return an empty List
    }
    else {
        if (Z->ref_data_type == IFT_REF_DATA_FILESET) {
            PyObject *list = PyList_New(Z->nsamples);

            iftFileSet *ref_data = (iftFileSet *) Z->ref_data;

            for (int i = 0; i < ref_data->n; i++)
                PyList_SetItem(list, i, PyString_FromString(ref_data->files[i]->path));
            return list;
        }
        if (Z->ref_data_type == IFT_REF_DATA_CSV) {
            iftCSV *ref_data = (iftCSV *) Z->ref_data;

            PyObject *csv = PyList_New(ref_data->nrows);

            for (int r = 0; r < ref_data->nrows; r++) {
                PyObject *row = PyList_New(ref_data->ncols);
                
                for (int c = 0; c < ref_data->ncols; c++)
                    PyList_SetItem(row, c, PyString_FromString(ref_data->data[r][c]));
                PyList_SetItem(csv, r, row);
            }
            
            return csv;
        }
        else return PyList_New(0); // return an empty List
    }
}


void SetRefData(PyObject *list) {

    /// Example of list: ref_data = [['00000001_000001.png'], ['00000001_000002.png'], ['00000002_000001.png'], ['00000002_000002.png']]

    iftDataSet *Z = ($self);

//    if (Z->ref_data != NULL)
//        SWIG_Error(0, "Existing Reference Data is not NULL");
 
    if (!PyList_CheckExact(list)){
        SWIG_Error(0, "Input \"list\" must be a list.");
        return;
    }

    long list_size = PyList_Size(list);

    if (list_size != Z->nsamples) {
        char error[512];
        sprintf(error, "Size of the reference data (list of paths) is != from the number of samples: %ld != %ld", list_size, Z->nsamples);
        SWIG_Error(0, error);
        return;
    }
    
    
    PyObject *item0 = PyList_GetItem(list, 0);

    if (!PyList_CheckExact(item0)){
        SWIG_Error(0, "Input \"item0\" must be a list.");
        return;
    }

    long n_cols = PyList_Size(item0);

    if (n_cols == 1) {
        Z->ref_data_type = IFT_REF_DATA_FILESET;
        
        iftFileSet *ref_data = iftCreateFileSet(Z->nsamples);
        
        for (int i = 0; i < Z->nsamples; i++) {
            PyObject *row = PyList_GetItem(list, i);

            if (!PyList_CheckExact(row)){
               SWIG_Error(0, "Input \"row\" must be a list.");
               return;
            }

            char *path = PyUnicodeToString(PyList_GetItem(row, 0));

            ref_data->files[i] = iftCreateFile(path);
            Z->sample[i].id = i;
        }
        Z->ref_data = ref_data;
    }
    else {
        Z->ref_data_type = IFT_REF_DATA_CSV;

        iftCSV *ref_data = iftCreateCSV(Z->nsamples, n_cols);

        for (int r = 0; r < Z->nsamples; r++) {
            PyObject *row = PyList_GetItem(list, r);
            
            for (int c = 0; c < n_cols; c++) {
                strcpy(ref_data->data[r][c], PyUnicodeToString(PyList_GetItem(row, c)));
                Z->sample[r].id = r;
            }
        }
        
        Z->ref_data = ref_data;
    }
}


void SetTrueLabels(PyObject *in_npy) {
    char error[512];
    iftDataSet *Z = ($self);

    PyArrayObject *npy = NULL;
    if (in_npy != Py_None) {
        if (PyArray_TYPE(in_npy) != NPY_INT32) {
            SWIG_Error(0, "Numpy array must be INT32");
            return;
        }

        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);

        int n_dims = array_numdims(npy);
        if (n_dims != 1) {
            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
            SWIG_Error(0, error);
            return;
        }

        npy_intp *shape = array_dimensions(npy); // (nsamples,)
        if (shape[0] != Z->nsamples) {
            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
            SWIG_Error(0, error);
            return;
        }

        int *ptr = (int *) array_data(npy);
        for (int s = 0; s < Z->nsamples; s++)
            Z->sample[s].truelabel = ptr[s];

        Z->nclasses = iftCountNumberOfClassesDataSet(Z);
    }

}


void SetId(PyObject *in_npy) {
    char error[512];
    iftDataSet *Z = ($self);

    PyArrayObject *npy = NULL;
    if (in_npy != Py_None) {
        if (PyArray_TYPE(in_npy) != NPY_INT32) {
            SWIG_Error(0, "Numpy array must be INT32");
            return;
        }

        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);

        int n_dims = array_numdims(npy);
        if (n_dims != 1) {
            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
            SWIG_Error(0, error);
            return;
        }

        npy_intp *shape = array_dimensions(npy); // (nsamples,)
        if (shape[0] != Z->nsamples) {
            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
            SWIG_Error(0, error);
            return;
        }

        int *ptr = (int *) array_data(npy);
        for (int s = 0; s < Z->nsamples; s++)
            Z->sample[s].id = ptr[s];
   
    }
}


void SetLabels(PyObject *in_npy) {
    char error[512];
    iftDataSet *Z = ($self);

    PyArrayObject *npy = NULL;
    if (in_npy != Py_None) {
        if (PyArray_TYPE(in_npy) != NPY_INT32) {
            SWIG_Error(0, "Numpy array must be INT32");
            return;
        }

        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);

        int n_dims = array_numdims(npy);
        if (n_dims != 1) {
            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
            SWIG_Error(0, error);
            return;
        }

        npy_intp *shape = array_dimensions(npy); // (nsamples,)
        if (shape[0] != Z->nsamples) {
            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
            SWIG_Error(0, error);
            return;
        }

        int *ptr = (int *) array_data(npy);
        for (int s = 0; s < Z->nsamples; s++)
            Z->sample[s].label = ptr[s];

	int tmp = iftCountNumberOfClassesDataSet(Z);
	if (tmp > Z->nclasses)
	    Z->nclasses = tmp;     
    }
}

void SetWeight(PyObject *input)
{
    char error[512];
    PyArrayObject *ary = obj_to_array_no_conversion(input, NPY_FLOAT);

    if (ary == NULL) return;

    if (!require_contiguous(ary))
        SWIG_Error(0, "Input numpy array is not contiguous");

    iftDataSet *Z = ($self);
    float* ptr = array_data(ary);


    int n_dims = array_numdims(ary);
    if (n_dims != 1) {
            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
            SWIG_Error(0, error);
            return;
        }

    npy_intp *shape = array_dimensions(ary); // (nsamples,)
    if (shape[0] != Z->nsamples) {
            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
            SWIG_Error(0, error);
            return;
    }

    for (int s = 0; s < Z->nsamples; s++)
        Z->sample[s].weight = ptr[s];

}

void SetGroups(PyObject *in_npy) {
    char error[512];
    iftDataSet *Z = ($self);

    PyArrayObject *npy = NULL;
    if (in_npy != Py_None) {
        if (PyArray_TYPE(in_npy) != NPY_INT32) {
            SWIG_Error(0, "Numpy array must be INT32");
            return;
        }

        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);

        int n_dims = array_numdims(npy);
        if (n_dims != 1) {
            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
            SWIG_Error(0, error);
            return;
        }

        npy_intp *shape = array_dimensions(npy); // (nsamples,)
        if (shape[0] != Z->nsamples) {
            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
            SWIG_Error(0, error);
            return;
        }

        int *ptr = (int *) array_data(npy);
        for (int s = 0; s < Z->nsamples; s++)
            Z->sample[s].group = ptr[s];

        Z->ngroups = iftCountNumberOfGroupsDataSet(Z);
    }
}

void SetGroupsFromImage(const iftImage *img, bool increment)
{
    iftDataSet *Z = ($self);
    if (Z->nsamples != img->n)
        return;

    for (int i = 0; i < img->n; i++)
        Z->sample[i].group = img->val[i] + increment;

    Z->ngroups = iftCountNumberOfGroupsDataSet(Z);
}


void SetNClasses(int n){
  
    iftDataSet *Z = ($self);
    Z->nclasses = n;
}


void SetData(PyObject *input)
{
    PyArrayObject *ary = obj_to_array_no_conversion(input, NPY_FLOAT);

    if (ary == NULL) return;

    if (!require_contiguous(ary))
        SWIG_Error(0, "Input numpy array is not contiguous");

    iftDataSet *Z = ($self);
    float* ptr = array_data(ary);

    int n_dims = array_numdims(ary);
    if (n_dims != 2)
        SWIG_Error(0, "Input must be a 2-D numpy array");

    int ary_xsize = array_size(ary, 1);
    int ary_ysize = array_size(ary, 0);

    if (Z->nsamples != ary_ysize || Z->nfeats != ary_xsize)
        SWIG_Error(0, "Input must be the same size as the previously set data");

    for (int i = 0; i < Z->nsamples; i++) {
        for (int j = 0; j < Z->nfeats; j++) {
            Z->sample[i].feat[j] = ptr[i * Z->nfeats + j];
        }
    }
}

void SetStatus(PyObject *in_npy) {
    char error[512];
    iftDataSet *Z = ($self);

    PyArrayObject *npy = NULL;
    if (in_npy != Py_None) {
        if (PyArray_TYPE(in_npy) != NPY_INT32) {
            SWIG_Error(0, "Numpy array must be INT32");
            return;
        }

        npy = obj_to_array_no_conversion(in_npy, NPY_INT32);

        int n_dims = array_numdims(npy);
        if (n_dims != 1) {
            sprintf(error, "Numpy array should be 1-dimensional. Found: %d dimensions", n_dims);
            SWIG_Error(0, error);
            return;
        }

        npy_intp *shape = array_dimensions(npy); // (nsamples,)
        if (shape[0] != Z->nsamples) {
            sprintf(error, "Numpy array size is different from the number of samples of the Dataset. %d != %d", shape[0], Z->nsamples);
            SWIG_Error(0, error);
            return;
        }

        int *ptr = (int *) array_data(npy);
        for (int s = 0; s < Z->nsamples; s++)
            Z->sample[s].status = ptr[s];
   
    }
}

void SetNTrainSamples(int n){
  
    iftDataSet *Z = ($self);
    Z->ntrainsamples = n;
}


