PyObject* AsList() {
    iftHist *hist = ($self);

    PyObject *list = PyList_New(hist->nbins);
    
    for(int i = 0; i < hist->nbins; i++){
        double val = (double) hist->val[i];
        PyList_SetItem(list, i, PyFloat_FromDouble(val));
    }

    return list;
}

void FromList(PyObject* list) {
    iftHist* hist = ($self);

    free(hist->val);
    hist->nbins = PyList_Size(list);
    hist->val = iftAlloc(hist->nbins, sizeof(float));

    for(int i = 0; i < hist->nbins; i++){
        float val = (float) PyFloat_AsDouble(PyList_GetItem(list, i));
        hist->val[i] = val;
    }
}

PyObject* AsNumPy() {
    iftHist* hist = ($self);

    npy_intp dims[1] = {hist->nbins};
    PyArrayObject* ref = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT32, (float*) hist->val);
    PyArrayObject* npy = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    PyArray_CopyInto(npy, ref);

    return PyArray_Return(npy);
}

void FromNumPy(PyObject* input) {

    PyArrayObject* ary = obj_to_array_no_conversion(input, NPY_FLOAT32);

    if (ary == NULL) return;

    if(PyArray_TYPE(ary) != NPY_FLOAT32){
        SWIG_Error(12, "Input must be a numpy double array");
    }

    iftHist* hist = ($self);

    require_dimensions(ary, 1);

    hist->nbins = array_size(ary, 0);
    hist->val = iftAlloc(hist->nbins, sizeof(float));

    float* data = (float*) array_data(ary);
    memcpy(hist->val, data, hist->nbins * sizeof(float));
}

