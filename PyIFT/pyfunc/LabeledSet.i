//! swig(newobject, stable)
iftLabeledSet *CreateLabeledSetFromCoord(PyObject *_index, PyObject *_labels)
{
    PyArrayObject* index = obj_to_array_no_conversion(_index, NPY_INT32);
    PyArrayObject* labels = obj_to_array_no_conversion(_labels, NPY_INT32);

    if (!index || !labels) return NULL;

    if (array_numdims(index) != 1 || array_numdims(labels) != 1) {
        SWIG_Error(0, "Input must be a 1-dimensional array");
        return NULL;
    }

    if (!require_contiguous(labels) || !require_contiguous(index)) {
        SWIG_Error(0, "Inputs array is not contiguous");
        return NULL;
    }

    if (array_size(index, 0) != array_size(labels, 0)) {
        SWIG_Error(0, "Inputs must have the same length");
        return NULL;
    }

    int *d_index = array_data(index);
    int *d_labels = array_data(labels);

    iftLabeledSet *s = NULL;

    int size = array_size(index, 0);
    for (int i = 0; i < size; i++) {
        iftInsertLabeledSet(&s, d_index[i], d_labels[i]);
    }

    return s;
}
