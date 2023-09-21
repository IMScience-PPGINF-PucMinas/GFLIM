//! swig(newobject, stable)
iftVoxel PyTupleToVoxel(PyObject *voxel_tuple) {
    int n_coords = PyTuple_Size(voxel_tuple);
    if (n_coords != 3) {
        char error[512];
        sprintf(error, "Number of Coordinates into Python Tuple to converting into a voxel is %d != 3", n_coords);
        SWIG_Error(0, error);
    }

    iftVoxel v;
    v.x = PyInt_AsLong(PyTuple_GetItem(voxel_tuple, 0));
    v.y = PyInt_AsLong(PyTuple_GetItem(voxel_tuple, 1));
    v.z = PyInt_AsLong(PyTuple_GetItem(voxel_tuple, 2));

    return v;
}

