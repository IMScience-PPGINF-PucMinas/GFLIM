// bb_tuple = ((x0,y0,z0), (x1,y1,z1))
//! swig(newobject, stable)
iftBoundingBox PyTupleToBoundingBox(PyObject *bb_tuple) {
    iftBoundingBox bb = {{0, 0, 0}, {0, 0, 0}};

    int n_voxels = PyTuple_Size(bb_tuple);

    if (n_voxels != 2) {
        char error[512];
        sprintf(error, "Number of Voxels into Python Tuple for creating a Bounding Box %d is != 2", n_voxels);
        SWIG_Error(0, error);
        return bb;
    }


    PyObject *begin = PyTuple_GetItem(bb_tuple, 0);
    int n_coords = PyTuple_Size(begin);
    if (n_coords != 3) {
        char error[512];
        sprintf(error, "Number of Coordinates into Python Tuple for creating the beginning voxel of the Bounding Box %d is != 3", n_coords);
        SWIG_Error(0, error);
        return bb;
    }

    PyObject *end = PyTuple_GetItem(bb_tuple, 1);
    n_coords = PyTuple_Size(end);
    if (n_coords != 3) {
        char error[512];
        sprintf(error, "Number of Coordinates into Python Tuple for creating the ending voxel of the Bounding Box %d is != 3", n_coords);
        SWIG_Error(0, error);
        return bb;
    }

    bb.begin.x = PyInt_AsLong(PyTuple_GetItem(begin, 0));
    bb.begin.y = PyInt_AsLong(PyTuple_GetItem(begin, 1));
    bb.begin.z = PyInt_AsLong(PyTuple_GetItem(begin, 2));
    bb.end.x = PyInt_AsLong(PyTuple_GetItem(end, 0));
    bb.end.y = PyInt_AsLong(PyTuple_GetItem(end, 1));
    bb.end.z = PyInt_AsLong(PyTuple_GetItem(end, 2));

    return bb;
}