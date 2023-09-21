PyObject* __getitem__(int i){
    PyObject *displacement = PyTuple_New(3);

    PyTuple_SetItem(displacement, 0, PyInt_FromLong(($self)->dx[i]));
    PyTuple_SetItem(displacement, 1, PyInt_FromLong(($self)->dy[i]));
    PyTuple_SetItem(displacement, 2, PyInt_FromLong(($self)->dz[i]));

    return displacement;
}