void __setitem__(int i, int value){
    if(i >= 0 && i < 3){
        ($self)->val[i] = value;
    }
}

void SetValues(PyObject *tuple){
    iftColor *color = ($self);

    color->val[0] = (int)PyInt_AsLong(PyTuple_GetItem(tuple, 0));
    color->val[1] = (int)PyInt_AsLong(PyTuple_GetItem(tuple, 1));
    color->val[2] = (int)PyInt_AsLong(PyTuple_GetItem(tuple, 2));
}

PyObject* GetValues(){
    iftColor *color = ($self);
    PyObject *tuple = PyTuple_New(3);

    PyTuple_SetItem(tuple, 0, PyInt_FromLong(color->val[0]));
    PyTuple_SetItem(tuple, 1, PyInt_FromLong(color->val[1]));
    PyTuple_SetItem(tuple, 2, PyInt_FromLong(color->val[2]));

    return tuple;
}