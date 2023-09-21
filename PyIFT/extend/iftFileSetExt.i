PyObject *GetPaths(void) {
    iftFileSet *fset = ($self);

    if (fset == NULL) {
        return PyList_New(0); // return an empty List
    }
    else {
        PyObject *list = PyList_New(fset->n);

        for (int i = 0; i < fset->n; i++)
            PyList_SetItem(list, i, PyString_FromString(fset->files[i]->path));
        return list;
    }
}
