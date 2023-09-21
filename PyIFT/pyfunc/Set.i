
//! swig(newobject)
iftSet *CreateSetFromList(PyObject *input)
{
    if (!PyList_CheckExact(input))
        SWIG_Error(0, "Input \"input\" must be a list.");

    iftSet *s = NULL;
    for (int i = 0; i < PyList_Size(input); i++) {
        int p = PyLong_AsLong(PyList_GetItem(input, i));
        iftInsertSet(&s, p);
    }
    return s;
}