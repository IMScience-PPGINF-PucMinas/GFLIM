PyObject *AsDict(void)
{
    iftCplGraph *graph = ($self);

    PyObject *dict = PyDict_New();
    for (int i = 0; i < graph->nnodes; i++)
    {
        PyObject *list = PyList_New(1);
        PyList_SetItem(list, 0, PyLong_FromLong((long) graph->node[i].pred));
        PyDict_SetItem(dict, PyLong_FromLong((long) i), list);
    }

    return dict;
}
