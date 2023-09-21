
PyObject *AsDict(void)
{
    iftKnnGraph *graph = ($self);

    PyObject *dict = PyDict_New();
    for (int i = 0; i < graph->nnodes; i++)
    {
        PyObject *list = PyList_New(graph->kmax);
        int j = 0;
        for (iftAdjSet *adj = graph->node[i].adj; adj != NULL; adj = adj->next)
        {
            PyList_SetItem(list, j, PyInt_FromLong((long) adj->node));
            ++j;
        }
        PyDict_SetItem(dict, PyInt_FromLong((long) i), list);
    }

    return dict;
}
