
iftImage *GetLabel(void)
{
    return iftCopyImage(($self)->label);
}


iftImage *GetRoot(void)
{
    return iftCopyImage(($self)->root);
}


iftImage *GetPred(void)
{
    return iftCopyImage(($self)->pred);
}


iftImage *GetOrder(void)
{
    return iftCopyImage(($self)->order);
}


iftImage *GetDelay(void)
{
    return iftCopyImage(($self)->delay);
}


iftFImage *GetCost(bool sqrt = false)
{
    iftFImage *aux = ($self)->cost;
    iftFImage *copy = iftCreateFImage(aux->xsize, aux->ysize, aux->zsize);

    for (int i = 0; i < aux->n; i++)
        copy->val[i] = aux->val[i];

    if (sqrt) {
       for (int i = 0; i < copy->n; i++) {
           copy->val[i] = sqrtf(copy->val[i]);
       }
    }
    return copy;
}


PyObject *GetSets(void)
{
    iftDTForest *forest = ($self);

    PyObject *dict = PyDict_New();
    for (int i = 0; i < forest->label->n; i++) {
        if (forest->dyn_trees[i]) {
            PyObject *val = SWIG_NewPointerObj(SWIG_as_voidptr(forest->dyn_trees[i]->begin),
                                               SWIGTYPE_p_iftSet, 0);
            PyDict_SetItem(dict, PyInt_FromLong(i), val);
        }
    }

    return dict;
}


